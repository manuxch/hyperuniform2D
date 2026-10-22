#include "thread_pool.hpp"
#include "fourier.hpp"

#include <algorithm>

#if defined(__linux__)
#include <sched.h>
#endif

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define HYPER2D_PAUSE() _mm_pause()
#else
#define HYPER2D_PAUSE() ((void)0)
#endif

void SpinBarrier::arrive_and_wait() noexcept {
    const std::uint64_t gen = generation_.load(std::memory_order_acquire);
    if (waiting_.fetch_add(1, std::memory_order_acq_rel) + 1 == count_) {
        waiting_.store(0, std::memory_order_relaxed);
        generation_.fetch_add(1, std::memory_order_release);  // releases everyone
        return;
    }
    for (unsigned spins = 0; generation_.load(std::memory_order_acquire) == gen; ++spins) {
        if (spins < 4096) HYPER2D_PAUSE();
        else std::this_thread::yield();  // oversubscribed machine: do not starve the other threads
    }
}

int available_cpus() {
#if defined(__linux__)
    // CPUs this process may run on. Under Slurm (cgroups / task affinity) this is the allocation, whereas
    // std::thread::hardware_concurrency() reports every CPU of the node.
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int n = CPU_COUNT(&set);
        if (n > 0) return n;
    }
#endif
    const auto n = static_cast<int>(std::thread::hardware_concurrency());
    return n > 0 ? n : 1;
}

ThreadPool::ThreadPool(int num_threads) {
    // Automatic: at most 4 threads. Measured gains saturate around 2-4 threads (the proposal is short and
    // memory-bound), and more spinning threads, especially on SMT siblings, made it slower.
    if (num_threads <= 0) num_threads = std::min(available_cpus(), 4);
    thread_count = num_threads;

    if (thread_count > 1) {
        // The calling thread takes chunk 0, so only thread_count - 1 workers are started: with spinning
        // barriers, one thread more than cores would be catastrophic.
        start_barrier = std::make_unique<SpinBarrier>(thread_count);
        end_barrier = std::make_unique<SpinBarrier>(thread_count);

        partials.resize(static_cast<std::size_t>(thread_count));

        for (int i = 1; i < thread_count; ++i) {
            workers.emplace_back(&ThreadPool::worker_loop, this, i);
        }
    }
}

ThreadPool::~ThreadPool() {
    if (thread_count > 1) {
        terminate_pool.store(true, std::memory_order_relaxed);
        start_barrier->arrive_and_wait(); // Wake up workers to terminate
        for (auto& worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }
}

double ThreadPool::parallel_reduce(std::size_t n, RangeFn fn, void* ctx) {
    if (thread_count <= 1) return fn(ctx, 0, n);

    // Set shared state for workers (published by the barrier)
    current_fn = fn;
    current_ctx = ctx;
    current_n = n;

    start_barrier->arrive_and_wait(); // Signal workers to start
    run_chunk(0);                     // The calling thread does chunk 0
    end_barrier->arrive_and_wait();   // Wait for workers to finish

    double total = 0.0;
    for (const auto& p : partials) total += p.value;
    return total;
}

double ThreadPool::parallel_reduce_delta_E(HyperuniformConfig& config,
                                           int xi, int yi, int xj, int yj) {
    auto range = [&](std::size_t begin, std::size_t end) {
        return delta_energy_range(config, xi, yi, xj, yj, begin, end);
    };
    return parallel_reduce(config.modes.size(), range);
}

void ThreadPool::run_chunk(int thread_id) {
    const auto threads = static_cast<std::size_t>(thread_count);
    const auto id = static_cast<std::size_t>(thread_id);
    const std::size_t chunk = (current_n + threads - 1) / threads;
    const std::size_t begin = std::min(id * chunk, current_n);
    const std::size_t end = std::min(begin + chunk, current_n);
    partials[id].value = current_fn(current_ctx, begin, end);
}

void ThreadPool::worker_loop(int thread_id) {
    while (true) {
        start_barrier->arrive_and_wait(); // Wait for work
        if (terminate_pool.load(std::memory_order_relaxed)) break;
        run_chunk(thread_id);
        end_barrier->arrive_and_wait(); // Signal work done
    }
}
