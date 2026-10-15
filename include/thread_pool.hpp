#pragma once

#include "config.hpp"

#include <atomic>
#include <cstdint>
#include <memory>
#include <thread>
#include <vector>

/// CPUs available to this process (affinity mask on Linux, e.g. the Slurm allocation), at least 1.
[[nodiscard]] int available_cpus();

/// Reusable barrier that spins instead of sleeping. Each annealing proposal crosses two barriers and
/// lasts only microseconds, so the futex wake-up of std::barrier (tens of µs) dominated the cost.
/// Spinning keeps every worker core busy for the whole run: use the pool only for large masks.
class SpinBarrier {
public:
    explicit SpinBarrier(int count) : count_(count), waiting_(0), generation_(0) {}
    void arrive_and_wait() noexcept;

private:
    const int count_;
    alignas(64) std::atomic<int> waiting_;
    alignas(64) std::atomic<std::uint64_t> generation_;
};

/// Lightweight thread pool optimized for synchronized parallel loops
class ThreadPool {
public:
    /// Creates the pool with the specified number of threads.
    /// If num_threads == 0, it uses min(available_cpus(), 4).
    explicit ThreadPool(int num_threads);
    ~ThreadPool();

    /// Not copyable/movable
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /// Splits the modes [0, |modes|) among the threads; each one fills its part of d_re / d_im
    /// (delta_energy_range). Returns the total energy change of moving (xi, yi) -> (xj, yj).
    double parallel_reduce_delta_E(HyperuniformConfig& config,
                                   int xi, int yi, int xj, int yj);

    [[nodiscard]] int get_thread_count() const noexcept { return thread_count; }

private:
    void worker_loop(int thread_id);
    void run_chunk(int thread_id);  ///< delta_energy_range on this thread's share of the modes

    int thread_count;
    std::vector<std::thread> workers;

    // Synchronization barriers
    std::unique_ptr<SpinBarrier> start_barrier;
    std::unique_ptr<SpinBarrier> end_barrier;

    // Shared state for the current job
    std::atomic<bool> terminate_pool{false};
    HyperuniformConfig* current_config = nullptr;
    int current_xi = 0, current_yi = 0;
    int current_xj = 0, current_yj = 0;

    // Partial sums calculated by each thread, one cache line each (no false sharing)
    struct alignas(64) Partial { double value = 0.0; };
    std::vector<Partial> partial_deltaE;
};
