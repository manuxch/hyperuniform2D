#include "thread_pool.hpp"
#include "constants.hpp"

#include <numeric>

ThreadPool::ThreadPool(int num_threads) {
    if (num_threads <= 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;
    }
    thread_count = num_threads;

    if (thread_count > 1) {
        // barrier for main thread + worker threads
        start_barrier = std::make_unique<std::barrier<>>(thread_count + 1);
        end_barrier = std::make_unique<std::barrier<>>(thread_count + 1);

        partial_deltaE.resize(thread_count, 0.0);

        for (int i = 0; i < thread_count; ++i) {
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

double ThreadPool::parallel_reduce_delta_E(HyperuniformConfig& config, 
                                           int xi, int yi, int xj, int yj) {
    // If sequential
    if (thread_count <= 1) {
        const int N = config.N;
        const auto num_modes = config.modes.size();
        
        for (std::size_t m = 0; m < num_modes; ++m) {
            const auto [kx, ky] = config.modes[m];
            const double theta_j = -2.0 * PI * (kx * xj + ky * yj) / N;
            const double theta_i = -2.0 * PI * (kx * xi + ky * yi) / N;
            config.delta_buffer[m] = std::polar(1.0, theta_j) - std::polar(1.0, theta_i);
        }

        double deltaE = 0.0;
        for (std::size_t m = 0; m < num_modes; ++m) {
            deltaE += 2.0 * std::real(std::conj(config.FT[m]) * config.delta_buffer[m])
                    + std::norm(config.delta_buffer[m]);
        }
        return deltaE;
    }

    // Set shared state for workers
    current_config = &config;
    current_xi = xi; current_yi = yi;
    current_xj = xj; current_yj = yj;

    // Signal workers to start
    start_barrier->arrive_and_wait();

    // Wait for workers to finish
    end_barrier->arrive_and_wait();

    // Sum partial results
    return std::accumulate(partial_deltaE.begin(), partial_deltaE.end(), 0.0);
}

void ThreadPool::worker_loop(int thread_id) {
    while (true) {
        start_barrier->arrive_and_wait(); // Wait for work
        if (terminate_pool.load(std::memory_order_relaxed)) break;

        const int N = current_config->N;
        const auto num_modes = current_config->modes.size();
        
        // Determine chunk for this thread
        const std::size_t chunk_size = (num_modes + thread_count - 1) / thread_count;
        const std::size_t start_idx = thread_id * chunk_size;
        const std::size_t end_idx = std::min(start_idx + chunk_size, num_modes);

        double local_deltaE = 0.0;

        for (std::size_t m = start_idx; m < end_idx; ++m) {
            const auto [kx, ky] = current_config->modes[m];
            const double theta_j = -2.0 * PI * (kx * current_xj + ky * current_yj) / N;
            const double theta_i = -2.0 * PI * (kx * current_xi + ky * current_yi) / N;
            current_config->delta_buffer[m] = std::polar(1.0, theta_j) - std::polar(1.0, theta_i);
        }

        for (std::size_t m = start_idx; m < end_idx; ++m) {
            local_deltaE += 2.0 * std::real(std::conj(current_config->FT[m]) * current_config->delta_buffer[m])
                          + std::norm(current_config->delta_buffer[m]);
        }

        partial_deltaE[thread_id] = local_deltaE;

        end_barrier->arrive_and_wait(); // Signal work done
    }
}
