#pragma once

#include "config.hpp"

#include <atomic>
#include <barrier>
#include <functional>
#include <memory>
#include <thread>
#include <vector>

/// Lightweight thread pool optimized for synchronized parallel loops
class ThreadPool {
public:
    /// Creates the pool with the specified number of threads. 
    /// If num_threads == 0, it auto-detects hardware concurrency.
    explicit ThreadPool(int num_threads);
    ~ThreadPool();

    /// Not copyable/movable
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /// Dispatches a loop over modes [0, total_modes) to the worker threads.
    /// The threads will evaluate the delta_buffer and partial deltaE.
    /// Returns the total summed deltaE.
    double parallel_reduce_delta_E(HyperuniformConfig& config, 
                                   int xi, int yi, int xj, int yj);

    [[nodiscard]] int get_thread_count() const noexcept { return thread_count; }

private:
    void worker_loop(int thread_id);

    int thread_count;
    std::vector<std::thread> workers;

    // Synchronization barriers
    std::unique_ptr<std::barrier<>> start_barrier;
    std::unique_ptr<std::barrier<>> end_barrier;

    // Shared state for the current job
    std::atomic<bool> terminate_pool{false};
    HyperuniformConfig* current_config = nullptr;
    int current_xi = 0, current_yi = 0;
    int current_xj = 0, current_yj = 0;

    // Partial sums calculated by each thread
    std::vector<double> partial_deltaE;
};
