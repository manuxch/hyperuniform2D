#include "config.hpp"

#include <algorithm>
#include <format>
#include <iostream>

void initialize_grid(HyperuniformConfig& config, int N, double density, int seed) {
    config.N = N;
    config.rng.seed(static_cast<std::mt19937::result_type>(seed));

    const int total = N * N;
    config.M = static_cast<int>(density * total);

    config.grid.assign(total, 0);
    std::fill(config.grid.begin(), config.grid.begin() + config.M, 1);
    std::shuffle(config.grid.begin(), config.grid.end(), config.rng);

    // Initialize uniform distribution of sites (member of struct — not recreated per call)
    config.site_dist = std::uniform_int_distribution<int>(0, total - 1);
}

void build_mask(HyperuniformConfig& config, double K) {
    config.K = K;
    const int    N     = config.N;
    const int    limit = N / 2;
    const double K2    = K * K; // Compare squares: avoids sqrt() in the loop

    config.modes.clear();
    for (int kx = -limit; kx < limit; ++kx) {
        for (int ky = -limit; ky < limit; ++ky) {
            if (kx == 0 && ky == 0) continue;
            if (static_cast<double>(kx * kx + ky * ky) <= K2) {
                config.modes.push_back({kx, ky});
            }
        }
    }

    config.FT.assign(config.modes.size(), Complex{0.0, 0.0});

    // Pre-allocate buffer for deltas — reused in each call to attempt_swap()
    config.delta_buffer.resize(config.modes.size());

    std::cout << std::format("Modes in mask: {}\n", config.modes.size());
}
