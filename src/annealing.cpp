#include "annealing.hpp"
#include "fourier.hpp"
#include "constants.hpp"

#include <cmath>
#include <complex>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>

std::optional<double> attempt_swap(HyperuniformConfig& config, double T) {
    const int N = config.N;

    // Distributions stored in config — not recreated per call
    const int i = config.site_dist(config.rng);
    const int j = config.site_dist(config.rng);

    if (i == j)                                      return std::nullopt;
    if (config.grid[i] != 1 || config.grid[j] != 0) return std::nullopt;

    const int xi = i % N, yi = i / N;
    const int xj = j % N, yj = j / N;
    const auto num_modes = config.modes.size();

    // Use pre-allocated delta_buffer — no new/delete per call
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

    // Metropolis criterion
    const bool accepted = (deltaE < 0.0)
        || (T > 0.0 && config.prob_dist(config.rng) < std::exp(-deltaE / T));

    if (!accepted) return std::nullopt;

    config.grid[i] = 0;
    config.grid[j] = 1;
    for (std::size_t m = 0; m < num_modes; ++m) {
        config.FT[m] += config.delta_buffer[m];
    }
    return deltaE; // Negative for energetically favorable swaps
}

void simulated_annealing(HyperuniformConfig& config, const Parameters& p) {
    double       T            = p.T_initial;
    const int    steps        = p.steps_per_temp;
    const double cooling_rate = p.cooling_rate;
    const double T_final      = p.T_final;

    long long attempts = 0;
    long long accepted = 0;

    // Incremental accumulator: avoids recalculating energy O(|modes|) at each T step
    double E_current = energy(config);

    // Open monitoring file
    std::ofstream fann(p.annealing_file);
    if (!fann) {
        std::cerr << std::format("Error: could not create '{}'\n", p.annealing_file);
    } else {
        // Header: columns aligned with data format
        fann << std::format("{:<16} {:<20} {:<14}\n", "# T", "E", "acc/rate");
    }

    std::cout << std::format("Starting simulated annealing → '{}'\n", p.annealing_file);
    while (T > T_final) {
        for (int step = 0; step < steps; ++step) {
            ++attempts;
            if (auto dE = attempt_swap(config, T)) {
                ++accepted;
                E_current += *dE;
            }
        }
        if (fann) {
            const double rate = static_cast<double>(accepted) / attempts;
            fann << std::format(
                "{:<16.8f} {:<20.8f} {:<14.8f}\n",
                T, E_current, rate);
        }
        T *= cooling_rate;
    }
    std::cout << std::format("Annealing finished. Monitoring saved to '{}'\n",
                             p.annealing_file);
}

