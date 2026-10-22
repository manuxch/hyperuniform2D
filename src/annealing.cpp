#include "annealing.hpp"
#include "fourier.hpp"
#include "constants.hpp"

#include "thread_pool.hpp"

#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>

namespace {

struct Proposal {
    std::size_t a = 0, b = 0;  ///< Positions in config.occupied / config.empty
    double dE = 0.0;
    bool valid = true;         ///< false: the move breaks the hard-core exclusion (dE not computed)
};

/// Moves a random occupied site to an empty site: with probability local_fraction a random empty site
/// within Chebyshev distance local_radius (a few draws; a site with no empty neighbour falls back to a
/// global move), otherwise any empty site. Fills d_re / d_im and returns the energy change.
Proposal propose(HyperuniformConfig& config, ThreadPool* pool, double local_fraction, int local_radius) {
    const int Nx = config.Nx, Ny = config.Ny;
    Proposal p;
    p.a = std::uniform_int_distribution<std::size_t>(0, config.occupied.size() - 1)(config.rng);
    const auto i = config.occupied[p.a];
    const int xi = static_cast<int>(i % static_cast<std::uint32_t>(Nx)), yi = static_cast<int>(i / static_cast<std::uint32_t>(Nx));
    bool found = false;
    if (local_fraction > 0.0 && config.prob_dist(config.rng) < local_fraction) {
        std::uniform_int_distribution<int> step(-local_radius, local_radius);
        for (int tries = 0; tries < 8 && !found; ++tries) {
            const int dx = step(config.rng), dy = step(config.rng);
            if (dx == 0 && dy == 0) continue;
            const auto j = static_cast<std::uint32_t>(periodic_idx(xi + dx, yi + dy, Nx, Ny));
            if (config.grid[j] == 0) { p.b = config.list_pos[j]; found = true; }
        }
    }
    if (!found) p.b = std::uniform_int_distribution<std::size_t>(0, config.empty.size() - 1)(config.rng);
    const auto j = config.empty[p.b];
    if (!config.excl.empty() && hardcore_blocked(config, j, i)) {
        p.valid = false;
        ++config.hc_rejections;
        return p;
    }
    const int xj = static_cast<int>(j % static_cast<std::uint32_t>(Nx)), yj = static_cast<int>(j / static_cast<std::uint32_t>(Nx));
    p.dE = pool ? pool->parallel_reduce_delta_E(config, xi, yi, xj, yj)
                : delta_energy_range(config, xi, yi, xj, yj, 0, config.modes.size());
    return p;
}

void accept(HyperuniformConfig& config, const Proposal& p) {
    const auto i = config.occupied[p.a];
    const auto j = config.empty[p.b];
    config.occupied[p.a] = j;
    config.empty[p.b] = i;
    config.grid[i] = 0;
    config.grid[j] = 1;
    config.list_pos[j] = static_cast<std::uint32_t>(p.a);
    config.list_pos[i] = static_cast<std::uint32_t>(p.b);
    for (std::size_t m = 0; m < config.modes.size(); ++m) {
        config.ft_re[m] += config.d_re[m];
        config.ft_im[m] += config.d_im[m];
    }
}

}  // namespace

std::optional<double> attempt_swap(HyperuniformConfig& config, double T, ThreadPool* pool,
                                   double local_fraction, int local_radius) {
    const Proposal p = propose(config, pool, local_fraction, local_radius);
    if (!p.valid) return std::nullopt;

    // Metropolis criterion
    const bool accepted = (p.dE < 0.0)
        || (T > 0.0 && config.prob_dist(config.rng) < std::exp(-p.dE / T));

    if (!accepted) return std::nullopt;
    accept(config, p);
    return p.dE; // Negative for energetically favorable swaps
}

AnnealReport simulated_annealing(HyperuniformConfig& config, const Parameters& p) {
    AnnealReport rep;
    const auto num_modes = static_cast<long long>(config.modes.size());

    std::unique_ptr<ThreadPool> pool;
    if (p.num_threads != 1 && num_modes >= p.parallel_min_modes) {
        pool = std::make_unique<ThreadPool>(p.num_threads);
        if (pool->get_thread_count() <= 1) pool.reset();
    }
    if (pool)
        std::cout << std::format("Initialized thread pool with {} threads.\n", pool->get_thread_count());
    else
        std::cout << std::format("Using sequential mode (num_threads = 1 or fewer than {} half-plane modes).\n",
                                 p.parallel_min_modes);

    // Incremental accumulator: avoids recalculating the energy O(|modes|) after every swap
    calculate_FT(config);
    double E_current = weighted_energy(config);  // the annealed quantity
    rep.E_initial = energy(config);
    rep.Ew_initial = E_current;

    long long since_resync = 0;
    const auto resync = [&] {
        calculate_FT(config);
        const double exact = weighted_energy(config);
        rep.max_drift = std::max(rep.max_drift, std::abs(E_current - exact));
        E_current = exact;
        since_resync = 0;
    };

    // Automatic initial temperature: a typical uphill swap from the random start is accepted with prob. 1/2.
    double T = p.T_initial;
    if (T <= 0.0) {
        double sum = 0.0;
        int n_up = 0;
        for (int s = 0; s < 400; ++s) {
            const Proposal q = propose(config, pool.get(), 0.0, 1);  // global moves, evaluated, not applied
            if (q.valid && q.dE > 0.0) { sum += q.dE; ++n_up; }
        }
        T = n_up > 0 ? (sum / n_up) / std::numbers::ln2 : 1.0;
        if (T <= p.T_final) T = 10.0 * p.T_final;
        std::cout << std::format("Automatic T_initial = {:.6g}\n", T);
    }
    rep.T_initial = T;

    // Open monitoring file
    std::ofstream fann(p.annealing_file);
    if (!fann) {
        std::cerr << std::format("Error: could not create '{}'\n", p.annealing_file);
    } else {
        // Header: columns aligned with data format
        fann << std::format("{:<16} {:<20} {:<14} {:<14}\n", "# T", config.weight_power != 0.0 ? "E_w" : "E", "acc/rate", "acc/level");
    }

    std::cout << std::format("Starting simulated annealing → '{}'\n", p.annealing_file);
    while (T > p.T_final && E_current > p.energy_tol) {
        long long level_accepted = 0;
        for (int step = 0; step < p.steps_per_temp; ++step) {
            ++rep.attempts;
            if (auto dE = attempt_swap(config, T, pool.get(), p.local_fraction, p.local_radius)) {
                ++rep.accepted;
                ++level_accepted;
                E_current += *dE;
                if (++since_resync >= p.resync_interval) resync();
            }
        }
        ++rep.levels;
        if (fann) {
            fann << std::format("{:<16.8g} {:<20.10g} {:<14.8f} {:<14.8f}\n", T, E_current,
                                static_cast<double>(rep.accepted) / static_cast<double>(rep.attempts),
                                static_cast<double>(level_accepted) / p.steps_per_temp);
        }
        rep.T_last = T;
        T *= p.cooling_rate;
    }

    // Zero-temperature quench: only strictly downhill swaps.
    for (long long s = 0; s < p.quench_steps && E_current > p.energy_tol; ++s) {
        ++rep.attempts;
        if (auto dE = attempt_swap(config, 0.0, pool.get(), p.local_fraction, p.local_radius)) {
            ++rep.accepted;
            E_current += *dE;
            if (++since_resync >= p.resync_interval) resync();
        }
    }

    resync();  // exact final energy, and the last drift check
    rep.Ew_final = E_current;
    rep.E_final = energy(config);
    rep.reached_tol = rep.Ew_final <= p.energy_tol;
    std::cout << std::format("Annealing finished: {} levels, {} proposals, acceptance {:.4f}. Monitoring saved to '{}'\n",
                             rep.levels, rep.attempts,
                             rep.attempts ? static_cast<double>(rep.accepted) / static_cast<double>(rep.attempts) : 0.0,
                             p.annealing_file);
    return rep;
}
