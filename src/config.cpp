#include "config.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

int occupied_count(int N, double density) noexcept {
    return static_cast<int>(std::llround(density * static_cast<double>(N) * static_cast<double>(N)));
}

void initialize_grid(HyperuniformConfig& config, int N, double density, unsigned seed) {
    config.N = N;
    config.rng.seed(seed);

    const int total = N * N;
    config.M = occupied_count(N, density);
    if (config.M < 1 || config.M > total - 1)
        throw std::invalid_argument(std::format(
            "density {} gives {} occupied sites out of {}: need at least one occupied and one empty site",
            density, config.M, total));

    // Random start: the first M sites of a shuffled list are occupied.
    std::vector<std::uint32_t> sites(static_cast<std::size_t>(total));
    std::iota(sites.begin(), sites.end(), 0u);
    std::shuffle(sites.begin(), sites.end(), config.rng);

    config.occupied.assign(sites.begin(), sites.begin() + config.M);
    config.empty.assign(sites.begin() + config.M, sites.end());
    config.grid.assign(static_cast<std::size_t>(total), 0);
    config.list_pos.assign(static_cast<std::size_t>(total), 0);
    for (std::size_t i = 0; i < config.occupied.size(); ++i) {
        config.grid[config.occupied[i]] = 1;
        config.list_pos[config.occupied[i]] = static_cast<std::uint32_t>(i);
    }
    for (std::size_t i = 0; i < config.empty.size(); ++i)
        config.list_pos[config.empty[i]] = static_cast<std::uint32_t>(i);
}

double chi_of(std::size_t n_modes, std::size_t n_points) noexcept {
    if (n_points < 2) return 0.0;
    return static_cast<double>(n_modes) / static_cast<double>(n_points - 1);
}

std::size_t chi_points(int N, int M) noexcept {
    const auto total = static_cast<std::size_t>(N) * static_cast<std::size_t>(N);
    const auto m = static_cast<std::size_t>(M);
    return std::min(m, total - m);
}

namespace {

/// Largest kx, |ky| in the mask: strictly below the Nyquist index N/2, where k and -k alias.
int mode_limit(int N) { return N / 2 - 1; }

}  // namespace

MaskChoice choose_mask_for_chi(int N, std::size_t n_points, double chi) {
    if (!(chi > 0.0)) throw std::invalid_argument(std::format("chi must be positive, got {}", chi));
    if (n_points < 2) throw std::invalid_argument("too few points to define chi");

    const int lim = mode_limit(N);
    const int max_k2 = lim * lim;  // largest disc fully inside the Nyquist box
    std::vector<std::size_t> per_shell(static_cast<std::size_t>(max_k2) + 1, 0);
    for (int kx = 0; kx <= lim; ++kx)
        for (int ky = -lim; ky <= lim; ++ky) {
            if (kx == 0 && ky <= 0) continue;
            const int k2 = kx * kx + ky * ky;
            if (k2 <= max_k2) ++per_shell[static_cast<std::size_t>(k2)];
        }

    const double target = chi * static_cast<double>(n_points - 1);
    MaskChoice best;
    double best_err = std::numeric_limits<double>::infinity();
    std::size_t cum = 0;
    for (int k2 = 1; k2 <= max_k2; ++k2) {
        if (per_shell[static_cast<std::size_t>(k2)] == 0) continue;
        cum += per_shell[static_cast<std::size_t>(k2)];
        const double err = std::abs(static_cast<double>(cum) - target);
        if (err < best_err) { best_err = err; best = {k2, cum, chi_of(cum, n_points)}; }
        if (static_cast<double>(cum) > target) break;  // the error can only grow from here
    }
    if (best.n_modes == 0)
        throw std::invalid_argument(std::format("chi = {} constrains no mode for N = {}", chi, N));
    if (static_cast<double>(cum) < target)
        std::cerr << std::format("Warning: chi = {} needs more modes than fit below the Nyquist limit; "
                                 "using the largest mask (chi = {:.4f})\n", chi, best.chi);
    return best;
}

void build_mask(HyperuniformConfig& config, double K2) {
    const int N = config.N;
    const int lim = mode_limit(N);
    if (K2 > static_cast<double>(lim * lim))
        throw std::invalid_argument(std::format(
            "K = {:.4g} reaches the Nyquist limit: K must be at most N/2 - 1 = {}", std::sqrt(K2), lim));
    config.K2 = K2;
    config.K = std::sqrt(std::max(K2, 0.0));

    config.modes.clear();
    for (int kx = 0; kx <= lim; ++kx)
        for (int ky = -lim; ky <= lim; ++ky) {
            if (kx == 0 && ky <= 0) continue;  // half plane; k = 0 excluded
            if (static_cast<double>(kx * kx + ky * ky) <= K2) config.modes.push_back({kx, ky});
        }
    if (config.modes.empty())
        throw std::invalid_argument(std::format("K = {:.4g} contains no mode (K must be >= 1)", config.K));

    const int half = N / 2;
    config.mt.resize(static_cast<std::size_t>(N) * static_cast<std::size_t>(N));
    for (int k = -half; k < half; ++k)
        for (int x = 0; x < N; ++x)
            config.mt[static_cast<std::size_t>((k + half) * N + x)] = static_cast<std::int32_t>(((k * x) % N + N) % N);
    config.w_re.resize(static_cast<std::size_t>(N));
    config.w_im.resize(static_cast<std::size_t>(N));
    for (int j = 0; j < N; ++j) {
        const double theta = 2.0 * PI * j / N;
        config.w_re[static_cast<std::size_t>(j)] = std::cos(theta);
        config.w_im[static_cast<std::size_t>(j)] = -std::sin(theta);
    }

    const std::size_t n = config.modes.size();
    config.kx_off.resize(n);
    config.ky_off.resize(n);
    for (std::size_t m = 0; m < n; ++m) {
        config.kx_off[m] = (config.modes[m].kx + half) * N;
        config.ky_off[m] = (config.modes[m].ky + half) * N;
    }
    config.ft_re.assign(n, 0.0);
    config.ft_im.assign(n, 0.0);
    config.d_re.assign(n, 0.0);
    config.d_im.assign(n, 0.0);

    std::cout << std::format("Modes in mask: {} half-plane ({} counting both ±k)\n", n, 2 * n);
}
