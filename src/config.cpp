#include "config.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

int occupied_count(int Nx, int Ny, double density) noexcept {
    return static_cast<int>(std::llround(density * static_cast<double>(Nx) * static_cast<double>(Ny)));
}

void initialize_grid(HyperuniformConfig& config, int Nx, int Ny, double density, unsigned seed,
                     double min_distance) {
    config.Nx = Nx;
    config.Ny = Ny;
    config.rng.seed(seed);

    const int total = Nx * Ny;
    config.M = occupied_count(Nx, Ny, density);
    if (config.M < 1 || config.M > total - 1)
        throw std::invalid_argument(std::format(
            "density {} gives {} occupied sites out of {}: need at least one occupied and one empty site",
            density, config.M, total));

    // Random start: the first M sites of a shuffled list are occupied.
    std::vector<std::uint32_t> sites(static_cast<std::size_t>(total));
    std::iota(sites.begin(), sites.end(), 0u);
    std::shuffle(sites.begin(), sites.end(), config.rng);

    config.grid.assign(static_cast<std::size_t>(total), 0);
    config.list_pos.assign(static_cast<std::size_t>(total), 0);
    config.min_distance = min_distance;
    config.excl.clear();
    config.hc_rejections = 0;
    if (min_distance > 0.0) {
        const int r = static_cast<int>(std::ceil(min_distance));
        if (2 * r >= std::min(Nx, Ny))
            throw std::invalid_argument(std::format("min_distance = {} must be below min(Nx, Ny)/2", min_distance));
        for (int dy = -r; dy <= r; ++dy)
            for (int dx = -r; dx <= r; ++dx) {
                const int d2 = dx * dx + dy * dy;
                if (d2 > 0 && static_cast<double>(d2) < min_distance * min_distance) config.excl.push_back({dx, dy});
            }
        // Random sequential addition in the shuffled order: keep a site if it respects the exclusion.
        std::vector<std::uint32_t> taken, rest;
        taken.reserve(static_cast<std::size_t>(config.M));
        for (const std::uint32_t s : sites) {
            if (static_cast<int>(taken.size()) < config.M && !hardcore_blocked(config, s, s)) {
                config.grid[s] = 1;
                taken.push_back(s);
            } else {
                rest.push_back(s);
            }
        }
        if (static_cast<int>(taken.size()) < config.M)
            throw std::invalid_argument(std::format(
                "min_distance = {} leaves room for only {} of {} occupied sites in a random start",
                min_distance, taken.size(), config.M));
        std::copy(taken.begin(), taken.end(), sites.begin());
        std::copy(rest.begin(), rest.end(), sites.begin() + config.M);
        std::fill(config.grid.begin(), config.grid.end(), std::uint8_t{0});
    }
    config.occupied.assign(sites.begin(), sites.begin() + config.M);
    config.empty.assign(sites.begin() + config.M, sites.end());
    for (std::size_t i = 0; i < config.occupied.size(); ++i) {
        config.grid[config.occupied[i]] = 1;
        config.list_pos[config.occupied[i]] = static_cast<std::uint32_t>(i);
    }
    for (std::size_t i = 0; i < config.empty.size(); ++i)
        config.list_pos[config.empty[i]] = static_cast<std::uint32_t>(i);
}

double min_pair_distance(const HyperuniformConfig& c, double r_max) {
    const int r = static_cast<int>(std::ceil(r_max));
    double best2 = r_max * r_max;
    for (const std::uint32_t s : c.occupied) {
        const int x = static_cast<int>(s % static_cast<std::uint32_t>(c.Nx)), y = static_cast<int>(s / static_cast<std::uint32_t>(c.Nx));
        for (int dy = -r; dy <= r; ++dy)
            for (int dx = -r; dx <= r; ++dx) {
                const int d2 = dx * dx + dy * dy;
                if (d2 == 0 || static_cast<double>(d2) >= best2) continue;
                if (c.grid[static_cast<std::size_t>(periodic_idx(x + dx, y + dy, c.Nx, c.Ny))]) best2 = d2;
            }
    }
    return std::sqrt(best2);
}

double chi_of(std::size_t n_modes, std::size_t n_points) noexcept {
    if (n_points < 2) return 0.0;
    return static_cast<double>(n_modes) / static_cast<double>(n_points - 1);
}

double chi_point_pattern(std::size_t n_modes, std::size_t n_points) noexcept {
    if (n_points < 2) return 0.0;
    return static_cast<double>(n_modes) / (2.0 * static_cast<double>(n_points - 1));
}

std::size_t chi_points(int Nx, int Ny, int M) noexcept {
    const auto total = static_cast<std::size_t>(Nx) * static_cast<std::size_t>(Ny);
    const auto m = static_cast<std::size_t>(M);
    return std::min(m, total - m);
}

namespace {

/// Largest |kx| (or |ky|) in the mask: strictly below the Nyquist index N/2, where k and -k alias.
int mode_limit(int N) { return N / 2 - 1; }

/// Largest q whose circle |k|² <= q fits inside the Nyquist box on both axes.
std::int64_t q_limit(int Nx, int Ny) {
    return std::min(mode_key(mode_limit(Nx), 0, Nx, Ny), mode_key(0, mode_limit(Ny), Nx, Ny));
}

/// floor(K2 Nx Ny): the largest q inside |k|² <= K2 (K2 in units of (2π/sqrt(Nx Ny))²).
std::int64_t q_of_K2(double K2, int Nx, int Ny) {
    const long double q = static_cast<long double>(K2) * static_cast<long double>(Nx) * static_cast<long double>(Ny);
    return q < 0 ? -1 : static_cast<std::int64_t>(std::floor(q));
}

double K2_of_q(std::int64_t q, int Nx, int Ny) {
    return static_cast<double>(q) / (static_cast<double>(Nx) * static_cast<double>(Ny));
}

/// Calls f(kx, ky, q) for every half-plane mode of the Nyquist box (kx > 0, or kx = 0 and ky > 0), in the
/// order kx = 0..limx, ky = -limy..limy (the order of the stored modes).
template <class F>
void for_each_box_mode(int Nx, int Ny, F&& f) {
    const int limx = mode_limit(Nx), limy = mode_limit(Ny);
    for (int kx = 0; kx <= limx; ++kx)
        for (int ky = -limy; ky <= limy; ++ky) {
            if (kx == 0 && ky <= 0) continue;  // half plane; k = 0 excluded
            f(kx, ky, mode_key(kx, ky, Nx, Ny));
        }
}

}  // namespace

double max_mask_K2(int Nx, int Ny) noexcept {
    return K2_of_q(q_limit(Nx, Ny), Nx, Ny);
}

MaskChoice choose_mask_for_chi(int Nx, int Ny, std::size_t n_points, double chi) {
    if (!(chi > 0.0)) throw std::invalid_argument(std::format("chi must be positive, got {}", chi));
    if (n_points < 2) throw std::invalid_argument("too few points to define chi");

    // Shells of the circle in physical k: modes with the same key q.
    const std::int64_t qlim = q_limit(Nx, Ny);
    std::vector<std::int64_t> keys;
    for_each_box_mode(Nx, Ny, [&](int, int, std::int64_t q) { if (q <= qlim) keys.push_back(q); });
    std::sort(keys.begin(), keys.end());

    const double target = chi * static_cast<double>(n_points - 1);
    MaskChoice best;
    double best_err = std::numeric_limits<double>::infinity();
    std::size_t cum = 0;
    for (std::size_t i = 0; i < keys.size();) {
        const std::int64_t q = keys[i];
        while (i < keys.size() && keys[i] == q) { ++cum; ++i; }
        const double err = std::abs(static_cast<double>(cum) - target);
        if (err < best_err) { best_err = err; best = {K2_of_q(q, Nx, Ny), q, cum, chi_of(cum, n_points)}; }
        if (static_cast<double>(cum) > target) break;  // the error can only grow from here
    }
    if (best.n_modes == 0)
        throw std::invalid_argument(std::format("chi = {} constrains no mode for {} x {}", chi, Nx, Ny));
    if (static_cast<double>(cum) < target)
        std::cerr << std::format("Warning: chi = {} needs more modes than fit below the Nyquist limit; "
                                 "using the largest mask (chi = {:.4f})\n", chi, best.chi);
    return best;
}

void build_mask(HyperuniformConfig& config, double K2, double weight_power) {
    const int Nx = config.Nx, Ny = config.Ny;
    const std::int64_t qlim = q_limit(Nx, Ny);
    const std::int64_t q_max = q_of_K2(K2, Nx, Ny);
    if (q_max > qlim) {
        if (config.square())
            throw std::invalid_argument(std::format(
                "K = {:.4g} reaches the Nyquist limit: K must be at most N/2 - 1 = {}", std::sqrt(K2), mode_limit(Nx)));
        throw std::invalid_argument(std::format(
            "K = {:.4g} (units of 2pi/sqrt(Nx Ny)) reaches the Nyquist limit of the {} x {} lattice: "
            "K must be at most {:.6g}", std::sqrt(K2), Nx, Ny, std::sqrt(max_mask_K2(Nx, Ny))));
    }
    const std::int64_t L = std::lcm(static_cast<std::int64_t>(Nx), static_cast<std::int64_t>(Ny));
    if (L > (std::int64_t{1} << 24))
        throw std::invalid_argument(std::format(
            "lcm(Nx, Ny) = {} is too large for the phase table (at most 2^24): choose Nx and Ny with a larger "
            "common factor", L));
    config.L = static_cast<int>(L);
    config.K2 = K2;
    config.K = std::sqrt(std::max(K2, 0.0));
    config.q_max = q_max;

    config.modes.clear();
    config.mode_k2.clear();
    config.q_next = 0;
    for_each_box_mode(Nx, Ny, [&](int kx, int ky, std::int64_t q) {
        if (q <= q_max) {
            config.modes.push_back({kx, ky});
            config.mode_k2.push_back(K2_of_q(q, Nx, Ny));
        } else if (q <= qlim && (config.q_next == 0 || q < config.q_next)) {
            config.q_next = q;
        }
    });
    if (config.modes.empty())
        throw std::invalid_argument(std::format("K = {:.4g} contains no mode (K must be >= 1)", config.K));

    // Axis tables: (k x L/N) mod L. For Nx = Ny = N, L = N and these are the 2.x values.
    const auto axis_table = [L](int N) {
        const int half = N / 2;
        const std::int64_t step = L / N;
        std::vector<std::int32_t> t(static_cast<std::size_t>(N) * static_cast<std::size_t>(N));
        for (int k = -half; k < half; ++k)
            for (int x = 0; x < N; ++x)
                t[static_cast<std::size_t>((k + half) * N + x)] =
                    static_cast<std::int32_t>(((static_cast<std::int64_t>(k) * x * step) % L + L) % L);
        return t;
    };
    config.mtx = axis_table(Nx);
    config.mty = axis_table(Ny);
    config.w_re.resize(static_cast<std::size_t>(L));
    config.w_im.resize(static_cast<std::size_t>(L));
    for (int j = 0; j < config.L; ++j) {
        const double theta = 2.0 * PI * j / config.L;
        config.w_re[static_cast<std::size_t>(j)] = std::cos(theta);
        config.w_im[static_cast<std::size_t>(j)] = -std::sin(theta);
    }

    const std::size_t n = config.modes.size();
    config.kx_off.resize(n);
    config.ky_off.resize(n);
    for (std::size_t m = 0; m < n; ++m) {
        config.kx_off[m] = (config.modes[m].kx + Nx / 2) * Nx;
        config.ky_off[m] = (config.modes[m].ky + Ny / 2) * Ny;
    }
    config.weight_power = weight_power;
    config.weight.assign(n, 1.0);
    if (weight_power != 0.0) {
        double sum = 0.0;
        for (std::size_t m = 0; m < n; ++m) {
            config.weight[m] = std::pow(config.mode_k2[m], -0.5 * weight_power);
            sum += config.weight[m];
        }
        const double mean = sum / static_cast<double>(n);
        for (double& w : config.weight) w /= mean;
    }
    config.ft_re.assign(n, 0.0);
    config.ft_im.assign(n, 0.0);
    config.d_re.assign(n, 0.0);
    config.d_im.assign(n, 0.0);

    std::cout << std::format("Modes in mask: {} half-plane ({} counting both ±k)\n", n, 2 * n);
}
