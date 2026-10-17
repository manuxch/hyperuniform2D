#pragma once

#include "constants.hpp"

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

/// Wave vector index pair, in units of 2π/N
struct ModeK { int kx, ky; };

/// Complete state of the hyperuniform configuration.
///
/// Only the half plane of k is stored (kx > 0, or kx = 0 and ky > 0): for a real occupation
/// rho_hat(-k) = conj(rho_hat(k)), so the other half carries no information. The energy still
/// counts both members of each ±k pair, E = 2 Σ_half |rho_hat(k)|², i.e. the sum over the full mask.
///
/// Phase factors are table look-ups, with no trigonometry in the annealing loop:
///   exp(-i 2π (kx x + ky y)/N) = w[(kx x + ky y) mod N],   (k x) mod N = mt[(k + N/2) N + x].
struct HyperuniformConfig {
    int    N  = 0;    ///< Grid side length
    int    M  = 0;    ///< Number of occupied sites
    double K2 = 0.0;  ///< Squared mask radius, units of (2π/N)²
    double K  = 0.0;  ///< Mask radius, sqrt(K2)

    std::vector<std::uint8_t>  grid;      ///< Grid occupation 0/1, index x + N y
    std::vector<std::uint32_t> occupied;  ///< Sites with grid = 1
    std::vector<std::uint32_t> empty;     ///< Sites with grid = 0
    std::vector<std::uint32_t> list_pos;  ///< Position of each site in its list (occupied or empty)

    std::vector<ModeK>        modes;      ///< Half-plane wave vectors in the mask
    std::vector<std::int32_t> mt;         ///< (k x) mod N for k in [-N/2, N/2), x in [0, N)
    std::vector<std::int32_t> kx_off;     ///< (kx + N/2) N for each mode: offset into mt
    std::vector<std::int32_t> ky_off;     ///< (ky + N/2) N for each mode
    std::vector<double>       w_re, w_im; ///< exp(-i 2π j / N), j in [0, N)

    std::vector<double> weight;           ///< w(k) of each mode in the annealing energy (all 1 if p = 0)
    double weight_power = 0.0;            ///< p in w(k) ∝ |k|^(-p)
    std::vector<double> ft_re, ft_im;     ///< rho_hat(k) for each mode
    std::vector<double> d_re, d_im;       ///< Reusable buffer: change of rho_hat for the proposed swap

    std::mt19937 rng;
    std::uniform_real_distribution<double> prob_dist{0.0, 1.0};
};

/// Index of the phase of mode m at site (x, y) in w_re / w_im.
[[nodiscard]] inline std::size_t phase_index(const HyperuniformConfig& c, std::size_t m, int x, int y) noexcept {
    std::int32_t j = c.mt[static_cast<std::size_t>(c.kx_off[m] + x)] + c.mt[static_cast<std::size_t>(c.ky_off[m] + y)];
    if (j >= c.N) j -= c.N;
    return static_cast<std::size_t>(j);
}

/// Direct index — only for coordinates already in [0, N). No modulo.
[[nodiscard]] inline int flat_idx(int x, int y, int N) noexcept {
    return x + N * y;
}

/// Index with periodic boundary conditions — for arbitrary coordinates.
[[nodiscard]] inline int periodic_idx(int x, int y, int N) noexcept {
    x = (x % N + N) % N;
    y = (y % N + N) % N;
    return x + N * y;
}

/// Number of occupied sites for a target density: round(density N²), not truncation
/// (0.29 * 100² = 2899.9999... in floating point).
[[nodiscard]] int occupied_count(int N, double density) noexcept;

/// Initializes the grid with M = occupied_count(N, density) randomly distributed particles.
/// Throws std::invalid_argument unless there is at least one occupied and one empty site.
void initialize_grid(HyperuniformConfig& config, int N, double density, unsigned seed);

/// Stealthiness chi = (half-plane constrained modes) / (n - 1), with n = min(M, N² - M).
/// It is the ratio of independent real constraints (2 per ±k pair) to real degrees of freedom
/// (2 per point minus the global translation): the χ of Torquato and Stillinger for d = 2.
/// The minority species is used because for k ≠ 0 the two species have the same |rho_hat|.
[[nodiscard]] double chi_of(std::size_t n_modes, std::size_t n_points) noexcept;

/// n = min(M, N² - M), the number of points that enters chi.
[[nodiscard]] std::size_t chi_points(int N, int M) noexcept;

struct MaskChoice {
    int         K2 = 0;       ///< Integer squared radius, units of (2π/N)²
    std::size_t n_modes = 0;  ///< Half-plane modes inside
    double      chi = 0.0;    ///< Achieved chi (discrete: a whole shell kx² + ky² = const enters at once)
};

/// The integer K2 whose chi is closest to `chi` for n_points points on an N x N lattice.
/// Throws std::invalid_argument if chi is not positive or no shell fits.
[[nodiscard]] MaskChoice choose_mask_for_chi(int N, std::size_t n_points, double chi);

/// Builds the half-plane mask 0 < |k|² <= K2, the phase tables, the buffers and the mode weights
/// w(k) = |k|^(-weight_power) / mean over the mask of |k|^(-weight_power) (mean weight 1, so the energy
/// scale of a random start, and thus the temperature scale, does not depend on weight_power).
/// Throws std::invalid_argument if the mask is empty or K reaches N/2 - 1 (Nyquist limit).
void build_mask(HyperuniformConfig& config, double K2, double weight_power = 0.0);
