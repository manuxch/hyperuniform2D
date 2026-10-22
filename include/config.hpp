#pragma once

#include "constants.hpp"

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

/// Wave vector index pair: k = 2π (kx / Nx, ky / Ny) in units of the inverse lattice spacing
struct ModeK { int kx, ky; };

/// Complete state of the hyperuniform configuration on a periodic Nx x Ny lattice (spacing 1).
///
/// Mask: a circle in physical k-space, 0 < |k| <= K, never a circle in index space (with Nx != Ny that
/// would be an ellipse in k, i.e. a directionally hyperuniform pattern). Each mode carries the exact
/// integer key q = kx² Ny² + ky² Nx², proportional to |k|²: |k|² = (2π)² q / (Nx Ny)². Masks and shells
/// are selected on q, so no shell is split by rounding.
///
/// K and K2 are expressed in units of 2π/sqrt(Nx Ny), so |k|² = K2_mode (2π)² / (Nx Ny) with
/// K2_mode = q / (Nx Ny). For a square lattice that is the unit 2π/N of version 2.x.
///
/// Only the half plane of k is stored (kx > 0, or kx = 0 and ky > 0): for a real occupation
/// rho_hat(-k) = conj(rho_hat(k)), so the other half carries no information. The energy still
/// counts both members of each ±k pair, E = 2 Σ_half |rho_hat(k)|², i.e. the sum over the full mask.
///
/// Phase factors are table look-ups, with no trigonometry in the annealing loop. With L = lcm(Nx, Ny),
///   exp(-i 2π (kx x / Nx + ky y / Ny)) = w[(kx x L/Nx + ky y L/Ny) mod L],
/// and each axis term comes from its own table mtx / mty. For Nx = Ny = N, L = N and the tables are
/// those of version 2.x, so results are bit-for-bit identical.
struct HyperuniformConfig {
    int    Nx = 0, Ny = 0;  ///< Lattice sides
    int    M  = 0;          ///< Number of occupied sites
    int    L  = 0;          ///< lcm(Nx, Ny): length of the phase table
    double K2 = 0.0;        ///< Squared mask radius, units of (2π/sqrt(Nx Ny))²
    double K  = 0.0;        ///< Mask radius, sqrt(K2)
    std::int64_t q_max = 0; ///< Mask: modes with q <= q_max
    std::int64_t q_next = 0;///< Smallest q outside the mask (0 if none fits below the Nyquist limit)

    std::vector<std::uint8_t>  grid;      ///< Grid occupation 0/1, index x + Nx y
    std::vector<std::uint32_t> occupied;  ///< Sites with grid = 1
    std::vector<std::uint32_t> empty;     ///< Sites with grid = 0
    std::vector<std::uint32_t> list_pos;  ///< Position of each site in its list (occupied or empty)

    std::vector<ModeK>        modes;      ///< Half-plane wave vectors in the mask
    std::vector<double>       mode_k2;    ///< |k|² of each mode, units of (2π/sqrt(Nx Ny))²
    std::vector<std::int32_t> mtx, mty;   ///< (k x L/Nx) mod L for k in [-Nx/2, Nx/2), x in [0, Nx); same for y
    std::vector<std::int32_t> kx_off;     ///< (kx + Nx/2) Nx for each mode: offset into mtx
    std::vector<std::int32_t> ky_off;     ///< (ky + Ny/2) Ny for each mode: offset into mty
    std::vector<double>       w_re, w_im; ///< exp(-i 2π j / L), j in [0, L)

    std::vector<double> weight;           ///< w(k) of each mode in the annealing energy (all 1 if p = 0)
    double weight_power = 0.0;            ///< p in w(k) ∝ |k|^(-p)
    std::vector<double> ft_re, ft_im;     ///< rho_hat(k) for each mode
    std::vector<double> d_re, d_im;       ///< Reusable buffer: change of rho_hat for the proposed swap

    /// Hard-core exclusion: no two occupied sites closer than min_distance (lattice spacings). excl holds
    /// the offsets (dx, dy) with 0 < dx² + dy² < min_distance²; empty = no exclusion.
    double min_distance = 0.0;
    std::vector<std::pair<int, int>> excl;
    long long hc_rejections = 0;          ///< Proposals rejected by the exclusion

    std::mt19937 rng;
    std::uniform_real_distribution<double> prob_dist{0.0, 1.0};

    [[nodiscard]] bool square() const noexcept { return Nx == Ny; }
    [[nodiscard]] std::size_t sites() const noexcept {
        return static_cast<std::size_t>(Nx) * static_cast<std::size_t>(Ny);
    }
};

/// Index of the phase of mode m at site (x, y) in w_re / w_im.
[[nodiscard]] inline std::size_t phase_index(const HyperuniformConfig& c, std::size_t m, int x, int y) noexcept {
    std::int32_t j = c.mtx[static_cast<std::size_t>(c.kx_off[m] + x)] + c.mty[static_cast<std::size_t>(c.ky_off[m] + y)];
    if (j >= c.L) j -= c.L;
    return static_cast<std::size_t>(j);
}

/// Direct index — only for coordinates already in [0, Nx) x [0, Ny). No modulo.
[[nodiscard]] inline int flat_idx(int x, int y, int Nx) noexcept {
    return x + Nx * y;
}

/// Index with periodic boundary conditions — for arbitrary coordinates.
[[nodiscard]] inline int periodic_idx(int x, int y, int Nx, int Ny) noexcept {
    x = (x % Nx + Nx) % Nx;
    y = (y % Ny + Ny) % Ny;
    return x + Nx * y;
}

/// Exact key q = kx² Ny² + ky² Nx² of a mode, proportional to |k|².
[[nodiscard]] inline std::int64_t mode_key(int kx, int ky, int Nx, int Ny) noexcept {
    const auto x = static_cast<std::int64_t>(kx), y = static_cast<std::int64_t>(ky);
    const auto nx = static_cast<std::int64_t>(Nx), ny = static_cast<std::int64_t>(Ny);
    return x * x * ny * ny + y * y * nx * nx;
}

/// True if an occupied site other than `ignore` lies closer than min_distance to site j.
[[nodiscard]] inline bool hardcore_blocked(const HyperuniformConfig& c, std::uint32_t j, std::uint32_t ignore) noexcept {
    const int x = static_cast<int>(j % static_cast<std::uint32_t>(c.Nx)), y = static_cast<int>(j / static_cast<std::uint32_t>(c.Nx));
    for (const auto& [dx, dy] : c.excl) {
        const auto s = static_cast<std::uint32_t>(periodic_idx(x + dx, y + dy, c.Nx, c.Ny));
        if (s != ignore && c.grid[s]) return true;
    }
    return false;
}

/// Smallest distance between two occupied sites (periodic, lattice spacings), searched up to `r_max`;
/// returns r_max if no pair is closer.
[[nodiscard]] double min_pair_distance(const HyperuniformConfig& c, double r_max);

/// Number of occupied sites for a target density: round(density Nx Ny), not truncation
/// (0.29 * 100² = 2899.9999... in floating point).
[[nodiscard]] int occupied_count(int Nx, int Ny, double density) noexcept;
[[nodiscard]] inline int occupied_count(int N, double density) noexcept { return occupied_count(N, N, density); }

/// Initializes the grid with M = occupied_count(Nx, Ny, density) randomly distributed particles.
/// With min_distance > 0 the start is a random sequential addition that respects the exclusion
/// (no two occupied sites closer than min_distance lattice spacings), and the exclusion is kept by the
/// annealing. Throws std::invalid_argument unless there is at least one occupied and one empty site, or
/// if the exclusion does not leave room for M sites.
void initialize_grid(HyperuniformConfig& config, int Nx, int Ny, double density, unsigned seed,
                     double min_distance = 0.0);
inline void initialize_grid(HyperuniformConfig& config, int N, double density, unsigned seed) {
    initialize_grid(config, N, N, density, seed);
}

/// Stealthiness chi = (half-plane constrained modes) / (n - 1), with n = min(M, Nx Ny - M): the
/// binary_lattice convention of the obstacle file format (hyperuniform2D 2.x). It is NOT the point-pattern
/// chi = M(K) / [2 (N - 1)] of Torquato et al.: for density < 1/2 and the same K it is about twice that.
/// The minority species is used because for k ≠ 0 the two species have the same |rho_hat|.
[[nodiscard]] double chi_of(std::size_t n_modes, std::size_t n_points) noexcept;

/// Point-pattern chi = M(K) / [2 (n - 1)] (d = 2), with n the number of points (the occupied sites).
[[nodiscard]] double chi_point_pattern(std::size_t n_modes, std::size_t n_points) noexcept;

/// n = min(M, Nx Ny - M), the number of points that enters the binary_lattice chi.
[[nodiscard]] std::size_t chi_points(int Nx, int Ny, int M) noexcept;
[[nodiscard]] inline std::size_t chi_points(int N, int M) noexcept { return chi_points(N, N, M); }

struct MaskChoice {
    double       K2 = 0.0;    ///< Squared radius of the outermost included shell, units of (2π/sqrt(Nx Ny))²
    std::int64_t q_max = 0;   ///< Key of that shell
    std::size_t  n_modes = 0; ///< Half-plane modes inside
    double       chi = 0.0;   ///< Achieved binary_lattice chi (discrete: a whole shell enters at once)
};

/// The shell whose binary_lattice chi is closest to `chi` for n_points points on an Nx x Ny lattice.
/// Throws std::invalid_argument if chi is not positive or no shell fits.
[[nodiscard]] MaskChoice choose_mask_for_chi(int Nx, int Ny, std::size_t n_points, double chi);
[[nodiscard]] inline MaskChoice choose_mask_for_chi(int N, std::size_t n_points, double chi) {
    return choose_mask_for_chi(N, N, n_points, chi);
}

/// Largest K2 (units of (2π/sqrt(Nx Ny))²) whose circle fits below the Nyquist limit on both axes,
/// |kx| <= Nx/2 - 1 and |ky| <= Ny/2 - 1.
[[nodiscard]] double max_mask_K2(int Nx, int Ny) noexcept;

/// Builds the half-plane mask 0 < |k|² <= K2 (a circle in physical k), the phase tables, the buffers and
/// the mode weights w(k) = |k|^(-weight_power) / mean over the mask of |k|^(-weight_power) (mean weight 1,
/// so the energy scale of a random start, and thus the temperature scale, does not depend on weight_power).
/// Throws std::invalid_argument if the mask is empty, the circle crosses the Nyquist limit, or
/// lcm(Nx, Ny) is too large for the phase table.
void build_mask(HyperuniformConfig& config, double K2, double weight_power = 0.0);
