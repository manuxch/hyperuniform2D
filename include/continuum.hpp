#pragma once

#include "config.hpp"

#include <complex>
#include <cstddef>
#include <random>
#include <vector>

class ThreadPool;

/// Point pattern in a periodic box [0, Lx) x [0, Ly) (lengths in grain diameters, d = 1).
///
/// Positions are stored interleaved, r = (x0, y0, x1, y1, ...): that vector is what the minimizer works on.
/// During a minimization they may leave the box; every energy is periodic in both directions, and positions
/// are folded back only for output.
///
/// Mask: the wave vectors k = 2π (nx/Lx, ny/Ly) with 0 < |k| <= K, a circle in physical k (an ellipse in
/// the indices when Lx != Ly). As on the lattice, only the half plane (nx > 0, or nx = 0 and ny > 0) is
/// stored, and E = 2 Σ_half w(k) |rho_hat(k)|² counts both members of each ±k pair.
///
/// Phases: exp(-i k·r) = ex[nx] ey[ny], with ex[n] = exp(-i 2π n x / Lx) obtained by recurrence along each
/// axis (one complex product per step, recomputed exactly every 32 steps), so there is no trigonometry
/// per mode.
struct ContinuumConfig {
    double      Lx = 0.0, Ly = 0.0;
    std::size_t N = 0;               ///< Number of points
    std::vector<double> r;           ///< Interleaved positions, size 2N

    double K = 0.0;                  ///< Mask radius, d^-1 (2π included)
    std::vector<ModeK>  modes;       ///< Half-plane index pairs (nx, ny)
    std::vector<double> kx, ky;      ///< Physical wave vector of each mode
    std::vector<double> weight;      ///< w(k), mean 1 (all 1 if weight_power = 0)
    double weight_power = 0.0;
    int nx_max = 0, ny_max = 0;      ///< Largest |nx|, |ny| in the mask

    /// The modes in storage order form runs of equal nx and consecutive ny (a circle cuts each column
    /// nx = const in one interval), so exp(-i k·r) = ex[nx] ey[ny] with ex fixed along a run. The kernels
    /// loop over runs with contiguous real arrays, which the compiler vectorizes.
    struct ModeRun { int nx = 0, ny_lo = 0; std::size_t begin = 0, count = 0; };
    std::vector<ModeRun> runs;

    std::vector<double> ft_re, ft_im;  ///< rho_hat(k) of each mode
    std::vector<double> part_re, part_im;  ///< Scratch: partial rho_hat of each point group (continuum_FT)

    double sigma = 0.0;              ///< Hard-core distance (0 = no core term)
    double core_weight = 1.0;        ///< E_core = core_weight Σ_{i<j, r<sigma} (1 - r/sigma)²

    std::mt19937 rng;

    [[nodiscard]] double density() const noexcept { return static_cast<double>(N) / (Lx * Ly); }
};

/// exp(-i 2π n u / L) for n = 0..n_max, by recurrence (exact every 32 steps). out must have n_max + 1 entries.
void axis_phases(double u, double L, int n_max, std::complex<double>* out) noexcept;

/// N uniform random points in [0, Lx) x [0, Ly). With sigma > 0, random sequential addition: no two points
/// (periodic distance) closer than sigma. Throws std::invalid_argument if that does not succeed within
/// 1000 N attempts or the arguments are invalid.
void continuum_init_random(ContinuumConfig& c, double Lx, double Ly, std::size_t N, unsigned seed, double sigma = 0.0);

/// Number of half-plane modes with 0 < |k| <= K.
[[nodiscard]] std::size_t continuum_mode_count(double Lx, double Ly, double K);

struct ContinuumMaskChoice {
    double      K = 0.0;        ///< Radius written midway between the last included shell and the next one
    std::size_t n_modes = 0;    ///< M(K), half space
    double      chi = 0.0;      ///< Achieved point-pattern chi = M(K) / [2 (N - 1)]
};

/// The shell whose point-pattern chi is closest to `chi`. Shells are the distinct values of |k|²; K is set
/// midway to the next shell. Throws std::invalid_argument if chi is not positive or N < 2.
[[nodiscard]] ContinuumMaskChoice choose_continuum_mask_for_chi(double Lx, double Ly, std::size_t N, double chi);

/// K for a given kappa = K / (2π sqrt(rho)), rho = N / (Lx Ly).
[[nodiscard]] double continuum_K_of_kappa(double Lx, double Ly, std::size_t N, double kappa) noexcept;

/// Builds the mask 0 < |k| <= K (half plane) and the weights w ∝ |k|^(-weight_power), mean 1.
/// Throws std::invalid_argument if the mask is empty.
void build_continuum_mask(ContinuumConfig& c, double K, double weight_power = 0.0);

/// rho_hat(k) of every mode from the positions (O(N |modes|)); fills ft_re / ft_im. The points are split
/// into a fixed number of groups (at most 32, independent of the thread count); each group sums into its
/// own buffer, the groups run in parallel on `pool` (if given) and are added in order, so the result is
/// the same for any thread count.
void continuum_FT(ContinuumConfig& c, ThreadPool* pool = nullptr);

/// Contribution of the points [p_begin, p_end) to rho_hat of every mode, accumulated into re / im
/// (arrays of |modes| entries).
void continuum_FT_points(const ContinuumConfig& c, std::size_t p_begin, std::size_t p_end, double* re, double* im);

/// E = 2 Σ_half |rho_hat|² and E_w = 2 Σ_half w |rho_hat|², from the stored rho_hat.
[[nodiscard]] double continuum_energy(const ContinuumConfig& c);
[[nodiscard]] double continuum_weighted_energy(const ContinuumConfig& c);

/// <S>_mask = E / (2 |modes| N): 1 for a Poisson pattern, 0 for a stealthy one.
[[nodiscard]] double continuum_residual(const ContinuumConfig& c, double E);

/// Mean S(k) / 1 over the mask modes with |k| <= k_max (d^-1).
[[nodiscard]] double continuum_band_residual(const ContinuumConfig& c, double k_max);

/// Gradient of E_w with respect to the points [p_begin, p_end), from the stored rho_hat, written into
/// g[2 p], g[2 p + 1] (g has size 2N):  dE_w/dr_j = 4 Σ_half w k Im(conj(rho_hat(k)) exp(-i k·r_j)).
/// Points are processed in blocks with the block as the inner loop, so the loop vectorizes without a
/// reduction across lanes; each point's sum always runs in the same order, whatever the split.
void continuum_gradient_range(const ContinuumConfig& c, std::size_t p_begin, std::size_t p_end, std::vector<double>& g);

/// Gradient of E_w for every point (resizes g to 2N), split by points on `pool` if given.
void continuum_gradient(const ContinuumConfig& c, std::vector<double>& g, ThreadPool* pool = nullptr);

/// Recomputes rho_hat and returns E_w; fills g (resized to 2N) with its gradient.
double continuum_spectral_energy_grad(ContinuumConfig& c, std::vector<double>& g, ThreadPool* pool = nullptr);

/// Hard-core energy E_core = core_weight Σ_{i<j, r_ij < sigma} (1 - r_ij/sigma)² (periodic distances),
/// O(N) with a cell list; adds its gradient to g (size 2N) if g is not null. 0 if sigma = 0.
double continuum_core_energy(const ContinuumConfig& c, std::vector<double>* g);

/// Smallest periodic distance between two points (both axes periodic); +inf for N < 2.
[[nodiscard]] double continuum_min_distance(const ContinuumConfig& c);

/// K for the obstacle file: midway between the outermost shell of the mask and the next shell, so that a
/// reader that rebuilds the circle |k| <= K gets the same modes despite rounding.
[[nodiscard]] double continuum_mask_K_midpoint(const ContinuumConfig& c);

/// Independent references for the tests: rho_hat over the full circle |k| <= K (all n in Z², by direct
/// trigonometry), E over the full mask, and the hard-core energy by brute force.
[[nodiscard]] double continuum_energy_direct(const ContinuumConfig& c);
[[nodiscard]] double continuum_core_energy_direct(const ContinuumConfig& c);
