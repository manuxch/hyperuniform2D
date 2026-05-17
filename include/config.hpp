#pragma once

#include "constants.hpp"

#include <vector>
#include <random>

/// Wave vector index pair
struct ModeK { int kx, ky; };

/// Complete state of the hyperuniform configuration.
struct HyperuniformConfig {
    int    N = 0;   ///< Grid side length
    int    M = 0;   ///< Number of particles
    double K = 0.0; ///< Exclusion radius in k-space

    std::vector<int>     grid;         ///< Grid occupation 0/1
    std::vector<Complex> FT;           ///< Fourier transform for each mode
    std::vector<ModeK>   modes;        ///< List of wave vectors in the mask
    std::vector<Complex> delta_buffer; ///< Reusable buffer for swap attempts

    std::mt19937                           rng;
    std::uniform_int_distribution<int>     site_dist; ///< Initialized in initialize_grid()
    std::uniform_real_distribution<double> prob_dist{0.0, 1.0};
};

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

/// Initializes the grid with M randomly distributed particles.
void initialize_grid(HyperuniformConfig& config, int N, double density, int seed);

/// Builds the mask of k-modes with |k| <= K. Reserves delta_buffer.
void build_mask(HyperuniformConfig& config, double K);
