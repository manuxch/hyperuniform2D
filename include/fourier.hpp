#pragma once

#include "config.hpp"

#include <cstddef>
#include <string>

/// Calculates the Fourier sums of the mask from scratch (O(M·|modes|), table look-ups).
void calculate_FT(HyperuniformConfig& config);

/// Energy E = Σ_{0<|k|<=K} |rho_hat(k)|² over the full mask (both members of each ±k pair),
/// from the stored Fourier sums.
[[nodiscard]] double energy(const HyperuniformConfig& config);

/// Mean S(k) = |rho_hat(k)|² / M over the mask, divided by its uncorrelated lattice-gas value 1 - rho:
/// 1 for a random start, 0 for a perfectly stealthy pattern. On a lattice the annealing usually stops at
/// a small residual (a few per cent at chi ~ 0.2-0.3), so this is the practical quality measure.
[[nodiscard]] double mask_residual(const HyperuniformConfig& config, double E);

/// Change of rho_hat over modes [begin, end) when the particle at (xi, yi) moves to the empty site (xj, yj).
/// Fills d_re / d_im in that range and returns that range's contribution to the energy change.
[[nodiscard]] double delta_energy_range(HyperuniformConfig& config, int xi, int yi, int xj, int yj,
                                        std::size_t begin, std::size_t end);

/// Independent reference: E over the full mask (both ±k) by direct trigonometry, no tables.
/// For tests and verification only: O(M·|modes|) with sin/cos.
[[nodiscard]] double energy_direct(const HyperuniformConfig& config);

/// Writes the radially averaged structure factor S(k) = |rho_hat(k)|² / M of the final configuration
/// over the whole Brillouin zone (k ≠ 0), in bins of width 1 in |k| (units of 2π/N).
void analyze_structure_factor(const HyperuniformConfig& config, const std::string& filename);
