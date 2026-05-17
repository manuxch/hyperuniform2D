#pragma once

#include "config.hpp"

/// Calculates the complete Fourier transform from scratch (O(N²·|modes|)).
void calculate_FT(HyperuniformConfig& config);

/// Calculates the energy E = Σ |FT[m]|² over all modes in the mask.
[[nodiscard]] double energy(const HyperuniformConfig& config);

/// Writes S(k) = |FT(kx,0)|²/N² for kx = 1..N/2-1 to the specified file.
void analyze_structure_factor(const HyperuniformConfig& config,
                              const std::string& filename);
