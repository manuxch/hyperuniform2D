#pragma once

#include "config.hpp"
#include "parametros.hpp"

#include <optional>

class ThreadPool;

/// Attempts to swap an occupied site with an empty one using the Metropolis criterion.
/// Returns deltaE if the swap was accepted, std::nullopt if rejected.
/// Uses config.delta_buffer (pre-allocated) — no dynamic allocation per call.
[[nodiscard]] std::optional<double> attempt_swap(HyperuniformConfig& config, double T, ThreadPool* pool);

/// Executes the simulated annealing. Maintains energy with an incremental accumulator.
void simulated_annealing(HyperuniformConfig& config, const Parameters& p);
