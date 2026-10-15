#pragma once

#include "config.hpp"
#include "parametros.hpp"

#include <optional>

class ThreadPool;

/// Proposes moving a random occupied site to an empty site and applies the Metropolis criterion
/// (T = 0: only downhill moves). Every call is a genuine proposal: the sites come from the occupied and
/// empty lists. With probability local_fraction the empty site is a neighbour within Chebyshev distance
/// local_radius (local moves), otherwise any empty site. Returns deltaE if accepted, std::nullopt if not.
/// Uses config.d_re / d_im (pre-allocated) — no dynamic allocation per call.
[[nodiscard]] std::optional<double> attempt_swap(HyperuniformConfig& config, double T, ThreadPool* pool,
                                                 double local_fraction = 0.0, int local_radius = 1);

struct AnnealReport {
    double    T_initial = 0.0;        ///< Temperature actually used (automatic if the parameter was <= 0)
    double    T_last = 0.0;           ///< Temperature of the last annealing level
    long long levels = 0;             ///< Temperature levels run
    long long attempts = 0;           ///< Proposals, annealing + quench
    long long accepted = 0;
    double    E_initial = 0.0;        ///< Exact values, recomputed from scratch
    double    E_final = 0.0;
    double    max_drift = 0.0;        ///< Largest |incremental E - exact E| seen at resyncs
    bool      reached_tol = false;    ///< E_final <= energy_tol
};

/// Executes the simulated annealing (and the optional zero-temperature quench). The energy is kept with
/// an incremental accumulator and resynchronised from scratch every p.resync_interval accepted swaps.
AnnealReport simulated_annealing(HyperuniformConfig& config, const Parameters& p);
