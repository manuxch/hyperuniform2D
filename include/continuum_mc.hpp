#pragma once

#include "continuum.hpp"
#include "parametros.hpp"

#include <ostream>

class ThreadPool;

struct ContinuumMcReport {
    double    T_initial = 0.0, T_last = 0.0;
    long long levels = 0, attempts = 0, accepted = 0;
    double    f_initial = 0.0, f_final = 0.0;  ///< Objective: E_w / (2 M N) + E_core
    double    step_final = 0.0;                ///< Displacement amplitude at the end (d)
    double    max_drift = 0.0;                 ///< Largest |incremental - exact| objective at the resyncs
};

/// Simulated annealing of a continuum pattern with single-point moves: a random point moves by a uniform
/// displacement in [-step, step]², accepted with the Metropolis probability of the change of the objective
/// f = E_w / (2 M N) + E_core. rho_hat is updated incrementally (O(|modes|) per move, split among the threads
/// of `pool`) and recomputed every p.resync_interval accepted moves; the hard core is evaluated locally with
/// a cell list, as a constraint (moves that create or deepen an overlap are rejected). Schedule as on the
/// lattice: T_initial (<= 0: automatic), cooling_rate, steps_per_temp per level down to T_final, then
/// quench_steps at T = 0; the step adapts after each level towards an acceptance of 0.3-0.5. Stops early when f <= p.residual_tol. Leaves c.ft_re / ft_im consistent with c.r.
ContinuumMcReport continuum_anneal(ContinuumConfig& c, const Parameters& p, ThreadPool* pool, std::ostream* log);
