#pragma once

#include <functional>
#include <limits>
#include <string>
#include <vector>

/// Limited-memory BFGS with a strong-Wolfe line search (Nocedal & Wright, Numerical Optimization,
/// 2nd ed., algorithms 7.4, 3.5 and 3.6). Self-contained, no external dependency; the interface lets
/// another implementation replace it without touching the callers.

/// f(x, g): returns the objective at x and writes its gradient into g (resized by the callee if needed).
using Objective = std::function<double(const std::vector<double>& x, std::vector<double>& g)>;

/// Called after every accepted step with the iteration number, f and the largest |g_i|.
/// Returning true stops the minimization (reason "callback").
using Progress = std::function<bool(long long iteration, double f, double g_max)>;

struct LbfgsOptions {
    int       memory = 10;              ///< Correction pairs kept
    long long max_iterations = 100000;
    double    grad_tol = 0.0;           ///< Stop when max |g_i| <= grad_tol
    double    f_target = -std::numeric_limits<double>::infinity();  ///< Stop when f <= f_target
    int       max_line_search = 40;     ///< Function evaluations per line search
    double    c1 = 1e-4, c2 = 0.9;      ///< Wolfe constants
};

struct LbfgsResult {
    long long   iterations = 0;
    long long   evaluations = 0;
    double      f = 0.0;
    double      g_max = 0.0;
    std::string reason;  ///< "f_target", "grad_tol", "callback", "max_iterations", "line_search", "no_progress"
};

/// Minimizes f from x (updated in place to the best point found).
LbfgsResult lbfgs_minimize(const Objective& f, std::vector<double>& x, const LbfgsOptions& opt,
                           const Progress& progress = {});
