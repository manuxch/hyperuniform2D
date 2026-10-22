#include "minimize.hpp"

#include <algorithm>
#include <cmath>
#include <deque>

namespace {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double max_abs(const std::vector<double>& a) {
    double m = 0.0;
    for (const double v : a) m = std::max(m, std::abs(v));
    return m;
}

/// One evaluated point of the line search: step alpha, phi(alpha) = f(x + alpha d), phi'(alpha) = g·d.
struct Trial {
    double alpha = 0.0, phi = 0.0, dphi = 0.0;
    std::vector<double> x, g;
};

/// Minimizer of the cubic through (a, fa, da) and (b, fb, db); NaN if it does not exist.
double cubic_min(double a, double fa, double da, double b, double fb, double db) {
    const double d1 = da + db - 3.0 * (fa - fb) / (a - b);
    const double disc = d1 * d1 - da * db;
    if (disc < 0.0) return std::numeric_limits<double>::quiet_NaN();
    const double d2 = std::copysign(std::sqrt(disc), b - a);
    return b - (b - a) * (db + d2 - d1) / (db - da + 2.0 * d2);
}

class LineSearch {
public:
    LineSearch(const Objective& f, const std::vector<double>& x0, double f0, const std::vector<double>& d,
               double dphi0, const LbfgsOptions& opt, long long& evaluations)
        : f_(f), x0_(x0), d_(d), phi0_(f0), dphi0_(dphi0), opt_(opt), evals_(evaluations) {}

    /// Strong-Wolfe step starting at alpha1; returns false if none was found within max_line_search
    /// evaluations (the best trial with sufficient decrease, if any, is still returned in `out`).
    bool run(double alpha1, Trial& out) {
        Trial prev;
        prev.alpha = 0.0;
        prev.phi = phi0_;
        prev.dphi = dphi0_;
        double alpha = alpha1;
        for (int i = 0; i < opt_.max_line_search; ++i) {
            Trial t = eval(alpha);
            if (!std::isfinite(t.phi) || t.phi > phi0_ + opt_.c1 * alpha * dphi0_ || (i > 0 && t.phi >= prev.phi))
                return zoom(prev, t, out);
            if (std::abs(t.dphi) <= -opt_.c2 * dphi0_) { out = std::move(t); return true; }
            if (t.dphi >= 0.0) return zoom(t, prev, out);
            prev = std::move(t);
            alpha *= 2.0;
        }
        out = std::move(prev);
        return false;
    }

private:
    Trial eval(double alpha) {
        Trial t;
        t.alpha = alpha;
        t.x.resize(x0_.size());
        for (std::size_t i = 0; i < x0_.size(); ++i) t.x[i] = x0_[i] + alpha * d_[i];
        t.phi = f_(t.x, t.g);
        t.dphi = dot(t.g, d_);
        ++evals_;
        ++used_;
        return t;
    }

    /// Algorithm 3.6: lo satisfies sufficient decrease and has the lowest phi so far.
    bool zoom(Trial lo, Trial hi, Trial& out) {
        while (used_ < opt_.max_line_search) {
            const double a = lo.alpha, b = hi.alpha;
            double alpha = std::isfinite(hi.phi) ? cubic_min(a, lo.phi, lo.dphi, b, hi.phi, hi.dphi)
                                                 : std::numeric_limits<double>::quiet_NaN();
            const double lo_b = std::min(a, b), hi_b = std::max(a, b), w = hi_b - lo_b;
            if (!(alpha > lo_b + 0.1 * w && alpha < hi_b - 0.1 * w)) alpha = 0.5 * (a + b);  // safeguard: bisect
            if (w <= 1e-16 * std::max(1.0, hi_b)) break;
            Trial t = eval(alpha);
            if (!std::isfinite(t.phi) || t.phi > phi0_ + opt_.c1 * alpha * dphi0_ || t.phi >= lo.phi) {
                hi = std::move(t);
            } else {
                if (std::abs(t.dphi) <= -opt_.c2 * dphi0_) { out = std::move(t); return true; }
                if (t.dphi * (hi.alpha - lo.alpha) >= 0.0) hi = std::move(lo);
                lo = std::move(t);
            }
        }
        out = std::move(lo);  // best point with sufficient decrease (alpha = 0 if none)
        return false;
    }

    const Objective& f_;
    const std::vector<double>& x0_;
    const std::vector<double>& d_;
    double phi0_, dphi0_;
    const LbfgsOptions& opt_;
    long long& evals_;
    int used_ = 0;
};

}  // namespace

LbfgsResult lbfgs_minimize(const Objective& f, std::vector<double>& x, const LbfgsOptions& opt,
                           const Progress& progress) {
    LbfgsResult res;
    std::vector<double> g;
    double fx = f(x, g);
    res.evaluations = 1;
    struct Pair { std::vector<double> s, y; double rho; };
    std::deque<Pair> mem;
    std::vector<double> d(x.size()), alpha_buf;

    for (;;) {
        res.f = fx;
        res.g_max = max_abs(g);
        if (fx <= opt.f_target) { res.reason = "f_target"; break; }
        if (res.g_max <= opt.grad_tol) { res.reason = "grad_tol"; break; }
        if (res.iterations >= opt.max_iterations) { res.reason = "max_iterations"; break; }

        // Two-loop recursion: d = -H g.
        for (std::size_t i = 0; i < x.size(); ++i) d[i] = -g[i];
        alpha_buf.assign(mem.size(), 0.0);
        for (std::size_t k = mem.size(); k-- > 0;) {
            alpha_buf[k] = mem[k].rho * dot(mem[k].s, d);
            for (std::size_t i = 0; i < d.size(); ++i) d[i] -= alpha_buf[k] * mem[k].y[i];
        }
        if (!mem.empty()) {
            const auto& last = mem.back();
            const double gamma = dot(last.s, last.y) / dot(last.y, last.y);
            for (double& v : d) v *= gamma;
        }
        for (std::size_t k = 0; k < mem.size(); ++k) {
            const double beta = mem[k].rho * dot(mem[k].y, d);
            for (std::size_t i = 0; i < d.size(); ++i) d[i] += (alpha_buf[k] - beta) * mem[k].s[i];
        }
        double dphi0 = dot(g, d);
        if (!(dphi0 < 0.0)) {  // not a descent direction: restart from steepest descent
            mem.clear();
            for (std::size_t i = 0; i < x.size(); ++i) d[i] = -g[i];
            dphi0 = dot(g, d);
        }
        // Without curvature pairs the direction is -g: start with a step of unit length. With pairs the
        // quasi-Newton step is already scaled (gamma), so start at 1.
        const double alpha1 = mem.empty() ? 1.0 / std::max(std::sqrt(dot(g, g)), 1e-300) : 1.0;

        Trial step;
        LineSearch ls(f, x, fx, d, dphi0, opt, res.evaluations);
        const bool ok = ls.run(alpha1, step);
        if (!ok && !(step.alpha > 0.0 && step.phi < fx)) {
            if (!mem.empty()) { mem.clear(); continue; }  // retry once from steepest descent
            res.reason = "line_search";
            break;
        }
        // Accept the step (strong Wolfe, or at least a decrease).
        Pair p;
        p.s.resize(x.size());
        p.y.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) { p.s[i] = step.x[i] - x[i]; p.y[i] = step.g[i] - g[i]; }
        const double sy = dot(p.s, p.y);
        const double f_old = fx;
        x = std::move(step.x);
        g = std::move(step.g);
        fx = step.phi;
        ++res.iterations;
        if (sy > 1e-12 * std::sqrt(dot(p.s, p.s) * dot(p.y, p.y))) {  // curvature condition: keep the pair
            p.rho = 1.0 / sy;
            mem.push_back(std::move(p));
            if (static_cast<int>(mem.size()) > opt.memory) mem.pop_front();
        }
        if (progress && progress(res.iterations, fx, max_abs(g))) { res.f = fx; res.g_max = max_abs(g); res.reason = "callback"; break; }
        if (!(fx < f_old)) { res.f = fx; res.g_max = max_abs(g); res.reason = "no_progress"; break; }
    }
    return res;
}
