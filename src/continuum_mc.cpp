#include "continuum_mc.hpp"
#include "thread_pool.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <vector>

namespace {

using cplx = std::complex<double>;

double fold(double u, double L) noexcept {
    double f = u - L * std::floor(u / L);
    if (f >= L || f < 0.0) f = 0.0;
    return f;
}

double min_image(double d, double L) noexcept { return d - L * std::round(d / L); }

/// Cell list that follows the moves (cells of side >= sigma; direct search if fewer than 3 per axis).
class MovingCells {
public:
    MovingCells(const ContinuumConfig& c) : c_(c) {
        if (!(c.sigma > 0.0)) return;
        ncx_ = static_cast<int>(c.Lx / c.sigma);
        ncy_ = static_cast<int>(c.Ly / c.sigma);
        direct_ = ncx_ < 3 || ncy_ < 3;
        if (direct_) return;
        hx_ = c.Lx / ncx_;
        hy_ = c.Ly / ncy_;
        bins_.assign(static_cast<std::size_t>(ncx_) * static_cast<std::size_t>(ncy_), {});
        cell_.resize(c.N);
        for (std::size_t i = 0; i < c.N; ++i) {
            cell_[i] = cell_of(c.r[2 * i], c.r[2 * i + 1]);
            bins_[cell_[i]].push_back(i);
        }
    }

    /// Hard-core energy of point i placed at (x, y) with all the others (i itself excluded).
    [[nodiscard]] double energy_at(std::size_t i, double x, double y) const {
        if (!(c_.sigma > 0.0)) return 0.0;
        double E = 0.0;
        const auto pair = [&](std::size_t j) {
            if (j == i) return;
            const double dx = min_image(x - c_.r[2 * j], c_.Lx), dy = min_image(y - c_.r[2 * j + 1], c_.Ly);
            const double d2 = dx * dx + dy * dy;
            if (d2 >= c_.sigma * c_.sigma) return;
            const double u = 1.0 - std::sqrt(d2) / c_.sigma;
            E += c_.core_weight * u * u;
        };
        if (direct_) {
            for (std::size_t j = 0; j < c_.N; ++j) pair(j);
            return E;
        }
        const std::size_t cc = cell_of(x, y);
        const int cx = static_cast<int>(cc % static_cast<std::size_t>(ncx_)), cy = static_cast<int>(cc / static_cast<std::size_t>(ncx_));
        for (int oy = -1; oy <= 1; ++oy)
            for (int ox = -1; ox <= 1; ++ox) {
                const int bx = (cx + ox + ncx_) % ncx_, by = (cy + oy + ncy_) % ncy_;
                for (const std::size_t j : bins_[static_cast<std::size_t>(bx + ncx_ * by)]) pair(j);
            }
        return E;
    }

    /// Point i has moved to its current position in c.r.
    void moved(std::size_t i) {
        if (!(c_.sigma > 0.0) || direct_) return;
        const std::size_t nc = cell_of(c_.r[2 * i], c_.r[2 * i + 1]);
        if (nc == cell_[i]) return;
        auto& old = bins_[cell_[i]];
        old.erase(std::find(old.begin(), old.end(), i));
        bins_[nc].push_back(i);
        cell_[i] = nc;
    }

private:
    [[nodiscard]] std::size_t cell_of(double x, double y) const noexcept {
        const int cx = std::min(ncx_ - 1, static_cast<int>(fold(x, c_.Lx) / hx_));
        const int cy = std::min(ncy_ - 1, static_cast<int>(fold(y, c_.Ly) / hy_));
        return static_cast<std::size_t>(cx + ncx_ * cy);
    }

    const ContinuumConfig& c_;
    int ncx_ = 0, ncy_ = 0;
    double hx_ = 0.0, hy_ = 0.0;
    bool direct_ = false;
    std::vector<std::vector<std::size_t>> bins_;
    std::vector<std::size_t> cell_;
};

struct Mover {
    ContinuumConfig& c;
    ThreadPool* pool;
    MovingCells cells;
    double norm;                       ///< 1 / (2 M N)
    std::vector<double> d_re, d_im;    ///< Change of rho_hat of the current proposal

    Mover(ContinuumConfig& cc, ThreadPool* pp)
        : c(cc), pool(pp), cells(cc),
          norm(1.0 / (2.0 * static_cast<double>(cc.modes.size()) * static_cast<double>(cc.N))),
          d_re(cc.modes.size()), d_im(cc.modes.size()) {}

    /// Change of the objective if point i moves to (x, y); fills d_re / d_im.
    double delta(std::size_t i, double x, double y) {
        const double xo = c.r[2 * i], yo = c.r[2 * i + 1];
        auto range = [&](std::size_t b, std::size_t e) {
            std::vector<cplx> exo(static_cast<std::size_t>(c.nx_max) + 1), eyo(static_cast<std::size_t>(c.ny_max) + 1);
            std::vector<cplx> exn(exo.size()), eyn(eyo.size());
            axis_phases(xo, c.Lx, c.nx_max, exo.data());
            axis_phases(yo, c.Ly, c.ny_max, eyo.data());
            axis_phases(x, c.Lx, c.nx_max, exn.data());
            axis_phases(y, c.Ly, c.ny_max, eyn.data());
            double dE = 0.0;
            for (std::size_t m = b; m < e; ++m) {
                const ModeK& k = c.modes[m];
                const cplx po = exo[static_cast<std::size_t>(k.kx)] * (k.ky >= 0 ? eyo[static_cast<std::size_t>(k.ky)] : std::conj(eyo[static_cast<std::size_t>(-k.ky)]));
                const cplx pn = exn[static_cast<std::size_t>(k.kx)] * (k.ky >= 0 ? eyn[static_cast<std::size_t>(k.ky)] : std::conj(eyn[static_cast<std::size_t>(-k.ky)]));
                const double dr = pn.real() - po.real(), di = pn.imag() - po.imag();
                d_re[m] = dr;
                d_im[m] = di;
                dE += c.weight[m] * (2.0 * (c.ft_re[m] * dr + c.ft_im[m] * di) + dr * dr + di * di);
            }
            return dE;
        };
        // Hard core as a constraint: a move that creates or deepens an overlap is rejected outright (+inf),
        // so a start without overlaps never gets one. An overlap already present (init_file) may only shrink.
        const double core_new = cells.energy_at(i, x, y), core_old = cells.energy_at(i, xo, yo);
        if (core_new > 0.0 && core_new >= core_old) return std::numeric_limits<double>::infinity();
        const double dEw = 2.0 * (pool ? pool->parallel_reduce(c.modes.size(), range) : range(0, c.modes.size()));
        return norm * dEw + (core_new - core_old);
    }

    void accept(std::size_t i, double x, double y) {
        c.r[2 * i] = x;
        c.r[2 * i + 1] = y;
        for (std::size_t m = 0; m < c.modes.size(); ++m) {
            c.ft_re[m] += d_re[m];
            c.ft_im[m] += d_im[m];
        }
        cells.moved(i);
    }

    [[nodiscard]] double exact() {
        continuum_FT(c);
        return norm * continuum_weighted_energy(c) + continuum_core_energy(c, nullptr);
    }
};

}  // namespace

ContinuumMcReport continuum_anneal(ContinuumConfig& c, const Parameters& p, ThreadPool* pool, std::ostream* log) {
    ContinuumMcReport rep;
    Mover mv(c, pool);
    double f = mv.exact();
    rep.f_initial = f;
    const double a = 1.0 / std::sqrt(c.density());
    double step = p.mc_step > 0.0 ? p.mc_step : 0.1 * a;
    std::uniform_int_distribution<std::size_t> pick(0, c.N - 1);
    std::uniform_real_distribution<double> unit(-1.0, 1.0), prob(0.0, 1.0);

    double T = p.T_initial;
    if (T <= 0.0) {  // a typical uphill move from the start is accepted with probability 1/2
        double sum = 0.0;
        int n_up = 0;
        for (int s = 0; s < 400; ++s) {
            const std::size_t i = pick(c.rng);
            const double dF = mv.delta(i, c.r[2 * i] + step * unit(c.rng), c.r[2 * i + 1] + step * unit(c.rng));
            if (dF > 0.0 && std::isfinite(dF)) { sum += dF; ++n_up; }
        }
        T = n_up > 0 ? (sum / n_up) / std::numbers::ln2 : 1.0;
        if (T <= p.T_final) T = 10.0 * p.T_final;
        std::cout << std::format("Automatic T_initial = {:.6g}\n", T);
    }
    rep.T_initial = T;
    if (log) *log << std::format("{:<16} {:<22} {:<12} {:<14}\n", "# T", "objective", "acc/level", "step");

    long long since_resync = 0;
    const auto resync = [&] {
        const double ex = mv.exact();
        rep.max_drift = std::max(rep.max_drift, std::abs(ex - f));
        f = ex;
        since_resync = 0;
    };
    const auto attempt = [&](double temp) {
        ++rep.attempts;
        const std::size_t i = pick(c.rng);
        const double x = c.r[2 * i] + step * unit(c.rng), y = c.r[2 * i + 1] + step * unit(c.rng);
        const double dF = mv.delta(i, x, y);
        if (dF < 0.0 || (temp > 0.0 && prob(c.rng) < std::exp(-dF / temp))) {
            mv.accept(i, x, y);
            f += dF;
            ++rep.accepted;
            if (++since_resync >= p.resync_interval) resync();
            return true;
        }
        return false;
    };

    while (T > p.T_final && f > p.residual_tol) {
        long long acc = 0;
        for (int s = 0; s < p.steps_per_temp; ++s) acc += attempt(T);
        const double rate = static_cast<double>(acc) / p.steps_per_temp;
        ++rep.levels;
        if (log) *log << std::format("{:<16.8g} {:<22.12e} {:<12.6f} {:<14.6g}\n", T, f, rate, step);
        if (rate > 0.5) step = std::min(a, step * 1.2);
        else if (rate < 0.3) step = std::max(1e-9 * a, step * 0.8);
        rep.T_last = T;
        T *= p.cooling_rate;
    }
    for (long long s = 0; s < p.quench_steps && f > p.residual_tol; ++s) attempt(0.0);
    resync();
    rep.f_final = f;
    rep.step_final = step;
    return rep;
}
