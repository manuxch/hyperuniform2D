#include "continuum.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

using cplx = std::complex<double>;

/// u folded into [0, L).
double fold(double u, double L) noexcept {
    double f = u - L * std::floor(u / L);
    if (f >= L) f = 0.0;  // u slightly below a multiple of L can round up to L
    return f;
}

/// Minimum-image difference along a periodic axis of length L.
double min_image(double d, double L) noexcept {
    return d - L * std::round(d / L);
}

/// Phase of mode m at a point, from its axis tables (ey holds non-negative ny; negative ny by conjugation).
inline cplx mode_phase(const ModeK& m, const cplx* ex, const cplx* ey) noexcept {
    const cplx a = ex[m.kx];
    const cplx b = m.ky >= 0 ? ey[m.ky] : std::conj(ey[-m.ky]);
    return a * b;
}

/// Cell list over the folded positions with cells of side >= h on both axes (periodic). Returns false if
/// there would be fewer than 3 cells on an axis (the caller then uses a direct search).
struct Cells {
    int ncx = 0, ncy = 0;
    double hx = 0.0, hy = 0.0;
    std::vector<std::vector<std::size_t>> bins;

    bool build(const ContinuumConfig& c, double h) {
        ncx = static_cast<int>(c.Lx / h);
        ncy = static_cast<int>(c.Ly / h);
        if (ncx < 3 || ncy < 3) return false;
        hx = c.Lx / ncx;
        hy = c.Ly / ncy;
        bins.assign(static_cast<std::size_t>(ncx) * static_cast<std::size_t>(ncy), {});
        for (std::size_t i = 0; i < c.N; ++i) bins[cell(c, i)].push_back(i);
        return true;
    }
    [[nodiscard]] std::pair<int, int> xy(const ContinuumConfig& c, std::size_t i) const noexcept {
        const int cx = std::min(ncx - 1, static_cast<int>(fold(c.r[2 * i], c.Lx) / hx));
        const int cy = std::min(ncy - 1, static_cast<int>(fold(c.r[2 * i + 1], c.Ly) / hy));
        return {cx, cy};
    }
    [[nodiscard]] std::size_t cell(const ContinuumConfig& c, std::size_t i) const noexcept {
        const auto [cx, cy] = xy(c, i);
        return static_cast<std::size_t>(cx + ncx * cy);
    }
    /// Calls f(j) for every point j in the 3 x 3 neighbourhood of point i's cell (each cell once).
    template <class F>
    void neighbours(const ContinuumConfig& c, std::size_t i, F&& f) const {
        const auto [cx, cy] = xy(c, i);
        for (int oy = -1; oy <= 1; ++oy)
            for (int ox = -1; ox <= 1; ++ox) {
                const int x = (cx + ox + ncx) % ncx, y = (cy + oy + ncy) % ncy;
                for (const std::size_t j : bins[static_cast<std::size_t>(x + ncx * y)]) f(j);
            }
    }
};

/// Shells of |k|² below K_lim, sorted: (k2, number of half-plane modes).
std::vector<std::pair<double, std::size_t>> shells_below(double Lx, double Ly, double K_lim) {
    const int nxl = static_cast<int>(K_lim * Lx / (2.0 * PI)) + 1, nyl = static_cast<int>(K_lim * Ly / (2.0 * PI)) + 1;
    std::vector<double> k2s;
    for (int nx = 0; nx <= nxl; ++nx)
        for (int ny = -nyl; ny <= nyl; ++ny) {
            if (nx == 0 && ny <= 0) continue;
            const double kx = 2.0 * PI * nx / Lx, ky = 2.0 * PI * ny / Ly;
            const double k2 = kx * kx + ky * ky;
            if (k2 <= K_lim * K_lim) k2s.push_back(k2);
        }
    std::sort(k2s.begin(), k2s.end());
    std::vector<std::pair<double, std::size_t>> shells;
    for (const double k2 : k2s) {
        if (!shells.empty() && k2 - shells.back().first <= 1e-12 * k2) ++shells.back().second;
        else shells.push_back({k2, 1});
    }
    return shells;
}

}  // namespace

void axis_phases(double u, double L, int n_max, std::complex<double>* out) noexcept {
    const double t = -2.0 * PI * fold(u, L) / L;
    const cplx base = std::polar(1.0, t);
    out[0] = 1.0;
    for (int n = 1; n <= n_max; ++n)
        out[n] = (n % 32 == 0) ? std::polar(1.0, t * n) : out[n - 1] * base;
}

void continuum_init_random(ContinuumConfig& c, double Lx, double Ly, std::size_t N, unsigned seed, double sigma) {
    if (!(Lx > 0.0 && Ly > 0.0) || N < 2)
        throw std::invalid_argument(std::format("need Lx, Ly > 0 and N >= 2 (got {}, {}, {})", Lx, Ly, N));
    c.Lx = Lx;
    c.Ly = Ly;
    c.N = N;
    c.sigma = sigma;
    c.rng.seed(seed);
    c.r.assign(2 * N, 0.0);
    std::uniform_real_distribution<double> ux(0.0, Lx), uy(0.0, Ly);
    if (!(sigma > 0.0)) {
        for (std::size_t i = 0; i < N; ++i) {
            c.r[2 * i] = ux(c.rng);
            c.r[2 * i + 1] = uy(c.rng);
        }
        return;
    }
    // Random sequential addition on a grid of cells of side >= sigma.
    const int ncx = std::max(1, static_cast<int>(Lx / sigma)), ncy = std::max(1, static_cast<int>(Ly / sigma));
    const double hx = Lx / ncx, hy = Ly / ncy;
    std::vector<std::vector<std::size_t>> bins(static_cast<std::size_t>(ncx) * static_cast<std::size_t>(ncy));
    std::size_t placed = 0, attempts = 0;
    while (placed < N) {
        if (++attempts > 1000 * N)
            throw std::invalid_argument(std::format(
                "random sequential addition with sigma = {} placed only {} of {} points", sigma, placed, N));
        const double x = ux(c.rng), y = uy(c.rng);
        const int cx = std::min(ncx - 1, static_cast<int>(x / hx)), cy = std::min(ncy - 1, static_cast<int>(y / hy));
        bool ok = true;
        for (int oy = -1; oy <= 1 && ok; ++oy)
            for (int ox = -1; ox <= 1 && ok; ++ox) {
                const int bx = (cx + ox + ncx) % ncx, by = (cy + oy + ncy) % ncy;
                for (const std::size_t j : bins[static_cast<std::size_t>(bx + ncx * by)]) {
                    const double dx = min_image(x - c.r[2 * j], Lx), dy = min_image(y - c.r[2 * j + 1], Ly);
                    if (dx * dx + dy * dy < sigma * sigma) { ok = false; break; }
                }
            }
        if (!ok) continue;
        c.r[2 * placed] = x;
        c.r[2 * placed + 1] = y;
        bins[static_cast<std::size_t>(cx + ncx * cy)].push_back(placed);
        ++placed;
    }
}

std::size_t continuum_mode_count(double Lx, double Ly, double K) {
    std::size_t n = 0;
    for (const auto& [k2, count] : shells_below(Lx, Ly, K)) n += count;
    return n;
}

ContinuumMaskChoice choose_continuum_mask_for_chi(double Lx, double Ly, std::size_t N, double chi) {
    if (!(chi > 0.0)) throw std::invalid_argument(std::format("chi must be positive, got {}", chi));
    if (N < 2) throw std::invalid_argument("too few points to define chi");
    const double target = chi * 2.0 * static_cast<double>(N - 1);
    // Half-plane modes below K: about K² Lx Ly / (8π).
    double K_lim = std::sqrt(8.0 * PI * target / (Lx * Ly)) * 1.5 + 2.0 * PI * (1.0 / Lx + 1.0 / Ly);
    while (true) {
        const auto shells = shells_below(Lx, Ly, K_lim);
        std::size_t cum = 0, best_n = 0, best_i = 0;
        double best_err = std::numeric_limits<double>::infinity();
        bool passed = false;
        for (std::size_t i = 0; i < shells.size(); ++i) {
            cum += shells[i].second;
            const double err = std::abs(static_cast<double>(cum) - target);
            if (err < best_err) { best_err = err; best_n = cum; best_i = i; }
            if (static_cast<double>(cum) > target) { passed = true; break; }
        }
        if (passed && best_i + 1 < shells.size()) {
            const double K = 0.5 * (std::sqrt(shells[best_i].first) + std::sqrt(shells[best_i + 1].first));
            return {K, best_n, static_cast<double>(best_n) / (2.0 * static_cast<double>(N - 1))};
        }
        K_lim *= 2.0;
    }
}

double continuum_K_of_kappa(double Lx, double Ly, std::size_t N, double kappa) noexcept {
    return 2.0 * PI * kappa * std::sqrt(static_cast<double>(N) / (Lx * Ly));
}

void build_continuum_mask(ContinuumConfig& c, double K, double weight_power) {
    c.K = K;
    c.modes.clear();
    c.kx.clear();
    c.ky.clear();
    const int nxl = static_cast<int>(K * c.Lx / (2.0 * PI)) + 1, nyl = static_cast<int>(K * c.Ly / (2.0 * PI)) + 1;
    c.nx_max = c.ny_max = 0;
    for (int nx = 0; nx <= nxl; ++nx)
        for (int ny = -nyl; ny <= nyl; ++ny) {
            if (nx == 0 && ny <= 0) continue;
            const double kx = 2.0 * PI * nx / c.Lx, ky = 2.0 * PI * ny / c.Ly;
            if (kx * kx + ky * ky > K * K) continue;
            c.modes.push_back({nx, ny});
            c.kx.push_back(kx);
            c.ky.push_back(ky);
            c.nx_max = std::max(c.nx_max, nx);
            c.ny_max = std::max(c.ny_max, std::abs(ny));
        }
    if (c.modes.empty())
        throw std::invalid_argument(std::format("K = {:.6g} contains no mode of the {} x {} box", K, c.Lx, c.Ly));
    const std::size_t n = c.modes.size();
    c.weight_power = weight_power;
    c.weight.assign(n, 1.0);
    if (weight_power != 0.0) {
        double sum = 0.0;
        for (std::size_t m = 0; m < n; ++m) {
            c.weight[m] = std::pow(c.kx[m] * c.kx[m] + c.ky[m] * c.ky[m], -0.5 * weight_power);
            sum += c.weight[m];
        }
        for (double& w : c.weight) w /= sum / static_cast<double>(n);
    }
    c.ft_re.assign(n, 0.0);
    c.ft_im.assign(n, 0.0);
}

void continuum_FT_range(const ContinuumConfig& c, std::size_t m_begin, std::size_t m_end,
                        std::vector<double>& re, std::vector<double>& im) {
    std::vector<cplx> ex(static_cast<std::size_t>(c.nx_max) + 1), ey(static_cast<std::size_t>(c.ny_max) + 1);
    for (std::size_t p = 0; p < c.N; ++p) {
        axis_phases(c.r[2 * p], c.Lx, c.nx_max, ex.data());
        axis_phases(c.r[2 * p + 1], c.Ly, c.ny_max, ey.data());
        for (std::size_t m = m_begin; m < m_end; ++m) {
            const cplx e = mode_phase(c.modes[m], ex.data(), ey.data());
            re[m] += e.real();
            im[m] += e.imag();
        }
    }
}

void continuum_FT(ContinuumConfig& c) {
    std::fill(c.ft_re.begin(), c.ft_re.end(), 0.0);
    std::fill(c.ft_im.begin(), c.ft_im.end(), 0.0);
    continuum_FT_range(c, 0, c.modes.size(), c.ft_re, c.ft_im);
}

double continuum_energy(const ContinuumConfig& c) {
    double E = 0.0;
    for (std::size_t m = 0; m < c.modes.size(); ++m) E += c.ft_re[m] * c.ft_re[m] + c.ft_im[m] * c.ft_im[m];
    return 2.0 * E;
}

double continuum_weighted_energy(const ContinuumConfig& c) {
    double E = 0.0;
    for (std::size_t m = 0; m < c.modes.size(); ++m)
        E += c.weight[m] * (c.ft_re[m] * c.ft_re[m] + c.ft_im[m] * c.ft_im[m]);
    return 2.0 * E;
}

double continuum_residual(const ContinuumConfig& c, double E) {
    return E / (2.0 * static_cast<double>(c.modes.size()) * static_cast<double>(c.N));
}

double continuum_band_residual(const ContinuumConfig& c, double k_max) {
    double sum = 0.0;
    std::size_t n = 0;
    for (std::size_t m = 0; m < c.modes.size(); ++m) {
        if (c.kx[m] * c.kx[m] + c.ky[m] * c.ky[m] > k_max * k_max) continue;
        sum += c.ft_re[m] * c.ft_re[m] + c.ft_im[m] * c.ft_im[m];
        ++n;
    }
    return n ? sum / static_cast<double>(n) / static_cast<double>(c.N) : 0.0;
}

void continuum_gradient_range(const ContinuumConfig& c, std::size_t p_begin, std::size_t p_end, std::vector<double>& g) {
    std::vector<cplx> ex(static_cast<std::size_t>(c.nx_max) + 1), ey(static_cast<std::size_t>(c.ny_max) + 1);
    for (std::size_t p = p_begin; p < p_end; ++p) {
        axis_phases(c.r[2 * p], c.Lx, c.nx_max, ex.data());
        axis_phases(c.r[2 * p + 1], c.Ly, c.ny_max, ey.data());
        double gx = 0.0, gy = 0.0;
        for (std::size_t m = 0; m < c.modes.size(); ++m) {
            const cplx e = mode_phase(c.modes[m], ex.data(), ey.data());
            const double s = c.weight[m] * (c.ft_re[m] * e.imag() - c.ft_im[m] * e.real());  // w Im(conj(F) e)
            gx += c.kx[m] * s;
            gy += c.ky[m] * s;
        }
        g[2 * p] = 4.0 * gx;
        g[2 * p + 1] = 4.0 * gy;
    }
}

double continuum_spectral_energy_grad(ContinuumConfig& c, std::vector<double>& g) {
    continuum_FT(c);
    g.assign(2 * c.N, 0.0);
    continuum_gradient_range(c, 0, c.N, g);
    return continuum_weighted_energy(c);
}

namespace {

/// Adds the pair (i, j) at separation (dx, dy) to the hard-core energy and gradient.
inline double core_pair(const ContinuumConfig& c, std::size_t i, std::size_t j, double dx, double dy,
                        std::vector<double>* g) noexcept {
    const double d2 = dx * dx + dy * dy;
    if (d2 >= c.sigma * c.sigma) return 0.0;
    const double d = std::sqrt(d2);
    const double u = 1.0 - d / c.sigma;
    if (g && d > 0.0) {
        const double f = -2.0 * c.core_weight * u / (c.sigma * d);  // dE/dd / d
        (*g)[2 * i] += f * dx;
        (*g)[2 * i + 1] += f * dy;
        (*g)[2 * j] -= f * dx;
        (*g)[2 * j + 1] -= f * dy;
    }
    return c.core_weight * u * u;
}

}  // namespace

double continuum_core_energy(const ContinuumConfig& c, std::vector<double>* g) {
    if (!(c.sigma > 0.0)) return 0.0;
    Cells cells;
    double E = 0.0;
    if (!cells.build(c, c.sigma)) {  // box narrower than 3 sigma on an axis: direct search
        for (std::size_t i = 0; i < c.N; ++i)
            for (std::size_t j = i + 1; j < c.N; ++j)
                E += core_pair(c, i, j, min_image(c.r[2 * i] - c.r[2 * j], c.Lx),
                               min_image(c.r[2 * i + 1] - c.r[2 * j + 1], c.Ly), g);
        return E;
    }
    for (std::size_t i = 0; i < c.N; ++i)
        cells.neighbours(c, i, [&](std::size_t j) {
            if (j <= i) return;
            E += core_pair(c, i, j, min_image(c.r[2 * i] - c.r[2 * j], c.Lx),
                           min_image(c.r[2 * i + 1] - c.r[2 * j + 1], c.Ly), g);
        });
    return E;
}

double continuum_min_distance(const ContinuumConfig& c) {
    if (c.N < 2) return std::numeric_limits<double>::infinity();
    double h = 0.75 * std::sqrt(c.Lx * c.Ly / static_cast<double>(c.N));
    while (true) {
        Cells cells;
        if (!cells.build(c, h)) {
            double best2 = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < c.N; ++i)
                for (std::size_t j = i + 1; j < c.N; ++j) {
                    const double dx = min_image(c.r[2 * i] - c.r[2 * j], c.Lx), dy = min_image(c.r[2 * i + 1] - c.r[2 * j + 1], c.Ly);
                    best2 = std::min(best2, dx * dx + dy * dy);
                }
            return std::sqrt(best2);
        }
        double best2 = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < c.N; ++i)
            cells.neighbours(c, i, [&](std::size_t j) {
                if (j == i) return;
                const double dx = min_image(c.r[2 * i] - c.r[2 * j], c.Lx), dy = min_image(c.r[2 * i + 1] - c.r[2 * j + 1], c.Ly);
                best2 = std::min(best2, dx * dx + dy * dy);
            });
        if (best2 <= std::min(cells.hx, cells.hy) * std::min(cells.hx, cells.hy)) return std::sqrt(best2);
        h *= 2.0;
    }
}

double continuum_mask_K_midpoint(const ContinuumConfig& c) {
    double k2_in = 0.0;
    for (std::size_t m = 0; m < c.modes.size(); ++m) k2_in = std::max(k2_in, c.kx[m] * c.kx[m] + c.ky[m] * c.ky[m]);
    const double step = 2.0 * PI / std::min(c.Lx, c.Ly);
    for (double lim = std::sqrt(k2_in) + 2.0 * step;; lim += 4.0 * step)
        for (const auto& [k2, count] : shells_below(c.Lx, c.Ly, lim))
            if (k2 > k2_in * (1.0 + 1e-12)) return 0.5 * (std::sqrt(k2_in) + std::sqrt(k2));
}

double continuum_energy_direct(const ContinuumConfig& c) {
    const int nxl = static_cast<int>(c.K * c.Lx / (2.0 * PI)) + 1, nyl = static_cast<int>(c.K * c.Ly / (2.0 * PI)) + 1;
    double E = 0.0;
    for (int nx = -nxl; nx <= nxl; ++nx)
        for (int ny = -nyl; ny <= nyl; ++ny) {
            if (nx == 0 && ny == 0) continue;
            const double kx = 2.0 * PI * nx / c.Lx, ky = 2.0 * PI * ny / c.Ly;
            if (kx * kx + ky * ky > c.K * c.K) continue;
            cplx sum{0.0, 0.0};
            for (std::size_t p = 0; p < c.N; ++p) sum += std::polar(1.0, -(kx * c.r[2 * p] + ky * c.r[2 * p + 1]));
            E += std::norm(sum);
        }
    return E;
}

double continuum_core_energy_direct(const ContinuumConfig& c) {
    if (!(c.sigma > 0.0)) return 0.0;
    double E = 0.0;
    for (std::size_t i = 0; i < c.N; ++i)
        for (std::size_t j = i + 1; j < c.N; ++j)
            E += core_pair(c, i, j, min_image(c.r[2 * i] - c.r[2 * j], c.Lx),
                           min_image(c.r[2 * i + 1] - c.r[2 * j + 1], c.Ly), nullptr);
    return E;
}
