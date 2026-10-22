#include "continuum.hpp"
#include "constants.hpp"
#include "thread_pool.hpp"

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
    c.runs.clear();
    for (std::size_t m = 0; m < n; ++m) {
        const ModeK& k = c.modes[m];
        if (!c.runs.empty() && c.runs.back().nx == k.kx &&
            c.runs.back().ny_lo + static_cast<int>(c.runs.back().count) == k.ky)
            ++c.runs.back().count;
        else
            c.runs.push_back({k.kx, k.ky, m, 1});
    }
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

namespace {

/// Axis tables of one point as real arrays: ex for nx = 0..nx_max, ey for ny = -ny_max..ny_max at
/// index ny + ny_max (negative ny by conjugation).
struct AxisTables {
    std::vector<cplx> tx, ty;
    std::vector<double> exr, exi, eyr, eyi;
    explicit AxisTables(const ContinuumConfig& c)
        : tx(static_cast<std::size_t>(c.nx_max) + 1), ty(static_cast<std::size_t>(c.ny_max) + 1),
          exr(tx.size()), exi(tx.size()), eyr(2 * ty.size() - 1), eyi(2 * ty.size() - 1) {}
    void fill(const ContinuumConfig& c, double x, double y) {
        axis_phases(x, c.Lx, c.nx_max, tx.data());
        axis_phases(y, c.Ly, c.ny_max, ty.data());
        for (std::size_t n = 0; n < tx.size(); ++n) { exr[n] = tx[n].real(); exi[n] = tx[n].imag(); }
        const auto o = static_cast<std::size_t>(c.ny_max);
        for (std::size_t n = 0; n < ty.size(); ++n) {
            eyr[o + n] = ty[n].real(); eyi[o + n] = ty[n].imag();
            eyr[o - n] = ty[n].real(); eyi[o - n] = -ty[n].imag();
        }
    }
};

}  // namespace

void continuum_FT_points(const ContinuumConfig& c, std::size_t p_begin, std::size_t p_end, double* re, double* im) {
    AxisTables t(c);
    double* __restrict out_r = re;
    double* __restrict out_i = im;
    for (std::size_t p = p_begin; p < p_end; ++p) {
        t.fill(c, c.r[2 * p], c.r[2 * p + 1]);
        for (const auto& run : c.runs) {
            const double xr = t.exr[static_cast<std::size_t>(run.nx)], xi = t.exi[static_cast<std::size_t>(run.nx)];
            // ey index of mode m: ny_lo + (m - begin) + ny_max
            const double* __restrict yr = t.eyr.data() + (run.ny_lo + c.ny_max) - static_cast<std::ptrdiff_t>(run.begin);
            const double* __restrict yi = t.eyi.data() + (run.ny_lo + c.ny_max) - static_cast<std::ptrdiff_t>(run.begin);
            const std::size_t e = run.begin + run.count;
            for (std::size_t m = run.begin; m < e; ++m) {
                out_r[m] += xr * yr[m] - xi * yi[m];
                out_i[m] += xr * yi[m] + xi * yr[m];
            }
        }
    }
}

void continuum_FT(ContinuumConfig& c, ThreadPool* pool) {
    const std::size_t M = c.modes.size();
    const std::size_t groups = std::min<std::size_t>(32, std::max<std::size_t>(1, c.N));
    c.part_re.assign(groups * M, 0.0);
    c.part_im.assign(groups * M, 0.0);
    auto by_group = [&c, M, groups](std::size_t b, std::size_t e) {
        for (std::size_t q = b; q < e; ++q)
            continuum_FT_points(c, q * c.N / groups, (q + 1) * c.N / groups, c.part_re.data() + q * M, c.part_im.data() + q * M);
        return 0.0;
    };
    auto add = [&c, M, groups](std::size_t b, std::size_t e) {
        for (std::size_t m = b; m < e; ++m) {
            double r = 0.0, i = 0.0;
            for (std::size_t q = 0; q < groups; ++q) { r += c.part_re[q * M + m]; i += c.part_im[q * M + m]; }
            c.ft_re[m] = r;
            c.ft_im[m] = i;
        }
        return 0.0;
    };
    if (pool) {
        (void)pool->parallel_reduce(groups, by_group);
        (void)pool->parallel_reduce(M, add);
    } else {
        (void)by_group(0, groups);
        (void)add(0, M);
    }
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
    // w Im(conj(F) e) = Im(A e) with A = w conj(F); along a run e = ex ey, so
    //   Σ_run kx Im(A e) = kx Im(ex S),  Σ_run ky Im(A e) = Im(ex T),  S = Σ A ey,  T = Σ ky A ey.
    // Blocks of B points: tables stored [index][point], and the point is the inner loop.
    constexpr std::size_t B = 32;
    const std::size_t n = c.modes.size();
    std::vector<double> ar(n), ai(n);
    for (std::size_t m = 0; m < n; ++m) { ar[m] = c.weight[m] * c.ft_re[m]; ai[m] = -c.weight[m] * c.ft_im[m]; }
    const std::size_t nxs = static_cast<std::size_t>(c.nx_max) + 1, nys = 2 * static_cast<std::size_t>(c.ny_max) + 1;
    std::vector<double> exr(nxs * B), exi(nxs * B), eyr(nys * B), eyi(nys * B);
    std::vector<double> Sr(B), Si(B), Tr(B), Ti(B), gx(B), gy(B);
    AxisTables t(c);
    for (std::size_t p0 = p_begin; p0 < p_end; p0 += B) {
        const std::size_t nb = std::min(B, p_end - p0);
        for (std::size_t q = 0; q < nb; ++q) {
            t.fill(c, c.r[2 * (p0 + q)], c.r[2 * (p0 + q) + 1]);
            for (std::size_t k = 0; k < nxs; ++k) { exr[k * B + q] = t.exr[k]; exi[k * B + q] = t.exi[k]; }
            for (std::size_t k = 0; k < nys; ++k) { eyr[k * B + q] = t.eyr[k]; eyi[k * B + q] = t.eyi[k]; }
        }
        std::fill(gx.begin(), gx.end(), 0.0);
        std::fill(gy.begin(), gy.end(), 0.0);
        double* __restrict sr = Sr.data();
        double* __restrict si = Si.data();
        double* __restrict tr = Tr.data();
        double* __restrict ti = Ti.data();
        for (const auto& run : c.runs) {
            std::fill(Sr.begin(), Sr.end(), 0.0);
            std::fill(Si.begin(), Si.end(), 0.0);
            std::fill(Tr.begin(), Tr.end(), 0.0);
            std::fill(Ti.begin(), Ti.end(), 0.0);
            for (std::size_t j = 0; j < run.count; ++j) {
                const std::size_t m = run.begin + j;
                const double a_r = ar[m], a_i = ai[m], kym = c.ky[m];
                const std::size_t row = static_cast<std::size_t>(run.ny_lo + static_cast<int>(j) + c.ny_max) * B;
                const double* __restrict yr = eyr.data() + row;
                const double* __restrict yi = eyi.data() + row;
                for (std::size_t q = 0; q < B; ++q) {
                    const double qr = a_r * yr[q] - a_i * yi[q];
                    const double qi = a_r * yi[q] + a_i * yr[q];
                    sr[q] += qr; si[q] += qi;
                    tr[q] += kym * qr; ti[q] += kym * qi;
                }
            }
            const double kx = c.kx[run.begin];
            const double* __restrict xr = exr.data() + static_cast<std::size_t>(run.nx) * B;
            const double* __restrict xi = exi.data() + static_cast<std::size_t>(run.nx) * B;
            for (std::size_t q = 0; q < B; ++q) {
                gx[q] += kx * (xr[q] * si[q] + xi[q] * sr[q]);
                gy[q] += xr[q] * ti[q] + xi[q] * tr[q];
            }
        }
        for (std::size_t q = 0; q < nb; ++q) {
            g[2 * (p0 + q)] = 4.0 * gx[q];
            g[2 * (p0 + q) + 1] = 4.0 * gy[q];
        }
    }
}

void continuum_gradient(const ContinuumConfig& c, std::vector<double>& g, ThreadPool* pool) {
    g.assign(2 * c.N, 0.0);
    auto job = [&c, &g](std::size_t b, std::size_t e) { continuum_gradient_range(c, b, e, g); return 0.0; };
    if (pool) (void)pool->parallel_reduce(c.N, job);
    else (void)job(0, c.N);
}

double continuum_spectral_energy_grad(ContinuumConfig& c, std::vector<double>& g, ThreadPool* pool) {
    continuum_FT(c, pool);
    continuum_gradient(c, g, pool);
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
