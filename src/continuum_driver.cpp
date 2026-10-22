#include "continuum_driver.hpp"
#include "constants.hpp"
#include "continuum.hpp"
#include "continuum_mc.hpp"
#include "minimize.hpp"
#include "obs.hpp"
#include "thread_pool.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

/// rho_hat of every mode, the modes split among the threads (each thread writes its own range, and each
/// mode sums the points in the same order, so the result does not depend on the thread count).
void parallel_FT(ContinuumConfig& c, ThreadPool* pool) {
    std::fill(c.ft_re.begin(), c.ft_re.end(), 0.0);
    std::fill(c.ft_im.begin(), c.ft_im.end(), 0.0);
    if (!pool) { continuum_FT_range(c, 0, c.modes.size(), c.ft_re, c.ft_im); return; }
    auto job = [&c](std::size_t b, std::size_t e) { continuum_FT_range(c, b, e, c.ft_re, c.ft_im); return 0.0; };
    (void)pool->parallel_reduce(c.modes.size(), job);
}

/// Gradient of E_w, the points split among the threads.
void parallel_gradient(const ContinuumConfig& c, std::vector<double>& g, ThreadPool* pool) {
    g.assign(2 * c.N, 0.0);
    if (!pool) { continuum_gradient_range(c, 0, c.N, g); return; }
    auto job = [&c, &g](std::size_t b, std::size_t e) { continuum_gradient_range(c, b, e, g); return 0.0; };
    (void)pool->parallel_reduce(c.N, job);
}

double fold(double u, double L) {
    double f = u - L * std::floor(u / L);
    if (f >= L || f < 0.0) f = 0.0;
    return f;
}

struct RadialBins { std::vector<double> k, S; std::vector<long long> count; };

/// Radially averaged S(k) = |rho_hat|² / N over the half plane up to k_max (S(-k) = S(k)), bins of width
/// 2π/sqrt(Lx Ly) centred on its multiples.
RadialBins radial_sk(const ContinuumConfig& c, double k_max, ThreadPool* pool) {
    ContinuumConfig s;
    s.Lx = c.Lx; s.Ly = c.Ly; s.N = c.N; s.r = c.r;
    build_continuum_mask(s, k_max, 0.0);
    parallel_FT(s, pool);
    const double dk = 2.0 * PI / std::sqrt(c.Lx * c.Ly);
    const auto nb = static_cast<std::size_t>(k_max / dk) + 2;
    std::vector<double> ss(nb, 0.0), sk(nb, 0.0);
    std::vector<long long> cnt(nb, 0);
    for (std::size_t m = 0; m < s.modes.size(); ++m) {
        const double k = std::sqrt(s.kx[m] * s.kx[m] + s.ky[m] * s.ky[m]);
        const auto b = static_cast<std::size_t>(std::lround(k / dk));
        ss[b] += (s.ft_re[m] * s.ft_re[m] + s.ft_im[m] * s.ft_im[m]) / static_cast<double>(c.N);
        sk[b] += k;
        ++cnt[b];
    }
    RadialBins out;
    for (std::size_t b = 0; b < nb; ++b) {
        if (cnt[b] == 0) continue;
        out.k.push_back(sk[b] / static_cast<double>(cnt[b]));
        out.S.push_back(ss[b] / static_cast<double>(cnt[b]));
        out.count.push_back(2 * cnt[b]);  // both members of each ±k pair
    }
    return out;
}

}  // namespace

void run_continuum(const Parameters& p, unsigned seed, int arr) {
    const bool several = p.realizations > 1;
    const std::size_t N = p.continuum_points();
    const double sigma = p.continuum_core_distance();

    ContinuumConfig c;
    if (p.init_file.empty()) {
        continuum_init_random(c, p.Lx, p.Ly, N, seed, sigma);
    } else {
        const ObsData init = read_obs(p.init_file);
        if (init.points.size() != N)
            throw std::runtime_error(std::format("init_file has {} points, expected {}", init.points.size(), N));
        c.Lx = p.Lx;
        c.Ly = p.Ly;
        c.N = N;
        c.sigma = sigma;
        c.rng.seed(seed);
        c.r.resize(2 * N);
        for (std::size_t i = 0; i < N; ++i) {
            c.r[2 * i] = init.points[i].x;
            c.r[2 * i + 1] = init.points[i].y;
        }
        std::cout << std::format("Start from '{}' (minimum distance {:.6g})\n", p.init_file, continuum_min_distance(c));
    }
    c.core_weight = p.core_weight;

    double K = 0.0;
    if (p.kappa > 0.0) {
        K = continuum_K_of_kappa(p.Lx, p.Ly, N, p.kappa);
    } else {
        K = choose_continuum_mask_for_chi(p.Lx, p.Ly, N, p.chi).K;
    }
    build_continuum_mask(c, K, p.weight_power);
    const std::size_t M = c.modes.size();
    const double chi = chi_point_pattern(M, N);
    const double kappa = K / (2.0 * PI * std::sqrt(c.density()));
    std::cout << std::format("Continuum: {} points in {} x {} (rho_O = {:.9g}), mean spacing a = {:.6g}\n", N, p.Lx,
                             p.Ly, c.density(), 1.0 / std::sqrt(c.density()));
    std::cout << std::format("Mask: K = {:.9g} d^-1, {} half-plane modes, index extents {} and {}; kappa = {:.6f}, "
                             "chi point_pattern = {:.6f}{}\n", K, M, c.nx_max, c.ny_max, kappa, chi,
                             p.kappa > 0.0 ? "" : std::format(" (target {})", p.chi));
    if (chi > 0.5)
        std::cerr << std::format("Warning: chi point_pattern = {:.4f} > 1/2: stealthy disordered ground states are not "
                                 "expected; the minimization may stall and the pattern tends to crystallize.\n", chi);
    if (sigma > 0.0)
        std::cout << std::format("Hard core: no two points closer than {:.9g} (2R = {}){}\n", sigma, 2.0 * p.R,
                                 p.init_file.empty() ? ", random sequential addition start" : "");

    std::unique_ptr<ThreadPool> pool;
    if (p.num_threads != 1 && static_cast<double>(N) * static_cast<double>(M) >= 1e6) {
        pool = std::make_unique<ThreadPool>(p.num_threads);
        if (pool->get_thread_count() <= 1) pool.reset();
    }
    std::cout << (pool ? std::format("Using {} threads.\n", pool->get_thread_count()) : std::string("Sequential.\n"));

    // Objective: weighted spectral residual E_w / (2 M N) (about 1 for a Poisson start) plus the hard core.
    const double norm = 1.0 / (2.0 * static_cast<double>(M) * static_cast<double>(N));
    std::vector<double> gs;
    const Objective f = [&](const std::vector<double>& x, std::vector<double>& g) {
        c.r = x;
        parallel_FT(c, pool.get());
        const double Ew = continuum_weighted_energy(c);
        parallel_gradient(c, gs, pool.get());
        g.resize(x.size());
        for (std::size_t i = 0; i < g.size(); ++i) g[i] = norm * gs[i];
        return norm * Ew + continuum_core_energy(c, &g);
    };

    continuum_FT(c);
    const double res0 = continuum_residual(c, continuum_energy(c));
    std::cout << std::format("Initial residual <S>_mask = {:.6g}\n", res0);

    const std::string log_path = path_with_arr(p.annealing_file, arr, several);
    std::ofstream log(log_path);
    std::string method;

    if (p.minimizer != "lbfgs") {  // Monte Carlo (alone, or before L-BFGS)
        const auto tm = std::chrono::steady_clock::now();
        if (log) log << "# Monte Carlo annealing, single-point moves\n";
        const ContinuumMcReport mc = continuum_anneal(c, p, pool.get(), log ? &log : nullptr);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tm).count();
        const double acc = mc.attempts ? static_cast<double>(mc.accepted) / static_cast<double>(mc.attempts) : 0.0;
        std::cout << std::format("Monte Carlo: {} levels, {} proposals, acceptance {:.3f}, {:.1f} s; objective {:.4g} -> {:.4g}, "
                                 "final step {:.3g}, max drift {:.2g}\n", mc.levels, mc.attempts, acc, sec,
                                 mc.f_initial, mc.f_final, mc.step_final, mc.max_drift);
        method = std::format("Monte Carlo {} proposals (acceptance {:.3f})", mc.attempts, acc);
    }

    std::vector<double> x = c.r;
    LbfgsResult lr;
    lr.reason = "not run";
    if (p.minimizer != "anneal") {
        if (log) log << std::format("# L-BFGS\n{:<12} {:<24} {:<24}\n", "# iteration", "objective", "max|gradient|");
        LbfgsOptions opt;
        opt.memory = p.lbfgs_memory;
        opt.max_iterations = p.max_iterations;
        opt.grad_tol = p.grad_tol;
        opt.f_target = p.residual_tol;
        const auto t0 = std::chrono::steady_clock::now();
        const Progress progress = [&](long long it, double fx, double gmax) {
            if (log) log << std::format("{:<12} {:<24.10e} {:<24.10e}\n", it, fx, gmax);
            return false;
        };
        // The objective is the weighted residual; the target is on the unweighted one. If the weighted target
        // is reached first, continue from there with a lower target (at most a few rounds).
        long long total_it = 0, total_ev = 0;
        for (int round = 0; round < 6; ++round) {
            lr = lbfgs_minimize(f, x, opt, progress);
            total_it += lr.iterations;
            total_ev += lr.evaluations;
            if (lr.reason != "f_target") break;
            c.r = x;
            continuum_FT(c);
            const double res = continuum_residual(c, continuum_energy(c));
            if (res <= p.residual_tol || continuum_core_energy(c, nullptr) > 0.0) break;
            opt.f_target = lr.f * 0.5 * p.residual_tol / res;
            opt.max_iterations = std::max<long long>(1, opt.max_iterations - lr.iterations);
        }
        lr.iterations = total_it;
        lr.evaluations = total_ev;
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << std::format("L-BFGS: {} iterations, {} evaluations, {:.1f} s, stopped by {}; objective {:.4g}\n",
                                 lr.iterations, lr.evaluations, seconds, lr.reason, lr.f);
        method += std::format("{}L-BFGS {} iterations ({})", method.empty() ? "" : " + ", lr.iterations, lr.reason);
    }

    // Positions folded into the box and rounded to what the file will hold.
    c.r = x;
    for (std::size_t i = 0; i < N; ++i) {
        c.r[2 * i] = parse_real(format_real(fold(c.r[2 * i], p.Lx)));
        c.r[2 * i + 1] = parse_real(format_real(fold(c.r[2 * i + 1], p.Ly)));
    }
    continuum_FT(c);
    const double residual = continuum_residual(c, continuum_energy(c));
    const double low = continuum_band_residual(c, K / 4.0);
    const double dmin_torus = continuum_min_distance(c);
    std::cout << std::format("Residual <S>_mask = {:.4g} (written coordinates); |k| <= K/4: {:.4g}\n", residual, low);
    std::cout << std::format("Minimum distance (periodic in x and y): {:.9g} = {:.4f} a\n", dmin_torus,
                             dmin_torus * std::sqrt(c.density()));
    if (residual > p.residual_tol)
        std::cerr << std::format("Warning: residual {:.4g} above residual_tol = {:g} ({})\n", residual, p.residual_tol, method);

    const double k_max = p.sk_kmax > 0.0 ? p.sk_kmax : 3.0 * K;
    const RadialBins sk = radial_sk(c, k_max, pool.get());
    {
        const std::string path = path_with_arr(p.sk_file, arr, several);
        std::ofstream out(path);
        if (!out) {
            std::cerr << std::format("Error: could not create '{}'\n", path);
        } else {
            out << std::format("# Radially averaged structure factor S(k) = |rho_hat(k)|^2 / N, N = {} points, "
                               "box {} x {}\n# mask radius K = {:.9g} d^-1; Poisson: S = 1\n", N, p.Lx, p.Ly, K);
            out << std::format("{:<14} {:<20} {:<8}\n", "# k (d^-1)", "S(k)", "modes");
            for (std::size_t i = 0; i < sk.k.size(); ++i)
                out << std::format("{:<14.8f} {:<20.10e} {:<8d}\n", sk.k[i], sk.S[i], sk.count[i]);
            std::cout << std::format("Structure factor saved to '{}'\n", path);
        }
    }

    if (p.obs_file.empty()) return;
    std::vector<ObsPoint> pts(N);
    for (std::size_t i = 0; i < N; ++i) pts[i] = {c.r[2 * i], c.r[2 * i + 1]};
    ObsMeta meta;
    meta.generator_version = HYPER2D_VERSION;
    meta.mode = "continuum";
    meta.cls = "stealthy";
    meta.Lx = p.Lx;
    meta.Ly = p.Ly;
    meta.boundary = p.boundary;
    meta.R = p.R;
    meta.seed = seed;
    meta.K = continuum_mask_K_midpoint(c);
    meta.n_constrained_modes = M;
    meta.chi_convention = "point_pattern";
    meta.chi = chi;
    meta.residual_S = residual;
    meta.notes = std::format("continuum, {}, weight_power {}, residual_low_k {:.4g}, hard core {:.9g}{}", method,
                             p.weight_power, low, sigma, p.init_file.empty() ? "" : ", start " + p.init_file);
    const ObsResult obs = build_obs(meta, pts);
    const std::string obs_path = p.obs_file == "auto"
        ? obs_auto_name(meta.cls, obs.rho_O, parse_real(format_real(meta.Lx)), parse_real(format_real(meta.Ly)), arr)
        : path_with_arr(p.obs_file, arr, several);
    {
        std::ofstream fo(obs_path, std::ios::binary);
        if (!fo || !(fo << obs.text)) throw std::runtime_error(std::format("could not write '{}'", obs_path));
    }
    if (obs.d_min < 2.0 * p.R)
        std::cerr << std::format("Warning: d_min = {} < 2R = {}: the file violates the no-overlap rule\n", obs.d_min, 2.0 * p.R);
    std::cout << std::format("Obstacle file saved to '{}' (n = {}, rho_O = {}, d_min = {}, sha256_body {})\n", obs_path,
                             obs.n, format_real(obs.rho_O), format_real(obs.d_min), obs.sha256_body);
    std::vector<double> k_d = sk.k;
    const std::string sk_path = path_with_extension(obs_path, ".sk");
    write_sk_companion(sk_path, {
        {"format", "1"}, {"generator", "hyperuniform2D"}, {"generator_version", HYPER2D_VERSION},
        {"obstacles_file", obs_path.substr(obs_path.find_last_of('/') + 1)}, {"sha256_body", obs.sha256_body},
        {"columns", "k S"}, {"k_units", "d^-1"},
        {"notes", std::format("radially averaged S(k) = |rho_hat(k)|^2 / n for 0 < |k| <= {:.9g}, bins of width "
                              "2pi/sqrt(Lx Ly)", k_max)},
    }, k_d, sk.S);
    std::cout << std::format("Structure factor companion saved to '{}'\n", sk_path);
}
