// Unit tests for hyperuniform2D. No framework: each check prints its result and the exit code
// counts failures (run with ctest, or ./build/test_hyper2D).
#include "annealing.hpp"
#include "config.hpp"
#include "constants.hpp"
#include "continuum.hpp"
#include "continuum_mc.hpp"
#include "fourier.hpp"
#include "minimize.hpp"
#include "obs.hpp"
#include "parametros.hpp"
#include "sha256.hpp"
#include "thread_pool.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
    if (!ok) ++failures;
}

bool close(double a, double b, double rel) {
    return std::abs(a - b) <= rel * std::max({1.0, std::abs(a), std::abs(b)});
}

HyperuniformConfig make(int N, double density, double K2, unsigned seed) {
    HyperuniformConfig c;
    initialize_grid(c, N, density, seed);
    build_mask(c, K2);
    calculate_FT(c);
    return c;
}

bool lists_consistent(const HyperuniformConfig& c) {
    if (static_cast<int>(c.occupied.size()) != c.M) return false;
    if (c.occupied.size() + c.empty.size() != c.grid.size()) return false;
    for (std::size_t i = 0; i < c.occupied.size(); ++i)
        if (c.grid[c.occupied[i]] != 1 || c.list_pos[c.occupied[i]] != i) return false;
    for (std::size_t i = 0; i < c.empty.size(); ++i)
        if (c.grid[c.empty[i]] != 0 || c.list_pos[c.empty[i]] != i) return false;
    return true;
}

std::filesystem::path tmp_file(const std::string& name) {
    return std::filesystem::temp_directory_path() / std::format("hyper2D_test_{}", name);
}

}  // namespace

int main() {
    std::cout << "occupied count and chi\n";
    check(occupied_count(100, 0.29) == 2900, "round, not truncate: 0.29 * 100^2 = 2900 (truncation gave 2899)");
    check(occupied_count(64, 0.4) == 1638, "0.4 * 64^2 = 1638.4 -> 1638");
    check(chi_points(10, 30) == 30 && chi_points(10, 70) == 30, "chi uses the minority species");
    check(close(chi_of(50, 101), 0.5, 1e-15), "chi = modes / (n - 1)");
    {
        const MaskChoice mc = choose_mask_for_chi(64, 1001, 0.3);
        HyperuniformConfig c;
        initialize_grid(c, 64, 0.25, 1);
        build_mask(c, mc.K2);
        check(c.modes.size() == mc.n_modes, "choose_mask_for_chi: reported modes = modes in the built mask");
        check(std::abs(static_cast<double>(mc.n_modes) - 300.0) <= 12.0, std::format(
            "choose_mask_for_chi: {} modes for target 300 (shell granularity)", mc.n_modes));
        bool threw = false;
        try { (void)choose_mask_for_chi(64, 1001, 0.0); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "chi = 0 rejected by choose_mask_for_chi");
        threw = false;
        try { build_mask(c, 31.0 * 31.0 + 1.0); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "K beyond N/2 - 1 rejected (Nyquist)");
    }

    std::cout << "energy\n";
    {
        auto c = make(32, 0.3, 30.0, 7);
        const double E = energy(c), Ed = energy_direct(c);
        check(close(E, Ed, 1e-10), std::format("half-plane tables = full-mask trigonometry ({:.10g} vs {:.10g})", E, Ed));

        // complement: rho_hat changes sign for k != 0, so E is the same
        HyperuniformConfig comp;
        initialize_grid(comp, 32, 0.7, 1);
        comp.occupied = c.empty;
        comp.empty = c.occupied;
        comp.M = static_cast<int>(comp.occupied.size());
        build_mask(comp, 30.0);
        calculate_FT(comp);
        check(close(energy(comp), E, 1e-10), "energy invariant under complement");
    }

    std::cout << "incremental updates\n";
    {
        auto c = make(40, 0.35, 50.0, 11);
        double E = energy(c);
        int accepted = 0;
        for (int s = 0; s < 20000; ++s)
            if (auto dE = attempt_swap(c, 1e9, nullptr)) { E += *dE; ++accepted; }
        check(accepted > 15000, std::format("every call is a real proposal: {} of 20000 accepted at T = 1e9", accepted));
        check(lists_consistent(c), "grid, occupied/empty lists and positions consistent after swaps");
        const double Ed = energy_direct(c);
        check(close(E, Ed, 1e-8), std::format("accumulated E = direct E after {} swaps ({:.10g} vs {:.10g})", accepted, E, Ed));
        calculate_FT(c);
        check(close(energy(c), Ed, 1e-10), "incremental rho_hat = recomputed rho_hat");
    }

    std::cout << "weighted energy\n";
    {
        HyperuniformConfig c;
        initialize_grid(c, 40, 0.3, 21);
        build_mask(c, 60.0, 2.0);
        calculate_FT(c);
        double wsum = 0.0;
        for (double w : c.weight) wsum += w;
        check(close(wsum / static_cast<double>(c.weight.size()), 1.0, 1e-12), "weights have mean 1");
        bool decreasing = true;
        for (std::size_t m = 0; m < c.modes.size(); ++m) {
            const int k2 = c.modes[m].kx * c.modes[m].kx + c.modes[m].ky * c.modes[m].ky;
            decreasing = decreasing && close(c.weight[m] * k2, c.weight[0] * (c.modes[0].kx * c.modes[0].kx + c.modes[0].ky * c.modes[0].ky), 1e-12);
        }
        check(decreasing, "w(k) proportional to |k|^-2 for weight_power = 2");
        double Ew = weighted_energy(c);
        for (int s = 0; s < 5000; ++s)
            if (auto dE = attempt_swap(c, 1e9, nullptr, 0.5, 1)) Ew += *dE;
        calculate_FT(c);
        check(close(Ew, weighted_energy(c), 1e-8), "accumulated weighted dE = recomputed weighted energy");
        check(close(energy(c), energy_direct(c), 1e-10), "unweighted energy still exact with weights on");
        auto u = make(40, 0.3, 60.0, 21);
        check(close(weighted_energy(u), energy(u), 1e-14), "weight_power = 0: weighted energy = energy");
    }

    std::cout << "thread pool\n";
    {
        auto a = make(64, 0.4, 400.0, 5);
        auto b = a;
        ThreadPool pool(4);
        const double d_seq = delta_energy_range(a, 3, 5, 40, 17, 0, a.modes.size());
        const double d_par = pool.parallel_reduce_delta_E(b, 3, 5, 40, 17);
        check(close(d_seq, d_par, 1e-12), std::format("parallel dE = sequential dE ({} modes)", a.modes.size()));
        bool same = true;
        for (std::size_t m = 0; m < a.modes.size(); ++m)
            same = same && a.d_re[m] == b.d_re[m] && a.d_im[m] == b.d_im[m];
        check(same, "parallel delta buffer = sequential delta buffer");
    }
    {
        ThreadPool pool(4);
        auto sum_range = [](std::size_t b, std::size_t e) {
            double s = 0.0;
            for (std::size_t i = b; i < e; ++i) s += static_cast<double>(i);
            return s;
        };
        bool ok = true;
        for (std::size_t n : {0u, 1u, 3u, 4u, 5u, 1001u})
            ok = ok && pool.parallel_reduce(n, sum_range) == static_cast<double>(n * (n - (n ? 1 : 0)) / 2);
        check(ok, "generic parallel_reduce: sum over [0, n) for n = 0, 1, 3, 4, 5, 1001");
        std::vector<int> hits(1000, 0);
        auto mark = [&](std::size_t b, std::size_t e) {
            for (std::size_t i = b; i < e; ++i) ++hits[i];
            return 0.0;
        };
        (void)pool.parallel_reduce(hits.size(), mark);
        check(std::all_of(hits.begin(), hits.end(), [](int h) { return h == 1; }),
              "generic parallel_reduce: every index visited exactly once");
    }

    std::cout << "threaded annealing is deterministic\n";
    {
        Parameters p;
        p.N = 48; p.density = 0.35; p.T_initial = 50.0; p.T_final = 1.0; p.cooling_rate = 0.8;
        p.steps_per_temp = 500; p.parallel_min_modes = 1;
        p.annealing_file = tmp_file("annealing_thr.dat").string();
        HyperuniformConfig seq, par;
        for (auto* c : {&seq, &par}) {
            initialize_grid(*c, p.N, p.density, 9);
            build_mask(*c, 200.0);
        }
        p.num_threads = 1;
        const AnnealReport r1 = simulated_annealing(seq, p);
        p.num_threads = 4;
        const AnnealReport r4 = simulated_annealing(par, p);
        // The threads sum dE in a different order, so in general the results agree only up to rounding
        // (a borderline Metropolis decision could flip); for this case they are identical.
        check(seq.grid == par.grid && r1.accepted == r4.accepted,
              "4 threads give the same configuration as 1 thread (same seed)");
        std::filesystem::remove(p.annealing_file);
    }

    std::cout << "annealing\n";
    {
        Parameters p;
        p.N = 32; p.density = 0.3; p.chi = 0.2; p.T_initial = 0.0; p.T_final = 1e-3; p.cooling_rate = 0.97;
        p.steps_per_temp = 300; p.quench_steps = 20000; p.num_threads = 1;
        p.annealing_file = tmp_file("annealing.dat").string();
        HyperuniformConfig c;
        initialize_grid(c, p.N, p.density, 3);
        build_mask(c, choose_mask_for_chi(p.N, chi_points(p.N, c.M), p.chi).K2);
        const AnnealReport rep = simulated_annealing(c, p);
        const double res = mask_residual(c, rep.E_final);
        check(close(mask_residual(c, rep.E_initial), 1.0, 0.25), std::format(
            "random start: residual S ratio ~ 1 ({:.3f})", mask_residual(c, rep.E_initial)));
        check(res < 0.05, std::format("N = 32, chi = {:.3f}: residual S ratio {:.4f} < 5 %",
                                      chi_of(c.modes.size(), chi_points(p.N, c.M)), res));
        check(close(rep.E_final, energy_direct(c), 1e-8) || rep.E_final < 1e-9, "reported E_final is exact");
        check(lists_consistent(c) && c.M == occupied_count(32, 0.3), "particle number conserved");

        // S(k) sum rule: sum over k != 0 of |rho_hat|^2 = N^2 M - M^2 (Parseval), i.e. the mode-weighted mean
        // of S = (N^2 - M) / (N^2 - 1)
        const auto sk = tmp_file("sk.dat");
        analyze_structure_factor(c, sk.string());
        std::ifstream in(sk);
        std::string line;
        double wsum = 0.0, n = 0.0;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            double k, s, cnt;
            ss >> k >> s >> cnt;
            wsum += s * cnt;
            n += cnt;
        }
        const double N2 = 32.0 * 32.0;
        check(n == N2 - 1, "S(k) covers every k != 0 once");
        check(close(wsum / n, (N2 - c.M) / (N2 - 1), 1e-10), "S(k) satisfies the Parseval sum rule");
        std::filesystem::remove(sk);
        std::filesystem::remove(p.annealing_file);
    }

    std::cout << "rectangular lattice\n";
    {
        // Mask = circle in physical k: (kx/Nx)² + (ky/Ny)² <= K2 / (Nx Ny), checked in plain doubles against
        // every mode of the Nyquist box, independently of the integer key used by build_mask.
        for (const auto [Nx, Ny] : {std::pair{64, 24}, std::pair{45, 30}, std::pair{24, 96}}) {
            HyperuniformConfig c;
            initialize_grid(c, Nx, Ny, 0.3, 4);
            const double K2 = 0.3 * max_mask_K2(Nx, Ny) + 0.37;
            build_mask(c, K2);
            const double r2 = K2 / (static_cast<double>(Nx) * Ny);
            std::set<std::pair<int, int>> in;
            for (const auto& m : c.modes) in.insert({m.kx, m.ky});
            bool circle = true;
            int max_kx = 0, max_ky = 0;
            for (int kx = 0; kx <= Nx / 2 - 1; ++kx)
                for (int ky = -(Ny / 2 - 1); ky <= Ny / 2 - 1; ++ky) {
                    if (kx == 0 && ky <= 0) continue;
                    const double k2 = (static_cast<double>(kx) / Nx) * (static_cast<double>(kx) / Nx)
                                    + (static_cast<double>(ky) / Ny) * (static_cast<double>(ky) / Ny);
                    const bool inside = k2 <= r2 * (1.0 + 1e-12);
                    const bool clearly = k2 <= r2 * (1.0 - 1e-12) || k2 > r2 * (1.0 + 1e-12);
                    if (clearly && inside != (in.count({kx, ky}) == 1)) circle = false;
                    if (in.count({kx, ky})) { max_kx = std::max(max_kx, kx); max_ky = std::max(max_ky, std::abs(ky)); }
                }
            check(circle, std::format("{} x {}: mask = circle in physical k ({} modes)", Nx, Ny, c.modes.size()));
            bool sym = true;
            for (const auto& m : c.modes)
                if (m.kx > 0 && !in.count({m.kx, -m.ky})) sym = false;
            check(sym, std::format("{} x {}: mask symmetric under ky -> -ky", Nx, Ny));
            // An index-space circle would give max |kx| = max |ky|; the physical circle gives extents
            // proportional to Nx and Ny (equal within one index step in physical units).
            const double ext_x = static_cast<double>(max_kx) / Nx, ext_y = static_cast<double>(max_ky) / Ny;
            check(std::abs(ext_x - ext_y) <= 1.0 / std::min(Nx, Ny) && max_kx != max_ky,
                  std::format("{} x {}: physical extents |kx|/Nx = {:.4f}, |ky|/Ny = {:.4f} (index extents {} and {})",
                              Nx, Ny, ext_x, ext_y, max_kx, max_ky));

            calculate_FT(c);
            const double E = energy(c), Ed = energy_direct(c);
            check(close(E, Ed, 1e-10), std::format("{} x {}: tables = direct trigonometry ({:.10g} vs {:.10g})", Nx, Ny, E, Ed));
            double Ei = E;
            int acc = 0;
            for (int st = 0; st < 5000; ++st)
                if (auto dE = attempt_swap(c, 1e9, nullptr, 0.5, 1)) { Ei += *dE; ++acc; }
            check(lists_consistent(c) && close(Ei, energy_direct(c), 1e-8),
                  std::format("{} x {}: accumulated E = direct E after {} swaps", Nx, Ny, acc));

            const auto sk = tmp_file("sk_rect.dat");
            analyze_structure_factor(c, sk.string());
            std::ifstream ins(sk);
            std::string line;
            double wsum = 0.0, n = 0.0;
            while (std::getline(ins, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream ss(line);
                double k, sv, cnt;
                ss >> k >> sv >> cnt;
                wsum += sv * cnt;
                n += cnt;
            }
            const double sites = static_cast<double>(Nx) * Ny;
            check(n == sites - 1 && close(wsum / n, (sites - c.M) / (sites - 1), 1e-10),
                  std::format("{} x {}: S(k) covers every k != 0 once and satisfies Parseval", Nx, Ny));
            std::filesystem::remove(sk);
        }
        {
            const MaskChoice mc = choose_mask_for_chi(64, 24, chi_points(64, 24, occupied_count(64, 24, 0.3)), 0.4);
            HyperuniformConfig c;
            initialize_grid(c, 64, 24, 0.3, 1);
            build_mask(c, mc.K2);
            check(c.modes.size() == mc.n_modes && c.q_max == mc.q_max && c.q_next > c.q_max,
                  "rectangular choose_mask_for_chi: reported modes = built mask; next shell outside");
            check(close(chi_point_pattern(c.modes.size(), static_cast<std::size_t>(c.M)),
                        static_cast<double>(c.modes.size()) / (2.0 * (c.M - 1)), 1e-15),
                  "chi point_pattern = M(K) / [2 (N - 1)]");
        }
        {
            bool threw = false;
            HyperuniformConfig c;
            initialize_grid(c, 64, 24, 0.3, 1);
            try { build_mask(c, max_mask_K2(64, 24) * 1.01); } catch (const std::invalid_argument&) { threw = true; }
            check(threw, "rectangular: circle beyond the Nyquist limit of the short axis rejected");
            HyperuniformConfig big;
            initialize_grid(big, 4099, 4097, 0.3, 1);
            threw = false;
            try { build_mask(big, 4.0); } catch (const std::invalid_argument&) { threw = true; }
            check(threw, "lcm(Nx, Ny) > 2^24 rejected (4099 x 4097)");
        }
        {
            auto a = HyperuniformConfig{};
            initialize_grid(a, 96, 40, 0.35, 5);
            build_mask(a, 300.0);
            calculate_FT(a);
            auto b = a;
            ThreadPool pool(4);
            const double d_seq = delta_energy_range(a, 3, 5, 70, 37, 0, a.modes.size());
            const double d_par = pool.parallel_reduce_delta_E(b, 3, 5, 70, 37);
            check(close(d_seq, d_par, 1e-12), "rectangular: parallel dE = sequential dE");
        }
    }

    std::cout << "hard-core exclusion\n";
    {
        HyperuniformConfig c;
        initialize_grid(c, 80, 60, 0.05, 3, 2.5);
        check(lists_consistent(c) && c.M == occupied_count(80, 60, 0.05), "exclusion: random start has M sites, lists consistent");
        check(min_pair_distance(c, 10.0) >= 2.5, std::format("exclusion: random start respects min_distance 2.5 (min {:.3f})",
                                                             min_pair_distance(c, 10.0)));
        build_mask(c, 40.0);
        calculate_FT(c);
        double E = energy(c);
        for (int st = 0; st < 20000; ++st)
            if (auto dE = attempt_swap(c, 1e9, nullptr, 0.5, 2)) E += *dE;
        check(c.hc_rejections > 0 && min_pair_distance(c, 10.0) >= 2.5,
              std::format("exclusion kept by the annealing ({} proposals rejected)", c.hc_rejections));
        check(lists_consistent(c) && close(E, energy_direct(c), 1e-8), "exclusion: accumulated E = direct E");
        HyperuniformConfig a, b;
        initialize_grid(a, 40, 0.3, 9);
        initialize_grid(b, 40, 40, 0.3, 9, 0.0);
        check(a.occupied == b.occupied && b.excl.empty(), "min_distance = 0: same start as without exclusion");
        bool threw = false;
        try { HyperuniformConfig d; initialize_grid(d, 40, 40, 0.3, 1, 3.0); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "exclusion with no room for M sites rejected");
    }

    std::cout << "obstacle file\n";
    {
        check(sha256_hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "SHA-256 of the empty string");
        check(sha256_hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "SHA-256 of 'abc'");
        check(sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
              == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1", "SHA-256 of the 448-bit NIST message");
        check(sha256_hex(std::string(1000000, 'a')) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
              "SHA-256 of one million 'a'");
        // Padding boundaries, expected digests from Python's hashlib
        const std::size_t lens[] = {55, 56, 63, 64, 65, 119, 120};
        const char* digests[] = {"d5e285683cd4efc02d021a5c62014694958901005d6f71e89e0989fac77e4072", "04c26261370ee7541549d16dee320c723e3fd14671e66a099afe0a377c16888e", "75220b47218278e656f2013bb8f0c455a25eaf01e86c64924e9d48d89776d6f2", "7ce100971f64e7001e8fe5a51973ecdfe1ced42befe7ee8d5fd6219506b5393c", "9537c5fdf120482f7d58d25e9ed583f52c02b4e304ea814db1633ad565aed7e9", "000b48d4edf0fa7bee3c6236ecd2785baa5db4eeb8bb54341b029e0d9fa5fb0c", "13f05a0b594787f5ecd315edc96141bd3243203d1b7d4f0836f37308b276ba98"};
        bool all = true;
        for (std::size_t i = 0; i < 7; ++i) all = all && sha256_hex(std::string(lens[i], 'x')) == digests[i];
        check(all, "SHA-256 padding boundaries (55, 56, 63, 64, 65, 119, 120 bytes) = hashlib");

        bool rt = true;
        std::mt19937 g(5);
        std::uniform_real_distribution<double> u(0.0, 1000.0);
        for (int i = 0; i < 10000; ++i) { const double v = u(g); rt = rt && parse_real(format_real(v)) == v; }
        check(rt && format_real(-0.0) == "0" && format_real(2.0) == "2", "format_real: exact round trip, -0 -> 0, 2 -> \"2\"");

        // d_min against brute force, periodic in x only and with walls
        std::vector<ObsPoint> pts;
        std::uniform_real_distribution<double> ux(0.0, 37.0), uy(0.0, 11.0);
        for (int i = 0; i < 600; ++i) pts.push_back({ux(g), uy(g)});
        pts.push_back({0.0005, 5.0});
        pts.push_back({36.9995, 5.0005});  // the closest pair, but only through the periodic x boundary
        double per_best = 0.0, wall_best = 0.0;
        for (const bool per : {true, false}) {
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < pts.size(); ++i)
                for (std::size_t j = i + 1; j < pts.size(); ++j) {
                    double dx = std::abs(pts[i].x - pts[j].x);
                    if (per) dx = std::min(dx, 37.0 - dx);
                    best = std::min(best, std::hypot(dx, pts[i].y - pts[j].y));
                }
            check(close(min_center_distance(pts, 37.0, 11.0, per), best, 1e-14),
                  std::format("d_min by cell list = brute force ({}, {:.6g})", per ? "periodic_x" : "walls_x", best));
            (per ? per_best : wall_best) = best;
        }
        check(per_best < 0.002 && wall_best > per_best, "d_min: the periodic_x pair is found only with periodic_x");
        std::vector<ObsPoint> sparse{{1.0, 1.0}, {90.0, 40.0}};
        check(close(min_center_distance(sparse, 100.0, 50.0, true), std::hypot(11.0, 39.0), 1e-15),
              "d_min for a very sparse pattern (cells grow until exact)");

        ObsMeta meta;
        meta.generator_version = "test";
        meta.mode = "lattice";
        meta.Lx = 4.5; meta.Ly = 3.0; meta.R = 0.1; meta.seed = 3; meta.created = "2026-01-01";
        meta.K = 2.0; meta.n_constrained_modes = 4; meta.chi_convention = "binary_lattice"; meta.chi = 0.5;
        meta.lattice_constant = 1.5;
        const ObsResult r = build_obs(meta, {{3.0, 1.5}, {0.0, 1.5}, {1.5, 0.0}});
        const std::string body = "1.5 0\n0 1.5\n3 1.5\n";
        check(r.text.ends_with("# columns: x y\n" + body) && r.sha256_body == sha256_hex(body),
              "build_obs: body ordered by y then x, checksum of the body bytes");
        check(r.n == 3 && r.rho_O == 3.0 / (4.5 * 3.0) && close(r.d_min, 1.5, 1e-15)
              && r.text.find("# rho_O: " + format_real(3.0 / 13.5) + "\n") != std::string::npos,
              "build_obs: n_obstacles, rho_O and d_min (x minimum image) from the written values");
        bool threw = false;
        try { (void)build_obs(meta, {{4.5, 0.0}}); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "build_obs: x = Lx rejected (domain is [0, Lx))");
        {
            const auto path = tmp_file("rt.obs");
            (void)write_obs(path.string(), meta, {{3.0, 1.5}, {0.0, 1.5}, {1.5, 0.0}});
            const ObsData d = read_obs(path.string());
            check(d.points.size() == 3 && d.points[0].x == 1.5 && d.points[2].x == 3.0 && d.number("Lx") == 4.5
                  && *d.get("chi_convention") == "binary_lattice",
                  "read_obs: header and body read back (ordered body)");
            std::string text = r.text;
            text.replace(text.rfind("1.5 0"), 5, "1.6 0");
            { std::ofstream o(path, std::ios::binary); o << text; }
            bool bad_sha = false;
            try { (void)read_obs(path.string()); } catch (const std::runtime_error&) { bad_sha = true; }
            check(bad_sha, "read_obs aborts on a sha256_body mismatch");
            std::filesystem::remove(path);
        }
        threw = false;
        ObsMeta bad = meta;
        bad.n_constrained_modes.reset();
        try { (void)build_obs(bad, {{0.0, 0.0}}); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "build_obs: class stealthy without n_constrained_modes rejected");
    }

    std::cout << "L-BFGS\n";
    {
        // Ill-conditioned quadratic, condition number 1e4, minimum at x_i = i
        const std::size_t n = 50;
        const Objective quad = [n](const std::vector<double>& x, std::vector<double>& g) {
            g.assign(n, 0.0);
            double f = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double a = std::pow(1e4, static_cast<double>(i) / (n - 1)), d = x[i] - static_cast<double>(i);
                f += 0.5 * a * d * d;
                g[i] = a * d;
            }
            return f;
        };
        std::vector<double> x(n, 0.0);
        LbfgsOptions opt;
        opt.grad_tol = 1e-10;
        LbfgsResult r = lbfgs_minimize(quad, x, opt);
        double err = 0.0;
        for (std::size_t i = 0; i < n; ++i) err = std::max(err, std::abs(x[i] - static_cast<double>(i)));
        check(r.reason == "grad_tol" && err < 1e-9, std::format("L-BFGS: quadratic, condition 1e4 ({} iterations, {} evaluations, error {:.1e})",
                                                                r.iterations, r.evaluations, err));
        // Rosenbrock in 2 and 100 dimensions, minimum at (1, ..., 1)
        for (const std::size_t dim : {2u, 100u}) {
            const Objective rosen = [dim](const std::vector<double>& y, std::vector<double>& g) {
                g.assign(dim, 0.0);
                double f = 0.0;
                for (std::size_t i = 0; i + 1 < dim; ++i) {
                    const double a = y[i + 1] - y[i] * y[i], b = 1.0 - y[i];
                    f += 100.0 * a * a + b * b;
                    g[i] += -400.0 * y[i] * a - 2.0 * b;
                    g[i + 1] += 200.0 * a;
                }
                return f;
            };
            std::vector<double> y(dim, -1.2);
            for (std::size_t i = 1; i < dim; i += 2) y[i] = 1.0;
            LbfgsOptions o;
            o.grad_tol = 1e-9;
            o.max_iterations = 10000;
            const LbfgsResult rr = lbfgs_minimize(rosen, y, o);
            double e = 0.0;
            for (const double v : y) e = std::max(e, std::abs(v - 1.0));
            check(rr.reason == "grad_tol" && e < 1e-7, std::format("L-BFGS: Rosenbrock, {} dimensions ({} iterations, error {:.1e})",
                                                                   dim, rr.iterations, e));
        }
        std::vector<double> z(n, 0.0);
        long long calls = 0;
        const LbfgsResult rc = lbfgs_minimize(quad, z, opt, [&](long long it, double, double) { ++calls; return it >= 5; });
        check(rc.reason == "callback" && rc.iterations == 5 && calls == 5, "L-BFGS: the progress callback stops the run");
        LbfgsOptions ft;
        ft.f_target = 1e-3;
        std::vector<double> w(n, 0.0);
        const LbfgsResult rf = lbfgs_minimize(quad, w, ft);
        check(rf.reason == "f_target" && rf.f <= 1e-3, "L-BFGS: stops at f_target");
    }

    std::cout << "continuum\n";
    {
        // Axis phases by recurrence against a long double reference (position folded in long double), up to
        // n = 600, also for positions outside the box. Any double evaluation of the angle 2π n x / L has an
        // error of order n 2π eps; the recurrence must stay within a few times that and be no worse than
        // std::polar in double.
        std::mt19937 g(17);
        std::uniform_real_distribution<double> u(-50.0, 150.0);
        double worst = 0.0, worst_polar = 0.0;
        std::vector<std::complex<double>> ph(601);
        const long double pi_l = 3.141592653589793238462643383279502884L;
        for (int t = 0; t < 500; ++t) {
            const double x = u(g), L = 97.3;
            axis_phases(x, L, 600, ph.data());
            const long double xf = static_cast<long double>(x) - static_cast<long double>(L) * std::floor(static_cast<long double>(x) / L);
            for (int n = 0; n <= 600; ++n) {
                const long double a = -2.0L * pi_l * n * xf / L;
                const std::complex<long double> ref(std::cos(a), std::sin(a));
                const auto& v = ph[static_cast<std::size_t>(n)];
                worst = std::max(worst, static_cast<double>(std::abs(std::complex<long double>(v.real(), v.imag()) - ref)));
                const auto pd = std::polar(1.0, -2.0 * PI * n * x / L);
                worst_polar = std::max(worst_polar, static_cast<double>(std::abs(std::complex<long double>(pd.real(), pd.imag()) - ref)));
            }
        }
        const double bound = 4.0 * 600 * 2.0 * PI * std::numeric_limits<double>::epsilon();
        check(worst < bound && worst <= worst_polar,
              std::format("axis phases by recurrence: max error {:.2e} up to n = 600 (bound {:.2e}; std::polar {:.2e})",
                          worst, bound, worst_polar));

        const double Lx = 60.0, Ly = 17.0;  // rectangular box, aspect ratio about 3.5
        ContinuumConfig c;
        continuum_init_random(c, Lx, Ly, 300, 5);
        const double K = continuum_K_of_kappa(Lx, Ly, 300, 0.5);
        build_continuum_mask(c, K, 0.0);
        bool circle = true;
        for (std::size_t m = 0; m < c.modes.size(); ++m) {
            const double kx = 2.0 * PI * c.modes[m].kx / Lx, ky = 2.0 * PI * c.modes[m].ky / Ly;
            circle = circle && kx * kx + ky * ky <= K * K;
        }
        check(circle && c.modes.size() == continuum_mode_count(Lx, Ly, K) && c.nx_max != c.ny_max,
              std::format("continuum mask = circle in physical k ({} modes, index extents {} and {})",
                          c.modes.size(), c.nx_max, c.ny_max));
        check(close(0.5 / std::sqrt(PI / 4.0 * 0.5 * 0.5), 1.0 / std::sqrt(PI / 4.0), 1e-15) &&
              std::abs(chi_point_pattern(c.modes.size(), 300) - PI * 0.25 / 4.0) < 0.03,
              std::format("kappa = 0.5 gives chi_point ~ pi kappa^2 / 4 ({:.4f} vs {:.4f})",
                          chi_point_pattern(c.modes.size(), 300), PI * 0.25 / 4.0));

        continuum_FT(c);
        const double E = continuum_energy(c), Ed = continuum_energy_direct(c);
        check(close(E, Ed, 1e-10), std::format("continuum E: recurrence tables = direct trigonometry ({:.10g} vs {:.10g})", E, Ed));
        check(std::abs(continuum_residual(c, E) - 1.0) < 0.2,
              std::format("Poisson start: residual <S>_mask ~ 1 ({:.3f})", continuum_residual(c, E)));

        // Gradient against central finite differences (weighted spectral energy)
        build_continuum_mask(c, K, 2.0);
        std::vector<double> grad;
        (void)continuum_spectral_energy_grad(c, grad);
        double worst_rel = 0.0;
        for (const std::size_t q : {0u, 1u, 77u, 300u, 599u}) {
            const double h = 1e-6, x0 = c.r[q];
            c.r[q] = x0 + h; continuum_FT(c); const double Ep = continuum_weighted_energy(c);
            c.r[q] = x0 - h; continuum_FT(c); const double Em = continuum_weighted_energy(c);
            c.r[q] = x0;
            const double fd = (Ep - Em) / (2.0 * h);
            worst_rel = std::max(worst_rel, std::abs(fd - grad[q]) / std::max(1.0, std::abs(fd)));
        }
        check(worst_rel < 1e-5, std::format("spectral gradient = finite differences (worst relative error {:.2e})", worst_rel));
        double sx = 0.0, sy = 0.0, scale = 0.0;
        for (std::size_t p = 0; p < c.N; ++p) { sx += grad[2 * p]; sy += grad[2 * p + 1]; scale += std::abs(grad[2 * p]); }
        check(std::abs(sx) < 1e-9 * scale && std::abs(sy) < 1e-9 * scale, "spectral gradient sums to zero (translation invariance)");
        double wmean = 0.0;
        for (double w : c.weight) wmean += w;
        check(close(wmean / static_cast<double>(c.weight.size()), 1.0, 1e-12), "continuum weights have mean 1");

        // Shell choice for a point-pattern chi; K midway, so recounting gives the same modes
        const ContinuumMaskChoice mc = choose_continuum_mask_for_chi(Lx, Ly, 300, 0.3);
        check(mc.n_modes == continuum_mode_count(Lx, Ly, mc.K) && std::abs(mc.chi - 0.3) < 0.02,
              std::format("choose_continuum_mask_for_chi: chi {:.4f} for target 0.3, {} modes, recount at K agrees",
                          mc.chi, mc.n_modes));

        // Hard core: energy and gradient
        ContinuumConfig h;
        continuum_init_random(h, 20.0, 7.0, 120, 8);
        h.sigma = 1.2;
        h.core_weight = 3.0;
        const double Ec = continuum_core_energy(h, nullptr), Ecd = continuum_core_energy_direct(h);
        check(Ec > 0.0 && close(Ec, Ecd, 1e-12), std::format("core energy: cell list = brute force ({:.6g})", Ec));
        std::vector<double> gc(2 * h.N, 0.0);
        (void)continuum_core_energy(h, &gc);
        double worst_core = 0.0, csum = 0.0;
        for (std::size_t q = 0; q < 2 * h.N; ++q) {
            const double hh = 1e-7, x0 = h.r[q];
            h.r[q] = x0 + hh; const double Ep = continuum_core_energy(h, nullptr);
            h.r[q] = x0 - hh; const double Em = continuum_core_energy(h, nullptr);
            h.r[q] = x0;
            worst_core = std::max(worst_core, std::abs((Ep - Em) / (2.0 * hh) - gc[q]));
            if (q % 2 == 0) csum += gc[q];
        }
        check(worst_core < 1e-5 && std::abs(csum) < 1e-10, std::format("core gradient = finite differences (worst {:.2e}), sums to 0", worst_core));
        ContinuumConfig rsa;
        continuum_init_random(rsa, 30.0, 10.0, 150, 3, 1.0);
        check(continuum_min_distance(rsa) >= 1.0 && continuum_core_energy(rsa, nullptr) == 0.0,
              std::format("continuum random sequential addition respects sigma (min distance {:.4f})", continuum_min_distance(rsa)));
        bool threw = false;
        try { ContinuumConfig full; continuum_init_random(full, 5.0, 5.0, 100, 1, 1.0); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "continuum: impossible random sequential addition rejected");
    }

    std::cout << "continuum kernels and threads\n";
    {
        ContinuumConfig c;
        continuum_init_random(c, 70.0, 20.0, 1400, 9);
        build_continuum_mask(c, continuum_K_of_kappa(70.0, 20.0, 1400, 0.65), 2.0);
        bool runs_ok = !c.runs.empty();
        std::size_t covered = 0;
        for (const auto& r : c.runs) {
            for (std::size_t j = 0; j < r.count; ++j)
                runs_ok = runs_ok && c.modes[r.begin + j].kx == r.nx && c.modes[r.begin + j].ky == r.ny_lo + static_cast<int>(j);
            runs_ok = runs_ok && r.begin == covered;
            covered += r.count;
        }
        check(runs_ok && covered == c.modes.size(), std::format("mode runs cover the mask in order ({} runs)", c.runs.size()));
        std::vector<double> g1, g4;
        const double E1 = continuum_spectral_energy_grad(c, g1, nullptr);
        const auto ft1 = c.ft_re;
        ThreadPool pool(4);
        const double E4 = continuum_spectral_energy_grad(c, g4, &pool);
        check(E1 == E4 && ft1 == c.ft_re && g1 == g4, "rho_hat and gradient with 4 threads = 1 thread, bit for bit");
        check(close(continuum_energy(c), continuum_energy_direct(c), 1e-10), "grouped rho_hat = direct trigonometry");
        // Gradient by blocks against a plain per-mode sum with std::polar
        double worst = 0.0;
        for (const std::size_t p : {0u, 33u, 700u, 1399u}) {
            double gx = 0.0, gy = 0.0;
            for (std::size_t m = 0; m < c.modes.size(); ++m) {
                const auto e = std::polar(1.0, -(c.kx[m] * c.r[2 * p] + c.ky[m] * c.r[2 * p + 1]));
                const double s = c.weight[m] * (c.ft_re[m] * e.imag() - c.ft_im[m] * e.real());
                gx += c.kx[m] * s;
                gy += c.ky[m] * s;
            }
            worst = std::max({worst, std::abs(4.0 * gx - g1[2 * p]) / (1.0 + std::abs(4.0 * gx)),
                              std::abs(4.0 * gy - g1[2 * p + 1]) / (1.0 + std::abs(4.0 * gy))});
        }
        check(worst < 1e-10, std::format("blocked gradient = per-mode sum with std::polar (worst relative {:.1e})", worst));
    }

    std::cout << "continuum Monte Carlo\n";
    {
        ContinuumConfig c;
        continuum_init_random(c, 30.0, 12.0, 360, 4, 0.25);
        c.core_weight = 1.0;
        build_continuum_mask(c, continuum_K_of_kappa(30.0, 12.0, 360, 0.6), 2.0);
        Parameters p;
        p.mode = "continuum";
        p.T_initial = 0.0; p.T_final = 1e-5; p.cooling_rate = 0.8; p.steps_per_temp = 2000; p.quench_steps = 2000;
        p.resync_interval = 1000000000;  // no resync before the end: the drift measures the incremental updates
        p.residual_tol = 0.0;
        const ContinuumMcReport r = continuum_anneal(c, p, nullptr, nullptr);
        check(r.accepted > 1000 && r.max_drift < 1e-10,
              std::format("MC: incremental objective = exact after {} accepted moves (drift {:.1e})", r.accepted, r.max_drift));
        check(r.f_final < 0.05 * r.f_initial, std::format("MC lowers the objective ({:.3g} -> {:.3g})", r.f_initial, r.f_final));
        check(continuum_min_distance(c) >= 0.25 && continuum_core_energy(c, nullptr) == 0.0,
              std::format("MC keeps the hard core as a constraint (min distance {:.4f})", continuum_min_distance(c)));
        ContinuumConfig fresh = c;
        continuum_FT(fresh);
        bool same = fresh.ft_re == c.ft_re && fresh.ft_im == c.ft_im;
        check(same, "MC leaves rho_hat consistent with the positions");
        ThreadPool pool(4);
        ContinuumConfig a, b;
        for (auto* q : {&a, &b}) {
            continuum_init_random(*q, 30.0, 12.0, 360, 4, 0.25);
            build_continuum_mask(*q, continuum_K_of_kappa(30.0, 12.0, 360, 0.6), 2.0);
        }
        p.T_final = 1e-3;
        (void)continuum_anneal(a, p, nullptr, nullptr);
        (void)continuum_anneal(b, p, &pool, nullptr);
        check(a.r == b.r, "MC with 4 threads = MC with 1 thread (same seed)");
    }

    std::cout << "parameters\n";
    {
        const auto f = tmp_file("params.in");
        { std::ofstream o(f); o << "N: 64abc\ndensity: 0.3   # trailing comment\nchi: 0.25\n"; }
        Parameters p;
        check(!p.read(f.string()), "malformed number reported (N: 64abc)");
        check(p.density == 0.3 && p.chi == 0.25, "trailing comments stripped");
        Parameters q;
        q.density = 1.2; q.cooling_rate = 1.0;
        check(q.validate().size() == 2, "validate: density and cooling_rate out of range");
        Parameters r;
        r.mode = "continuum";
        check(r.validate().size() == 2, "validate: mode = continuum without box or mask rejected");
        r.Lx = 40.0; r.Ly = 10.0; r.rho_O = 1.0; r.kappa = 0.6;
        check(r.validate().empty() && r.continuum_points() == 400, "validate: mode = continuum with Lx, Ly, rho_O and kappa accepted");
        r.K_given = true;
        check(r.validate().size() == 1, "validate: K rejected in mode continuum");
        r.K_given = false;
        r.R = 0.1;
        check(close(r.continuum_core_distance(), 0.2 * (1.0 + 1e-6), 1e-15), "continuum hard core = 2R (1 + 1e-6)");
        r.mode = "lattice";
        check(r.validate().empty(), "validate: mode = lattice accepted (default parameters valid)");
        Parameters u;
        u.Nx = 64; u.Ny = 24;
        check(u.validate().size() == 1, "validate: rectangular lattice with only K rejected");
        u.kappa = 0.6;
        check(u.validate().empty() && u.nx() == 64 && u.ny() == 24, "validate: rectangular lattice with kappa accepted");
        u.Nx = 3;
        check(!u.validate().empty(), "validate: Nx = 3 rejected");
        Parameters w;
        w.Lx = 160.0; w.Ly = 40.0; w.rho_O = 0.2; w.density = 0.02;
        (void)w.derive_lattice();
        check(w.nx() == 506 && w.ny() == 127 && w.nx() * w.lattice_constant == 160.0 && close(w.density, 0.2 * w.lattice_constant * w.lattice_constant, 1e-15),
              std::format("derive_lattice: Lx kept exact, {} x {} sites, a = {:.6f}", w.nx(), w.ny(), w.lattice_constant));
        w.R = 0.2;
        check(close(w.effective_min_distance(), 0.4 / w.lattice_constant, 1e-8) && w.effective_min_distance() > 0.4 / w.lattice_constant,
              "R implies min_distance slightly above 2R / lattice_constant");
        w.obs_file = "auto"; w.kappa = 0.6;
        check(w.validate().empty(), "validate: physical box with R and obs_file accepted");
        w.R = 0.0;
        check(w.validate().size() == 1, "validate: obs_file without R rejected");
        std::filesystem::remove(f);
    }

    std::cout << (failures ? std::format("{} FAILED\n", failures) : std::string("all tests passed\n"));
    return failures ? 1 : 0;
}
