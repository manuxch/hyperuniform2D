// Unit tests for hyperuniform2D. No framework: each check prints its result and the exit code
// counts failures (run with ctest, or ./build/test_hyper2D).
#include "annealing.hpp"
#include "config.hpp"
#include "fourier.hpp"
#include "parametros.hpp"
#include "thread_pool.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

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
        std::filesystem::remove(f);
    }

    std::cout << (failures ? std::format("{} FAILED\n", failures) : std::string("all tests passed\n"));
    return failures ? 1 : 0;
}
