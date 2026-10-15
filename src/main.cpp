#include "annealing.hpp"
#include "config.hpp"
#include "fourier.hpp"
#include "io.hpp"
#include "parametros.hpp"

#include <cmath>
#include <exception>
#include <format>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    const std::string parameter_file = (argc > 1) ? argv[1] : "params.in";

    Parameters p;
    if (!p.read(parameter_file)) {
        std::cerr << "Could not read the parameter file.\n";
        return 1;
    }
    if (const auto errors = p.validate(); !errors.empty()) {
        for (const auto& e : errors) std::cerr << "Error: " << e << '\n';
        return 1;
    }
    if (p.chi > 0.0 && p.K_given)
        std::cerr << "Warning: both chi and K given; K is derived from chi and the K value is ignored.\n";
    p.print();

    try {
        HyperuniformConfig config;
        initialize_grid(config, p.N, p.density, p.seed);
        const std::size_t n_chi = chi_points(p.N, config.M);

        double K2 = p.K * p.K;
        if (p.chi > 0.0) {
            const MaskChoice mc = choose_mask_for_chi(p.N, n_chi, p.chi);
            K2 = mc.K2;
        }
        build_mask(config, K2);
        const double chi = chi_of(config.modes.size(), n_chi);
        std::cout << std::format("Occupied sites: {} (density {:.6f}); chi uses n = {} points\n",
                                 config.M, static_cast<double>(config.M) / (static_cast<double>(p.N) * p.N), n_chi);
        std::cout << std::format("Mask: K = {:.6g} (K² = {:.6g}), chi = {:.6f}{}\n", config.K, config.K2, chi,
                                 p.chi > 0.0 ? std::format(" (target {})", p.chi) : "");
        if (chi > 0.5)
            std::cerr << std::format("Warning: chi = {:.4f} > 1/2: stealthy disordered ground states are not expected; "
                                     "the annealing may not reach E = 0 and the pattern tends to crystallize.\n", chi);

        const AnnealReport rep = simulated_annealing(config, p);
        std::cout << std::format("Initial energy: {:.6g}\n", rep.E_initial);
        const double residual = mask_residual(config, rep.E_final);
        std::cout << std::format("Final energy: {:.6g}{}, max incremental drift {:.3g}\n", rep.E_final,
                                 rep.reached_tol ? " (<= energy_tol)" : "", rep.max_drift);
        std::cout << std::format("Residual <S>_mask / (1 - rho) = {:.4g} (1 = random start, 0 = perfectly stealthy)\n",
                                 residual);

        save_PGM(config, p.output_file, {
            "hyperuniform2D: stealthy hyperuniform pattern, 1 = occupied",
            std::format("N {} density {:.6f} occupied {}", p.N, static_cast<double>(config.M) / (static_cast<double>(p.N) * p.N), config.M),
            std::format("K {:.6g} K2 {:.6g} modes_half_plane {} chi {:.6f} chi_target {}", config.K, config.K2,
                        config.modes.size(), chi, p.chi),
            std::format("seed {} E_initial {:.6g} E_final {:.6g} residual_S_ratio {:.4g}", p.seed, rep.E_initial,
                        rep.E_final, residual),
        });
        analyze_structure_factor(config, p.sk_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
