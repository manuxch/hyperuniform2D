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
    const std::string arg = (argc > 1) ? argv[1] : "";
    if (arg == "--version" || arg == "-v") {
        std::cout << "hyper2D " HYPER2D_VERSION "\n";
        return 0;
    }
    if (arg == "--help" || arg == "-h") {
        std::cout << "Usage: hyper2D [parameter_file]   (default: params.in)\n"
                     "       hyper2D --version\n";
        return 0;
    }
    const std::string parameter_file = arg.empty() ? "params.in" : arg;
    std::cout << "hyper2D " HYPER2D_VERSION "\n";

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
        build_mask(config, K2, p.weight_power);
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
        if (p.weight_power != 0.0)
            std::cout << std::format("Weighted energy (w ∝ |k|^-{}): {:.6g} -> {:.6g}\n", p.weight_power,
                                     rep.Ew_initial, rep.Ew_final);
        const double residual = mask_residual(config, rep.E_final);
        std::cout << std::format("Final energy: {:.6g}{}, max incremental drift {:.3g}\n", rep.E_final,
                                 rep.reached_tol ? " (<= energy_tol)" : "", rep.max_drift);
        const double low = band_residual(config, config.K / 4.0);
        std::cout << std::format("Residual <S>_mask / (1 - rho) = {:.4g} (1 = random start, 0 = perfectly stealthy)\n",
                                 residual);
        std::cout << std::format("Residual of the longest wavelengths, |k| <= K/4 = {:.4g}: {:.4g}\n", config.K / 4.0, low);

        save_PGM(config, p.output_file, {
            "hyperuniform2D " HYPER2D_VERSION ": stealthy hyperuniform pattern, 1 = occupied",
            std::format("N {} density {:.6f} occupied {}", p.N, static_cast<double>(config.M) / (static_cast<double>(p.N) * p.N), config.M),
            std::format("K {:.6g} K2 {:.6g} modes_half_plane {} chi {:.6f} chi_target {}", config.K, config.K2,
                        config.modes.size(), chi, p.chi),
            std::format("seed {} E_initial {:.6g} E_final {:.6g} residual_S_ratio {:.4g} residual_low_k {:.4g} weight_power {}",
                        p.seed, rep.E_initial, rep.E_final, residual, low, p.weight_power),
        });
        analyze_structure_factor(config, p.sk_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
