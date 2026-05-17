#include "annealing.hpp"
#include "config.hpp"
#include "fourier.hpp"
#include "io.hpp"
#include "parametros.hpp"

#include <format>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    const std::string parameter_file = (argc > 1) ? argv[1] : "params.in";

    Parameters p;
    if (!p.read(parameter_file)) {
        std::cerr << "Could not read parameters. Default values will be used.\n";
    }
    p.print();

    HyperuniformConfig config;
    initialize_grid(config, p.N, p.density, p.seed);
    build_mask(config, p.K);
    calculate_FT(config);

    std::cout << std::format("Initial energy: {:.4f}\n", energy(config));

    simulated_annealing(config, p);

    std::cout << std::format("Final energy: {:.4f}\n", energy(config));

    save_PGM(config, p.output_file);
    analyze_structure_factor(config, p.sk_file);

    return 0;
}
