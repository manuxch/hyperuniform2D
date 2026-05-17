#include "fourier.hpp"
#include "constants.hpp"

#include <complex>
#include <format>
#include <fstream>
#include <iostream>

void calculate_FT(HyperuniformConfig& config) {
    const int N = config.N;
    const auto num_modes = config.modes.size();
    std::fill(config.FT.begin(), config.FT.end(), Complex{0.0, 0.0});

    for (int x = 0; x < N; ++x) {
        for (int y = 0; y < N; ++y) {
            // flat_idx: no modulo (x,y are already in [0,N))
            if (config.grid[flat_idx(x, y, N)] == 0) continue;
            for (std::size_t m = 0; m < num_modes; ++m) {
                const auto [kx, ky] = config.modes[m]; // C++17 structured binding
                const double theta  = -2.0 * PI * (kx * x + ky * y) / N;
                config.FT[m] += std::polar(1.0, theta); // std::polar instead of Complex(cos,sin)
            }
        }
    }
}

double energy(const HyperuniformConfig& config) {
    double E = 0.0;
    for (const auto& val : config.FT) {
        E += std::norm(val);
    }
    return E;
}

void analyze_structure_factor(const HyperuniformConfig& config,
                              const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << std::format("Error: could not create '{}'\n", filename);
        return;
    }

    const int N      = config.N;
    const int limit  = N / 2;

    // File header
    fout << "# Structure factor S(k) for ky = 0\n";
    fout << std::format("{:<8} {:<20}\n", "# kx", "S(k)");

    for (int kx = 1; kx < limit; ++kx) {
        Complex sum{0.0, 0.0};
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                if (config.grid[flat_idx(x, y, N)] == 1) {
                    const double theta = -2.0 * PI * kx * x / N;
                    sum += std::polar(1.0, theta);
                }
            }
        }
        const double S = std::norm(sum) / static_cast<double>(N * N);
        fout << std::format("{:<8d} {:<20.10e}\n", kx, S);
    }
    std::cout << std::format("Structure factor saved to '{}'\n", filename);
}
