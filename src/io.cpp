#include "io.hpp"

#include <format>
#include <fstream>
#include <iostream>

void save_PGM(const HyperuniformConfig& config, const std::string& filename,
              const std::vector<std::string>& comments) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << std::format("Error: could not open '{}' for writing.\n", filename);
        return;
    }

    // Accumulate everything in a string before dumping it to disk (one write call)
    const auto N = static_cast<std::size_t>(config.N);
    std::string buf = "P2\n";
    for (const auto& c : comments) buf += "# " + c + '\n';
    buf += std::format("{} {}\n1\n", config.N, config.N);
    buf.reserve(buf.size() + 2 * N * N);
    for (std::size_t y = 0; y < N; ++y) {
        for (std::size_t x = 0; x < N; ++x) {
            buf += config.grid[x + N * y] ? '1' : '0';
            buf += (x + 1 < N) ? ' ' : '\n';
        }
    }
    fout << buf;
    // fout is closed automatically when going out of scope (RAII)

    std::cout << std::format("Configuration saved to '{}'\n", filename);
}
