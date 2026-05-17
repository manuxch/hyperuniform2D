#include "io.hpp"

#include <format>
#include <fstream>
#include <iostream>
#include <sstream>

void save_PGM(const HyperuniformConfig& config, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << std::format("Error: could not open '{}' for writing.\n", filename);
        return;
    }

    // Accumulate everything in an internal string buffer before dumping to disk
    // (reduces the number of calls to the OS)
    std::ostringstream buf;
    buf << "P2\n" << config.N << ' ' << config.N << "\n1\n";

    for (int y = 0; y < config.N; ++y) {
        for (int x = 0; x < config.N; ++x) {
            buf << config.grid[flat_idx(x, y, config.N)];
            if (x + 1 < config.N) buf << ' ';
        }
        buf << '\n';
    }
    fout << buf.str();
    // fout is closed automatically when going out of scope (RAII)

    std::cout << std::format("Configuration saved to '{}'\n", filename);
}
