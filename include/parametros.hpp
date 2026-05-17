#pragma once

#include <string>

/// Simulation parameters, read from a text file.
struct Parameters {
    int         N                 = 64;
    double      density           = 0.4;
    double      K                 = 4.0;
    double      T_initial         = 10.0;
    double      T_final           = 0.01;
    double      cooling_rate      = 0.999;
    int         steps_per_temp    = 1000;
    int         seed              = 42;
    int         num_threads       = 0;  ///< 0 means auto-detect, 1 means sequential
    std::string output_file       = "hyperuniform.pgm";
    std::string annealing_file    = "annealing.dat";  ///< T/E/rate monitoring per step
    std::string sk_file           = "sk.dat";         ///< Structure factor S(k)

    /// Reads the parameter file. Returns false if it couldn't be opened.
    [[nodiscard]] bool read(const std::string &filename);

    /// Prints the active parameters to stdout.
    void print() const;
};
