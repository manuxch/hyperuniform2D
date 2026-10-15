#pragma once

#include <string>
#include <vector>

/// Simulation parameters, read from a text file.
struct Parameters {
    int         N                 = 64;
    double      density           = 0.4;
    double      K                 = 4.0;    ///< Mask radius (units of 2π/N); ignored if chi > 0
    double      chi               = 0.0;    ///< Stealthiness; > 0 derives K from it, 0 uses K
    double      T_initial         = 10.0;   ///< <= 0: automatic (typical uphill swap accepted with prob. 1/2)
    double      T_final           = 0.01;
    double      cooling_rate      = 0.999;
    int         steps_per_temp    = 1000;   ///< Swap proposals (occupied <-> empty) per temperature level
    double      local_fraction    = 0.5;    ///< Fraction of proposals that move a point to a nearby empty site
    int         local_radius      = 1;      ///< "Nearby": Chebyshev distance <= local_radius
    long long   quench_steps      = 0;      ///< Zero-temperature proposals after the annealing (0 = none)
    double      energy_tol        = 1e-6;   ///< Stop as soon as E <= energy_tol (stealthy reached)
    long long   resync_interval   = 1000000; ///< Recompute rho_hat from scratch every so many accepted swaps
    unsigned    seed              = 1729;
    int         num_threads       = 0;      ///< 0: min(available CPUs, 4); 1: sequential
    int         parallel_min_modes = 2000;  ///< Fewer half-plane modes than this: sequential
    std::string output_file       = "hyperuniform.pgm";
    std::string annealing_file    = "annealing.dat";  ///< T/E/rate monitoring per temperature level
    std::string sk_file           = "sk.dat";         ///< Radially averaged structure factor S(k)

    bool K_given = false;  ///< K appeared in the parameter file

    /// Reads the parameter file. Returns false if it couldn't be opened or a value couldn't be parsed.
    [[nodiscard]] bool read(const std::string &filename);

    /// Checks ranges; returns one message per problem (empty if all valid).
    [[nodiscard]] std::vector<std::string> validate() const;

    /// Prints the active parameters to stdout.
    void print() const;
};
