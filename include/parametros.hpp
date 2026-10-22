#pragma once

#include <cstddef>
#include <string>
#include <vector>

/// Simulation parameters, read from a text file.
struct Parameters {
    std::string mode              = "lattice";  ///< "lattice" (sites of a periodic grid); "continuum" is planned
    int         N                 = 64;     ///< Side of a square lattice; Nx, Ny override it
    int         Nx                = 0;      ///< Lattice width in sites (0: N)
    int         Ny                = 0;      ///< Lattice height in sites (0: N)
    double      density           = 0.4;
    double      K                 = 4.0;    ///< Mask radius (units of 2π/N), square lattice only; ignored if kappa or chi > 0
    double      chi               = 0.0;    ///< Stealthiness (binary_lattice convention); > 0 derives K from it
    double      min_distance      = 0.0;    ///< Hard-core exclusion between occupied sites, lattice spacings (0 = none)
    double      kappa             = 0.0;    ///< kappa = K a / 2π, a = mean obstacle spacing; > 0 sets K (wins over chi and K)
    double      T_initial         = 10.0;   ///< <= 0: automatic (typical uphill swap accepted with prob. 1/2)
    double      T_final           = 0.01;
    double      cooling_rate      = 0.999;
    int         steps_per_temp    = 1000;   ///< Swap proposals (occupied <-> empty) per temperature level
    double      local_fraction    = 0.5;    ///< Fraction of proposals that move a point to a nearby empty site
    int         local_radius      = 1;      ///< "Nearby": Chebyshev distance <= local_radius
    long long   quench_steps      = 0;      ///< Zero-temperature proposals after the annealing (0 = none)
    double      weight_power      = 0.0;    ///< Annealing energy Σ w(k)|rho_hat|², w ∝ |k|^(-weight_power); 0 = standard
    double      energy_tol        = 1e-6;   ///< Stop as soon as E <= energy_tol (stealthy reached)
    long long   resync_interval   = 1000000; ///< Recompute rho_hat from scratch every so many accepted swaps
    unsigned    seed              = 1729;
    int         num_threads       = 0;      ///< 0: min(available CPUs, 4); 1: sequential
    int         parallel_min_modes = 2000;  ///< Fewer half-plane modes than this: sequential
    // Obstacle file (format of hyperGalton2D, version 1); lengths in grain diameters d
    std::string obs_file          = "";     ///< Output .obs; "auto" = <class>_rho.._Lx.._Ly.._arrNNN.obs; empty = none
    double      lattice_constant  = 1.0;    ///< Lattice spacing in d
    double      R                 = 0.0;    ///< Obstacle radius in d (> 0 for obs_file); implies min_distance >= 2R / lattice_constant
    std::string boundary          = "periodic_x";  ///< Declared in the .obs: periodic_x or walls_x
    double      Lx = 0.0, Ly = 0.0, rho_O = 0.0;   ///< If all > 0: Nx, Ny, lattice_constant derived from them and density
    int         arr_index         = 0;      ///< Realization index of the first file
    int         realizations      = 1;      ///< Independent realizations, seeds seed, seed + 1, ...
    // Continuum mode (mode: continuum): Lx, Ly, rho_O give the box and N = round(rho_O Lx Ly) points
    std::string minimizer         = "lbfgs";  ///< lbfgs, anneal (Monte Carlo) or anneal+lbfgs
    std::string init_file         = "";     ///< Continuum: start from this obstacle file (its box is used)
    double      mc_step           = 0.0;    ///< Continuum Monte Carlo: initial displacement (d); 0 = 0.1 mean spacing
    int         lbfgs_memory      = 10;
    long long   max_iterations    = 100000;
    double      grad_tol          = 0.0;    ///< Stop when the largest gradient component <= grad_tol (0 = off)
    double      residual_tol      = 1e-20;  ///< Stop when <S>_mask <= residual_tol and no pair is closer than the core
    double      core_weight       = 1.0;    ///< Hard-core term relative to the spectral residual
    double      sk_kmax           = 0.0;    ///< S(k) output up to this |k| (d^-1); 0 = 3 K
    std::string output_file       = "hyperuniform.pgm";
    std::string annealing_file    = "annealing.dat";  ///< T/E/rate monitoring per temperature level
    std::string sk_file           = "sk.dat";         ///< Radially averaged structure factor S(k)

    bool K_given = false;  ///< K appeared in the parameter file

    /// With Lx, Ly, rho_O > 0: lattice_constant = sqrt(density / rho_O), Nx = round(Lx / a), Ny = round(Ly / a).
    /// Returns a description of the derived values (empty if nothing was derived).
    std::string derive_lattice();

    /// Hard-core distance actually imposed, in lattice spacings: max(min_distance, 2R / lattice_constant).
    [[nodiscard]] double effective_min_distance() const noexcept;

    /// Continuum: hard-core distance in d, max(min_distance, 2R) (1 + 1e-6), so that a residual overlap of the
    /// minimizer cannot bring a pair below 2R. 0 = no core term.
    [[nodiscard]] double continuum_core_distance() const noexcept;

    /// Continuum: number of points round(rho_O Lx Ly).
    [[nodiscard]] std::size_t continuum_points() const noexcept;

    [[nodiscard]] int nx() const noexcept { return Nx > 0 ? Nx : N; }
    [[nodiscard]] int ny() const noexcept { return Ny > 0 ? Ny : N; }

    /// Reads the parameter file. Returns false if it couldn't be opened or a value couldn't be parsed.
    [[nodiscard]] bool read(const std::string &filename);

    /// Checks ranges; returns one message per problem (empty if all valid).
    [[nodiscard]] std::vector<std::string> validate() const;

    /// Checks of mode continuum (called by validate).
    [[nodiscard]] std::vector<std::string> validate_continuum() const;

    /// Prints the active parameters to stdout.
    void print() const;
};
