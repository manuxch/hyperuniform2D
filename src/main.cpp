#include "annealing.hpp"
#include "config.hpp"
#include "continuum_driver.hpp"
#include "fourier.hpp"
#include "io.hpp"
#include "obs.hpp"
#include "parametros.hpp"
#include "constants.hpp"

#include <cmath>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <format>
#include <iostream>
#include <string>
#include <vector>

namespace {

/// |k| of the shell with key q, in d^-1: |k| = 2π sqrt(q) / (Nx Ny a).
double k_of_q(std::int64_t q, int Nx, int Ny, double a) {
    return 2.0 * PI * std::sqrt(static_cast<double>(q)) / (static_cast<double>(Nx) * Ny * a);
}

void run_realization(const Parameters& p, unsigned seed, int arr) {
    const bool several = p.realizations > 1;
    HyperuniformConfig config;
    const int Nx = p.nx(), Ny = p.ny();
    const double sigma = p.effective_min_distance();
    initialize_grid(config, Nx, Ny, p.density, seed, sigma);
    const std::size_t n_chi = chi_points(Nx, Ny, config.M);

    // K in units of 2π/sqrt(Nx Ny) (2π/N if square). kappa = K_phys a / 2π with a = sqrt(Nx Ny / M) the mean
    // spacing of the occupied sites, so K = kappa sqrt(M).
    double K2 = p.K * p.K;
    if (p.kappa > 0.0) {
        K2 = p.kappa * p.kappa * static_cast<double>(config.M);
    } else if (p.chi > 0.0) {
        const MaskChoice mc = choose_mask_for_chi(Nx, Ny, n_chi, p.chi);
        K2 = mc.K2;
    }
    build_mask(config, K2, p.weight_power);
    const double chi = chi_of(config.modes.size(), n_chi);
    const double chi_pp = chi_point_pattern(config.modes.size(), static_cast<std::size_t>(config.M));
    const double kappa = config.K / std::sqrt(static_cast<double>(config.M));
    std::cout << std::format("Occupied sites: {} (density {:.6f}); chi uses n = {} points\n",
                             config.M, static_cast<double>(config.M) / (static_cast<double>(Nx) * Ny), n_chi);
    std::cout << std::format("Mask: K = {:.6g} (K² = {:.6g}{}), chi = {:.6f}{}\n", config.K, config.K2,
                             config.square() ? "" : ", units of 2pi/sqrt(Nx Ny)", chi,
                             p.kappa > 0.0 ? "" : p.chi > 0.0 ? std::format(" (target {})", p.chi) : "");
    std::cout << std::format("      kappa = K a / 2pi = {:.6f}; chi binary_lattice = {:.6f}, "
                             "chi point_pattern M(K)/[2(M-1)] = {:.6f}\n", kappa, chi, chi_pp);
    if (chi_pp > 0.5)
        std::cerr << std::format("Warning: chi point_pattern = {:.4f} > 1/2: stealthy disordered ground states are "
                                 "not expected; the pattern tends to crystallize.\n", chi_pp);
    else if (chi > 0.5)
        std::cerr << std::format("Warning: chi binary_lattice = {:.4f} > 1/2: on the lattice the annealing leaves a "
                                 "residual that grows with K and with the occupation density.\n", chi);

    Parameters run = p;
    run.annealing_file = path_with_arr(p.annealing_file, arr, several);
    const AnnealReport rep = simulated_annealing(config, run);
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
    const double a_mean = std::sqrt(static_cast<double>(Nx) * Ny / config.M);
    const double dmin = min_pair_distance(config, 2.0 * a_mean);
    std::cout << std::format("Minimum distance between occupied sites: {:.4g} lattice spacings = {:.4f} a "
                             "(a = {:.4g}, mean spacing){}\n", dmin, dmin / a_mean, a_mean,
                             sigma > 0.0 ? std::format("; exclusion {:.6g} rejected {} proposals",
                                                       sigma, config.hc_rejections) : "");

    std::vector<std::string> header{
        "hyperuniform2D " HYPER2D_VERSION ": stealthy hyperuniform pattern, 1 = occupied",
    };
    if (config.square()) {
        header.push_back(std::format("N {} density {:.6f} occupied {}", Nx,
                                     static_cast<double>(config.M) / (static_cast<double>(Nx) * Ny), config.M));
        header.push_back(std::format("K {:.6g} K2 {:.6g} modes_half_plane {} chi {:.6f} chi_target {}", config.K,
                                     config.K2, config.modes.size(), chi, p.chi));
    } else {
        header.push_back(std::format("Nx {} Ny {} density {:.6f} occupied {}", Nx, Ny,
                                     static_cast<double>(config.M) / (static_cast<double>(Nx) * Ny), config.M));
        header.push_back(std::format("K {:.6g} K2 {:.6g} K_units 2pi/sqrt(NxNy) kappa {:.6f} modes_half_plane {} "
                                     "chi {:.6f} chi_convention binary_lattice chi_point_pattern {:.6f} chi_target {}",
                                     config.K, config.K2, kappa, config.modes.size(), chi, chi_pp, p.chi));
    }
    header.push_back(std::format("seed {} E_initial {:.6g} E_final {:.6g} residual_S_ratio {:.4g} residual_low_k {:.4g} weight_power {}",
                                 seed, rep.E_initial, rep.E_final, residual, low, p.weight_power));
    save_PGM(config, path_with_arr(p.output_file, arr, several), header);
    const RadialSk sk = radial_structure_factor(config);
    analyze_structure_factor(config, path_with_arr(p.sk_file, arr, several), sk);

    if (p.obs_file.empty()) return;

    // Obstacle file: site (i, j) at (i a, j a), lengths in grain diameters.
    const double a = p.lattice_constant;
    std::vector<ObsPoint> pts;
    pts.reserve(config.occupied.size());
    for (const std::uint32_t s : config.occupied)
        pts.push_back({static_cast<double>(s % static_cast<std::uint32_t>(Nx)) * a,
                       static_cast<double>(s / static_cast<std::uint32_t>(Nx)) * a});
    // K written midway between the outermost included shell and the first excluded one, so a reader that
    // rebuilds the mask from K, Lx and Ly gets the same modes despite rounding.
    std::int64_t q_in = 0;
    for (const auto& m : config.modes) q_in = std::max(q_in, mode_key(m.kx, m.ky, Nx, Ny));
    const double K_in = k_of_q(q_in, Nx, Ny, a);
    const double K_obs = config.q_next > 0 ? 0.5 * (K_in + k_of_q(config.q_next, Nx, Ny, a)) : K_in * (1.0 + 1e-9);

    ObsMeta meta;
    meta.generator_version = HYPER2D_VERSION;
    meta.mode = "lattice";
    meta.cls = "stealthy";
    meta.Lx = Nx * a;
    meta.Ly = Ny * a;
    meta.boundary = p.boundary;
    meta.R = p.R;
    meta.seed = seed;
    meta.K = K_obs;
    meta.n_constrained_modes = config.modes.size();
    meta.chi_convention = "binary_lattice";
    meta.chi = chi;
    meta.residual_S = residual;
    meta.lattice_constant = a;
    meta.notes = std::format("lattice {} x {}, occupation {:.6g}, weight_power {}, residual_low_k {:.4g}, "
                             "exclusion {:.6g} lattice spacings, chi_point_pattern {:.6f}",
                             Nx, Ny, static_cast<double>(config.M) / (static_cast<double>(Nx) * Ny), p.weight_power,
                             low, sigma, chi_pp);
    ObsResult obs = build_obs(meta, pts);
    const std::string obs_path = p.obs_file == "auto"
        ? obs_auto_name(meta.cls, obs.rho_O, parse_real(format_real(meta.Lx)), parse_real(format_real(meta.Ly)), arr)
        : path_with_arr(p.obs_file, arr, several);
    {
        std::ofstream f(obs_path, std::ios::binary);
        if (!f || !(f << obs.text)) throw std::runtime_error(std::format("could not write '{}'", obs_path));
    }
    if (obs.d_min < 2.0 * p.R)
        std::cerr << std::format("Warning: d_min = {} < 2R = {}: the file violates the no-overlap rule\n", obs.d_min, 2.0 * p.R);
    std::cout << std::format("Obstacle file saved to '{}' (n = {}, rho_O = {}, d_min = {}, sha256_body {})\n",
                             obs_path, obs.n, format_real(obs.rho_O), format_real(obs.d_min), obs.sha256_body);

    std::vector<double> k_d, S_d;
    const double k_unit = 2.0 * PI / (std::sqrt(static_cast<double>(Nx) * Ny) * a);
    for (std::size_t i = 0; i < sk.k.size(); ++i) { k_d.push_back(sk.k[i] * k_unit); S_d.push_back(sk.S[i]); }
    const std::string sk_path = path_with_extension(obs_path, ".sk");
    write_sk_companion(sk_path, {
        {"format", "1"}, {"generator", "hyperuniform2D"}, {"generator_version", HYPER2D_VERSION},
        {"obstacles_file", obs_path.substr(obs_path.find_last_of('/') + 1)}, {"sha256_body", obs.sha256_body},
        {"columns", "k S"}, {"k_units", "d^-1"},
        {"notes", "radially averaged S(k) = |rho_hat(k)|^2 / n over every lattice wave vector k != 0, bins of "
                  "width 2pi/sqrt(Lx Ly)"},
    }, k_d, S_d);
    std::cout << std::format("Structure factor companion saved to '{}'\n", sk_path);
}

}  // namespace

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
    if (const std::string derived = p.derive_lattice(); !derived.empty()) std::cout << derived << '\n';
    if (p.mode == "continuum" && !p.init_file.empty()) {
        // The box comes from the initial file; Lx, Ly, rho_O of the parameter file, if given, must agree.
        try {
            const ObsData init = read_obs(p.init_file);
            const double Lx = init.number("Lx"), Ly = init.number("Ly");
            const double rho = static_cast<double>(init.points.size()) / (Lx * Ly);
            if ((p.Lx > 0.0 && std::abs(p.Lx - Lx) > 1e-9 * Lx) || (p.Ly > 0.0 && std::abs(p.Ly - Ly) > 1e-9 * Ly) ||
                (p.rho_O > 0.0 && std::abs(p.rho_O - rho) > 1e-9 * rho)) {
                std::cerr << std::format("Error: Lx, Ly, rho_O = {}, {}, {} disagree with init_file '{}' ({}, {}, {})\n",
                                         p.Lx, p.Ly, p.rho_O, p.init_file, Lx, Ly, rho);
                return 1;
            }
            p.Lx = Lx;
            p.Ly = Ly;
            p.rho_O = rho;
            std::cout << std::format("Initial configuration '{}': {} points in {} x {}\n", p.init_file, init.points.size(), Lx, Ly);
        } catch (const std::exception& e) {
            std::cerr << "Error: init_file: " << e.what() << '\n';
            return 1;
        }
    }
    if (const auto errors = p.validate(); !errors.empty()) {
        for (const auto& e : errors) std::cerr << "Error: " << e << '\n';
        return 1;
    }
    if (p.kappa > 0.0 && (p.chi > 0.0 || p.K_given))
        std::cerr << "Warning: kappa given together with chi or K; the mask is set by kappa and the others are ignored.\n";
    else if (p.chi > 0.0 && p.K_given)
        std::cerr << "Warning: both chi and K given; K is derived from chi and the K value is ignored.\n";
    p.print();

    try {
        for (int r = 0; r < p.realizations; ++r) {
            const unsigned seed = p.seed + static_cast<unsigned>(r);
            const int arr = p.arr_index + r;
            if (p.realizations > 1) std::cout << std::format("\n=== Realization {} (seed {}) ===\n", arr, seed);
            if (p.mode == "continuum") run_continuum(p, seed, arr);
            else run_realization(p, seed, arr);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
