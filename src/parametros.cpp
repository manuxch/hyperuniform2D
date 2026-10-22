#include "parametros.hpp"

#include <charconv>
#include <cmath>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <algorithm>

namespace {

/// Removes leading and trailing whitespace from a string_view.
[[nodiscard]] std::string_view trim(std::string_view sv) noexcept {
  const auto start = sv.find_first_not_of(" \t\r");
  if (start == sv.npos)
    return {};
  const auto end = sv.find_last_not_of(" \t\r");
  return sv.substr(start, end - start + 1);
}

/// Parses the whole value as a number; throws std::invalid_argument otherwise
/// (std::stoi would accept "64abc" and abort the program on "abc").
template <typename T>
[[nodiscard]] T parse_number(std::string_view v) {
  T out{};
  const auto [ptr, ec] = std::from_chars(v.data(), v.data() + v.size(), out);
  if (ec != std::errc{} || ptr != v.data() + v.size())
    throw std::invalid_argument(std::format("'{}' is not a valid number", v));
  return out;
}

} // namespace

bool Parameters::read(const std::string &filename) {
  std::ifstream file(filename);
  if (!file) {
    std::cerr << std::format("Error: could not open '{}'\n", filename);
    return false;
  }

  // Parsing table: each key maps to a function that assigns its value.
  const std::unordered_map<std::string, std::function<void(std::string_view)>>
      parsers = {
          {"mode", [&](std::string_view v) { mode = std::string(v); }},
          {"N", [&](std::string_view v) { N = parse_number<int>(v); }},
          {"Nx", [&](std::string_view v) { Nx = parse_number<int>(v); }},
          {"Ny", [&](std::string_view v) { Ny = parse_number<int>(v); }},
          {"kappa", [&](std::string_view v) { kappa = parse_number<double>(v); }},
          {"min_distance", [&](std::string_view v) { min_distance = parse_number<double>(v); }},
          {"obs_file", [&](std::string_view v) { obs_file = std::string(v); }},
          {"minimizer", [&](std::string_view v) { minimizer = std::string(v); }},
          {"init_file", [&](std::string_view v) { init_file = std::string(v); }},
          {"mc_step", [&](std::string_view v) { mc_step = parse_number<double>(v); }},
          {"lbfgs_memory", [&](std::string_view v) { lbfgs_memory = parse_number<int>(v); }},
          {"max_iterations", [&](std::string_view v) { max_iterations = parse_number<long long>(v); }},
          {"grad_tol", [&](std::string_view v) { grad_tol = parse_number<double>(v); }},
          {"residual_tol", [&](std::string_view v) { residual_tol = parse_number<double>(v); }},
          {"core_weight", [&](std::string_view v) { core_weight = parse_number<double>(v); }},
          {"sk_kmax", [&](std::string_view v) { sk_kmax = parse_number<double>(v); }},
          {"lattice_constant", [&](std::string_view v) { lattice_constant = parse_number<double>(v); }},
          {"R", [&](std::string_view v) { R = parse_number<double>(v); }},
          {"boundary", [&](std::string_view v) { boundary = std::string(v); }},
          {"Lx", [&](std::string_view v) { Lx = parse_number<double>(v); }},
          {"Ly", [&](std::string_view v) { Ly = parse_number<double>(v); }},
          {"rho_O", [&](std::string_view v) { rho_O = parse_number<double>(v); }},
          {"arr_index", [&](std::string_view v) { arr_index = parse_number<int>(v); }},
          {"realizations", [&](std::string_view v) { realizations = parse_number<int>(v); }},
          {"density", [&](std::string_view v) { density = parse_number<double>(v); }},
          {"K", [&](std::string_view v) { K = parse_number<double>(v); K_given = true; }},
          {"chi", [&](std::string_view v) { chi = parse_number<double>(v); }},
          {"T_initial", [&](std::string_view v) { T_initial = parse_number<double>(v); }},
          {"T_final", [&](std::string_view v) { T_final = parse_number<double>(v); }},
          {"cooling_rate", [&](std::string_view v) { cooling_rate = parse_number<double>(v); }},
          {"steps_per_temp", [&](std::string_view v) { steps_per_temp = parse_number<int>(v); }},
          {"local_fraction", [&](std::string_view v) { local_fraction = parse_number<double>(v); }},
          {"local_radius", [&](std::string_view v) { local_radius = parse_number<int>(v); }},
          {"quench_steps", [&](std::string_view v) { quench_steps = parse_number<long long>(v); }},
          {"weight_power", [&](std::string_view v) { weight_power = parse_number<double>(v); }},
          {"energy_tol", [&](std::string_view v) { energy_tol = parse_number<double>(v); }},
          {"resync_interval", [&](std::string_view v) { resync_interval = parse_number<long long>(v); }},
          {"seed", [&](std::string_view v) { seed = parse_number<unsigned>(v); }},
          {"num_threads", [&](std::string_view v) { num_threads = parse_number<int>(v); }},
          {"parallel_min_modes", [&](std::string_view v) { parallel_min_modes = parse_number<int>(v); }},
          {"output_file", [&](std::string_view v) { output_file = std::string(v); }},
          {"annealing_file", [&](std::string_view v) { annealing_file = std::string(v); }},
          {"sk_file", [&](std::string_view v) { sk_file = std::string(v); }},
      };

  bool ok = true;
  std::string line;
  int line_no = 0;
  while (std::getline(file, line)) {
    ++line_no;
    auto sv = trim(std::string_view(line));
    if (sv.empty() || sv[0] == '#')
      continue;
    if (const auto hash = sv.find('#'); hash != std::string_view::npos)
      sv = trim(sv.substr(0, hash));  // trailing comment

    const auto pos = sv.find(':');
    if (pos == std::string_view::npos) {
      std::cerr << std::format(
          "Warning: malformed line {} (ignored): '{}'\n", line_no, line);
      continue;
    }

    const auto key = std::string(trim(sv.substr(0, pos)));
    const auto value = trim(sv.substr(pos + 1));

    if (const auto it = parsers.find(key); it != parsers.end()) {
      try {
        it->second(value);
      } catch (const std::exception &e) {
        std::cerr << std::format("Error: line {}, key '{}': {}\n", line_no, key, e.what());
        ok = false;
      }
    } else {
      std::cerr << std::format(
          "Warning: unknown key '{}' ignored.\n", key);
    }
  }
  return ok;
}

std::vector<std::string> Parameters::validate() const {
  std::vector<std::string> err;
  if (mode == "continuum") return validate_continuum();
  if (mode != "lattice")
    err.push_back(std::format("mode = '{}' must be 'lattice' or 'continuum'", mode));
  if (Nx == 0 || Ny == 0) {
    if (N < 4 || N > 46340)
      err.push_back(std::format("N = {} must be in [4, 46340]", N));
  }
  for (const auto& [name, v] : {std::pair{"Nx", Nx}, std::pair{"Ny", Ny}})
    if (v != 0 && (v < 4 || v > 46340))
      err.push_back(std::format("{} = {} must be in [4, 46340] (or 0 to use N)", name, v));
  if (static_cast<long long>(nx()) * ny() > 46340LL * 46340LL)
    err.push_back(std::format("Nx Ny = {} x {} has too many sites (at most 46340^2)", nx(), ny()));
  if (!(density > 0.0 && density < 1.0))
    err.push_back(std::format("density = {} must be in (0, 1)", density));
  if (chi < 0.0 || !std::isfinite(chi))
    err.push_back(std::format("chi = {} must be >= 0 (0 = use K)", chi));
  if (min_distance < 0.0 || !std::isfinite(min_distance) || (min_distance > 0.0 && 2.0 * std::ceil(min_distance) >= std::min(nx(), ny())))
    err.push_back(std::format("min_distance = {} must be >= 0 and below min(Nx, Ny)/2", min_distance));
  if (!(lattice_constant > 0.0) || !std::isfinite(lattice_constant))
    err.push_back(std::format("lattice_constant = {} must be > 0", lattice_constant));
  if (R < 0.0 || !std::isfinite(R))
    err.push_back(std::format("R = {} must be >= 0", R));
  if (!obs_file.empty() && !(R > 0.0))
    err.push_back("obs_file needs the obstacle radius R > 0");
  if (boundary != "periodic_x" && boundary != "walls_x")
    err.push_back(std::format("boundary = '{}' must be periodic_x or walls_x", boundary));
  {
    const int given = (Lx > 0.0) + (Ly > 0.0) + (rho_O > 0.0);
    if (given != 0 && given != 3)
      err.push_back("Lx, Ly and rho_O go together: give all three (> 0) or none");
    if (Lx < 0.0 || Ly < 0.0 || rho_O < 0.0)
      err.push_back("Lx, Ly and rho_O must be >= 0");
  }
  if (arr_index < 0)
    err.push_back(std::format("arr_index = {} must be >= 0", arr_index));
  if (realizations < 1)
    err.push_back(std::format("realizations = {} must be >= 1", realizations));
  if (const double s = effective_min_distance(); s > 0.0 && 2.0 * std::ceil(s) >= std::min(nx(), ny()))
    err.push_back(std::format("2R / lattice_constant = {} must be below min(Nx, Ny)/2", s));
  if (kappa < 0.0 || !std::isfinite(kappa))
    err.push_back(std::format("kappa = {} must be >= 0 (0 = use chi or K)", kappa));
  if (kappa == 0.0 && chi == 0.0) {
    if (nx() != ny())
      err.push_back(std::format("a {} x {} lattice needs kappa or chi: K (units of 2pi/N) is defined only for "
                                "square lattices", nx(), ny()));
    else if (!(K >= 1.0))
      err.push_back(std::format("K = {} must be >= 1 (the first shell |k| = 1)", K));
  }
  if (!(cooling_rate > 0.0 && cooling_rate < 1.0))
    err.push_back(std::format("cooling_rate = {} must be in (0, 1)", cooling_rate));
  if (!(T_final > 0.0))
    err.push_back(std::format("T_final = {} must be > 0", T_final));
  if (T_initial > 0.0 && T_initial <= T_final)
    err.push_back(std::format("T_initial = {} must be > T_final = {} (or <= 0 for automatic)", T_initial, T_final));
  if (steps_per_temp < 1)
    err.push_back(std::format("steps_per_temp = {} must be >= 1", steps_per_temp));
  if (!(local_fraction >= 0.0 && local_fraction <= 1.0))
    err.push_back(std::format("local_fraction = {} must be in [0, 1]", local_fraction));
  if (local_radius < 1 || local_radius >= std::min(nx(), ny()) / 2)
    err.push_back(std::format("local_radius = {} must be in [1, min(Nx, Ny)/2)", local_radius));
  if (quench_steps < 0)
    err.push_back(std::format("quench_steps = {} must be >= 0", quench_steps));
  if (!(weight_power >= 0.0 && weight_power <= 8.0))
    err.push_back(std::format("weight_power = {} must be in [0, 8]", weight_power));
  if (energy_tol < 0.0)
    err.push_back(std::format("energy_tol = {} must be >= 0", energy_tol));
  if (resync_interval < 1)
    err.push_back(std::format("resync_interval = {} must be >= 1", resync_interval));
  if (num_threads < 0)
    err.push_back(std::format("num_threads = {} must be >= 0", num_threads));
  return err;
}

std::vector<std::string> Parameters::validate_continuum() const {
  std::vector<std::string> err;
  if (!(Lx > 0.0 && Ly > 0.0 && rho_O > 0.0) || !std::isfinite(Lx * Ly * rho_O))
    err.push_back("mode continuum needs Lx, Ly and rho_O > 0");
  else if (continuum_points() < 2)
    err.push_back(std::format("rho_O Lx Ly = {} gives fewer than 2 points", rho_O * Lx * Ly));
  if (!(kappa > 0.0) && !(chi > 0.0))
    err.push_back("mode continuum needs kappa or chi (point-pattern convention) > 0");
  if (K_given)
    err.push_back("K is not used in mode continuum: give kappa or chi");
  if (kappa < 0.0 || chi < 0.0 || !std::isfinite(kappa) || !std::isfinite(chi))
    err.push_back("kappa and chi must be >= 0");
  if (minimizer != "lbfgs" && minimizer != "anneal" && minimizer != "anneal+lbfgs")
    err.push_back(std::format("minimizer = '{}' must be lbfgs, anneal or anneal+lbfgs", minimizer));
  if (minimizer != "lbfgs") {
    if (!(cooling_rate > 0.0 && cooling_rate < 1.0))
      err.push_back(std::format("cooling_rate = {} must be in (0, 1)", cooling_rate));
    if (!(T_final > 0.0))
      err.push_back(std::format("T_final = {} must be > 0", T_final));
    if (T_initial > 0.0 && T_initial <= T_final)
      err.push_back(std::format("T_initial = {} must be > T_final = {} (or <= 0 for automatic)", T_initial, T_final));
    if (steps_per_temp < 1 || quench_steps < 0 || resync_interval < 1)
      err.push_back("steps_per_temp and resync_interval must be >= 1, quench_steps >= 0");
    if (mc_step < 0.0)
      err.push_back(std::format("mc_step = {} must be >= 0", mc_step));
  }
  if (lbfgs_memory < 1 || lbfgs_memory > 100)
    err.push_back(std::format("lbfgs_memory = {} must be in [1, 100]", lbfgs_memory));
  if (max_iterations < 1)
    err.push_back(std::format("max_iterations = {} must be >= 1", max_iterations));
  if (grad_tol < 0.0 || residual_tol < 0.0)
    err.push_back("grad_tol and residual_tol must be >= 0");
  if (!(core_weight > 0.0))
    err.push_back(std::format("core_weight = {} must be > 0", core_weight));
  if (sk_kmax < 0.0)
    err.push_back(std::format("sk_kmax = {} must be >= 0", sk_kmax));
  if (!(weight_power >= 0.0 && weight_power <= 8.0))
    err.push_back(std::format("weight_power = {} must be in [0, 8]", weight_power));
  if (R < 0.0 || min_distance < 0.0)
    err.push_back("R and min_distance must be >= 0");
  if (!obs_file.empty() && !(R > 0.0))
    err.push_back("obs_file needs the obstacle radius R > 0");
  if (boundary != "periodic_x" && boundary != "walls_x")
    err.push_back(std::format("boundary = '{}' must be periodic_x or walls_x", boundary));
  if (arr_index < 0 || realizations < 1)
    err.push_back("arr_index must be >= 0 and realizations >= 1");
  if (num_threads < 0)
    err.push_back(std::format("num_threads = {} must be >= 0", num_threads));
  return err;
}

double Parameters::continuum_core_distance() const noexcept {
  const double s = std::max(min_distance, 2.0 * R);
  return s > 0.0 ? s * (1.0 + 1e-6) : 0.0;
}

std::size_t Parameters::continuum_points() const noexcept {
  const double n = rho_O * Lx * Ly;
  return n > 0.0 && std::isfinite(n) ? static_cast<std::size_t>(std::llround(n)) : 0;
}

std::string Parameters::derive_lattice() {
  if (mode != "lattice") return {};
  if (!(Lx > 0.0 && Ly > 0.0 && rho_O > 0.0 && density > 0.0)) return {};
  // Lx is kept exact (it is the swept variable of a campaign): Nx = round(Lx / a0) with a0 = sqrt(density / rho_O),
  // a = Lx / Nx; Ly is the closest multiple of a. The occupation fraction is then set so that the number of
  // occupied sites is round(rho_O Lx Ly_eff), i.e. the obstacle density stays rho_O.
  const double target_density = density;
  const double a0 = std::sqrt(target_density / rho_O);
  Nx = std::max(1, static_cast<int>(std::lround(Lx / a0)));
  lattice_constant = Lx / Nx;
  Ny = std::max(1, static_cast<int>(std::lround(Ly / lattice_constant)));
  density = rho_O * lattice_constant * lattice_constant;
  return std::format("Lattice derived from Lx = {}, Ly = {}, rho_O = {}, occupation {}: lattice_constant = {:.9g}, "
                     "{} x {} sites, Lx = {:.9g}, effective Ly = {:.9g}, occupation {:.6g}",
                     Lx, Ly, rho_O, target_density, lattice_constant, Nx, Ny, Nx * lattice_constant,
                     Ny * lattice_constant, density);
}

double Parameters::effective_min_distance() const noexcept {
  // Slightly above 2R / a, so that a pair at exactly that distance, once written with rounding, still
  // measures >= 2R.
  const double from_R = R > 0.0 && lattice_constant > 0.0 ? 2.0 * R / lattice_constant * (1.0 + 1e-9) : 0.0;
  return std::max(min_distance, from_R);
}

void Parameters::print() const {
  if (mode == "continuum") {
    std::cout << std::format(
        "Simulation parameters:\n"
        "  mode              = continuum\n"
        "  Lx x Ly           = {} x {}\n"
        "  rho_O             = {} ({} points)\n"
        "  kappa             = {}{}\n"
        "  chi               = {}{}\n"
        "  R                 = {} (hard core {:.9g})\n"
        "  weight_power      = {}\n"
        "  minimizer         = {}, memory {}, max_iterations {}, grad_tol {}, residual_tol {}\n"
        "  core_weight       = {}\n"
        "  seed              = {}, arr_index {}, realizations {}\n"
        "  num_threads       = {}\n"
        "  boundary          = {}\n"
        "  obs_file          = {}\n"
        "  init_file         = {}\n"
        "  annealing_file    = {} (minimization log)\n"
        "  sk_file           = {}\n",
        Lx, Ly, rho_O, continuum_points(), kappa, kappa > 0.0 ? "" : " (unused)",
        chi, kappa > 0.0 ? " (ignored: mask from kappa)" : " (point_pattern convention)",
        R, continuum_core_distance(), weight_power, minimizer, lbfgs_memory, max_iterations, grad_tol, residual_tol,
        core_weight, seed, arr_index, realizations, num_threads, boundary, obs_file.empty() ? "(none)" : obs_file,
        init_file.empty() ? "(random start)" : init_file, annealing_file, sk_file);
    if (minimizer != "lbfgs")
      std::cout << std::format("  Monte Carlo       : T_initial {}{}, T_final {}, cooling_rate {}, steps_per_temp {}, "
                               "quench_steps {}, mc_step {}{}\n", T_initial, T_initial <= 0.0 ? " (automatic)" : "",
                               T_final, cooling_rate, steps_per_temp, quench_steps, mc_step,
                               mc_step > 0.0 ? "" : " (0.1 mean spacing)");
    return;
  }
    std::cout << std::format(
        "Simulation parameters:\n"
        "  mode              = {}\n"
        "  Nx x Ny           = {} x {}\n"
        "  density           = {}\n"
        "  K                 = {}{}\n"
        "  chi               = {}{}\n"
        "  kappa             = {}{}\n"
        "  min_distance      = {}{}\n"
        "  lattice_constant  = {}\n"
        "  R                 = {}\n"
        "  boundary          = {}\n"
        "  obs_file          = {}\n"
        "  arr_index         = {}, realizations = {}\n"
        "  T_initial         = {}{}\n"
        "  T_final           = {}\n"
        "  cooling_rate      = {}\n"
        "  steps_per_temp    = {}\n"
        "  local_fraction    = {}\n"
        "  local_radius      = {}\n"
        "  quench_steps      = {}\n"
        "  weight_power      = {}\n"
        "  energy_tol        = {}\n"
        "  resync_interval   = {}\n"
        "  seed              = {}\n"
        "  num_threads       = {}\n"
        "  parallel_min_modes= {}\n"
        "  output_file       = {}\n"
        "  annealing_file    = {}\n"
        "  sk_file           = {}\n",
        mode, nx(), ny(), density, K, (kappa > 0.0 || chi > 0.0) ? " (ignored: derived from kappa or chi)" : "",
        chi, kappa > 0.0 ? " (ignored: mask from kappa)" : chi > 0.0 ? " (binary_lattice convention)" : " (unused: mask from K)",
        kappa, kappa > 0.0 ? "" : " (unused)",
        min_distance, effective_min_distance() > 0.0 ? std::format(" (imposed: {:.6g} lattice spacings)", effective_min_distance()) : " (no exclusion)",
        lattice_constant, R, boundary, obs_file.empty() ? "(none)" : obs_file, arr_index, realizations,
        T_initial, T_initial <= 0.0 ? " (automatic)" : "",
        T_final, cooling_rate, steps_per_temp, local_fraction, local_radius, quench_steps, weight_power, energy_tol, resync_interval,
        seed, num_threads, parallel_min_modes,
        output_file, annealing_file, sk_file);
}
