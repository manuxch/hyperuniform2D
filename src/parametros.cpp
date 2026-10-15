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
          {"N", [&](std::string_view v) { N = parse_number<int>(v); }},
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
  if (N < 4 || N > 46340)
    err.push_back(std::format("N = {} must be in [4, 46340]", N));
  if (!(density > 0.0 && density < 1.0))
    err.push_back(std::format("density = {} must be in (0, 1)", density));
  if (chi < 0.0 || !std::isfinite(chi))
    err.push_back(std::format("chi = {} must be >= 0 (0 = use K)", chi));
  if (chi == 0.0 && !(K >= 1.0))
    err.push_back(std::format("K = {} must be >= 1 (the first shell |k| = 1)", K));
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
  if (local_radius < 1 || local_radius >= N / 2)
    err.push_back(std::format("local_radius = {} must be in [1, N/2)", local_radius));
  if (quench_steps < 0)
    err.push_back(std::format("quench_steps = {} must be >= 0", quench_steps));
  if (energy_tol < 0.0)
    err.push_back(std::format("energy_tol = {} must be >= 0", energy_tol));
  if (resync_interval < 1)
    err.push_back(std::format("resync_interval = {} must be >= 1", resync_interval));
  if (num_threads < 0)
    err.push_back(std::format("num_threads = {} must be >= 0", num_threads));
  return err;
}

void Parameters::print() const {
    std::cout << std::format(
        "Simulation parameters:\n"
        "  N                 = {}\n"
        "  density           = {}\n"
        "  K                 = {}{}\n"
        "  chi               = {}{}\n"
        "  T_initial         = {}{}\n"
        "  T_final           = {}\n"
        "  cooling_rate      = {}\n"
        "  steps_per_temp    = {}\n"
        "  local_fraction    = {}\n"
        "  local_radius      = {}\n"
        "  quench_steps      = {}\n"
        "  energy_tol        = {}\n"
        "  resync_interval   = {}\n"
        "  seed              = {}\n"
        "  num_threads       = {}\n"
        "  parallel_min_modes= {}\n"
        "  output_file       = {}\n"
        "  annealing_file    = {}\n"
        "  sk_file           = {}\n",
        N, density, K, chi > 0.0 ? " (ignored: derived from chi)" : "",
        chi, chi > 0.0 ? "" : " (unused: mask from K)",
        T_initial, T_initial <= 0.0 ? " (automatic)" : "",
        T_final, cooling_rate, steps_per_temp, local_fraction, local_radius, quench_steps, energy_tol, resync_interval,
        seed, num_threads, parallel_min_modes,
        output_file, annealing_file, sk_file);
}
