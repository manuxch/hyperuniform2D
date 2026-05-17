#include "parametros.hpp"

#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>

namespace {

/// Removes leading and trailing whitespace from a string_view.
[[nodiscard]] std::string_view trim(std::string_view sv) noexcept {
  const auto start = sv.find_first_not_of(" \t");
  if (start == sv.npos)
    return {};
  const auto end = sv.find_last_not_of(" \t");
  return sv.substr(start, end - start + 1);
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
          {"N", [&](std::string_view v) { N = std::stoi(std::string(v)); }},
          {"density",
           [&](std::string_view v) { density = std::stod(std::string(v)); }},
          {"K", [&](std::string_view v) { K = std::stod(std::string(v)); }},
          {"T_initial",
           [&](std::string_view v) { T_initial = std::stod(std::string(v)); }},
          {"T_final",
           [&](std::string_view v) { T_final = std::stod(std::string(v)); }},
          {"cooling_rate",
           [&](std::string_view v) {
             cooling_rate = std::stod(std::string(v));
           }},
          {"steps_per_temp",
           [&](std::string_view v) {
             steps_per_temp = std::stoi(std::string(v));
           }},
          {"seed",
           [&](std::string_view v) { seed = std::stoi(std::string(v)); }},
          {"num_threads",
           [&](std::string_view v) { num_threads = std::stoi(std::string(v)); }},
          {"output_file",
           [&](std::string_view v) { output_file = std::string(v); }},
          {"annealing_file",
           [&](std::string_view v) { annealing_file = std::string(v); }},
          {"sk_file",
           [&](std::string_view v) { sk_file = std::string(v); }},
      };

  std::string line;
  while (std::getline(file, line)) {
    const auto sv = trim(std::string_view(line));
    if (sv.empty() || sv[0] == '#')
      continue;

    const auto pos = sv.find(':');
    if (pos == std::string_view::npos) {
      std::cerr << std::format(
          "Warning: malformed line (ignored): '{}'\n", line);
      continue;
    }

    const auto key = std::string(trim(sv.substr(0, pos)));
    const auto value = trim(sv.substr(pos + 1));

    if (const auto it = parsers.find(key); it != parsers.end()) {
      it->second(value);
    } else {
      std::cerr << std::format(
          "Warning: unknown key '{}' ignored.\n", key);
    }
  }
  return true;
}

void Parameters::print() const {
    std::cout << std::format(
        "Simulation parameters:\n"
        "  N                 = {}\n"
        "  density           = {}\n"
        "  K                 = {}\n"
        "  T_initial         = {}\n"
        "  T_final           = {}\n"
        "  cooling_rate      = {}\n"
        "  steps_per_temp    = {}\n"
        "  seed              = {}\n"
        "  num_threads       = {}\n"
        "  output_file       = {}\n"
        "  annealing_file    = {}\n"
        "  sk_file           = {}\n",
        N, density, K, T_initial, T_final, cooling_rate,
        steps_per_temp, seed, num_threads,
        output_file, annealing_file, sk_file);
}
