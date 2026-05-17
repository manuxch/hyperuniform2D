#pragma once

#include "config.hpp"

#include <string>

/// Saves the configuration in PGM format (P2, text).
/// Uses an internal string buffer to minimize I/O calls.
void save_PGM(const HyperuniformConfig& config, const std::string& filename);
