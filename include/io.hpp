#pragma once

#include "config.hpp"

#include <string>
#include <vector>

/// Saves the configuration in PGM format (P2, text), occupied = 1.
/// Each entry of `comments` is written as a "# ..." line after the magic number
/// (allowed by the PGM format; readers skip them). Uses an internal string buffer.
void save_PGM(const HyperuniformConfig& config, const std::string& filename,
              const std::vector<std::string>& comments = {});
