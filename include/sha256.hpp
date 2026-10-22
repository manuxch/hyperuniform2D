#pragma once

#include <string>
#include <string_view>

/// SHA-256 (FIPS 180-4) of a byte string, as 64 lowercase hexadecimal digits.
/// Used for sha256_body in the obstacle file format; no external dependency.
[[nodiscard]] std::string sha256_hex(std::string_view data);
