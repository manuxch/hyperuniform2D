#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

/// Writer for the obstacle file format of hyperGalton2D, version 1 (docs/formats/obstacles-format.pdf
/// of that project). Every value the header declares about the body (n_obstacles, rho_O, d_min,
/// sha256_body, kappa) is computed from the numbers exactly as they are written, so a reader that
/// recomputes them finds the same values.

struct ObsPoint { double x, y; };

/// Header values that the writer does not derive from the body.
struct ObsMeta {
    std::string generator_version;
    std::string mode;                  ///< "lattice" or "continuum"
    std::string cls = "stealthy";      ///< class
    double      Lx = 0.0, Ly = 0.0;    ///< Domain [0, Lx) x [0, Ly), reduced units (d = 1)
    std::string boundary = "periodic_x";
    double      R = 0.0;               ///< Obstacle radius
    unsigned    seed = 0;
    std::string created;               ///< ISO-8601 date; empty = today (UTC)

    // Hyperuniformity diagnostics (required for class stealthy)
    std::optional<double>      K;                    ///< Mask radius, d^-1, 2π included
    std::optional<std::size_t> n_constrained_modes;  ///< M(K), half space, k = 0 excluded
    std::string                chi_convention;       ///< "binary_lattice" or "point_pattern"
    std::optional<double>      chi;
    bool                       write_kappa = true;   ///< kappa = K / (2π sqrt(rho_O)), from the written values
    std::optional<double>      residual_S;
    std::optional<double>      lattice_constant;
    std::string                notes;                ///< One line
};

struct ObsResult {
    std::string text;         ///< The whole file
    std::string sha256_body;
    std::size_t n = 0;
    double      rho_O = 0.0;
    double      d_min = 0.0;
};

/// Shortest decimal representation that reads back as the same double (C locale); -0 is written as 0.
[[nodiscard]] std::string format_real(double v);

/// Parses a number written by format_real (C locale). Throws std::invalid_argument on failure.
[[nodiscard]] double parse_real(const std::string& s);

/// Smallest centre-to-centre distance. y is never periodic; x is periodic (minimum image) if periodic_x.
/// O(N) with a cell list. Returns +inf for fewer than two points.
[[nodiscard]] double min_center_distance(const std::vector<ObsPoint>& pts, double Lx, double Ly, bool periodic_x);

/// Builds the file: formats the coordinates, orders the body by y and then x (on the written values),
/// checks 0 <= x < Lx and 0 <= y < Ly, and computes the derived header keys. Throws std::invalid_argument
/// if a coordinate is outside the domain, the boundary is unknown, or a stealthy diagnostic is missing.
[[nodiscard]] ObsResult build_obs(const ObsMeta& meta, const std::vector<ObsPoint>& pts);

/// build_obs, then writes the text to `path`. Throws std::runtime_error if the file cannot be written.
ObsResult write_obs(const std::string& path, const ObsMeta& meta, const std::vector<ObsPoint>& pts);

/// A file read back: header keys in order, and the body points.
struct ObsData {
    std::vector<std::pair<std::string, std::string>> header;
    std::vector<ObsPoint> points;
    [[nodiscard]] const std::string* get(const std::string& key) const;  ///< nullptr if absent
    [[nodiscard]] double number(const std::string& key) const;           ///< throws if absent or not a number
};

/// Reads an obstacle file. Checks the format version, the two-column body, the record count and
/// sha256_body (a reader that finds a different digest must abort); throws std::runtime_error otherwise.
[[nodiscard]] ObsData read_obs(const std::string& path);

/// Companion file <name>.sk: header lines "# key: value" and two columns, k (d^-1) and S(k).
void write_sk_companion(const std::string& path, const std::vector<std::pair<std::string, std::string>>& header,
                        const std::vector<double>& k, const std::vector<double>& S);

/// Recommended file name: <class>_rho<value>_Lx<value>_Ly<value>_arr<index>.obs (index on three digits).
[[nodiscard]] std::string obs_auto_name(const std::string& cls, double rho_O, double Lx, double Ly, int arr_index);

/// "name.ext" -> "name_arrNNN.ext" when `several` is true; unchanged otherwise.
[[nodiscard]] std::string path_with_arr(const std::string& path, int arr, bool several);

/// Replaces (or adds) the extension of the file name: "dir/a.obs", ".sk" -> "dir/a.sk".
[[nodiscard]] std::string path_with_extension(const std::string& path, const std::string& ext);

/// Today's date in UTC, YYYY-MM-DD.
[[nodiscard]] std::string iso_date_today();
