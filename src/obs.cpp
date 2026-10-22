#include "obs.hpp"
#include "constants.hpp"
#include "sha256.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <iterator>
#include <string_view>
#include <limits>
#include <stdexcept>

std::string format_real(double v) {
    if (v == 0.0) v = 0.0;  // -0 -> 0
    return std::format("{}", v);
}

double parse_real(const std::string& s) {
    double v = 0.0;
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || ptr != s.data() + s.size())
        throw std::invalid_argument(std::format("'{}' is not a valid number", s));
    return v;
}

double min_center_distance(const std::vector<ObsPoint>& pts, double Lx, double Ly, bool periodic_x) {
    const std::size_t n = pts.size();
    if (n < 2) return std::numeric_limits<double>::infinity();
    // Cell list of side h: any pair closer than h lies in neighbouring cells. Start near the mean
    // spacing and double h until some pair closer than h is found (or one cell covers the domain).
    double h = 0.75 * std::sqrt(Lx * Ly / static_cast<double>(n));
    while (true) {
        const int ncx = std::max(1, static_cast<int>(Lx / h));
        const int ncy = std::max(1, static_cast<int>(Ly / h));
        const double hx = Lx / ncx, hy = Ly / ncy;
        std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(ncx) * static_cast<std::size_t>(ncy));
        const auto cell_of = [&](const ObsPoint& p) {
            const int cx = std::min(ncx - 1, static_cast<int>(p.x / hx));
            const int cy = std::min(ncy - 1, static_cast<int>(p.y / hy));
            return std::pair{cx, cy};
        };
        for (std::size_t i = 0; i < n; ++i) {
            const auto [cx, cy] = cell_of(pts[i]);
            cells[static_cast<std::size_t>(cx + ncx * cy)].push_back(i);
        }
        double best2 = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; ++i) {
            const auto [cx, cy] = cell_of(pts[i]);
            for (int oy = -1; oy <= 1; ++oy) {
                const int y = cy + oy;
                if (y < 0 || y >= ncy) continue;
                for (int ox = -1; ox <= 1; ++ox) {
                    int x = cx + ox;
                    if (x < 0 || x >= ncx) {
                        if (!periodic_x) continue;
                        x = (x + ncx) % ncx;
                    }
                    for (const std::size_t j : cells[static_cast<std::size_t>(x + ncx * y)]) {
                        if (j <= i) continue;
                        double dx = std::abs(pts[i].x - pts[j].x);
                        if (periodic_x) dx = std::min(dx, Lx - dx);
                        const double dy = pts[i].y - pts[j].y;
                        best2 = std::min(best2, dx * dx + dy * dy);
                    }
                }
            }
        }
        const double hmin = std::min(hx, hy);
        // Exact once a pair closer than the cell side is found, or when the cells cannot get larger
        // (with fewer than 3 cells along an axis every pair is already visited).
        if (best2 <= hmin * hmin || (ncx < 3 && ncy < 3)) return std::sqrt(best2);
        h *= 2.0;
    }
}

ObsResult build_obs(const ObsMeta& meta, const std::vector<ObsPoint>& pts) {
    if (meta.boundary != "periodic_x" && meta.boundary != "walls_x")
        throw std::invalid_argument(std::format("boundary '{}' must be periodic_x or walls_x", meta.boundary));
    if (meta.cls == "stealthy" && (!meta.K || !meta.n_constrained_modes || meta.chi_convention.empty()))
        throw std::invalid_argument("class stealthy needs K, n_constrained_modes and chi_convention");
    if (meta.chi && meta.chi_convention.empty())
        throw std::invalid_argument("chi must be accompanied by chi_convention");

    const std::string Lx_s = format_real(meta.Lx), Ly_s = format_real(meta.Ly);
    const double Lx = parse_real(Lx_s), Ly = parse_real(Ly_s);

    // Format, read back, and order on the written values.
    struct Row { double x, y; std::string sx, sy; };
    std::vector<Row> rows;
    rows.reserve(pts.size());
    for (const auto& p : pts) {
        Row r{0.0, 0.0, format_real(p.x), format_real(p.y)};
        r.x = parse_real(r.sx);
        r.y = parse_real(r.sy);
        if (!(r.x >= 0.0 && r.x < Lx && r.y >= 0.0 && r.y < Ly))
            throw std::invalid_argument(std::format("obstacle ({}, {}) outside [0, {}) x [0, {})", r.sx, r.sy, Lx_s, Ly_s));
        rows.push_back(std::move(r));
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.y != b.y ? a.y < b.y : a.x < b.x; });

    std::string body;
    body.reserve(rows.size() * 24);
    std::vector<ObsPoint> written;
    written.reserve(rows.size());
    for (const auto& r : rows) {
        body += r.sx;
        body += ' ';
        body += r.sy;
        body += '\n';
        written.push_back({r.x, r.y});
    }

    ObsResult res;
    res.n = rows.size();
    res.rho_O = static_cast<double>(res.n) / (Lx * Ly);
    res.d_min = min_center_distance(written, Lx, Ly, meta.boundary == "periodic_x");
    res.sha256_body = sha256_hex(body);

    std::string h;
    const auto kv = [&h](const std::string& k, const std::string& v) { h += "# " + k + ": " + v + '\n'; };
    kv("format", "1");
    kv("generator", "hyperuniform2D");
    kv("generator_version", meta.generator_version);
    kv("mode", meta.mode);
    kv("class", meta.cls);
    kv("Lx", Lx_s);
    kv("Ly", Ly_s);
    kv("boundary", meta.boundary);
    kv("n_obstacles", std::to_string(res.n));
    kv("rho_O", format_real(res.rho_O));
    kv("R", format_real(meta.R));
    kv("d_min", std::isfinite(res.d_min) ? format_real(res.d_min) : std::string("inf"));
    kv("seed", std::to_string(meta.seed));
    kv("created", meta.created.empty() ? iso_date_today() : meta.created);
    kv("sha256_body", res.sha256_body);
    if (meta.K) {
        const std::string K_s = format_real(*meta.K);
        kv("K", K_s);
        if (meta.n_constrained_modes) kv("n_constrained_modes", std::to_string(*meta.n_constrained_modes));
        if (!meta.chi_convention.empty()) kv("chi_convention", meta.chi_convention);
        if (meta.chi) kv("chi", format_real(*meta.chi));
        if (meta.write_kappa) kv("kappa", format_real(parse_real(K_s) / (2.0 * PI * std::sqrt(res.rho_O))));
    }
    if (meta.residual_S) kv("residual_S", format_real(*meta.residual_S));
    if (meta.lattice_constant) kv("lattice_constant", format_real(*meta.lattice_constant));
    kv("columns", "x y");
    if (!meta.notes.empty()) kv("notes", meta.notes);

    res.text = h + body;
    return res;
}

ObsResult write_obs(const std::string& path, const ObsMeta& meta, const std::vector<ObsPoint>& pts) {
    ObsResult res = build_obs(meta, pts);
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(std::format("could not create '{}'", path));
    f << res.text;
    if (!f) throw std::runtime_error(std::format("could not write '{}'", path));
    return res;
}

const std::string* ObsData::get(const std::string& key) const {
    for (const auto& [k, v] : header)
        if (k == key) return &v;
    return nullptr;
}

double ObsData::number(const std::string& key) const {
    const std::string* v = get(key);
    if (!v) throw std::runtime_error(std::format("obstacle file: missing key '{}'", key));
    try {
        return parse_real(*v);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error(std::format("obstacle file: key '{}' = '{}' is not a number", key, *v));
    }
}

ObsData read_obs(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(std::format("could not open '{}'", path));
    const std::string raw((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    ObsData d;
    std::size_t pos = 0;
    while (pos < raw.size() && raw[pos] == '#') {
        const std::size_t end = raw.find('\n', pos);
        const std::string line = raw.substr(pos + 1, (end == std::string::npos ? raw.size() : end) - pos - 1);
        if (const auto colon = line.find(':'); colon != std::string::npos) {
            const auto trim = [](std::string s) {
                const auto b = s.find_first_not_of(" \t\r"), e = s.find_last_not_of(" \t\r");
                return b == std::string::npos ? std::string{} : s.substr(b, e - b + 1);
            };
            d.header.emplace_back(trim(line.substr(0, colon)), trim(line.substr(colon + 1)));
        }
        if (end == std::string::npos) { pos = raw.size(); break; }
        pos = end + 1;
    }
    const std::string_view body(raw.data() + pos, raw.size() - pos);
    const std::string* fmt = d.get("format");
    if (!fmt || *fmt != "1") throw std::runtime_error(std::format("'{}': format must be 1", path));
    if (const std::string* cols = d.get("columns"); cols && *cols != "x y")
        throw std::runtime_error(std::format("'{}': columns '{}' not supported (x y)", path, *cols));
    const std::string* sha = d.get("sha256_body");
    if (!sha || sha256_hex(body) != *sha) throw std::runtime_error(std::format("'{}': sha256_body does not match the body", path));
    std::size_t p = 0;
    while (p < body.size()) {
        const std::size_t e = body.find('\n', p);
        const std::string_view line = body.substr(p, (e == std::string_view::npos ? body.size() : e) - p);
        p = e == std::string_view::npos ? body.size() : e + 1;
        if (line.empty()) continue;
        const auto sp = line.find(' ');
        if (sp == std::string_view::npos || line.find(' ', sp + 1) != std::string_view::npos)
            throw std::runtime_error(std::format("'{}': body line '{}' must have two columns", path, line));
        d.points.push_back({parse_real(std::string(line.substr(0, sp))), parse_real(std::string(line.substr(sp + 1)))});
    }
    if (d.points.size() != static_cast<std::size_t>(d.number("n_obstacles")))
        throw std::runtime_error(std::format("'{}': {} records, n_obstacles = {}", path, d.points.size(), *d.get("n_obstacles")));
    return d;
}

void write_sk_companion(const std::string& path, const std::vector<std::pair<std::string, std::string>>& header,
                        const std::vector<double>& k, const std::vector<double>& S) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error(std::format("could not create '{}'", path));
    std::string out;
    for (const auto& [key, value] : header) out += "# " + key + ": " + value + '\n';
    for (std::size_t i = 0; i < k.size() && i < S.size(); ++i) out += format_real(k[i]) + ' ' + format_real(S[i]) + '\n';
    f << out;
}

std::string obs_auto_name(const std::string& cls, double rho_O, double Lx, double Ly, int arr_index) {
    return std::format("{}_rho{:.6g}_Lx{:.6g}_Ly{:.6g}_arr{:03}.obs", cls, rho_O, Lx, Ly, arr_index);
}

std::string path_with_arr(const std::string& path, int arr, bool several) {
    if (!several) return path;
    const auto slash = path.find_last_of('/');
    const auto dot = path.find_last_of('.');
    const std::string suffix = std::format("_arr{:03}", arr);
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + suffix;
    return path.substr(0, dot) + suffix + path.substr(dot);
}

std::string path_with_extension(const std::string& path, const std::string& ext) {
    const auto slash = path.find_last_of('/');
    const auto dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + ext;
    return path.substr(0, dot) + ext;
}

std::string iso_date_today() {
    const auto today = std::chrono::floor<std::chrono::days>(std::chrono::system_clock::now());
    const std::chrono::year_month_day ymd{today};
    return std::format("{:04}-{:02}-{:02}", static_cast<int>(ymd.year()), static_cast<unsigned>(ymd.month()),
                       static_cast<unsigned>(ymd.day()));
}
