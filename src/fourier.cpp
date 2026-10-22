#include "fourier.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <format>
#include <fstream>
#include <iostream>

void calculate_FT(HyperuniformConfig& config) {
    const int N = config.Nx;
    const std::size_t num_modes = config.modes.size();
    std::fill(config.ft_re.begin(), config.ft_re.end(), 0.0);
    std::fill(config.ft_im.begin(), config.ft_im.end(), 0.0);

    for (const std::uint32_t s : config.occupied) {
        const int x = static_cast<int>(s % static_cast<std::uint32_t>(N));
        const int y = static_cast<int>(s / static_cast<std::uint32_t>(N));
        for (std::size_t m = 0; m < num_modes; ++m) {
            const auto j = phase_index(config, m, x, y);
            config.ft_re[m] += config.w_re[j];
            config.ft_im[m] += config.w_im[j];
        }
    }
}

double energy(const HyperuniformConfig& config) {
    double E = 0.0;
    for (std::size_t m = 0; m < config.ft_re.size(); ++m)
        E += config.ft_re[m] * config.ft_re[m] + config.ft_im[m] * config.ft_im[m];
    return 2.0 * E;  // each stored mode stands for the pair ±k
}

double weighted_energy(const HyperuniformConfig& config) {
    double E = 0.0;
    for (std::size_t m = 0; m < config.ft_re.size(); ++m)
        E += config.weight[m] * (config.ft_re[m] * config.ft_re[m] + config.ft_im[m] * config.ft_im[m]);
    return 2.0 * E;
}

double band_residual(const HyperuniformConfig& config, double k_max) {
    const double rho = config.M / (static_cast<double>(config.Nx) * config.Ny);
    const double k2max = k_max * k_max;
    double sum = 0.0;
    std::size_t n = 0;
    for (std::size_t m = 0; m < config.modes.size(); ++m) {
        if (config.mode_k2[m] > k2max) continue;
        sum += config.ft_re[m] * config.ft_re[m] + config.ft_im[m] * config.ft_im[m];
        ++n;
    }
    return n ? sum / static_cast<double>(n) / config.M / (1.0 - rho) : 0.0;
}

double mask_residual(const HyperuniformConfig& config, double E) {
    const double sites = static_cast<double>(config.Nx) * config.Ny;
    const double rho = config.M / sites;
    const double full_modes = 2.0 * static_cast<double>(config.modes.size());
    return E / (full_modes * config.M) / (1.0 - rho);
}

double delta_energy_range(HyperuniformConfig& config, int xi, int yi, int xj, int yj,
                          std::size_t begin, std::size_t end) {
    double dE = 0.0;
    for (std::size_t m = begin; m < end; ++m) {
        const auto ji = phase_index(config, m, xi, yi);
        const auto jj = phase_index(config, m, xj, yj);
        const double dr = config.w_re[jj] - config.w_re[ji];
        const double di = config.w_im[jj] - config.w_im[ji];
        config.d_re[m] = dr;
        config.d_im[m] = di;
        dE += config.weight[m] * (2.0 * (config.ft_re[m] * dr + config.ft_im[m] * di) + dr * dr + di * di);
    }
    return 2.0 * dE;  // ±k pair
}

double energy_direct(const HyperuniformConfig& config) {
    // Independent of the tables and of the stored mode list: the full circle |k| <= K in physical k,
    // with the Nyquist limit of each axis, by direct trigonometry.
    const int Nx = config.Nx, Ny = config.Ny;
    const int limx = Nx / 2 - 1, limy = Ny / 2 - 1;
    double E = 0.0;
    for (int kx = -limx; kx <= limx; ++kx)
        for (int ky = -limy; ky <= limy; ++ky) {
            if ((kx == 0 && ky == 0) || mode_key(kx, ky, Nx, Ny) > config.q_max) continue;
            Complex sum{0.0, 0.0};
            for (const std::uint32_t s : config.occupied) {
                const int x = static_cast<int>(s % static_cast<std::uint32_t>(Nx));
                const int y = static_cast<int>(s / static_cast<std::uint32_t>(Nx));
                sum += std::polar(1.0, -2.0 * PI * (static_cast<double>(kx) * x / Nx + static_cast<double>(ky) * y / Ny));
            }
            E += std::norm(sum);
        }
    return E;
}

RadialSk radial_structure_factor(const HyperuniformConfig& config) {
    const int Nx = config.Nx, Ny = config.Ny;
    const auto nx = static_cast<std::size_t>(Nx), ny = static_cast<std::size_t>(Ny);

    // Separable DFT with tables of roots of unity, O(Nx Ny (Nx + Ny)): first along x for every row y,
    // then along y.
    const auto roots = [](int N) {
        std::vector<Complex> w(static_cast<std::size_t>(N));
        for (std::size_t j = 0; j < w.size(); ++j) w[j] = std::polar(1.0, -2.0 * PI * static_cast<double>(j) / N);
        return w;
    };
    const std::vector<Complex> wx = roots(Nx), wy = Nx == Ny ? wx : roots(Ny);
    std::vector<Complex> rows(nx * ny, Complex{0.0, 0.0});  // rows[kx + Nx y] = Σ_x g(x, y) wx^(kx x)
    for (int y = 0; y < Ny; ++y)
        for (int x = 0; x < Nx; ++x) {
            if (config.grid[static_cast<std::size_t>(flat_idx(x, y, Nx))] == 0) continue;
            for (int kx = 0; kx < Nx; ++kx)
                rows[static_cast<std::size_t>(kx + Nx * y)] += wx[static_cast<std::size_t>((kx * x) % Nx)];
        }

    // Radial bins of width 1 centred on integers, |k| in units of 2π/sqrt(Nx Ny) (2π/N if square), from the
    // exact key q; k folded into (-N/2, N/2] on each axis.
    const double area = static_cast<double>(Nx) * static_cast<double>(Ny);
    const double k_far = std::sqrt(static_cast<double>(mode_key(Nx / 2, Ny / 2, Nx, Ny)) / area);
    const int n_bins = static_cast<int>(std::ceil(k_far)) + 2;
    std::vector<double> sum_s(static_cast<std::size_t>(n_bins), 0.0), sum_k(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<long long> count(static_cast<std::size_t>(n_bins), 0);
    for (int kx = 0; kx < Nx; ++kx)
        for (int ky = 0; ky < Ny; ++ky) {
            if (kx == 0 && ky == 0) continue;
            Complex f{0.0, 0.0};
            for (int y = 0; y < Ny; ++y)
                f += rows[static_cast<std::size_t>(kx + Nx * y)] * wy[static_cast<std::size_t>((ky * y) % Ny)];
            const int qx = kx > Nx / 2 ? kx - Nx : kx;
            const int qy = ky > Ny / 2 ? ky - Ny : ky;
            const double k = std::sqrt(static_cast<double>(mode_key(qx, qy, Nx, Ny)) / area);
            const auto b = static_cast<std::size_t>(std::lround(k));
            sum_s[b] += std::norm(f) / config.M;
            sum_k[b] += k;
            ++count[b];
        }

    RadialSk out;
    for (int b = 1; b < n_bins; ++b) {
        const auto i = static_cast<std::size_t>(b);
        if (count[i] == 0) continue;
        out.k.push_back(sum_k[i] / static_cast<double>(count[i]));
        out.S.push_back(sum_s[i] / static_cast<double>(count[i]));
        out.count.push_back(count[i]);
    }
    return out;
}

void analyze_structure_factor(const HyperuniformConfig& config, const std::string& filename) {
    analyze_structure_factor(config, filename, radial_structure_factor(config));
}

void analyze_structure_factor(const HyperuniformConfig& config, const std::string& filename, const RadialSk& sk) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << std::format("Error: could not create '{}'\n", filename);
        return;
    }
    const int Nx = config.Nx, Ny = config.Ny;
    const double area = static_cast<double>(Nx) * static_cast<double>(Ny);
    const double rho = static_cast<double>(config.M) / area;
    fout << std::format("# Radially averaged structure factor S(k) = |rho_hat(k)|^2 / M, M = {} occupied sites\n",
                        config.M);
    if (config.square())
        fout << std::format("# N = {}, mask radius K = {:.6g} (units of 2pi/N); uncorrelated lattice gas: S = 1 - rho = {:.6g}\n",
                            Nx, config.K, 1.0 - rho);
    else
        fout << std::format("# Nx = {}, Ny = {}, mask radius K = {:.6g} (units of 2pi/sqrt(Nx Ny), lattice spacing 1); "
                            "uncorrelated lattice gas: S = 1 - rho = {:.6g}\n", Nx, Ny, config.K, 1.0 - rho);
    fout << std::format("{:<12} {:<20} {:<8}\n", "# k", "S(k)", "modes");
    for (std::size_t i = 0; i < sk.k.size(); ++i)
        fout << std::format("{:<12.6f} {:<20.10e} {:<8d}\n", sk.k[i], sk.S[i], sk.count[i]);
    std::cout << std::format("Structure factor saved to '{}'\n", filename);
}
