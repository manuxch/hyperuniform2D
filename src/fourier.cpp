#include "fourier.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <format>
#include <fstream>
#include <iostream>

void calculate_FT(HyperuniformConfig& config) {
    const int N = config.N;
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
    const double rho = config.M / (static_cast<double>(config.N) * config.N);
    const double k2max = k_max * k_max;
    double sum = 0.0;
    std::size_t n = 0;
    for (std::size_t m = 0; m < config.modes.size(); ++m) {
        const auto [kx, ky] = config.modes[m];
        if (static_cast<double>(kx * kx + ky * ky) > k2max) continue;
        sum += config.ft_re[m] * config.ft_re[m] + config.ft_im[m] * config.ft_im[m];
        ++n;
    }
    return n ? sum / static_cast<double>(n) / config.M / (1.0 - rho) : 0.0;
}

double mask_residual(const HyperuniformConfig& config, double E) {
    const double sites = static_cast<double>(config.N) * config.N;
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
    const int N = config.N;
    const int lim = N / 2 - 1;
    double E = 0.0;
    for (int kx = -lim; kx <= lim; ++kx)
        for (int ky = -lim; ky <= lim; ++ky) {
            if ((kx == 0 && ky == 0) || static_cast<double>(kx * kx + ky * ky) > config.K2) continue;
            Complex sum{0.0, 0.0};
            for (const std::uint32_t s : config.occupied) {
                const int x = static_cast<int>(s % static_cast<std::uint32_t>(N));
                const int y = static_cast<int>(s / static_cast<std::uint32_t>(N));
                sum += std::polar(1.0, -2.0 * PI * (kx * x + ky * y) / N);
            }
            E += std::norm(sum);
        }
    return E;
}

void analyze_structure_factor(const HyperuniformConfig& config, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << std::format("Error: could not create '{}'\n", filename);
        return;
    }
    const int N = config.N;
    const auto n = static_cast<std::size_t>(N);

    // Separable DFT with a table of roots of unity, O(N³): first along x for every row y, then along y.
    std::vector<Complex> w(n);
    for (std::size_t j = 0; j < n; ++j) w[j] = std::polar(1.0, -2.0 * PI * static_cast<double>(j) / N);
    std::vector<Complex> rows(n * n, Complex{0.0, 0.0});  // rows[kx + N y] = Σ_x g(x, y) w^(kx x)
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            if (config.grid[static_cast<std::size_t>(flat_idx(x, y, N))] == 0) continue;
            for (int kx = 0; kx < N; ++kx)
                rows[static_cast<std::size_t>(kx + N * y)] += w[static_cast<std::size_t>((kx * x) % N)];
        }

    // Radial bins of width 1 centred on integers; k folded into (-N/2, N/2] on each axis.
    const int n_bins = static_cast<int>(std::ceil(std::sqrt(2.0) * N / 2.0)) + 1;
    std::vector<double> sum_s(static_cast<std::size_t>(n_bins), 0.0), sum_k(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<long long> count(static_cast<std::size_t>(n_bins), 0);
    for (int kx = 0; kx < N; ++kx)
        for (int ky = 0; ky < N; ++ky) {
            if (kx == 0 && ky == 0) continue;
            Complex f{0.0, 0.0};
            for (int y = 0; y < N; ++y)
                f += rows[static_cast<std::size_t>(kx + N * y)] * w[static_cast<std::size_t>((ky * y) % N)];
            const int qx = kx > N / 2 ? kx - N : kx;
            const int qy = ky > N / 2 ? ky - N : ky;
            const double k = std::sqrt(static_cast<double>(qx * qx + qy * qy));
            const auto b = static_cast<std::size_t>(std::lround(k));
            sum_s[b] += std::norm(f) / config.M;
            sum_k[b] += k;
            ++count[b];
        }

    const double rho = static_cast<double>(config.M) / (static_cast<double>(N) * N);
    fout << std::format("# Radially averaged structure factor S(k) = |rho_hat(k)|^2 / M, M = {} occupied sites\n",
                        config.M);
    fout << std::format("# N = {}, mask radius K = {:.6g} (units of 2pi/N); uncorrelated lattice gas: S = 1 - rho = {:.6g}\n",
                        N, config.K, 1.0 - rho);
    fout << std::format("{:<12} {:<20} {:<8}\n", "# k", "S(k)", "modes");
    for (int b = 1; b < n_bins; ++b) {
        const auto i = static_cast<std::size_t>(b);
        if (count[i] == 0) continue;
        fout << std::format("{:<12.6f} {:<20.10e} {:<8d}\n",
                            sum_k[i] / static_cast<double>(count[i]), sum_s[i] / static_cast<double>(count[i]), count[i]);
    }
    std::cout << std::format("Structure factor saved to '{}'\n", filename);
}
