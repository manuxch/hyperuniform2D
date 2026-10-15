#!/usr/bin/env python3
"""Hyperuniformity analysis of a binary pattern stored as a plain PGM (P2) image.

For the occupied sites (pixel == max value, i.e. 1 in the files written by hyper2D) it computes:

  1. the number variance sigma_N^2(R) of the count in a disk of radius R, over ALL N^2 lattice centres
     (exact: a periodic FFT convolution of the occupation with the disk, no random sampling);
  2. the radially averaged structure factor S(k) = |rho_hat(k)|^2 / M, k in units of 2 pi / N.

Poisson-like (uncorrelated) patterns: sigma_N^2 ~ R^2.  Hyperuniform: slower than R^2 (stealthy: ~ R).
The exponent is fitted on [--rmin, --rmax]; the default rmax = N/4 keeps the disk well below the box, where
the fixed total number of particles would suppress the variance artificially.
sigma_rho^2 = sigma_N^2 / (pi R^2)^2, so its exponent is exactly alpha - 4: it is shown, but it is not
an independent test.

Usage:
  python3 test_hyper2D.py image.pgm [--rmin 5] [--rmax N/4] [--save fig.pdf] [--no-show]
  python3 test_hyper2D.py image.pgm --method sample --samples 5000   # old random-window estimator
"""
import argparse
import shutil
import sys

import numpy as np


def read_pgm_p2(filename):
    """Reads a PGM file in text mode (P2) ignoring comments."""
    words = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                words.extend(line.split('#')[0].split())
    except FileNotFoundError:
        sys.exit(f"Error: the file '{filename}' does not exist.")
    if not words:
        raise ValueError("The file is empty.")
    if words[0] != 'P2':
        raise ValueError("The file is not in PGM P2 format (text mode).")
    width, height, max_val = int(words[1]), int(words[2]), int(words[3])
    pixels = np.array(words[4:4 + width * height], dtype=np.int32)
    if pixels.size != width * height:
        raise ValueError(f"Expected {width * height} pixels, found {pixels.size}.")
    return pixels.reshape((height, width)), width, height, max_val


def disk_kernel(L, R):
    """Indicator of the lattice sites within distance R of the origin, on a periodic L x L grid."""
    d = np.minimum(np.arange(L), L - np.arange(L))
    return ((d[:, None] ** 2 + d[None, :] ** 2) <= R * R).astype(float)


def number_variance_fft(occ, radii):
    """Exact variance of the disk counts over all lattice centres (periodic)."""
    L = occ.shape[0]
    F = np.fft.rfft2(occ.astype(float))
    var = []
    for R in radii:
        counts = np.fft.irfft2(F * np.fft.rfft2(disk_kernel(L, R)), s=occ.shape)
        var.append(np.var(np.rint(counts)))
    return np.array(var)


def number_variance_sample(points, L, radii, samples, rng):
    """Old estimator: random continuous window centres, periodic minimum-image distances."""
    var = []
    for R in radii:
        cx, cy = rng.uniform(0, L, samples), rng.uniform(0, L, samples)
        counts = np.empty(samples)
        for s in range(samples):
            dx = np.abs(points[:, 0] - cx[s]); dx = np.minimum(dx, L - dx)
            dy = np.abs(points[:, 1] - cy[s]); dy = np.minimum(dy, L - dy)
            counts[s] = np.count_nonzero(dx * dx + dy * dy <= R * R)
        var.append(np.var(counts))
    return np.array(var)


def structure_factor(occ):
    """Radially averaged S(k) = |rho_hat|^2 / M, bins of width 1 in |k| (units of 2 pi / L), k != 0."""
    L = occ.shape[0]
    S = np.abs(np.fft.fft2(occ.astype(float))) ** 2 / occ.sum()
    q = np.fft.fftfreq(L, 1.0 / L)
    k = np.sqrt(q[:, None] ** 2 + q[None, :] ** 2)
    b = np.rint(k).astype(int)
    mask = k > 0
    nb = b.max() + 1
    cnt = np.bincount(b[mask], minlength=nb)
    ssum = np.bincount(b[mask], weights=S[mask], minlength=nb)
    ksum = np.bincount(b[mask], weights=k[mask], minlength=nb)
    ok = cnt > 0
    return ksum[ok] / cnt[ok], ssum[ok] / cnt[ok]


def main():
    parser = argparse.ArgumentParser(description="Analyzes hyperuniformity: number variance exponent and S(k).")
    parser.add_argument("filename", help="Path to the PGM file in text mode (P2)")
    parser.add_argument("--rmin", type=float, default=5.0, help="smallest radius of the fit (default 5)")
    parser.add_argument("--rmax", type=float, default=None, help="largest radius of the fit (default N/4)")
    parser.add_argument("--n-radii", type=int, default=30, help="number of radii (default 30)")
    parser.add_argument("--method", choices=("fft", "sample"), default="fft",
                        help="fft: exact over all lattice centres (default); sample: random windows")
    parser.add_argument("--samples", type=int, default=1000, help="windows per radius for --method sample")
    parser.add_argument("--seed", type=int, default=12345, help="random seed for --method sample")
    parser.add_argument("--save", help="save the figure to this file")
    parser.add_argument("--no-show", action="store_true", help="do not open the interactive window")
    args = parser.parse_args()

    print(f"Reading file: {args.filename}...")
    pixels, width, height, max_val = read_pgm_p2(args.filename)
    if width != height:
        sys.exit("Error: the analysis assumes a square periodic lattice.")
    L = width
    occ = pixels == max_val  # occupied sites (1 in hyper2D output), whatever their fraction
    M = int(occ.sum())
    rho = M / (L * L)
    print(f"Image: {L}x{L}, {M} occupied sites, density {rho:.4f}")

    rmax = args.rmax if args.rmax is not None else L / 4.0
    if not (0 < args.rmin < rmax <= L / 2):
        sys.exit(f"Error: need 0 < rmin < rmax <= N/2 (got {args.rmin}, {rmax}).")
    radii = np.linspace(args.rmin, rmax, args.n_radii)

    print(f"Number variance ({args.method}) for R in [{args.rmin:g}, {rmax:g}]...")
    if args.method == "fft":
        var_N = number_variance_fft(occ, radii)
    else:
        pts = np.argwhere(occ)[:, [1, 0]].astype(float)
        var_N = number_variance_sample(pts, L, radii, args.samples, np.random.default_rng(args.seed))
    if np.any(var_N <= 0):
        sys.exit("Error: zero variance at some radius (perfectly regular pattern?): no log-log fit possible.")
    var_rho = var_N / (np.pi * radii ** 2) ** 2

    log_R = np.log(radii)
    alpha, c_N = np.polyfit(log_R, np.log(var_N), 1)
    beta = alpha - 4.0
    k, S = structure_factor(occ)
    low = k <= max(2.0, 0.05 * L)

    print("\n" + "=" * 64)
    print("                 EXPERIMENTAL FIT RESULTS")
    print("=" * 64)
    print(f"-> Exponent of var_N (alpha)   : {alpha:.3f}")
    print("   [Poisson = 2 | hyperuniform < 2 (stealthy ~ 1)]")
    print(f"-> Exponent of var_rho (beta)  : {beta:.3f}  (= alpha - 4 by definition)")
    print(f"-> <S(k)> for |k| <= {k[low].max():.1f} (2pi/N) : {S[low].mean():.4g}"
          f"   [uncorrelated lattice gas: 1 - rho = {1 - rho:.4f}]")
    print("=" * 64 + "\n")

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'lines.linewidth': 1.5, 'font.size': 10,
                                'text.usetex': shutil.which("latex") is not None})

    fit_N = np.exp(c_N) * radii ** alpha
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    ax = axs[0]
    ax.loglog(radii, var_N, 'o', color='tab:blue', label='measured')
    ax.loglog(radii, fit_N, '-', color='tab:blue', label=fr'fit, $\alpha = {alpha:.2f}$')
    ax.loglog(radii, radii ** 2 * var_N[0] / radii[0] ** 2, '--', color='gray', label=r'$R^2$ (Poisson)')
    ax.loglog(radii, radii * var_N[0] / radii[0], ':', color='tab:red', label=r'$R$ (stealthy)')
    ax.set_xlabel(r'$R$'); ax.set_ylabel(r'$\sigma_N^2$'); ax.legend(); ax.grid(True, which="both", ls="--", alpha=0.4)

    ax = axs[1]
    ax.loglog(radii, var_rho, 's', color='tab:purple', label='measured')
    ax.loglog(radii, fit_N / (np.pi * radii ** 2) ** 2, '-', color='tab:purple', label=fr'$\beta = \alpha - 4 = {beta:.2f}$')
    ax.loglog(radii, radii ** -2 * var_rho[0] / radii[0] ** -2, '--', color='gray', label=r'$R^{-2}$ (Poisson)')
    ax.loglog(radii, radii ** -3 * var_rho[0] / radii[0] ** -3, ':', color='tab:red', label=r'$R^{-3}$ (stealthy)')
    ax.set_xlabel(r'$R$'); ax.set_ylabel(r'$\sigma_\rho^2$'); ax.legend(); ax.grid(True, which="both", ls="--", alpha=0.4)

    ax = axs[2]
    ax.plot(k, S, '.-', color='tab:green')
    ax.axhline(1 - rho, color='gray', ls='--', label=r'$1-\rho$ (uncorrelated)')
    ax.set_xlabel(r'$|k|$ ($2\pi/N$)'); ax.set_ylabel(r'$S(k)$'); ax.set_xlim(0, L / 2)
    ax.legend(); ax.grid(True, ls="--", alpha=0.4)

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save)
        print(f"Figure saved to {args.save}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
