#!/usr/bin/env python3
"""Hyperuniformity analysis of an obstacle file (format of hyperGalton2D, version 1), lattice or continuum,
in any box Lx x Ly. Lengths in grain diameters; wave vectors in d^-1.

  1. rho_hat(k) on every wave vector k = 2 pi (nx/Lx, ny/Ly) with |k| <= k_max, as a matrix product of the
     axis phase tables, rho_hat(nx, ny) = sum_j ex[j, nx] ey[j, ny];
  2. S(k) = |rho_hat|^2 / N: radial average, and the largest S inside the mask |k| <= K by angular sector
     (isotropy: the mask is a circle in physical k);
  3. the number variance sigma_N^2(R) of the count in a disk of radius R:
       - spectral, exact for the pattern on the torus (both axes periodic, as hyperuniform2D generates it):
           sigma^2 = (N / A^2) sum_{k != 0} S(k) m(k)^2,   m(k) = 2 pi R J1(kR) / k,
         with the modes beyond k_max closed by Parseval, sum_{k != 0} m^2 = A pi R^2 - (pi R^2)^2, taking
         S = 1 there (its value for |k| -> infinity);
       - direct count over random window centres, as an independent check (and the only estimate when the
         pattern is not periodic in y);
     and the exponent alpha of sigma_N^2 ~ R^alpha (Poisson 2, stealthy hyperuniform 1);
  4. the pair correlation g(r).

  python3 scripts/analyze_obs.py file.obs [--save fig.pdf] [--no-show] [--write-alpha]
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obs_tools  # noqa: E402


def axis_tables(u, L, n_max, negative):
    """exp(-i 2 pi n u / L) for n = 0..n_max (or -n_max..n_max if negative), rows = points."""
    n = np.arange(-n_max, n_max + 1) if negative else np.arange(0, n_max + 1)
    return np.exp(-2j * np.pi * np.outer(u, n) / L), n


def rho_hat_grid(xy, Lx, Ly, k_max, chunk=4096):
    """rho_hat on the half plane nx >= 0 of the rectangle |nx| <= k_max Lx / 2 pi, |ny| <= k_max Ly / 2 pi.
    Returns rho (nx, ny), and the physical kx, ky grids."""
    nxm, nym = int(k_max * Lx / (2 * np.pi)) + 1, int(k_max * Ly / (2 * np.pi)) + 1
    rho = np.zeros((nxm + 1, 2 * nym + 1), dtype=complex)
    for b in range(0, len(xy), chunk):  # chunks keep the tables small
        ex, nx = axis_tables(xy[b:b + chunk, 0], Lx, nxm, False)
        ey, ny = axis_tables(xy[b:b + chunk, 1], Ly, nym, True)
        rho += ex.T @ ey
    kx = 2 * np.pi * nx[:, None] / Lx * np.ones_like(ny)[None, :]
    ky = 2 * np.pi * ny[None, :] / Ly * np.ones_like(nx)[:, None]
    return rho, kx, ky, nx, ny


def disk_ft(k, R):
    """Fourier transform of the indicator of a disk of radius R: 2 pi R J1(kR) / k (pi R^2 at k = 0)."""
    from scipy.special import j1
    k = np.asarray(k, dtype=float)
    out = np.full_like(k, np.pi * R * R)
    nz = k > 0
    out[nz] = 2 * np.pi * R * j1(k[nz] * R) / k[nz]
    return out


def number_variance_spectral(S_half, k_half, N, A, radii):
    """sigma_N^2(R) on the torus from S on the half plane (k and -k give the same S), with the Parseval tail."""
    var = []
    for R in radii:
        m2 = disk_ft(k_half, R) ** 2
        inside = 2.0 * np.sum(m2)                      # both members of each +-k pair
        total = A * np.pi * R * R - (np.pi * R * R) ** 2  # sum over every k != 0
        var.append(N / A ** 2 * (2.0 * np.sum(S_half * m2) + (total - inside)))
    return np.array(var)


def number_variance_count(xy, Lx, Ly, radii, samples, rng, periodic_y):
    """Variance of the counts in disks with random centres; x periodic, y periodic or with the centres kept
    at least R from the walls."""
    from scipy.spatial import cKDTree
    if periodic_y:
        tree = cKDTree(xy, boxsize=[Lx, Ly])
    else:
        tree = cKDTree(xy, boxsize=[Lx, 4.0 * Ly + 1.0])  # y effectively open
    var = []
    for R in radii:
        y0, y1 = (0.0, Ly) if periodic_y else (R, Ly - R)
        c = np.c_[rng.uniform(0, Lx, samples), rng.uniform(y0, y1, samples)]
        var.append(np.var(tree.query_ball_point(c, R, return_length=True)))
    return np.array(var)


def pair_correlation(xy, Lx, Ly, r_max, bins, periodic_y):
    from scipy.spatial import cKDTree
    tree = cKDTree(xy, boxsize=[Lx, Ly] if periodic_y else [Lx, 4.0 * Ly + 1.0])
    d = tree.sparse_distance_matrix(tree, r_max, output_type="ndarray")["v"]
    d = d[d > 0]
    edges = np.linspace(0, r_max, bins + 1)
    h, _ = np.histogram(d, bins=edges)  # ordered pairs
    r = 0.5 * (edges[1:] + edges[:-1])
    N = len(xy)
    ideal = N * (N - 1) / (Lx * Ly) * np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    return r, h / ideal


def write_alpha(path, alpha):
    """Adds or replaces alpha_numvar in the header; the body (and sha256_body) is not touched."""
    o = obs_tools.ObsFile(path)
    lines = o.raw[: len(o.raw) - len(o.body)].decode("utf-8").splitlines()
    new = f"# alpha_numvar: {alpha:.6g}"
    for i, line in enumerate(lines):
        if line.startswith("# alpha_numvar:"):
            lines[i] = new
            break
    else:
        at = next((i for i, l in enumerate(lines) if l.startswith("# columns:")), len(lines))
        lines.insert(at, new)
    with open(path, "wb") as f:
        f.write(("\n".join(lines) + "\n").encode("utf-8") + o.body)


def analyze(path, args):
    o = obs_tools.ObsFile(path)
    h = o.header
    Lx, Ly, xy = float(h["Lx"]), float(h["Ly"]), o.xy
    N, A = len(xy), float(h["Lx"]) * float(h["Ly"])
    rho = N / A
    a = rho ** -0.5
    if args.periodic == "auto":
        periodic_y = h.get("generator") == "hyperuniform2D"  # generated on the torus
    else:
        periodic_y = args.periodic == "xy"
    K = float(h["K"]) if "K" in h else None
    print(f"{path}: {N} obstacles in {Lx:g} x {Ly:g}, rho_O = {rho:.6g}, mean spacing a = {a:.6g}; "
          f"mode {h.get('mode')}, class {h.get('class')}, periodic {'x and y' if periodic_y else 'x only'}")

    # S(k)
    k_max = args.kmax if args.kmax else (4.0 * K if K else 8.0 * np.pi / a)
    rho_k, kx, ky, nx, ny = rho_hat_grid(xy, Lx, Ly, k_max)
    k = np.hypot(kx, ky)
    half = (nx[:, None] > 0) | ((nx[:, None] == 0) & (ny[None, :] > 0))
    sel = half & (k <= k_max)
    S = np.abs(rho_k[sel]) ** 2 / N
    ks, kxs, kys = k[sel], kx[sel], ky[sel]
    if "lattice_constant" in h and h.get("mode") == "lattice":
        k_nyq = np.pi / float(h["lattice_constant"])
        if k_max > k_nyq:
            print(f"  note: k_max = {k_max:.4g} exceeds the lattice Nyquist limit pi/a_lattice = {k_nyq:.4g}: "
                  "S(k) aliases above it")
    report = {}
    if K:
        inside = ks <= K
        theta = np.degrees(np.arctan2(np.abs(kys), np.abs(kxs)))
        mean_in = S[inside].mean()
        print(f"  mask |k| <= K = {K:.6g}: {inside.sum()} half-plane modes (header: {h.get('n_constrained_modes')}), "
              f"<S> = {mean_in:.3e}, max S = {S[inside].max():.3e}")
        for lo, hi in ((0, 30), (30, 60), (60, 90)):
            m = inside & (theta >= lo) & (theta < hi)
            if m.any():
                print(f"    sector {lo:2d}-{hi:2d} deg: {m.sum():6d} modes, <S> = {S[m].mean():.3e}, max S = {S[m].max():.3e}")
        report["S_mask_mean"] = mean_in

    dk = 2 * np.pi / math.sqrt(A)
    b = np.rint(ks / dk).astype(int)
    cnt = np.bincount(b)
    with np.errstate(invalid="ignore"):
        k_rad = np.bincount(b, weights=ks) / cnt
        S_rad = np.bincount(b, weights=S) / cnt
    ok = cnt > 0
    k_rad, S_rad = k_rad[ok], S_rad[ok]

    # Number variance
    rmin = args.rmin if args.rmin else a
    rmax = args.rmax if args.rmax else min(Lx, Ly) / 4.0
    if not (0 < rmin < rmax <= min(Lx, Ly) / 2):
        sys.exit(f"Error: need 0 < rmin < rmax <= min(Lx, Ly)/2 (got {rmin}, {rmax})")
    radii = np.linspace(rmin, rmax, args.n_radii)
    rng = np.random.default_rng(args.seed)
    var_count = number_variance_count(xy, Lx, Ly, radii, args.samples, rng, periodic_y)
    var_spec = None
    if periodic_y:
        var_spec = number_variance_spectral(S, ks, N, A, radii)
        alpha = np.polyfit(np.log(radii), np.log(var_spec), 1)[0]
        dev = np.max(np.abs(var_count / var_spec - 1))
        print(f"  number variance on R in [{rmin:.4g}, {rmax:.4g}]: alpha = {alpha:.3f} (spectral, exact); "
              f"direct count ({args.samples} windows) agrees within {100 * dev:.1f} %")
    else:
        alpha = np.polyfit(np.log(radii), np.log(var_count), 1)[0]
        print(f"  number variance on R in [{rmin:.4g}, {rmax:.4g}]: alpha = {alpha:.3f} (direct count, "
              f"{args.samples} windows away from the y walls)")
    alpha_count = np.polyfit(np.log(radii), np.log(var_count), 1)[0]
    print("    [Poisson 2 | hyperuniform < 2 | stealthy ~ 1]")
    report.update(alpha=alpha, alpha_count=alpha_count, radii=radii, var_spec=var_spec, var_count=var_count)

    r, g = pair_correlation(xy, Lx, Ly, args.gr_rmax * a, args.gr_bins, periodic_y)
    nn = np.nonzero(g > 0)[0]
    print(f"  g(r): first nonzero bin at r = {r[nn[0]] / a:.3f} a" if len(nn) else "  g(r): empty")
    report.update(k=k_rad, S=S_rad, r=r, g=g, a=a, K=K)

    if args.write_alpha:
        write_alpha(path, alpha)
        errors, _ = obs_tools.validate(path)
        print(f"  alpha_numvar = {alpha:.6g} written to the header; file {'still valid' if not errors else 'INVALID: ' + errors[0]}")

    if not args.no_show or args.save:
        plot(o, report, args)
    return report


def plot(o, rep, args):
    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Lx, Ly = float(o.header["Lx"]), float(o.header["Ly"])
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(2, 1, 1)
    ax.scatter(o.xy[:, 0], o.xy[:, 1], s=1.0, c="k")
    ax.set_xlim(0, Lx); ax.set_ylim(0, Ly); ax.set_aspect("equal"); ax.set_title(os.path.basename(o.path))
    ax = fig.add_subplot(2, 3, 4)
    radii = rep["radii"]
    if rep["var_spec"] is not None:
        ax.loglog(radii, rep["var_spec"], "-", label=f"spectral, alpha = {rep['alpha']:.2f}")
    ax.loglog(radii, rep["var_count"], "o", ms=3, label=f"direct count, alpha = {rep['alpha_count']:.2f}")
    v0 = rep["var_count"][0]
    ax.loglog(radii, v0 * (radii / radii[0]) ** 2, "--", c="gray", label="R^2 (Poisson)")
    ax.loglog(radii, v0 * radii / radii[0], ":", c="tab:red", label="R (stealthy)")
    ax.set_xlabel("R (d)"); ax.set_ylabel("sigma_N^2"); ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
    ax = fig.add_subplot(2, 3, 5)
    ax.plot(rep["k"], rep["S"], ".-", ms=3)
    if rep["K"]:
        ax.axvline(rep["K"], c="tab:red", ls=":", label="K")
    ax.axhline(1.0, c="gray", ls="--", label="Poisson")
    ax.set_xlabel("|k| (1/d)"); ax.set_ylabel("S(k)"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax = fig.add_subplot(2, 3, 6)
    ax.plot(rep["r"] / rep["a"], rep["g"], "-")
    ax.axhline(1.0, c="gray", ls="--")
    ax.set_xlabel("r / a"); ax.set_ylabel("g(r)"); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save)
        print(f"  figure saved to {args.save}")
    if not args.no_show:
        plt.show()


def build_parser():
    p = argparse.ArgumentParser(description="Hyperuniformity analysis of an obstacle file (.obs, format version 1).")
    p.add_argument("filename")
    p.add_argument("--kmax", type=float, default=None, help="largest |k| (1/d) for S(k) (default 4K, or 8 pi/a)")
    p.add_argument("--rmin", type=float, default=None, help="smallest R of the fit, in d (default a)")
    p.add_argument("--rmax", type=float, default=None, help="largest R of the fit, in d (default min(Lx, Ly)/4)")
    p.add_argument("--n-radii", type=int, default=24)
    p.add_argument("--samples", type=int, default=20000, help="window centres per radius for the direct count")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--periodic", choices=("auto", "xy", "x"), default="auto",
                   help="auto: x and y for files of hyperuniform2D (generated on the torus), else x only")
    p.add_argument("--gr-rmax", type=float, default=4.0, help="g(r) up to this many mean spacings")
    p.add_argument("--gr-bins", type=int, default=120)
    p.add_argument("--write-alpha", action="store_true", help="write alpha_numvar into the header (body untouched)")
    p.add_argument("--save", help="save the figure to this file")
    p.add_argument("--no-show", action="store_true", help="do not open the interactive window")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    analyze(args.filename, args)


if __name__ == "__main__":
    main()
