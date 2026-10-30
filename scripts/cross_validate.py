#!/usr/bin/env python3
"""Cross-validation of two families of obstacle files at equal K (kappa) and equal rho_O, e.g. continuum stealthy
patterns against their fine-lattice control. Each family is a set of realizations (seeds); every statistic is
averaged over them, with the standard error of the mean (SEM) as error bar, and the two families are compared
bin by bin through z = (mean_1 - mean_2) / sqrt(SEM_1^2 + SEM_2^2).

Statistics, in the units that pair the families (k / K and r / a, a = rho_O^(-1/2)):
  - S(k): radial average in bins of k / K; <S>_mask (point-pattern normalisation) and the header residual_S;
  - H = S(0) / S(k_peak) (Torquato 2018), S(0) from S0 + c k^2 on |k| <= K/8 (0 when the fit is negative), and
    the upper bound S(k < 0.05 K) / S(k_peak) from the lowest bin;
  - number variance sigma_N^2(R), spectral and exact on the torus (scripts/analyze_obs.py), its exponent alpha
    and the surface coefficient Lambda / phi^(1/2) = <sigma^2 / R> 2 / sqrt(pi rho);
  - nearest-neighbour distance: distribution function F(r_nn), mean, standard deviation, minimum;
  - the cumulative coordination Z(r) / (rho pi r^2), the integral of g(r) (histograms of r_nn and g(r) are not
    comparable bin by bin: lattice distances are discrete and give a comb).
Both families are taken as periodic in x and y (hyperuniform2D generates them on the torus).

  python3 scripts/cross_validate.py --a 'cont/*.obs' --b 'lat*/*.obs' --labels continuum lattice \\
          --out results/ [--title 'kappa = 0.618']
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_obs as ana  # noqa: E402
import obs_tools  # noqa: E402

KB = np.arange(0.0, 4.0 + 1e-9, 0.05)       # bins of k / K
RA = np.linspace(0.25, 8.5, 34)              # radii / a for sigma^2
# Distances on a lattice are discrete (sqrt(i^2 + j^2) lattice spacings), so histograms of r_nn and g(r) show a
# comb at any bin width comparable to the gaps between allowed distances; the comparison uses cumulative functions.
RNN = np.arange(0.05, 1.5 + 1e-9, 0.025)     # r_nn / a at which the distribution function F(r_nn) is evaluated
RZ = np.arange(0.1, 4.0 + 1e-9, 0.05)        # r / a for the cumulative coordination Z(r) / (rho pi r^2)


def one(path, check=True):
    o = obs_tools.ObsFile(path)
    h = o.header
    errors, _ = obs_tools.validate(path, quiet=True) if check else ([], [])
    if errors:
        sys.exit(f"{path}: invalid obstacle file: {errors[0]}")
    Lx, Ly, xy = float(h["Lx"]), float(h["Ly"]), o.xy
    N, A = len(xy), Lx * Ly
    rho = N / A
    a = rho ** -0.5
    K = float(h["K"])
    k_max = 4.0 * K
    rk, kx, ky, nx, ny = ana.rho_hat_grid(xy, Lx, Ly, k_max)
    k = np.hypot(kx, ky)
    half = (nx[:, None] > 0) | ((nx[:, None] == 0) & (ny[None, :] > 0))
    sel = half & (k <= k_max)
    S, ks = np.abs(rk[sel]) ** 2 / N, k[sel]
    x = ks / K
    b = np.digitize(x, KB) - 1
    ok = (b >= 0) & (b < len(KB) - 1)
    cnt = np.bincount(b[ok], minlength=len(KB) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        S_bin = np.bincount(b[ok], weights=S[ok], minlength=len(KB) - 1) / cnt
        k_bin = np.bincount(b[ok], weights=ks[ok], minlength=len(KB) - 1) / cnt
    good = cnt > 0
    H, S0, S0_se, S_peak, k_peak = ana.hyperuniformity_index(ks, S, k_bin[good], S_bin[good], K / 8.0)

    radii = RA * a
    var = ana.number_variance_spectral(S, ks, N, A, radii)
    fit = radii >= a
    alpha = np.polyfit(np.log(radii[fit]), np.log(var[fit]), 1)[0]
    lam_phi, _ = ana.surface_coefficient(radii[fit], var[fit], rho, a)

    tree = cKDTree(np.mod(xy, [Lx, Ly]), boxsize=[Lx, Ly])
    dnn = tree.query(np.mod(xy, [Lx, Ly]), k=2)[0][:, 1] / a
    nn_cdf = np.searchsorted(np.sort(dnn), RNN, side="right") / N
    dp = tree.sparse_distance_matrix(tree, RZ[-1] * a, output_type="ndarray")["v"]
    dp = np.sort(dp[dp > 0]) / a  # ordered pairs
    zcum = np.searchsorted(dp, RZ, side="right") / N / (np.pi * RZ ** 2)  # r in units of a: rho a^2 = 1

    return dict(file=os.path.basename(path), mode=h.get("mode"), N=N, Lx=Lx, Ly=Ly, rho=rho, a=a, K=K,
                kappa=float(h.get("kappa", "nan")), chi=float(h.get("chi", "nan")), chi_convention=h.get("chi_convention"),
                residual_S=float(h.get("residual_S", "nan")), S_mask=float(S[ks <= K].mean()),
                S_mask_max=float(S[ks <= K].max()), S_low=float(S[ks <= K / 4].mean()),
                H=H, H_up=float(S_bin[0]) / S_peak, S0=S0, S0_se=S0_se, S_peak=S_peak, kpeak_K=k_peak / K, alpha=alpha, lambda_phi=lam_phi,
                d_min_a=float(dnn.min()), nn_mean=float(dnn.mean()), nn_std=float(dnn.std()),
                frac_nn_below_0p3=float(np.mean(dnn < 0.3)),
                S_bin=S_bin, var=var, var_over_R=var / RA, nn_cdf=nn_cdf, zcum=zcum)


SCALARS = [("residual_S", "header residual_S", "{:.3e}"), ("S_mask", "<S>_mask (point pattern)", "{:.3e}"),
           ("S_mask_max", "max S in mask", "{:.3e}"), ("S_low", "<S> on |k| <= K/4", "{:.3e}"),
           ("S0", "S(0) (fit)", "{:.3e}"), ("H", "H = S(0)/S(k_peak)", "{:.3e}"),
           ("H_up", "S(k < 0.05 K) / S(k_peak)", "{:.3e}"), ("S_peak", "S(k_peak)", "{:.4f}"),
           ("kpeak_K", "k_peak / K", "{:.4f}"), ("alpha", "alpha", "{:.4f}"), ("lambda_phi", "Lambda/phi^1/2", "{:.4f}"),
           ("nn_mean", "<r_nn>/a", "{:.4f}"), ("nn_std", "std r_nn / a", "{:.4f}"), ("d_min_a", "d_min / a", "{:.4f}"),
           ("frac_nn_below_0p3", "P(r_nn < 0.3 a)", "{:.4f}")]
CURVES = [("S_bin", "S(k)", KB[:-1] + 0.025, "k/K"), ("var", "sigma^2(R)", RA, "R/a"),
          ("nn_cdf", "F(r_nn)", RNN, "r_nn/a"), ("zcum", "Z(r)/(rho pi r^2)", RZ, "r/a")]


def stats(rows, key):
    v = np.array([r[key] for r in rows], dtype=float)
    return np.nanmean(v, axis=0), np.nanstd(v, axis=0, ddof=1) / math.sqrt(len(rows))


def zscore(m1, s1, m2, s2):
    den = np.sqrt(s1 ** 2 + s2 ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (m1 - m2) / den, np.where(m1 == m2, 0.0, np.inf))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--a", required=True, nargs="+", help="globs of the first family")
    p.add_argument("--b", required=True, nargs="+", help="globs of the second family")
    p.add_argument("--labels", nargs=2, default=("A", "B"))
    p.add_argument("--out", required=True, help="output directory (summary.json, tables, figure)")
    p.add_argument("--title", default="")
    p.add_argument("--no-validate", action="store_true",
                   help="do not require valid obstacle files (derived controls, e.g. patterns rounded to a lattice)")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    fam = []
    for g in (args.a, args.b):
        files = sorted({f for pat in g for f in glob.glob(pat)})
        if len(files) < 2:
            sys.exit(f"need at least two realizations for {g}")
        fam.append([one(f, not args.no_validate) for f in files])
    la, lb = args.labels
    print(f"{la}: {len(fam[0])} realizations, {lb}: {len(fam[1])} realizations")
    for lab, rows in zip(args.labels, fam):
        r = rows[0]
        print(f"  {lab}: N = {r['N']}, box {r['Lx']:g} x {r['Ly']:g}, rho = {r['rho']:.6g}, a = {r['a']:.5g}, "
              f"K = {r['K']:.6g}, kappa = {r['kappa']:.4g}, chi = {r['chi']:.4g} ({r['chi_convention']})")

    summary = {"labels": list(args.labels), "n": [len(fam[0]), len(fam[1])], "scalars": {}, "curves": {}}
    lines = [f"| quantity | {la} | {lb} | z |", "|---|---|---|---|"]
    for key, name, fmt in SCALARS:
        m1, s1 = stats(fam[0], key)
        m2, s2 = stats(fam[1], key)
        z = float(zscore(m1, s1, m2, s2))
        summary["scalars"][key] = dict(a=[float(m1), float(s1)], b=[float(m2), float(s2)], z=z)
        lines.append(f"| {name} | {fmt.format(m1)} ± {fmt.format(s1)} | {fmt.format(m2)} ± {fmt.format(s2)} | {z:+.1f} |")
    table = "\n".join(lines)
    print(table)

    print("curves (bins with data in both families):")
    cur = {}
    for key, name, xs, xl in CURVES:
        m1, s1 = stats(fam[0], key)
        m2, s2 = stats(fam[1], key)
        z = zscore(m1, s1, m2, s2)
        ok = np.isfinite(m1) & np.isfinite(m2)
        if key == "S_bin":  # inside the mask the lattice residual is the known difference; compare outside it
            ok_out = ok & (xs > 1.0)
            zs = {"mask (k < K)": z[ok & (xs < 1.0)], "outside (k > K)": z[ok_out]}
        else:
            zs = {"all": z[ok]}
        cur[key] = dict(x=xs.tolist(), a=m1.tolist(), a_sem=s1.tolist(), b=m2.tolist(), b_sem=s2.tolist(), z=z.tolist())
        for part, zz in zs.items():
            zz = zz[np.isfinite(zz)]
            big = np.abs(zz) > 3
            print(f"  {name} [{part}]: {len(zz)} bins, max |z| = {np.max(np.abs(zz)):.1f}, |z| > 3 in {big.sum()}, "
                  f"|z| > 2 in {(np.abs(zz) > 2).sum()}")
            summary["curves"][f"{key} {part}"] = dict(bins=int(len(zz)), max_abs_z=float(np.max(np.abs(zz))),
                                                       n_z_gt_3=int(big.sum()), n_z_gt_2=int((np.abs(zz) > 2).sum()))
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    with open(os.path.join(args.out, "curves.json"), "w") as f:
        json.dump(cur, f)
    with open(os.path.join(args.out, "table.md"), "w") as f:
        f.write(table + "\n")
    with open(os.path.join(args.out, "realizations.csv"), "w") as f:
        keys = ["file", "mode"] + [k for k, _, _ in SCALARS]
        f.write(",".join(keys) + "\n")
        for rows in fam:
            for r in rows:
                f.write(",".join(str(r[k]) for k in keys) + "\n")
    plot(cur, args)


def plot(cur, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    la, lb = args.labels
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))

    def band(ax, key, scale=None, logy=False):
        c = cur[key]
        x = np.array(c["x"])
        for lab, m, s, col in ((la, c["a"], c["a_sem"], "tab:blue"), (lb, c["b"], c["b_sem"], "tab:orange")):
            m, s = np.array(m), np.array(s)
            if scale is not None:
                m, s = m / scale, s / scale
            ax.plot(x, m, "-", c=col, label=lab, lw=1.2)
            ax.fill_between(x, m - 2 * s, m + 2 * s, color=col, alpha=0.3, lw=0)
        if logy:
            ax.set_yscale("log")
        return x

    ax = axs[0, 0]
    band(ax, "S_bin")
    ax.axvline(1, c="gray", ls=":")
    ax.axhline(1, c="gray", ls="--", lw=0.8)
    ax.set_xlabel("k / K"); ax.set_ylabel("S(k)"); ax.legend(); ax.set_title("S(k), mean ± 2 SEM")
    ax = axs[0, 1]
    band(ax, "S_bin", logy=True)
    ax.axvline(1, c="gray", ls=":")
    ax.set_ylim(1e-24, 3); ax.set_xlim(0, 2)
    ax.set_xlabel("k / K"); ax.set_ylabel("S(k)"); ax.set_title("S(k), log scale")
    ax = axs[0, 2]
    c = cur["S_bin"]
    x, z = np.array(c["x"]), np.array(c["z"])
    m = (x > 1.0) & np.isfinite(z)
    ax.plot(x[m], z[m], "o", ms=3, c="k")
    ax.axhspan(-2, 2, color="gray", alpha=0.2)
    ax.set_xlabel("k / K"); ax.set_ylabel("z"); ax.set_title(f"z = ({la} - {lb}) / SEM, S(k) for k > K")
    ax = axs[1, 0]
    c = cur["var"]
    x = np.array(c["x"])
    for lab, mm, ss, col in ((la, c["a"], c["a_sem"], "tab:blue"), (lb, c["b"], c["b_sem"], "tab:orange")):
        mm, ss = np.array(mm) / x, np.array(ss) / x
        ax.plot(x, mm, "-", c=col, label=lab)
        ax.fill_between(x, mm - 2 * ss, mm + 2 * ss, color=col, alpha=0.3, lw=0)
    ax.set_xlabel("R / a"); ax.set_ylabel("sigma_N^2 / (R / a)"); ax.legend(); ax.set_title("number variance")
    ax = axs[1, 1]
    band(ax, "nn_cdf")
    ax.set_xlabel("r_nn / a"); ax.set_ylabel("F(r_nn)"); ax.legend(); ax.set_title("nearest-neighbour distance, distribution function")
    ax = axs[1, 2]
    band(ax, "zcum")
    ax.axhline(1, c="gray", ls="--", lw=0.8)
    ax.set_xlabel("r / a"); ax.set_ylabel("Z(r) / (rho pi r^2)"); ax.legend(); ax.set_title("cumulative coordination (integrated g)")
    for a_ in axs.flat:
        a_.grid(True, alpha=0.3)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out, f"cross_validation.{ext}"), dpi=130)
    print(f"figure saved to {os.path.join(args.out, 'cross_validation.png')}")


if __name__ == "__main__":
    main()
