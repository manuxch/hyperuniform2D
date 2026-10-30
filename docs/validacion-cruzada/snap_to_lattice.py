#!/usr/bin/env python3
"""Control of the lattice <-> continuum cross-validation: rounds each continuum obstacle file to the nearest site
of the fine lattice of the control (spacing Lx/Nx, Ly/Ny) and writes it as a new obstacle file, so that the effect
of the discretisation alone can be measured. The header is kept (K, kappa, chi of the continuum pattern); the
body, n_obstacles, rho_O, d_min, residual_S and sha256_body are recomputed, and duplicate sites are reported.

  python3 snap_to_lattice.py --nx 1012 --ny 242 --out DIR file.obs ...
"""
import argparse
import hashlib
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
import analyze_obs as ana  # noqa: E402
import obs_tools  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--nx", type=int, required=True)
p.add_argument("--ny", type=int, required=True)
p.add_argument("--out", required=True)
p.add_argument("files", nargs="+")
args = p.parse_args()
os.makedirs(args.out, exist_ok=True)
for path in args.files:
    o = obs_tools.ObsFile(path)
    h = o.header
    Lx, Ly = float(h["Lx"]), float(h["Ly"])
    ax, ay = Lx / args.nx, Ly / args.ny
    i = np.mod(np.rint(o.xy[:, 0] / ax).astype(int), args.nx)
    j = np.mod(np.rint(o.xy[:, 1] / ay).astype(int), args.ny)
    sites = np.unique(np.c_[j, i], axis=0)  # sorted by y, then x
    dup = len(o.xy) - len(sites)
    xy = np.c_[sites[:, 1] * ax, sites[:, 0] * ay]
    body = "".join(f"{float(x)!r} {float(y)!r}\n" for x, y in xy)
    N = len(xy)
    K = float(h["K"])
    rk, kx, ky, nx, ny = ana.rho_hat_grid(xy, Lx, Ly, K)
    k = np.hypot(kx, ky)
    half = (nx[:, None] > 0) | ((nx[:, None] == 0) & (ny[None, :] > 0))
    res = float((np.abs(rk[half & (k <= K)]) ** 2 / N).mean())
    new = {"n_obstacles": str(N), "rho_O": repr(N / (Lx * Ly)),
           "d_min": repr(obs_tools.min_distance(xy, Lx, Ly, True)), "residual_S": f"{res:.6g}",
           "sha256_body": hashlib.sha256(body.encode()).hexdigest()}
    lines = []
    for line in o.header_lines:
        key = line[1:].partition(":")[0].strip()
        lines.append(f"# {key}: {new[key]}" if key in new else line)
    lines.append(f"# notes: rounded to the {args.nx} x {args.ny} lattice (spacing {ax:.6g} x {ay:.6g}), {dup} duplicate sites merged")
    out = os.path.join(args.out, os.path.basename(path).replace(".obs", "_snap.obs"))
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n" + body)
    print(f"{out}: {N} obstacles ({dup} duplicates), residual_S = {res:.4g}")
