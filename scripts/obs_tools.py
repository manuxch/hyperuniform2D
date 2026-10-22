#!/usr/bin/env python3
"""Tools for obstacle files in the format of hyperGalton2D, version 1 (docs/formats/obstacles-format.pdf
of that project).

  validate FILE...   checks every rule a version-1 reader enforces, plus the file conventions, the
                     diagnostics required for class: stealthy and the mask they declare.
                     Exit code 1 if any file fails.
  topgm FILE         converts to a PGM image for viewing.
                     mode lattice: one pixel per lattice site (the PGM of hyper2D), unless --px-per-d is given;
                     mode continuum: disks of radius R drawn at --px-per-d pixels per grain diameter.

  python3 scripts/obs_tools.py validate stealthy_*.obs
  python3 scripts/obs_tools.py topgm map.obs -o map.pgm [--y-up] [--px-per-d 10] [--binary]

Needs numpy; scipy (cKDTree) makes d_min fast for large files, otherwise a chunked direct search is used.
"""
import argparse
import hashlib
import math
import sys

import numpy as np

REQUIRED = ["format", "generator", "generator_version", "mode", "class", "Lx", "Ly", "boundary",
            "n_obstacles", "rho_O", "R", "d_min", "seed", "created", "sha256_body"]
STEALTHY_REQUIRED = ["K", "n_constrained_modes", "chi_convention"]
ENUMS = {"mode": {"lattice", "continuum"},
         "class": {"lattice", "perturbed", "stealthy", "rsa", "poisson"},
         "boundary": {"periodic_x", "walls_x"},
         "chi_convention": {"point_pattern", "binary_lattice"}}


class ObsFile:
    """Header (ordered dict of strings), raw body bytes and coordinates."""

    def __init__(self, path):
        self.path = path
        with open(path, "rb") as f:
            self.raw = f.read()
        lines = self.raw.split(b"\n")
        offset = 0
        self.header = {}
        self.header_lines = []
        for line in lines:
            if not line.startswith(b"#"):
                break
            offset += len(line) + 1
            text = line[1:].decode("utf-8").strip()
            self.header_lines.append(line.decode("utf-8"))
            key, sep, value = text.partition(":")
            if sep:
                self.header[key.strip()] = value.strip()
        self.body = self.raw[offset:]
        rows = [r for r in self.body.decode("utf-8").split("\n") if r.strip()]
        self.columns = [r.split() for r in rows]
        self.xy = np.array([[float(c[0]), float(c[1])] for c in self.columns]) if rows else np.zeros((0, 2))

    def f(self, key):
        return float(self.header[key])


def min_distance(xy, Lx, Ly, periodic_x):
    """Smallest centre-to-centre distance; x minimum image if periodic_x, y never periodic."""
    n = len(xy)
    if n < 2:
        return math.inf
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        cKDTree = None
    if cKDTree is not None:
        pts = xy
        if periodic_x:  # images shifted by +-Lx; query each original point against all
            pts = np.vstack([xy, xy + [Lx, 0.0], xy - [Lx, 0.0]])
        d, _ = cKDTree(pts).query(xy, k=2)
        return float(d[:, 1].min())
    best = math.inf
    for i in range(0, n, 1000):
        dx = np.abs(xy[i:i + 1000, None, 0] - xy[None, :, 0])
        if periodic_x:
            dx = np.minimum(dx, Lx - dx)
        dy = xy[i:i + 1000, None, 1] - xy[None, :, 1]
        r2 = dx * dx + dy * dy
        idx = np.arange(i, min(i + 1000, n))
        r2[idx - i, idx] = np.inf
        best = min(best, float(r2.min()))
    return math.sqrt(best)


def count_mask_modes(K, Lx, Ly):
    """M(K): half-space wave vectors k = 2 pi (nx/Lx, ny/Ly), 0 < |k| <= K (a circle in physical k)."""
    nxmax, nymax = int(K * Lx / (2 * math.pi)) + 1, int(K * Ly / (2 * math.pi)) + 1
    nx = np.arange(0, nxmax + 1)[:, None]
    ny = np.arange(-nymax, nymax + 1)[None, :]
    k2 = (2 * math.pi) ** 2 * ((nx / Lx) ** 2 + (ny / Ly) ** 2)
    half = (nx > 0) | ((nx == 0) & (ny > 0))
    inside = half & (k2 <= K * K)
    # distance of the closest shell to the boundary, relative: a small value means rounding could flip a shell
    margin = float(np.min(np.abs(np.sqrt(k2[half]) - K)) / K)
    return int(inside.sum()), margin


def validate(path, quiet=False):
    errors, warnings = [], []
    try:
        o = ObsFile(path)
    except (OSError, UnicodeDecodeError, ValueError, IndexError) as e:
        return [f"cannot read: {e}"], []
    h = o.header

    # File conventions
    if b"\r" in o.raw:
        errors.append("CR characters: line endings must be LF")
    if any(line != line.rstrip() for line in o.raw.decode("utf-8").split("\n")):
        errors.append("trailing whitespace")
    if not o.raw.endswith(b"\n"):
        errors.append("the body does not end with a newline")

    # 1. format and required keys
    missing = [k for k in REQUIRED if k not in h]
    if h.get("class") == "stealthy":
        missing += [k for k in STEALTHY_REQUIRED if k not in h]
    if "chi" in h and "chi_convention" not in h and "chi_convention" not in missing:
        missing.append("chi_convention (required with chi)")
    if missing:
        errors.append("missing keys: " + ", ".join(missing))
    if h.get("format") != "1":
        errors.append(f"format = {h.get('format')!r}: only version 1 is supported")
    for key, allowed in ENUMS.items():
        if key in h and h[key] not in allowed:
            errors.append(f"{key} = {h[key]!r} not in {sorted(allowed)}")
    if missing:
        return errors, warnings

    cols = h.get("columns", "x y")
    if cols != "x y":
        errors.append(f"columns = {cols!r}: version 1 of the simulator requires a common radius (x y)")
    bad = [i for i, c in enumerate(o.columns) if len(c) != 2]
    if bad:
        errors.append(f"{len(bad)} body lines without exactly two columns (first: line {bad[0] + 1} of the body)")
        return errors, warnings

    Lx, Ly, R, rho, n = o.f("Lx"), o.f("Ly"), o.f("R"), o.f("rho_O"), int(h["n_obstacles"])
    N = len(o.xy)
    # 2. record count
    if N != n:
        errors.append(f"{N} body records, header says n_obstacles = {n}")
    # 3. checksum
    sha = hashlib.sha256(o.body).hexdigest()
    if sha != h["sha256_body"]:
        errors.append(f"sha256_body mismatch: header {h['sha256_body']}, body {sha}")
    # 4. density
    if abs(N / (Lx * Ly) - rho) >= 1e-9:
        errors.append(f"|N/(Lx Ly) - rho_O| = {abs(N / (Lx * Ly) - rho):.3e} >= 1e-9")
    # 5. domain
    x, y = o.xy[:, 0], o.xy[:, 1]
    outside = int(np.count_nonzero((x < 0) | (x >= Lx) | (y < 0) | (y >= Ly)))
    if outside:
        errors.append(f"{outside} obstacles outside [0, Lx) x [0, Ly)")
    # ordering: ascending y, then x
    order = np.lexsort((x, y))
    if not np.array_equal(order, np.arange(N)):
        errors.append("body not ordered by ascending y, then x")
    # 6. and 7. distances
    if not R > 0:
        errors.append(f"R = {R} must be > 0")
    dmin = min_distance(o.xy, Lx, Ly, h["boundary"] == "periodic_x")
    if not abs(dmin - o.f("d_min")) <= 1e-6:
        errors.append(f"measured d_min = {dmin:.12g}, header says {h['d_min']} (boundary {h['boundary']})")
    if dmin < 2 * R:
        errors.append(f"overlap: d_min = {dmin:.12g} < 2R = {2 * R:.12g}")
    # 8. lattice positions
    if h["mode"] == "lattice":
        if "lattice_constant" not in h:
            errors.append("mode lattice needs lattice_constant (rule 8)")
        else:
            a = o.f("lattice_constant")
            dev = np.max(np.abs(o.xy - a * np.rint(o.xy / a))) if N else 0.0
            if dev >= 1e-9:
                errors.append(f"coordinates are not multiples of lattice_constant (largest deviation {dev:.3e})")

    # Hyperuniformity diagnostics
    if "K" in h:
        K = o.f("K")
        if "n_constrained_modes" not in h:
            errors.append("K present without n_constrained_modes")
        else:
            m, margin = count_mask_modes(K, Lx, Ly)
            if m != int(h["n_constrained_modes"]):
                errors.append(f"the circle |k| <= K holds {m} half-space modes, header says {h['n_constrained_modes']}")
            elif margin < 1e-9:
                warnings.append(f"a shell lies within a relative {margin:.1e} of K: rounding could change the mask")
        if "kappa" in h:
            kappa = K / (2 * math.pi * math.sqrt(rho))
            if abs(kappa - o.f("kappa")) > 1e-9 * kappa:
                errors.append(f"kappa = {h['kappa']}, K / (2 pi sqrt(rho_O)) = {kappa:.12g}")
        if "chi" in h and "n_constrained_modes" in h:
            M = int(h["n_constrained_modes"])
            if h["chi_convention"] == "point_pattern":
                chi = M / (2 * (N - 1))
            else:
                a = o.f("lattice_constant") if "lattice_constant" in h else None
                if a is None:
                    chi = None
                    warnings.append("binary_lattice chi not checked: no lattice_constant")
                else:
                    sites = round(Lx / a) * round(Ly / a)
                    chi = M / (min(N, sites - N) - 1)
            if chi is not None and abs(chi - o.f("chi")) > 1e-9 * chi:
                errors.append(f"chi = {h['chi']}, {h['chi_convention']} formula gives {chi:.12g}")
    gap = (rho ** -0.5 - 2 * R - 1.0) / 1.0
    if gap < 0:
        warnings.append(f"free gap gamma = {gap:.4g} < 0 (extreme clogging regime)")
    return errors, warnings


def cmd_validate(args):
    failed = 0
    for path in args.files:
        errors, warnings = validate(path)
        status = "FAIL" if errors else "ok"
        print(f"{status:4}  {path}")
        for e in errors:
            print(f"      error:   {e}")
        for w in warnings:
            print(f"      warning: {w}")
        failed += bool(errors)
    return 1 if failed else 0


def cmd_topgm(args):
    o = ObsFile(args.file)
    h = o.header
    Lx, Ly = float(h["Lx"]), float(h["Ly"])
    lattice_pixels = h.get("mode") == "lattice" and args.px_per_d is None and "lattice_constant" in h
    if lattice_pixels:
        a = float(h["lattice_constant"])
        W, H = round(Lx / a), round(Ly / a)
        img = np.zeros((H, W), dtype=np.uint8)
        i = np.rint(o.xy[:, 0] / a).astype(int)
        j = np.rint(o.xy[:, 1] / a).astype(int)
        img[j, i] = 1
        what = f"one pixel per lattice site, lattice_constant {h['lattice_constant']}"
    else:
        ppd = args.px_per_d if args.px_per_d is not None else 10.0
        R = args.radius if args.radius is not None else float(h["R"])
        W, H = max(1, math.ceil(Lx * ppd)), max(1, math.ceil(Ly * ppd))
        img = np.zeros((H, W), dtype=np.uint8)
        r_px = R * ppd
        span = int(math.ceil(r_px)) + 1
        offs = np.arange(-span, span + 1)
        periodic = h.get("boundary") == "periodic_x"
        for x, y in o.xy:
            cx, cy = x * ppd, y * ppd
            px = np.floor(cx).astype(int) + offs
            py = np.floor(cy).astype(int) + offs
            inside = ((px[None, :] + 0.5 - cx) ** 2 + (py[:, None] + 0.5 - cy) ** 2) <= r_px * r_px
            ys, xs = np.nonzero(inside)
            gx, gy = px[xs], py[ys]
            if periodic:
                gx = gx % W
            keep = (gx >= 0) & (gx < W) & (gy >= 0) & (gy < H)
            img[gy[keep], gx[keep]] = 1
        what = f"disks of radius {R} at {ppd} pixels per d"
    if args.y_up:
        img = img[::-1]
    out = args.output or (args.file.rsplit(".", 1)[0] + ".pgm")
    comments = [f"obs_tools topgm: {what}; 1 = obstacle; "
                f"row 0 = {'top of the domain (y up)' if args.y_up else 'y = 0'}",
                f"Nx {W} Ny {H} occupied {int(img.sum()) if lattice_pixels else len(o.xy)}"]
    comments += [line.lstrip("#").strip() for line in o.header_lines]
    with open(out, "wb") as f:
        f.write(b"P5\n" if args.binary else b"P2\n")
        for c in comments:
            f.write(f"# {c}\n".encode("utf-8"))
        f.write(f"{W} {H}\n1\n".encode())
        if args.binary:
            f.write(img.tobytes())
        else:
            f.write("\n".join(" ".join("1" if v else "0" for v in row) for row in img).encode() + b"\n")
    print(f"{out}: {W} x {H} pixels, {what}")
    return 0


def main():
    parser = argparse.ArgumentParser(description="Tools for hyperGalton2D obstacle files (format version 1).")
    sub = parser.add_subparsers(dest="command", required=True)
    v = sub.add_parser("validate", help="check files against the format rules")
    v.add_argument("files", nargs="+")
    t = sub.add_parser("topgm", help="convert to a PGM image")
    t.add_argument("file")
    t.add_argument("-o", "--output", help="output PGM (default: same name, .pgm)")
    t.add_argument("--px-per-d", type=float, default=None,
                   help="pixels per grain diameter; draws disks (default 10 for continuum; lattice: one pixel per site)")
    t.add_argument("--radius", type=float, default=None, help="disk radius in d (default: R of the file)")
    t.add_argument("--y-up", action="store_true", help="first image row = top of the domain (y increasing upwards)")
    t.add_argument("--binary", action="store_true", help="binary PGM (P5) instead of plain text (P2)")
    args = parser.parse_args()
    sys.exit(cmd_validate(args) if args.command == "validate" else cmd_topgm(args))


if __name__ == "__main__":
    main()
