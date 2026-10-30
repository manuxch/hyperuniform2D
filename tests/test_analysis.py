#!/usr/bin/env python3
"""Tests of scripts/analyze_obs.py: spectral number variance against a direct count and against the Poisson
value, the Parseval tail, the index H, the surface-area coefficient of the triangular and square lattices, a stealthy
pattern written by hyper2D, and --write-alpha.

  python3 tests/test_analysis.py path/to/hyper2D path/to/scripts
"""
import hashlib
import io
import contextlib
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

failures = 0


def check(ok, what):
    global failures
    print(("  ok   " if ok else "  FAIL ") + what)
    failures += not ok


def write_points(path, xy, Lx, Ly):
    pts = sorted(((float(x), float(y)) for x, y in xy), key=lambda p: (p[1], p[0]))
    body = "".join(f"{x!r} {y!r}\n" for x, y in pts)
    header = [("format", "1"), ("generator", "test"), ("generator_version", "0"), ("mode", "continuum"),
              ("class", "poisson"), ("Lx", repr(Lx)), ("Ly", repr(Ly)), ("boundary", "periodic_x"),
              ("n_obstacles", str(len(xy))), ("rho_O", repr(len(xy) / (Lx * Ly))), ("R", "0.001"),
              ("d_min", "0"), ("seed", "0"), ("created", "2026-01-01"),
              ("sha256_body", hashlib.sha256(body.encode()).hexdigest())]
    with open(path, "w") as f:
        f.write("".join(f"# {k}: {v}\n" for k, v in header) + body)


def run(ana, argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rep = ana.analyze(argv[0], ana.build_parser().parse_args(argv + ["--no-show"]))
    return rep, out.getvalue()


def main():
    binary, scripts = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    sys.path.insert(0, scripts)
    import analyze_obs as ana
    import obs_tools

    with tempfile.TemporaryDirectory() as d:
        # Poisson: the spectral variance equals the direct count and the Poisson value rho pi R^2 (1 - pi R^2 / A)
        rng = np.random.default_rng(3)
        Lx, Ly, N = 60.0, 30.0, 1800
        poisson = os.path.join(d, "poisson.obs")
        write_points(poisson, np.c_[rng.uniform(0, Lx, N), rng.uniform(0, Ly, N)], Lx, Ly)
        rep, _ = run(ana, [poisson, "--periodic", "xy", "--rmin", "1", "--rmax", "3", "--kmax", "40"])
        R = rep["radii"]
        expected = N / (Lx * Ly) * np.pi * R ** 2 * (1 - np.pi * R ** 2 / (Lx * Ly))
        check(np.max(np.abs(rep["var_spec"] / rep["var_count"] - 1)) < 0.05,
              "Poisson: spectral sigma^2 = direct count within 5 %")
        check(np.max(np.abs(rep["var_spec"] / expected - 1)) < 0.12 and 1.8 < rep["alpha"] < 2.2,
              f"Poisson: sigma^2 = rho pi R^2 (1 - pi R^2/A) within 12 %, alpha = {rep['alpha']:.3f}")
        check(rep["H"] > 0.3, f"Poisson: H = S(0)/S(k_peak) = {rep['H']:.2f} is of order 1")

        # Surface-area coefficient of lattices (Torquato and Stillinger 2003), from the periodic replica: the
        # triangular lattice of the format example of hyperGalton2D (a = 2, 8 x 8 cells) and a square lattice
        tri = os.path.join(d, "tri.obs")
        write_points(tri, [(2.0 * i + (j % 2), math.sqrt(3.0) * j) for j in range(8) for i in range(8)], 16.0, 8 * math.sqrt(3.0))
        # k_max a >~ 400 and R up to ~1600 a are needed for 5 digits (Lambda(R) oscillates)
        rep, out = run(ana, [tri, "--periodic", "xy", "--kmax", "200", "--rmin", "2", "--rmax", "6",
                             "--surface", "3200", "--surface-radii", "32000"])
        check(abs(rep["surface_phi"] - 0.508347) < 3e-5, f"triangular lattice: Lambda/phi^(1/2) = {rep['surface_phi']:.6f} (0.508347)")
        sq = os.path.join(d, "square.obs")
        write_points(sq, [(2.0 * i, 2.0 * j) for j in range(8) for i in range(8)], 16.0, 16.0)
        rep, out = run(ana, [sq, "--periodic", "xy", "--kmax", "200", "--rmin", "2", "--rmax", "6",
                             "--surface", "3200", "--surface-radii", "32000"])
        check(abs(rep["surface_D"] - 0.4576) < 1e-4 and abs(rep["surface_phi"] - 0.516401) < 3e-5,
              f"square lattice: Lambda = {rep['surface_D']:.5f} with D = a (0.4576), Lambda/phi^(1/2) = {rep['surface_phi']:.6f} (0.516401)")

        # The Parseval tail makes the result nearly independent of k_max
        r1, _ = run(ana, [poisson, "--periodic", "xy", "--rmin", "1", "--rmax", "3", "--kmax", "10"])
        r2, _ = run(ana, [poisson, "--periodic", "xy", "--rmin", "1", "--rmax", "3", "--kmax", "30"])
        check(np.max(np.abs(r1["var_spec"] / r2["var_spec"] - 1)) < 0.03,
              "Parseval tail: k_max = 10 and 30 give sigma^2 within 3 %")

        # A stealthy pattern written by hyper2D
        with open(os.path.join(d, "s.in"), "w") as f:
            f.write("mode: continuum\nLx: 80\nLy: 20\nrho_O: 1\nkappa: 0.618\nR: 0.1\nweight_power: 2\n"
                    "num_threads: 2\nobs_file: s.obs\nannealing_file: m.dat\nsk_file: sk.dat\n")
        subprocess.run([binary, "s.in"], cwd=d, check=True, capture_output=True)
        s_obs = os.path.join(d, "s.obs")
        rep, out = run(ana, [s_obs])
        o = obs_tools.ObsFile(s_obs)
        check(f"half-plane modes (header: {o.header['n_constrained_modes']})" in out
              and f"{o.header['n_constrained_modes']} half-plane modes" in out,
              "stealthy: the circle |k| <= K holds n_constrained_modes modes")
        check(rep["S_mask_mean"] < 1e-18, f"stealthy: <S>_mask = {rep['S_mask_mean']:.2e}")
        check(rep["H"] < 1e-4, f"stealthy: H = {rep['H']:.1e} below 1e-4 (S(0) fit {rep['S0']:.1e})")
        check(0.85 < rep["alpha"] < 1.15 and np.max(np.abs(rep["var_spec"] / rep["var_count"] - 1)) < 0.06,
              f"stealthy: alpha = {rep['alpha']:.3f}, direct count agrees within 6 %")
        check(rep["g"][rep["r"] < 0.1].sum() == 0, "g(r) = 0 below the hard core (0.2 d)")

        # --write-alpha keeps the body and the checksum, and replaces an existing value
        body = o.body
        run(ana, [s_obs, "--write-alpha"])
        run(ana, [s_obs, "--write-alpha"])
        o2 = obs_tools.ObsFile(s_obs)
        errors, _ = obs_tools.validate(s_obs)
        n_alpha = sum(1 for l in o2.header_lines if l.startswith("# alpha_numvar:"))
        check(o2.body == body and not errors and n_alpha == 1 and abs(float(o2.header["alpha_numvar"]) - rep["alpha"]) < 1e-5,
              "--write-alpha: one alpha_numvar key, body and checksum unchanged, file valid")

        # test_hyper2D.py hands .obs files to analyze_obs.py
        r = subprocess.run([sys.executable, os.path.join(scripts, "test_hyper2D.py"), s_obs, "--no-show"],
                           capture_output=True, text=True)
        check(r.returncode == 0 and "number variance" in r.stdout, "test_hyper2D.py analyses .obs files")

    print(f"{failures} FAILED" if failures else "all analysis tests passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
