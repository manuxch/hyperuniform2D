#!/usr/bin/env python3
"""End-to-end test of the obstacle file output: hyper2D writes .obs files, scripts/obs_tools.py validates
them and converts them back to the PGM that hyper2D wrote.

  python3 tests/test_obs.py path/to/hyper2D path/to/scripts/obs_tools.py
"""
import glob
import hashlib
import importlib.util
import os
import subprocess
import sys
import tempfile

PARAMS = """Lx: 40
Ly: 20
rho_O: 0.2
density: 0.04
R: 0.2
kappa: 0.5
weight_power: 2
T_initial: 0
T_final: 0.01
cooling_rate: 0.8
steps_per_temp: 2000
quench_steps: 5000
seed: 21
num_threads: 1
realizations: 2
arr_index: 4
obs_file: auto
output_file: map.pgm
annealing_file: ann.dat
sk_file: sk.dat
"""

failures = 0


def check(ok, what):
    global failures
    print(("  ok   " if ok else "  FAIL ") + what)
    failures += not ok


def pgm_pixels(path):
    words = []
    with open(path) as f:
        for line in f:
            words.extend(line.partition("#")[0].split())
    w, h = int(words[1]), int(words[2])
    return w, h, words[4:4 + w * h]


def run(binary, workdir):
    with open(os.path.join(workdir, "p.in"), "w") as f:
        f.write(PARAMS)
    r = subprocess.run([binary, "p.in"], cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
    return r.returncode


def main():
    binary, tools_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    spec = importlib.util.spec_from_file_location("obs_tools", tools_path)
    tools = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tools)

    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        check(run(binary, d1) == 0, "hyper2D runs with obs_file: auto and two realizations")
        obs = sorted(glob.glob(os.path.join(d1, "*.obs")))
        check(len(obs) == 2 and obs[0].endswith("_arr004.obs") and obs[1].endswith("_arr005.obs"),
              f"two files named by the format recommendation: {[os.path.basename(p) for p in obs]}")
        check(all(os.path.exists(p[:-4] + ".sk") for p in obs), "a .sk companion next to each .obs")
        for p in obs:
            errors, warnings = tools.validate(p)
            check(not errors, f"{os.path.basename(p)} passes every reader rule {errors if errors else ''}")
        a, b = (tools.ObsFile(p) for p in obs)
        check(a.header["seed"] == "21" and b.header["seed"] == "22" and a.header["sha256_body"] != b.header["sha256_body"],
              "consecutive seeds give different realizations")
        check(hashlib.sha256(a.body).hexdigest() == a.header["sha256_body"], "sha256_body = hashlib of the body bytes")
        check(a.header["mode"] == "lattice" and a.header["chi_convention"] == "binary_lattice"
              and a.header["Lx"] == "40", "mode lattice, chi_convention binary_lattice, Lx kept exact")

        # topgm reproduces the pixels of the PGM written by hyper2D
        out = os.path.join(d1, "back.pgm")
        subprocess.run([sys.executable, tools_path, "topgm", obs[0], "-o", out], check=True, capture_output=True)
        check(pgm_pixels(out) == pgm_pixels(os.path.join(d1, "map_arr004.pgm")),
              "topgm of the .obs = PGM written by hyper2D, pixel for pixel")

        # Byte-reproducible body: same parameters, same body
        check(run(binary, d2) == 0, "second run with the same parameters")
        c = tools.ObsFile(sorted(glob.glob(os.path.join(d2, "*.obs")))[0])
        check(c.body == a.body, "the body is reproducible byte for byte")

        # The validator catches altered files
        raw = open(obs[0], "rb").read()
        cases = {
            "edited coordinate": raw.replace(a.body.split(b"\n")[0], a.body.split(b"\n")[0] + b"1", 1),
            "missing record": raw[: len(raw) - len(a.body.split(b"\n")[-2]) - 1],
            "radius larger than half the spacing": raw.replace(b"# R: 0.2\n", b"# R: 5\n"),
            "missing n_constrained_modes": b"".join(l + b"\n" for l in raw.split(b"\n")[:-1]
                                                    if not l.startswith(b"# n_constrained_modes")),
        }
        rows = a.body.split(b"\n")[:-1]
        swapped = b"\n".join([rows[1], rows[0]] + rows[2:]) + b"\n"
        cases["unordered body (checksum updated)"] = raw.replace(a.body, swapped).replace(
            a.header["sha256_body"].encode(), hashlib.sha256(swapped).hexdigest().encode())
        for name, data in cases.items():
            bad = os.path.join(d1, "bad.obs")
            with open(bad, "wb") as f:
                f.write(data)
            errors, _ = tools.validate(bad)
            check(bool(errors), f"validate rejects: {name} ({errors[0] if errors else 'accepted'})")

    # Continuum mode: stealthy to 1e-20, valid file, same body with 1 and 4 threads
    cont = """mode: continuum
Lx: 96
Ly: 24
rho_O: 1
kappa: 0.618
R: 0.1
weight_power: 2
seed: 5
obs_file: cont.obs
annealing_file: min.dat
sk_file: sk.dat
"""
    bodies = []
    for threads in (1, 4):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "c.in"), "w") as f:
                f.write(cont + f"num_threads: {threads}\n")
            r = subprocess.run([binary, "c.in"], cwd=d, capture_output=True, text=True)
            check(r.returncode == 0 and "stopped by f_target" in r.stdout
                  and (threads == 1 or f"Using {threads} threads" in r.stdout),
                  f"continuum, {threads} thread(s): L-BFGS reaches the residual target")
            p = os.path.join(d, "cont.obs")
            errors, _ = tools.validate(p)
            o = tools.ObsFile(p)
            check(not errors and o.header["mode"] == "continuum" and o.header["chi_convention"] == "point_pattern"
                  and float(o.header["residual_S"]) <= 1e-20,
                  f"continuum file valid, point_pattern chi, residual_S = {o.header.get('residual_S')} {errors if errors else ''}")
            check(os.path.exists(os.path.join(d, "cont.sk")), "continuum .sk companion written")
            bodies.append(o.body)
    check(bodies[0] == bodies[1], "continuum: 4 threads give the same body as 1 thread")

    # init_file and the Monte Carlo minimizers: restart the continuum file at a larger kappa
    with tempfile.TemporaryDirectory() as d:
        start = os.path.join(d, "start.obs")
        with open(start, "wb") as f:
            f.write(o.raw)
        for minimizer in ("anneal+lbfgs", "lbfgs"):
            with open(os.path.join(d, "r.in"), "w") as f:
                f.write(f"""mode: continuum
init_file: start.obs
kappa: 0.7
R: 0.1
weight_power: 2
minimizer: {minimizer}
T_final: 1e-5
cooling_rate: 0.8
steps_per_temp: 2000
num_threads: 1
obs_file: restart.obs
""")
            r = subprocess.run([binary, "r.in"], cwd=d, capture_output=True, text=True)
            errors, _ = tools.validate(os.path.join(d, "restart.obs"))
            check(r.returncode == 0 and "Start from 'start.obs'" in r.stdout and "stopped by f_target" in r.stdout
                  and not errors, f"init_file + {minimizer}: restart at kappa 0.7 reaches the target, valid file {errors if errors else ''}")
        with open(os.path.join(d, "x.in"), "w") as f:
            f.write("mode: continuum\ninit_file: start.obs\nLx: 50\nkappa: 0.6\n")
        r = subprocess.run([binary, "x.in"], cwd=d, capture_output=True, text=True)
        check(r.returncode != 0 and "disagree" in r.stderr, "init_file with a different Lx rejected")

    print(f"{failures} FAILED" if failures else "all obs tests passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
