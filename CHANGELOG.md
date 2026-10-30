# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Effective-hyperuniformity index** H = S(0)/S(k_peak) (Torquato 2018) in `scripts/analyze_obs.py`, printed next to the header `residual_S`. S(0) comes from a least-squares fit S0 + c k² on |k| ≤ K/8 (`--k-fit`; up to K/2 the fit is biased negative on a lattice control, whose S grows faster than k² in the mask); H uses max(S0, 0). H ≲ 1e-4 means effectively hyperuniform.
- **Surface-area coefficient** of the number variance, σ²_N ≈ Λ R/D (Torquato and Stillinger 2003). `analyze_obs.py` prints Λ/φ^{1/2} = ⟨σ²/R⟩ 2/√(πρ), which does not depend on D, and Λ with D = d_min. With `--surface R_MAX` it averages over R up to R_MAX in the periodic replica of the pattern, with the modes beyond k_max closed exactly, πρR²[J0² + J1²](k_max R). It reproduces the lattice values: triangular 0.508347, square 0.516401 (Λ = 0.4576 with D = a).
- **`scripts/cross_validate.py`.** Compares two families of obstacle files at equal K and ρ_O, each a set of realizations: S(k) in bins of k/K, ⟨S⟩_mask, H, σ²_N(R), α, Λ, the distribution function of nearest-neighbour distances and the cumulative coordination (integrated g(r); lattice distances are discrete, so histograms are not compared bin by bin), as means with their standard errors and a z score per bin. Writes a table, JSON summaries and a figure; `--no-validate` admits derived files.
- **Tests** of H (Poisson, stealthy) and of the surface coefficient of the triangular lattice of the hyperGalton2D format example and of the square lattice.
- **`docs/validacion-cruzada.md`**: lattice ↔ continuum cross-validation at equal κ (the pending item of hyperGalton2D milestone M0a), with 16 realizations per family at κ = 0.618 and 0.757, and two controls: continuum patterns rounded to the lattice (`docs/validacion-cruzada/snap_to_lattice.py`) and lattice patterns relaxed in the continuum. Both families are hyperuniform (H ≤ 1.3e-5). At χ_point = 0.45 they agree except for the lattice residual, 0.0114, above the accepted 0.011. At χ_point = 0.3 they differ outside the mask (lattice: S = 1.37 at 1.08 K, Λ 4 % larger), a difference of the preparation protocol that survives relaxing the lattice patterns to S = 0.

## [2.2.0] — 2026-10-22

Feature release. Existing parameter files run unchanged, and a square lattice with the same seed gives the
same pattern, `sk.dat` and `annealing.dat` as 2.1.0, bit for bit; only the version in the PGM header changes.
The new features produce the obstacle files of hyperGalton2D: rectangular boxes, a continuum mode and the
file format itself.

### Added

- **Rectangular lattices**, keys `Nx` and `Ny` (default `N`).
  - The mask is a circle in physical k, (k_x/N_x)² + (k_y/N_y)² ≤ (K/2π)², i.e. an ellipse in the integer indices. A circle in the indices would give a directionally hyperuniform pattern.
  - Modes and shells are selected with the exact integer key q = k_x²N_y² + k_y²N_x², so rounding never splits a shell.
  - The phase tables use L = lcm(N_x, N_y), at most 2²⁴; the circle must fit below the Nyquist limit of both axes.
  - K is expressed in units of 2π/√(N_x N_y), which is 2π/N on a square lattice.
  - The PGM is N_x wide and N_y high; on a rectangular lattice its header also records N_x, N_y, κ and χ in both conventions.
- **Key `kappa`.** κ = K a / 2π, with a the mean spacing of the occupied sites (or of the points), sets the mask directly. It wins over `chi` and `K`, and it is the variable that pairs patterns of different generators and box shapes. `K` is accepted only on square lattices.
- **Both χ conventions** in the output: `binary_lattice`, M_K/(n_min − 1), the lattice value used so far, and `point_pattern`, M_K/[2(N − 1)].
- **Key `min_distance`.** Hard-core exclusion between occupied sites, in lattice spacings: the random start is a random sequential addition that respects it, and the annealing rejects the moves that would break it. The program prints the smallest distance between occupied sites.
- **Obstacle files for hyperGalton2D** (format version 1), keys `obs_file`, `lattice_constant`, `R`, `boundary`, `Lx`, `Ly`, `rho_O`, `arr_index` and `realizations`.
  - Every required key and the `class: stealthy` diagnostics: `K` in d⁻¹ (halfway between the last shell of the mask and the next), `n_constrained_modes`, `chi_convention`, `chi`, `kappa`, `residual_S`, `lattice_constant`.
  - `n_obstacles`, `rho_O`, `d_min`, `kappa` and `sha256_body` are computed from the values as written (shortest round-trip representation); the body is ordered by y, then x, and reproducible byte for byte.
  - `R` implies a hard core of 2R, so d_min ≥ 2R.
  - With `Lx`, `Ly`, `rho_O` the lattice is derived from the occupation fraction, keeping Lx exact.
  - `obs_file: auto` uses the recommended file name; `realizations` writes a series with consecutive seeds, and the other outputs get `_arrNNN`.
  - A `.sk` companion with S(k), k in d⁻¹.
  - SHA-256 implemented in the program (no new dependency), tested against the NIST vectors.
- **Fine-lattice control.** At the K of a continuum pattern, the lattice residual is about 0.2 ρ (χ_point = 0.3) and 0.6 ρ (χ_point = 0.45), set by the occupation ρ rather than by the schedule. With ρ = 0.02 it stays below 0.011, the residual accepted for the hyperGalton2D campaign. `min_distance` up to 0.5 a does not change it.
- **Continuum mode** (`mode: continuum`): point patterns with real coordinates in a periodic Lx × Ly box, N = round(rho_O Lx Ly).
  - Circular mask in physical k, set by `kappa` or by a point-pattern `chi` (K midway to the next shell).
  - L-BFGS with a strong-Wolfe line search, implemented in the program (same iteration counts as SciPy's L-BFGS-B on the test problems), on the weighted residual plus a hard-core term; it stops on the unweighted residual (`residual_tol`, default 1e-20).
  - Phases by recurrence along each axis, with no trigonometry per mode. The mask is stored as runs of equal n_x, so the inner loops are contiguous real arrays; the gradient is factorized along each run and computed for blocks of 32 points; ρ̂ is summed in 32 fixed groups of points. One evaluation of the energy and its gradient costs 0.7–0.95 ns per point and mode with one thread and scales to 4.5× on 6 cores. The results are identical for any thread count. Standard C++20 and `std::thread` only.
  - `init_file` starts from an obstacle file, lattice or continuum, whose box is used; the reader checks the format, the record count and `sha256_body`.
  - `minimizer: anneal` and `anneal+lbfgs`: single-point Monte Carlo with incremental ρ̂, adaptive step (`mc_step`), the lattice schedule keys and the hard core as a constraint.
  - Keys `minimizer`, `init_file`, `mc_step`, `lbfgs_memory`, `max_iterations`, `grad_tol`, `residual_tol`, `core_weight`, `sk_kmax`.
  - Outputs: minimization log, S(k), and the obstacle file (`chi_convention: point_pattern`) with its `.sk`.
  - 4 900 points in 140 × 35, χ_point = 0.3, `weight_power: 2`, 6 threads: 1–4 s to a residual below 5e-21, α ≈ 1.
- **`scripts/obs_tools.py`.** `validate` (the eight reader rules, the file conventions, the stealthy keys and the mask they declare) and `topgm` (one pixel per lattice site, identical to the PGM of hyper2D, or disks of radius R).
- **`scripts/analyze_obs.py`.** Analysis of obstacle files in grain diameters: S(k) radial and by angular sector, the check of `n_constrained_modes`, the spectral number variance (exact on the torus, with the modes beyond k_max closed by Parseval), a direct count, α and g(r). `--write-alpha` adds `alpha_numvar` to the header without touching the body. `test_hyper2D.py` hands `.obs` files to it.
- **Examples.** `examples/continuum.in` and `examples/lattice_control.in`, paired by κ.
- **Tests.** Rectangular masks, the exclusion, SHA-256, the obstacle writer and reader, L-BFGS, the continuum kernels and threads, the Monte Carlo; end-to-end tests of the obstacle files (`tests/test_obs.py`) and of the analysis (`tests/test_analysis.py`), run by `ctest` when Python is found.
- **`docs/PLAN-continuo.md`**: the development plan of this version and its measurements.

### Changed

- **Thread pool.** It runs any range function (`parallel_reduce`), no longer only the lattice ΔE.
- **Analysis script.** Rectangular PGM images; |k| in units of 2π/√(N_x N_y). Square images give the same output as before.
- **`hyper2D --help`** lists the modes, the examples and the tools.

### Fixed

- **Definition of χ in the README.** The displayed equation was M_K/[2(n − 1)], but the program has always computed M_K/(n − 1) (the README text said so too). The equation is corrected; results do not change. This lattice χ is about twice the point-pattern χ of the continuum literature at the same K.
- **χ warning.** It used to say that χ > 1/2 means no disordered ground states, which holds for the point-pattern χ. It now warns for χ_point > 1/2, and separately when the lattice χ > 1/2, where the lattice residual becomes large.
- **README equations invisible in some setups.** The display equations used a `<picture>` with light and dark variants chosen by `prefers-color-scheme`. That follows the operating system, not the GitHub theme, so a dark system with a light GitHub page showed light text on white. They are now a single SVG each, dark text on a white background, readable with any theme.
- **README map figures.** The two map figures are shown at half size and side by side.

## [2.1.0] — 2026-10-16

### Added

- **Option `weight_power` (p).**
  - The annealing minimises E_w = Σ w(k)|ρ̂(k)|², with w ∝ |k|^(−p) normalised to mean 1. All weights are positive, so E_w = 0 exactly when E = 0.
  - The unavoidable lattice residual is pushed away from the longest wavelengths, which control the number variance at large R.
  - With ρ = 0.3 and χ = 0.3, p = 2 gives α = 0.98 at N = 400, against 1.69 with p = 0. The residual for |k| ≤ K/4 drops from 3.7 % to 0.11 %, and the mask average is unchanged.
  - Default 0, the standard unweighted energy.
- **New reported quantity.** The residual of the longest wavelengths (|k| ≤ K/4) is reported and written to the PGM header, together with `weight_power`.
- **Analysis script, `--species {auto,max,min}`.** In `auto` mode the occupied pixel value comes from the `occupied M` count in the hyper2D header, so images inverted afterwards (e.g. `magick -negate`) are analysed correctly.

### Changed

- **`params.in`.** It is now the documented example of the README: a 400 × 400 lattice, density 0.3, χ = 0.3 and `weight_power: 2`, with every key commented.
- **README.**
  - Rewritten around that example.
  - New figures: the generated map with an enlarged corner, and the validation plots of the analysis script. They replace the 1.0.0 figures.
  - Math no longer depends on GitHub's MathML rendering, which varies with the math fonts of each browser and gave uneven sizes and spacing.
    - Inline math is plain Unicode and HTML (e.g. *N*², σ<sub>*N*</sub>²).
    - Display equations are SVG images rendered with LaTeX, with light and dark theme variants, by `docs/eq/make_equations.py`. Their LaTeX source is kept as alt text.
- **Analysis script.**
  - The *S*(*k*) panel shows |*k*| ≤ *N*/2 with a fitting vertical range; the corners of the zone used to stretch it.
  - Larger fonts.

## [2.0.0] — 2026-10-15

Major release. Existing parameter files still run unchanged, but several outputs and the meaning of
`steps_per_temp` changed (see **Changed**), so results are not directly comparable with 1.0.0.

### Added

- **Stealthiness parameter `chi`.**
  - Definition: χ = M_K / [d(n − 1)], with d = 2. M_K is the number of independently constrained wave vectors (each ±k pair counted once) and n = min(M, N² − M).
  - With `chi > 0` the program derives the mask radius K: it picks the integer K² whose χ is closest to the target.
  - It reports the achieved χ, which is discrete because whole shells enter at once.
  - It warns for χ > 1/2.
  - `chi: 0` (default) keeps using `K`.
- **Annealing options:**
  - `local_fraction` and `local_radius`: proposals to a nearby empty site;
  - `quench_steps`: zero-temperature quench after the annealing;
  - `energy_tol`: early stop;
  - `resync_interval`: periodic exact recomputation of the Fourier sums;
  - `parallel_min_modes`: below this mask size the run is sequential;
  - `T_initial ≤ 0` selects an automatic initial temperature.
- **Quality report.** The program prints the residual ⟨S⟩_mask / (1 − ρ): 1 for the random start, 0 for a perfectly stealthy pattern. On a lattice, E = 0 is not reached at χ ≈ 0.2–0.3, and the residual is typically 1–2 % at χ = 0.3.
- **PGM provenance.** The PGM output carries comment lines with version, N, density, K, χ (achieved and target), seed, energies and residual.
- **`annealing.dat`.** New 4th column: acceptance rate of each temperature level.
- **Command line.** `hyper2D --version` and `hyper2D --help`.
- **Tests.** Unit tests (`ctest`) cover:
  - energy against a direct trigonometric evaluation;
  - incremental updates;
  - the thread pool and threaded determinism;
  - χ selection;
  - the Parseval sum rule of S(k);
  - parameter parsing.
- **CMake.** Option `HYPER2D_NATIVE` (default ON) to build without `-march=native`, e.g. for clusters.
- **Analysis script options:**
  - `--rmin`, `--rmax` and `--n-radii`;
  - `--method {fft,sample}` and `--seed`;
  - `--save` and `--no-show`.

### Changed

- **Swap proposals.**
  - The occupied and empty sites are drawn from their own lists, so every proposal is a genuine swap.
  - In 1.0.0 two random sites were drawn, and only a fraction ρ(1 − ρ) of the attempts were valid pairs.
  - For the same `steps_per_temp` the annealing now does about 1/[ρ(1 − ρ)] times more work (≈ 4× at ρ = 0.4). Divide `steps_per_temp` by that factor to keep the 1.0.0 effort.
- **Acceptance rate.** The rate in `annealing.dat` now counts genuine proposals only.
- **`sk.dat`.**
  - It is now the radially averaged S(k) = |ρ̂(k)|² / M over the whole Brillouin zone (all nonzero k, bins of width 1 in units of 2π/N), with columns k, S(k) and number of modes.
  - In 1.0.0 it was the k_y = 0 axis only, normalised by N².
- **Energy definition.**
  - Only the half plane of k is stored, because ρ̂(−k) = ρ̂(k)*.
  - The energy is unchanged, E = Σ over the full mask = 2 Σ over the half plane, so temperatures keep their meaning.
- **Threading.**
  - The pool uses spinning barriers and the calling thread does one share of the modes. With `std::barrier`, the wake-up latency made the threaded runs slower than the sequential ones.
  - `num_threads: 0` now means min(4, CPUs in the process affinity mask), i.e. the Slurm allocation rather than every core of the node.
- **Parameter parsing.**
  - Values are parsed strictly: `64abc` is rejected.
  - Trailing `#` comments are allowed.
  - Out-of-range values stop the program with a message naming the line.
- **Analysis script (`scripts/test_hyper2D.py`).**
  - The number variance is computed exactly over all lattice centres with an FFT convolution, instead of random windows.
  - The default fit range is 5 ≤ R ≤ N/4.
  - The script also plots S(k).
  - It reports that the exponent of σ²_ρ is α − 4 by definition, so it is not an independent test.
  - It uses LaTeX only if it is installed.
- **README.** Math rewritten in GitHub's robust syntax: inline math as ``$`...`$`` and display math in ```` ```math ```` blocks.

### Fixed

- **Number of occupied sites.** It is rounded instead of truncated: `density: 0.29` with `N: 100` gave 2899 sites instead of 2900.
- **Final energy.** The reported final energy is recomputed exactly. It used to be the incrementally accumulated value, subject to drift.
- **Masks beyond Nyquist.** `K ≥ N/2` is rejected. The mask then included the aliased k_x = −N/2 row without its partner.
- **Analysis script, occupied sites.** It took the minority pixel value as the points, so above density 0.5 it analysed the empty sites and reported the wrong density.
- **Analysis script, reproducibility.** Its random sampling had no seed.
- **README.**
  - It claimed that phase factors were precomputed, but they were recomputed with sin/cos at every attempt. They now come from tables.
  - Example script name corrected (`test_hypper2D.py`).
  - Default output name corrected (`hiperuniform.pgm`).

### Performance

- Cost per genuine proposal about 30× lower: N = 120 with 3456 modes, from ≈ 120 µs to ≈ 4 µs. This comes from three changes:
  - the half plane of k;
  - table look-ups instead of trigonometry;
  - no wasted attempts.
- The default `params.in` runs in 0.8 s instead of 4.1 s, with about 4× more genuine proposals and a lower final energy.

## [1.0.0] — 2026-07-31

First release.

- Stealthy hyperuniform binary patterns on an N × N periodic lattice, generated by Metropolis simulated annealing of E = Σ_{0<|k|≤K} |ρ̂(k)|².
- Parameter file in `key: value` format.
- Thread pool for the mode sums.
- Plain PGM (P2) output, annealing monitor and S(k) along k_y = 0.
- Analysis script for the number-variance exponent.

[Unreleased]: https://github.com/manuxch/hyperuniform2D/compare/v2.2.0...HEAD
[2.2.0]: https://github.com/manuxch/hyperuniform2D/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/manuxch/hyperuniform2D/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/manuxch/hyperuniform2D/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/manuxch/hyperuniform2D/releases/tag/v1.0.0
