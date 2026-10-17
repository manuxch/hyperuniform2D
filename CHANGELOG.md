# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

[2.1.0]: https://github.com/manuxch/hyperuniform2D/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/manuxch/hyperuniform2D/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/manuxch/hyperuniform2D/releases/tag/v1.0.0
