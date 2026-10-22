#!/usr/bin/env python3
"""Renders the display equations of README.md as SVG images: dark text on a white background.

GitHub draws README math as native MathML, whose look depends on the math fonts of each browser;
SVG images look the same everywhere. The white background keeps them readable whatever the theme:
a <picture> with a prefers-color-scheme source follows the operating system, not the GitHub theme,
and showed light text on a white page when the two differed.
Needs Matplotlib and a LaTeX installation (amsmath, amssymb).

  python3 docs/eq/make_equations.py            # all, from the repository root
  python3 docs/eq/make_equations.py chi fourier  # only these
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
                            "svg.fonttype": "path", "font.size": 15})

EQUATIONS = {
    "chi": r"\chi = \frac{M_K}{n_{\min}-1}, \qquad n_{\min} = \min(M,\, N_xN_y - M)",
    "weighted_energy": r"E_w = \sum_{0<|\mathbf{k}|\le K} w(\mathbf{k})\,\bigl|\hat\rho(\mathbf{k})\bigr|^2, "
                       r"\qquad w(\mathbf{k}) \propto |\mathbf{k}|^{-p}",
    "hyperuniform": r"\lim_{|\mathbf{k}|\to 0} S(\mathbf{k}) = 0",
    "fourier": r"\hat{\rho}(\mathbf{k}) = \sum_{\mathbf{r}} s(\mathbf{r})\, e^{-i\,\mathbf{k}\cdot\mathbf{r}}, "
               r"\qquad \mathbf{k} = 2\pi\Bigl(\frac{k_x}{N_x}, \frac{k_y}{N_y}\Bigr), \quad k_x, k_y \in \mathbb{Z}",
    "structure_factor": r"S(\mathbf{k}) = \frac{1}{M}\,\bigl|\hat{\rho}(\mathbf{k})\bigr|^2",
    "energy": r"E = \sum_{0<|\mathbf{k}|\le K} \bigl|\hat{\rho}(\mathbf{k})\bigr|^2",
}
COLOR = "#1f2328"  # GitHub's light-theme text colour

out = Path(__file__).resolve().parent
selected = sys.argv[1:] or list(EQUATIONS)
for name in selected:
    tex = EQUATIONS[name]
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, rf"$\displaystyle {tex}$", color=COLOR)
    fig.savefig(out / f"{name}.svg", facecolor="white", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(out / f"{name}.svg")
