#!/usr/bin/env python3
"""Renders the display equations of README.md as SVG, one file per GitHub theme (light, dark).

GitHub draws README math as native MathML, whose look depends on the math fonts of each browser;
SVG images look the same everywhere. Needs Matplotlib and a LaTeX installation (amsmath, amssymb).

  python3 docs/eq/make_equations.py      # from the repository root
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
                            "svg.fonttype": "path", "font.size": 15})

EQUATIONS = {
    "chi": r"\chi = \frac{M_K}{d\,(n-1)} = \frac{M_K}{2\,(n-1)}",
    "weighted_energy": r"E_w = \sum_{0<|\mathbf{k}|\le K} w(\mathbf{k})\,\bigl|\hat\rho(\mathbf{k})\bigr|^2, "
                       r"\qquad w(\mathbf{k}) \propto |\mathbf{k}|^{-p}",
    "hyperuniform": r"\lim_{|\mathbf{k}|\to 0} S(\mathbf{k}) = 0",
    "fourier": r"\hat{\rho}(\mathbf{k}) = \sum_{\mathbf{r}} s(\mathbf{r})\, e^{-i\,\mathbf{k}\cdot\mathbf{r}}, "
               r"\qquad \mathbf{k} = \frac{2\pi}{N}\,(k_x, k_y), \quad k_x, k_y \in \mathbb{Z}",
    "structure_factor": r"S(\mathbf{k}) = \frac{1}{M}\,\bigl|\hat{\rho}(\mathbf{k})\bigr|^2",
    "energy": r"E = \sum_{0<|\mathbf{k}|\le K} \bigl|\hat{\rho}(\mathbf{k})\bigr|^2",
}
THEMES = {"light": "#1f2328", "dark": "#e6edf3"}  # GitHub's text colours

out = Path(__file__).resolve().parent
for name, tex in EQUATIONS.items():
    for theme, color in THEMES.items():
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, rf"$\displaystyle {tex}$", color=color)
        fig.savefig(out / f"{name}_{theme}.svg", transparent=True, bbox_inches="tight", pad_inches=0.04)
        plt.close(fig)
        print(out / f"{name}_{theme}.svg")
