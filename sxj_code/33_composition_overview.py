"""Overview of composition changes, WT vs TNNC1-I4M: all cell types and all subclusters in one dot plot
(effect size and significance, in the style of a scCODA effect-size dot plot).

Input:  results_dir/genotype/composition_test.csv (12_genotype_prep.py), subcluster_composition_test.csv
        (32_subcluster_figures.py); both asin-sqrt OLS ~ sex + genotype per mouse, BH.
x = log2(mean fraction I4M / mean fraction WT); dot size = -log10(FDR); colour = direction; black ring = FDR <= 0.05.
Output: results_dir/genotype/styled/composition_overview.png
Usage: python 33_composition_overview.py <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
UP, DOWN = "#C0392B", "#2E6DA4"
CLIP = 2.0
LOWQ = ("Low quality", "Glia-like", "Minor")


def panel(ax, d, labels, groups=None):
    y = np.arange(len(d))
    lfc = np.log2(d.mean_I4M.clip(1e-6) / d.mean_WT.clip(1e-6))
    x = lfc.clip(-CLIP, CLIP)
    size = 30 + 110 * (-np.log10(d.fdr.clip(1e-4))).clip(0, 3)
    col = [UP if v > 0 else DOWN for v in lfc]
    alpha = [0.35 if str(n).startswith(LOWQ) else 1 for n in labels]
    for yi, xi, si, ci, ai, f in zip(y, x, size, col, alpha, d.fdr):
        ax.scatter(xi, yi, s=si, color=ci, alpha=ai, edgecolors="black" if f <= 0.05 else "none", lw=1.2, zorder=3)
    for yi, v in zip(y, lfc):
        if abs(v) > CLIP:
            ax.text(np.sign(v) * CLIP * 1.08, yi, f"{v:.1f}", va="center", ha="left" if v > 0 else "right", fontsize=7)
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(y, labels, fontsize=9)
    ax.set_ylim(len(d) - 0.5, -0.5)
    ax.set_xlim(-CLIP * 1.35, CLIP * 1.35)
    ax.set_xlabel("log2 (I4M / WT), mean fraction")
    if groups is not None:
        prev = None
        for i, g in enumerate(groups):
            if g != prev:
                if prev is not None:
                    ax.axhline(i - 0.5, color="#BBBBBB", lw=0.8)
                ax.text(CLIP * 1.33, i - 0.3, g, fontsize=9, weight="bold", ha="right", va="top", color="#444444")
                prev = g


def main(res_dir):
    rd = Path(res_dir) / "genotype"
    ct = pd.read_csv(rd / "composition_test.csv").sort_values("mean_WT", ascending=False).reset_index(drop=True)
    sb = pd.read_csv(rd / "subcluster_composition_test.csv")
    order = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
    sb["o"] = sb.cell_type.map({k: i for i, k in enumerate(order)})
    sb = sb.sort_values(["o", "mean_WT"], ascending=[True, False]).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 9.2), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.55})
    panel(axes[0], ct, [c.replace("_", " ") for c in ct.cell_type])
    axes[0].set_title("Step 1: all cell types\n(% of all nuclei)", fontsize=12)
    panel(axes[1], sb, list(sb.subcluster), list(sb.cell_type))
    axes[1].set_title("Step 2: subclusters within each cell type\n(% of that cell type)", fontsize=12)
    hs = [plt.scatter([], [], s=30 + 110 * v, color="#888888") for v in [1, 2, 3]]
    hs += [plt.scatter([], [], s=80, color=UP), plt.scatter([], [], s=80, color=DOWN),
           plt.scatter([], [], s=80, color="white", edgecolors="black", lw=1.2)]
    fig.legend(hs, ["FDR 0.1", "FDR 0.01", "FDR 0.001", "higher in I4M", "lower in I4M", "FDR <= 0.05"],
               loc="lower center", ncol=6, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.text(0.5, -0.045, "Faded: low-quality / doublet-like groups. Values beyond +/-2 are clipped and printed.",
             ha="center", fontsize=9, color="#555555")
    fig.savefig(rd / "styled" / "composition_overview.png", dpi=220, bbox_inches="tight")
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
