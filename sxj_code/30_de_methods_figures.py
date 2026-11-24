"""Figures for the DE method comparison (29_de_methods.R).

  styled/de_methods_counts.png   robust DEGs per method for the four main cell types, with the part shared with the
                                 main method (edgeR QL) shaded; plus the FDR <= 0.1 exploratory list and consensus
  styled/de_methods_upset.png    overlap of robust DEG sets between the four methods (UpSet style), CM / endocardial /
                                 fibroblast
  styled/de_methods_logfc.png    per-gene log2FC: edgeR QL vs limma-voom with quality weights (cardiomyocyte, fibroblast)
Usage: python 30_de_methods_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
METH = ["edgeR_QL", "edgeR_LRT", "voom", "voom_qw"]
ML = {"edgeR_QL": "edgeR QL\n(main)", "edgeR_LRT": "edgeR\nLRT", "voom": "limma-\nvoom", "voom_qw": "voom +\nquality wts",
      "consensus_3of4": "consensus\n(>=3 of 4)", "edgeR_QL_FDR0.05": "edgeR QL\nFDR<=0.05"}
COL = {"edgeR_QL": "#4C72B0", "edgeR_LRT": "#DD8452", "voom": "#55A868", "voom_qw": "#8172B3",
       "consensus_3of4": "#333333", "edgeR_QL_FDR0.05": "#9DB6D8"}


def main(work_dir, res_dir):
    wd, rd = Path(work_dir) / "genotype", Path(res_dir) / "genotype"
    out = rd / "styled"
    s = pd.read_csv(rd / "de_methods_summary.csv")
    order = METH + ["consensus_3of4", "edgeR_QL_FDR0.05"]

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.3))
    for ax, ct in zip(axes, CTS):
        d = s[s.cell_type == ct].set_index("method").reindex(order)
        x = np.arange(len(order))
        ax.bar(x, d.robust, color=[COL[m] for m in order], alpha=0.35)
        ax.bar(x, d.overlap_with_edgeR_QL, color=[COL[m] for m in order])
        for i, (r, o) in enumerate(zip(d.robust, d.overlap_with_edgeR_QL)):
            ax.text(i, r + max(d.robust.max() * 0.015, 0.3), f"{int(r)}", ha="center", fontsize=9)
        ax.set_xticks(x, [ML[m] for m in order], fontsize=8)
        ax.set_ylim(0, max(d.robust.max() * 1.12, 5))
        ax.set_title(ct, fontsize=12)
    axes[0].set_ylabel("Robust DEGs (FDR<=0.1, |log2FC|>=0.25,\nno technical / ambient / shared flag)")
    fig.text(0.5, -0.04, "Dark part: also robust in the main method (edgeR QL); light part: found only by that method",
             ha="center", fontsize=10, color="#444444")
    fig.tight_layout()
    fig.savefig(out / "de_methods_counts.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- UpSet style ----
    fig, axes = plt.subplots(2, 3, figsize=(16, 6.2), gridspec_kw={"height_ratios": [2.2, 1]})
    for j, ct in enumerate(["Cardiomyocyte", "Endocardial", "Fibroblast"]):
        d = pd.read_csv(wd / f"dem_{ct}.csv.gz")
        rb = d[[f"{m}_robust" for m in METH]].values
        pat = pd.Series([tuple(r) for r in rb])[rb.any(1)]
        cnt = pat.value_counts().head(10)
        ax, axm = axes[0, j], axes[1, j]
        x = np.arange(len(cnt))
        ax.bar(x, cnt.values, color=["#333333" if sum(p) >= 3 else "#999999" for p in cnt.index])
        for i, v in enumerate(cnt.values):
            ax.text(i, v + cnt.values.max() * 0.02, int(v), ha="center", fontsize=8)
        ax.set_xticks([]); ax.set_xlim(-0.6, len(cnt) - 0.4)
        ax.set_title(ct, fontsize=12)
        for i, p in enumerate(cnt.index):
            on = [k for k in range(4) if p[k]]
            axm.scatter([i] * 4, range(4), s=40, color="#DDDDDD", zorder=1)
            axm.scatter([i] * len(on), on, s=40, color="#333333", zorder=2)
            if len(on) > 1:
                axm.plot([i, i], [min(on), max(on)], color="#333333", lw=1.5, zorder=1)
        axm.set_yticks(range(4), [ML[m].replace("-\n", "-").replace("\n", " ") for m in METH], fontsize=9)
        axm.set_ylim(3.5, -0.5); axm.set_xlim(-0.6, len(cnt) - 0.4); axm.set_xticks([])
        for sp in axm.spines.values():
            sp.set_visible(False)
    axes[0, 0].set_ylabel("Robust DEGs in this\ncombination of methods")
    fig.text(0.5, 0.0, "Black bars: genes robust in >= 3 of 4 methods (consensus); each column = one combination",
             ha="center", fontsize=10, color="#444444")
    fig.tight_layout()
    fig.savefig(out / "de_methods_upset.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- logFC agreement ----
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, ct in zip(axes, ["Cardiomyocyte", "Fibroblast"]):
        d = pd.read_csv(wd / f"dem_{ct}.csv.gz")
        both = d.edgeR_QL_robust & d.voom_qw_robust
        only_e = d.edgeR_QL_robust & ~d.voom_qw_robust
        only_v = ~d.edgeR_QL_robust & d.voom_qw_robust
        ax.scatter(d.edgeR_QL_logFC, d.voom_qw_logFC, s=3, color="#CCCCCC", rasterized=True, label="not robust")
        for m, c, l in [(only_v, "#8172B3", "voom+qw only"), (only_e, "#4C72B0", "edgeR QL only"),
                        (both, "#C0392B", "both")]:
            ax.scatter(d.edgeR_QL_logFC[m], d.voom_qw_logFC[m], s=7, color=c, label=f"{l} ({int(m.sum())})")
        lim = np.nanpercentile(np.abs(np.r_[d.edgeR_QL_logFC, d.voom_qw_logFC]), 99.8)
        ax.plot([-lim, lim], [-lim, lim], color="#666666", lw=0.8, ls="--")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        r = np.corrcoef(d.edgeR_QL_logFC, d.voom_qw_logFC)[0, 1]
        ax.set_title(f"{ct} (r = {r:.2f})", fontsize=12)
        ax.set_xlabel("log2FC, edgeR QL"); ax.set_ylabel("log2FC, limma-voom + quality weights")
        ax.legend(fontsize=8, frameon=False, loc="upper left", markerscale=2)
    fig.tight_layout()
    fig.savefig(out / "de_methods_logfc.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
