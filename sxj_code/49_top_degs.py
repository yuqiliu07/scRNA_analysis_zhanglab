"""Unbiased top-DEG figures per major cell type: the most significant robust DEGs up and down in I4M.

  genotype/styled/top_deg_<type>.png   left = top robust DEGs higher in I4M, right = lower in I4M (ranked by P value,
                                       up to N each); per row: per-mouse z-score heatmap (A-D WT, E-H I4M) and the
                                       log2FC bar with FDR
Genes are chosen only by statistics (main analysis, edgeR QL ~ sex + genotype, robust = FDR <= 0.1, |log2FC| >= 0.25,
not flagged technical / ambient / shared), not by prior interest.
Usage: python 49_top_degs.py <work_dir> <results_dir> [N=15]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = {"Cardiomyocyte": "cardiomyocyte", "Endocardial": "endocardial", "Fibroblast": "fibroblast",
       "Endothelial": "endothelial"}
RED, BLUE = "#B2182B", "#2166AC"
CMAP = LinearSegmentedColormap.from_list("bwr2", ["#1A1AA6", "#7B7BD0", "#E8E8E8", "#F07070", "#E8141E"])


def fig_for(wd, ct, meta, out, n):
    de = pd.read_csv(wd / f"de_{ct}.csv.gz")
    pb = pd.read_csv(wd / f"pb_{ct}.csv.gz", index_col=0)
    order = list(meta.index[meta.genotype == "WT"]) + list(meta.index[meta.genotype == "I4M"])
    cpm = pb[order] / pb[order].sum() * 1e6
    rb = de[de.robust].sort_values(["PValue", "FDR"])
    up, dn = rb[rb.logFC > 0].head(n), rb[rb.logFC < 0].head(n)
    nr = max(len(up), len(dn))
    H = 0.27 * nr + 1.55
    fig = plt.figure(figsize=(12.3, H))
    bot, top = 0.62 / H, 1 - 0.85 / H
    for k, (d, lab, col) in enumerate([(up, "higher in I4M", RED), (dn, "lower in I4M", BLUE)]):
        x0 = 0.1 + k * 0.5
        hm = fig.add_axes([x0, top - (top - bot) * len(d) / nr, 0.16, (top - bot) * len(d) / nr])
        br = fig.add_axes([x0 + 0.17, top - (top - bot) * len(d) / nr, 0.2, (top - bot) * len(d) / nr])
        lc = np.log2(cpm.loc[d.gene] + 1)
        z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, 1), axis=0).clip(-2, 2)
        hm.imshow(z.values, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
        hm.set_yticks(range(len(d)), d.gene, fontsize=10, style="italic")
        hm.set_xticks(range(8), order, fontsize=9)
        hm.xaxis.tick_top()
        hm.axvline(3.5, color="white", lw=2.5)
        hm.text(1.5, -1.5, "WT", ha="center", fontsize=9.5, weight="bold")
        hm.text(5.5, -1.5, "I4M", ha="center", fontsize=9.5, weight="bold")
        for s in hm.spines.values():
            s.set_visible(False)
        hm.tick_params(length=0)
        y = np.arange(len(d))
        br.barh(y, d.logFC, color=col, height=0.7)
        lim = max(d.logFC.abs().max(), 0.5)
        for yy, r in zip(y, d.itertuples()):
            br.text(r.logFC + (0.05 if r.logFC > 0 else -0.05) * lim, yy, f"{r.logFC:+.2f}  (FDR {r.FDR:.1g})",
                    va="center", ha="left" if r.logFC > 0 else "right", fontsize=8.5, color="#333333")
        br.set_ylim(len(d) - 0.5, -0.5)
        br.set_xlim((0, lim * 2.1) if k == 0 else (-lim * 2.1, 0))
        br.axvline(0, color="#333333", lw=0.8)
        br.set_yticks([])
        br.spines["left"].set_visible(False)
        br.tick_params(axis="x", labelsize=8.5)
        br.set_xlabel("log2FC (I4M vs WT)", fontsize=9)
        fig.text(x0 - 0.08, 1 - 0.32 / H, f"Top {len(d)} {lab}" + (f" (of {int((rb.logFC > 0).sum()) if k == 0 else int((rb.logFC < 0).sum())} robust)"),
                 fontsize=12, weight="bold", color=col)
    fig.text(0.5, 0.012, f"{ct}: most significant robust DEGs (ranked by P value; FDR <= 0.1, |log2FC| >= 0.25, not "
             "flagged). Heatmap = per-mouse z-score of log2(CPM + 1); red = high, blue = low.", ha="center",
             fontsize=9, color="#555555")
    name = f"top_deg_{CTS[ct]}.png"
    fig.savefig(out / name, dpi=220)
    plt.close(fig)
    return up, dn


def main(work_dir, res_dir, n=15):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    meta = pd.read_csv(HERE / "sample_sheet.csv").set_index("sample")
    for ct in CTS:
        up, dn = fig_for(wd, ct, meta, out, int(n))
        print(ct, "UP:", " ".join(up.gene), "| DOWN:", " ".join(dn.gene))


if __name__ == "__main__":
    main(*sys.argv[1:4])
