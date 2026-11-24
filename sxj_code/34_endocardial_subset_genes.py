"""Endocardial cells: do the genes that define each subcluster change between WT and I4M? (pseudobulk heatmap)

Rows: subcluster-specific genes of the four biological endocardial subclusters (31/32): Wilcoxon vs the other
endocardial nuclei, padj < 0.05, detected in >= 25% of the subcluster and >= 1.5x more often than in the rest,
technical genes removed; top 8 by score per subcluster.
Columns: the 8 mice (WT A-D, I4M E-H); values: log2 CPM of the summed endocardial counts per mouse (pb_Endocardial),
z-scored per gene. Right-hand column: log2FC and significance of the pseudobulk DE (13_pseudobulk_de.R; * FDR <= 0.1
and |log2FC| >= 0.25).
Output: results_dir/genotype/styled/endocardial_subset_genes.png
Usage: python 34_endocardial_subset_genes.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial"})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
SUBS = [("0", "Endo Nrg1+"), ("1", "Endo Stox2+ Sulf1+"), ("2", "Endo Ptprb+ (EC-like)"), ("5", "Valve endocardium")]
TECH = r"^(?:mt-|Rp[ls]\d|Gm\d|ENSMUSG)|Rik$|^(?:Malat1|Cmss1|Lars2|Camk1d|Gphn|AY036118|Cdk8)$"
N = 8


def main(work_dir, res_dir):
    wd = Path(work_dir)
    rd = Path(res_dir) / "genotype"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    mk = pd.read_csv(wd / "subclusters" / "markers_Endocardial.csv.gz", dtype={"group": str})
    pb = pd.read_csv(wd / "genotype" / "pb_Endocardial.csv.gz", index_col=0)[list(meta.index)]
    de = pd.read_csv(wd / "genotype" / "de_Endocardial.csv.gz").set_index("gene")
    cpm = np.log2(pb / pb.sum() * 1e6 + 1)
    rows, groups = [], []
    for c, name in SUBS:
        g = mk[(mk.group == c) & (mk.pvals_adj < 0.05) & (mk.pct_nz_group >= 0.25)
               & (mk.pct_nz_group >= 1.5 * mk.pct_nz_reference) & ~mk.names.str.contains(TECH)]
        g = [x for x in g.sort_values("scores", ascending=False).names if x in de.index and x not in rows][:N]
        rows += g
        groups += [name] * len(g)
    Z = cpm.loc[rows]
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1).replace(0, 1), axis=0)
    fig, (ax, axd) = plt.subplots(1, 2, figsize=(7.6, 0.27 * len(rows) + 1.6), gridspec_kw={"width_ratios": [8, 1.3],
                                                                                         "wspace": 0.05})
    im = ax.imshow(Z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(8), [f"{s}\n{meta.loc[s, 'sex']}" for s in meta.index])
    for t, s in zip(ax.get_xticklabels(), meta.index):
        t.set_color(GCOL[meta.loc[s, "genotype"]])
    ax.set_yticks(range(len(rows)), rows, style="italic", fontsize=9)
    ax.axvline(3.5, color="white", lw=2)
    ax.xaxis.tick_top()
    ax.text(1.5, -2.7, "WT", ha="center", color=GCOL["WT"], fontsize=11, weight="bold")
    ax.text(5.5, -2.7, "I4M", ha="center", color=GCOL["I4M"], fontsize=11, weight="bold")
    prev = None
    for i, g in enumerate(groups):
        if g != prev:
            if prev is not None:
                ax.axhline(i - 0.5, color="black", lw=1.2)
                axd.axhline(i - 0.5, color="black", lw=1.2)
            n = groups.count(g)
            ax.text(-1.9, i + n / 2 - 0.5, g, rotation=90, ha="center", va="center", fontsize=9, weight="bold")
            prev = g
    lfc = de.loc[rows, "logFC"].values
    sig = ((de.loc[rows, "FDR"] <= 0.1) & (de.loc[rows, "logFC"].abs() >= 0.25)).values
    axd.barh(range(len(rows)), lfc, color=["#C0392B" if v > 0 else "#2E6DA4" for v in lfc], height=0.7)
    for i, (v, s_) in enumerate(zip(lfc, sig)):
        if s_:
            axd.text(v + (0.05 if v > 0 else -0.05), i, "*", va="center", ha="left" if v > 0 else "right", fontsize=11)
    axd.axvline(0, color="#666666", lw=0.8)
    axd.set_ylim(len(rows) - 0.5, -0.5)
    axd.set_yticks([])
    lim = max(0.6, np.abs(lfc).max() * 1.35)
    axd.set_xlim(-lim, lim)
    axd.set_title("log2FC\nI4M vs WT", fontsize=9)
    for sp in ["top", "right", "left"]:
        axd.spines[sp].set_visible(False)
    cax = fig.add_axes([0.93, 0.35, 0.015, 0.25])
    fig.colorbar(im, cax=cax).set_label("z-score (per gene)", fontsize=8)
    fig.savefig(rd / "styled" / "endocardial_subset_genes.png", dpi=220, bbox_inches="tight")
    out = pd.DataFrame({"subcluster": groups, "gene": rows, "logFC": lfc, "FDR": de.loc[rows, "FDR"].values,
                        "DE": sig})
    print(out.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
