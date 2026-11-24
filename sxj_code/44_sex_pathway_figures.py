"""Figures for the GSEA by sex (43_sex_pathways.R).

  styled/sex_gsea_scatter.png   NES in males vs females per cell type (sets significant in either sex, all 5 databases);
                                % of those sets changing in the same direction in both sexes
  styled/sex_gsea_hallmark.png  MSigDB Hallmark NES in males (M) and females (F) side by side for the 4 main cell types
                                (sets significant in at least one column; * padj <= 0.05)
Usage: python 44_sex_pathway_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
DBS = {"GO_BP": "#7F7F7F", "KEGG": "#E69F00", "Reactome": "#0072B2", "WikiPathways": "#009E73", "Hallmark": "#CC3311"}


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    S = pd.read_csv(Path(res_dir) / "genotype" / "sex_gsea_summary.csv")

    # ---- scatter: NES male vs female ----
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.0))
    for ax, ct in zip(axes, CTS):
        n_all = same = 0
        for db, col in DBS.items():
            m = pd.read_csv(wd / f"sexgsea_{db}_{ct}.csv.gz").dropna()
            m = m[(m.padj_male <= 0.05) | (m.padj_female <= 0.05)]
            ax.scatter(m.NES_male, m.NES_female, s=10, color=col, alpha=0.6, lw=0, label=db)
            n_all += len(m); same += int((np.sign(m.NES_male) == np.sign(m.NES_female)).sum())
        ax.axhline(0, color="#999999", lw=0.6); ax.axvline(0, color="#999999", lw=0.6)
        ax.plot([-3.5, 3.5], [-3.5, 3.5], color="#CCCCCC", lw=0.8, ls="--")
        ax.set_xlim(-3.5, 3.5); ax.set_ylim(-3.5, 3.5)
        ax.set_xlabel("NES in males (I4M vs WT)")
        r = S[(S.cell_type == ct) & (S.database == "Hallmark")].r_NES.iloc[0]
        ax.set_title(f"{ct}\n{n_all} sets; {100 * same / max(n_all, 1):.0f}% same direction", fontsize=11)
    axes[0].set_ylabel("NES in females (I4M vs WT)")
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=2, fontsize=9)
    fig.suptitle("GSEA by sex: pathways significant (padj <= 0.05) in at least one sex", fontsize=12, y=1.04)
    fig.savefig(out / "sex_gsea_scatter.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- Hallmark heatmap: M / F per cell type ----
    cols, N, P = [], [], []
    for ct in CTS:
        m = pd.read_csv(wd / f"sexgsea_Hallmark_{ct}.csv.gz").set_index("set")
        for s, lab in [("male", "M"), ("female", "F")]:
            cols.append(f"{ {'Cardiomyocyte': 'CM', 'Endocardial': 'Endo', 'Fibroblast': 'Fib', 'Endothelial': 'EC'}[ct]} {lab}")
            N.append(m[f"NES_{s}"]); P.append(m[f"padj_{s}"])
    N = pd.concat(N, axis=1); P = pd.concat(P, axis=1)
    N.columns = P.columns = cols
    keep = (P <= 0.05).any(axis=1)
    N, P = N[keep], P[keep]
    order = N.mean(axis=1).sort_values().index
    N, P = N.loc[order], P.loc[order]
    names = [s.replace("HALLMARK_", "").replace("_", " ").capitalize() for s in N.index]
    fig, ax = plt.subplots(figsize=(6.4, 0.27 * len(N) + 1.4))
    im = ax.imshow(N.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -3, 3), aspect="auto")
    for i in range(N.shape[0]):
        for j in range(N.shape[1]):
            if P.values[i, j] <= 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
    for x in [1.5, 3.5, 5.5]:
        ax.axvline(x, color="white", lw=3)
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right", fontsize=9.5)
    ax.set_yticks(range(len(N)), names, fontsize=9)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("NES (I4M vs WT)")
    ax.set_title("MSigDB Hallmark GSEA by sex (* padj <= 0.05)", fontsize=11)
    fig.savefig(out / "sex_gsea_hallmark.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(S[S.database == "Hallmark"].to_string(index=False), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
