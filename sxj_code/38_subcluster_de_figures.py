"""Figures for the within-subcluster WT vs I4M test (37_subcluster_de.R).

  subclusters/subcluster_de_counts.png     robust DEGs per subcluster (up / down), dark = also robust at the cell-type
                                           level, light = only in the subcluster; nuclei of the smallest mouse
  subclusters/subcluster_de_consistency.png  log2FC (I4M vs WT) of the top cell-type DEGs in each subcluster: does the
                                           whole cell type change, or one subcluster?
  subclusters/subcluster_ora.png           top ORA terms (KEGG / Reactome / WikiPathways) per subcluster comparison
Usage: python 38_subcluster_de_figures.py <work_dir> <results_dir>
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
UP, DN = "#C0392B", "#2E6DB4"
FDR = 0.1


def short(s):
    return s.replace(" (main)", "").replace(" (EC-like)", "")


def main(work_dir, res_dir):
    sd, od, gd = Path(work_dir) / "subclusters", Path(res_dir) / "subclusters", Path(work_dir) / "genotype"
    S = pd.read_csv(od / "subcluster_de_summary.csv", dtype={"cluster": str})

    # ---- counts ----
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), gridspec_kw={"width_ratios": [6, 3, 4, 6]})
    for ax, ct in zip(axes, CTS):
        d = S[S.cell_type == ct].reset_index(drop=True)
        y = np.arange(len(d))
        for i, r in d.iterrows():
            sd_ = pd.read_csv(sd / f"sde_{ct}_{r.cluster}.csv.gz")
            rb = sd_[sd_.robust]
            for sign, col, xs in [(1, UP, 1), (-1, DN, -1)]:
                g = rb[np.sign(rb.logFC) == sign]
                both = int(g.robust_celltype.sum()); only = len(g) - both
                ax.barh(i, xs * both, color=col, height=0.7)
                ax.barh(i, xs * only, left=xs * both, color=col, alpha=0.35, height=0.7)
                if len(g):
                    ax.text(xs * (len(g)) + xs * d.robust.max() * 0.02, i, str(len(g)), va="center",
                            ha="left" if xs > 0 else "right", fontsize=8.5)
        lim = max(d[["up", "down"]].max().max() * 1.3, 5)
        ax.set_xlim(-lim, lim); ax.axvline(0, color="#555555", lw=0.8)
        ax.set_yticks(y, [f"{short(n)}\n(min {m} nuclei/mouse)" for n, m in zip(d.subcluster, d.min_nuclei)], fontsize=8.5)
        ax.set_ylim(len(d) - 0.5, -0.5)
        ax.set_xticks(ax.get_xticks(), [str(abs(int(v))) for v in ax.get_xticks()], fontsize=8)
        ax.set_xlim(-lim, lim)
        ax.set_title(ct, fontsize=12)
        ax.set_xlabel("down  |  up  (robust DEGs)", fontsize=9)
    fig.text(0.5, -0.03, "Dark: also robust at the whole cell-type level; light: found only within the subcluster. "
             "FDR <= 0.1, |log2FC| >= 0.25, no technical / ambient / shared flag", ha="center", fontsize=10, color="#444444")
    fig.tight_layout()
    fig.savefig(od / "subcluster_de_counts.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- consistency heatmaps ----
    NTOP = {"Cardiomyocyte": 30, "Endocardial": 30, "Fibroblast": 30, "Endothelial": 17}
    fig, axes = plt.subplots(1, 4, figsize=(19, 8.6), gridspec_kw={"width_ratios": [7, 4, 5, 7], "wspace": 0.9})
    cmap = LinearSegmentedColormap.from_list("d", [DN, "#F7F7F7", UP])
    for ax, ct in zip(axes, CTS):
        de = pd.read_csv(gd / f"de_{ct}.csv.gz")
        de = de[de.robust & ~de.gene.str.match(r"^(Gm\d|.*Rik$)")].sort_values("PValue").head(NTOP[ct])
        de = de.sort_values("logFC", ascending=False)
        d = S[S.cell_type == ct]
        cols, M, F = ["cell type"], [de.logFC.values], [de.FDR.values]
        for r in d.itertuples():
            x = pd.read_csv(sd / f"sde_{ct}_{r.cluster}.csv.gz").set_index("gene").reindex(de.gene)
            cols.append(short(r.subcluster)); M.append(x.logFC.values); F.append(x.FDR.values)
        M, F = np.array(M).T, np.array(F).T
        im = ax.imshow(M, cmap=cmap, norm=TwoSlopeNorm(0, -1.5, 1.5), aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    ax.text(j, i, "n.t.", ha="center", va="center", fontsize=6, color="#999999")
                elif F[i, j] <= FDR:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=9)
        ax.axvline(0.5, color="black", lw=1.2)
        ax.set_xticks(range(len(cols)), cols, rotation=55, ha="right", fontsize=8.5)
        ax.set_yticks(range(len(de)), de.gene, fontsize=8, style="italic")
        ax.set_title(ct, fontsize=12)
        for sp in ax.spines.values():
            sp.set_visible(False)
    cb = fig.colorbar(im, ax=axes, fraction=0.012, pad=0.01)
    cb.set_label("log2FC, I4M vs WT")
    fig.text(0.5, -0.06, "Rows: top robust DEGs of the whole cell type (by P). First column: whole cell type; others: "
             "the same test inside each subcluster. * FDR <= 0.1; n.t. = not tested (too low in that subcluster)",
             ha="center", fontsize=10, color="#444444")
    fig.savefig(od / "subcluster_de_consistency.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- ORA ----
    O = pd.read_csv(od / "subcluster_ora_top.csv", dtype={"cluster": str})
    O = O[O.database.isin(["KEGG", "Reactome", "WikiPathways"])]
    ab = {"KEGG": "K", "Reactome": "R", "WikiPathways": "W"}

    def clean(s, db):
        if db != "KEGG":
            s = re.sub(r"^(REACTOME|WP)_", "", s).replace("_", " ").lower()
            s = s[0].upper() + s[1:]
        s = f"[{ab[db]}] {s}"
        return s if len(s) <= 58 else s[:56] + "..."

    O["term"] = [clean(s, db) for s, db in zip(O.set, O.database)]
    O["group"] = [f"{short(n)}\n{'up' if d == 'Up' else 'down'} (n={k})" for n, d, k in zip(O.subcluster, O.direction, O.n_input)]
    keep = O[O.padj <= 0.05].sort_values("p").groupby("group").head(4)
    terms = list(dict.fromkeys(keep.sort_values(["database", "p"]).term))
    groups = [g for g in dict.fromkeys(O.sort_values(["cell_type", "cluster", "direction"], key=lambda c: c.map(
        {c_: i for i, c_ in enumerate(CTS)}) if c.name == "cell_type" else c).group) if g in set(keep.group)]
    # read each term's value in every shown group from the full tables
    rows = []
    for g in groups:
        r0 = O[O.group == g].iloc[0]
        for db in ab:
            f = sd / f"sora_{r0.cell_type}_{r0.cluster}_{db}_{r0.direction}.csv.gz"
            t = pd.read_csv(f)
            t["term"] = [clean(s, db) for s in t.set]
            t = t[t.term.isin(terms) & (t.overlap > 0)]
            for r in t.itertuples():
                rows.append({"group": g, "term": r.term, "ratio": r.overlap / r0.n_input, "padj": r.padj})
    F = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(1.15 * len(groups) + 6, 0.3 * len(terms) + 2.2))
    cm = LinearSegmentedColormap.from_list("ora", ["#D7263D", "#F2A3A8", "#3B5BA9"])
    sc_ = ax.scatter([groups.index(g) for g in F.group], [terms.index(t) for t in F.term], s=(F.ratio * 1500).clip(6, 380),
                     c=np.log10(F.padj.clip(1e-6)), cmap=cm, vmin=-4, vmax=0,
                     edgecolors=["black" if p <= 0.05 else "#BBBBBB" for p in F.padj],
                     linewidths=[1.1 if p <= 0.05 else 0.4 for p in F.padj])
    ax.set_xticks(range(len(groups)), groups, fontsize=8, rotation=50, ha="right")
    ax.set_yticks(range(len(terms)), terms, fontsize=8.5)
    ax.set_ylim(len(terms) - 0.5, -0.5); ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.grid(color="#EEEEEE", lw=0.6); ax.set_axisbelow(True)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.03, pad=0.02)
    cb.set_ticks([-4, -3, -2, -1, 0], labels=["1e-4", "1e-3", "0.01", "0.1", "1"]); cb.set_label("p.adjust")
    hs = [plt.scatter([], [], s=v * 1500, c="#999999") for v in [0.02, 0.05, 0.1]]
    hs.append(plt.scatter([], [], s=80, c="white", edgecolors="black", linewidths=1.1))
    ax.legend(hs, ["0.02", "0.05", "0.10", "padj <= 0.05"], title="GeneRatio", loc="upper left", bbox_to_anchor=(1.12, 1),
              frameon=False, labelspacing=1.2)
    ax.set_title("Pathways over-represented among DEGs within subclusters ([K] KEGG, [R] Reactome, [W] WikiPathways)",
                 fontsize=11)
    fig.savefig(od / "subcluster_ora.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("groups", len(groups), "terms", len(terms), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
