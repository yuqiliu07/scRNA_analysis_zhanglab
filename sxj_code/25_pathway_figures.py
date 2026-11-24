"""Figures comparing pathway databases / methods (23_pathways.R, 24_pathways_robust.R).

  styled/pathway_db_comparison.png  GSEA: significant sets per database for four cell types, split into
                                    shared (same direction in >= 3 cell types, batch-like) and specific;
                                    plus ORA counts
  styled/pathway_themes.png         GSEA NES of selected sets (cardiomyocyte / endocardial themes and shared examples)
                                    across cardiomyocyte, endocardial, fibroblast, endothelial
Usage: python 25_pathway_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
DBS = ["Hallmark", "GO_BP", "KEGG", "Reactome", "WikiPathways"]
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
DBL = {"Hallmark": "Hallmark", "GO_BP": "GO BP", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WikiPathways"}
# (label, database, set name, theme)
THEMES = [
    ("BCAA degradation", "KEGG", "Valine, leucine and isoleucine degradation", "Cardiomyocyte"),
    ("BCAA catabolism", "Reactome", "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM", "Cardiomyocyte"),
    ("Unfolded protein response", "Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Cardiomyocyte"),
    ("HSF1-dependent transactivation", "Reactome", "REACTOME_HSF1_DEPENDENT_TRANSACTIVATION", "Cardiomyocyte"),
    ("Inflammatory response", "Hallmark", "HALLMARK_INFLAMMATORY_RESPONSE", "Cardiomyocyte"),
    ("Angiogenesis", "GO_BP", "angiogenesis (GO:0001525)", "Endocardial"),
    ("Response to wounding", "GO_BP", "response to wounding (GO:0009611)", "Endocardial"),
    ("TGF-beta signaling", "Hallmark", "HALLMARK_TGF_BETA_SIGNALING", "Endocardial"),
    ("Epithelial-mesenchymal transition", "Hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "Endocardial"),
    ("mRNA splicing", "Reactome", "REACTOME_MRNA_SPLICING", "Shared (batch-like)"),
    ("Spliceosome", "KEGG", "Spliceosome", "Shared (batch-like)"),
    ("MYC targets V1", "Hallmark", "HALLMARK_MYC_TARGETS_V1", "Shared (batch-like)"),
    ("mRNA processing", "WikiPathways", "WP_MRNA_PROCESSING", "Shared (batch-like)"),
]


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    rd = Path(res_dir) / "genotype"
    out = rd / "styled"
    rs = pd.read_csv(rd / "pathway_robust_summary.csv")
    ps = pd.read_csv(rd / "pathway_summary.csv")

    # ---- comparison bars ----
    fig, axes = plt.subplots(1, len(CTS) + 1, figsize=(19, 4.2), gridspec_kw={"width_ratios": [1] * len(CTS) + [1.2]})
    for ax, ct in zip(axes, CTS):
        d = rs[rs.cell_type == ct].set_index("database").reindex(DBS).fillna(0)
        x = np.arange(len(DBS))
        ax.bar(x, d["shared"], color="#BDBDBD", label="shared by >= 3 cell types (batch-like)")
        ax.bar(x, d["specific"], bottom=d["shared"], color="#E07A7A", label="specific to this cell type")
        for i, v in enumerate(d["sig"]):
            ax.text(i, v + 2, int(v), ha="center", fontsize=8)
        ax.set_xticks(x, [DBL[k] for k in DBS], rotation=40, ha="right")
        ax.set_ylim(0, d["sig"].max() * 1.15 + 3)
        ax.set_title(f"{ct}: GSEA", fontsize=11, pad=10)
    axes[0].set_ylabel("Significant pathways (padj <= 0.05)")
    ax = axes[-1]
    o = ps[ps.method == "ORA"].copy()
    o["grp"] = o.cell_type.str.replace("Cardiomyocyte", "CM").str.replace("Endocardial", "Endo").str.replace(
        "Fibroblast", "FB") + " " + o.direction.str.lower()
    grps = list(dict.fromkeys(o.grp))
    w = 0.15
    cols = ["#E8C27A", "#7FB3D5", "#A9DFBF", "#D7BDE2", "#F5B7B1"]
    for j, k in enumerate(DBS):
        v = [int(o[(o.grp == g) & (o.database == k)].n_sig.sum()) for g in grps]
        ax.bar(np.arange(len(grps)) + (j - 2) * w, v, w, color=cols[j], label=DBL[k])
    ax.set_xticks(range(len(grps)), grps, rotation=40, ha="right")
    ax.set_title("ORA on robust DEGs", fontsize=11)
    ax.legend(fontsize=7, frameon=False)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.42, -0.1))
    fig.tight_layout()
    fig.savefig(out / "pathway_db_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- theme dot plot ----
    rows = []
    for lab, db, s, th in THEMES:
        for ct in CTS:
            d = pd.read_csv(wd / f"pw_GSEA_{db}_{ct}.csv.gz")
            r = d[d.set == s]
            rows.append({"label": f"{lab} ({DBL[db]})", "theme": th, "cell_type": ct,
                         "NES": float(r.NES.iloc[0]) if len(r) else np.nan,
                         "padj": float(r.padj.iloc[0]) if len(r) else np.nan})
    t = pd.DataFrame(rows)
    labels = list(dict.fromkeys(t.label))
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    for i, lab in enumerate(labels):
        for j, ct in enumerate(CTS):
            r = t[(t.label == lab) & (t.cell_type == ct)].iloc[0]
            if np.isnan(r.NES):
                continue
            size = 30 + 90 * min(3, -np.log10(max(r.padj, 1e-3)))
            ax.scatter(j, i, s=size, c=[r.NES], cmap="RdBu_r", norm=norm, edgecolors="black" if r.padj <= 0.05 else "none",
                       lw=1.0)
    ax.set_xticks(range(len(CTS)), CTS, rotation=20)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xlim(-0.6, len(CTS) - 0.4)
    prev = None
    for i, lab in enumerate(labels):
        th = t[t.label == lab].theme.iloc[0]
        if th != prev:
            if prev is not None:
                ax.axhline(i - 0.5, color="#999999", lw=0.8)
            ax.text(len(CTS) - 0.35, i - 0.35, th, fontsize=9, weight="bold", va="top", color="#444444")
            prev = th
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.25)
    cb.set_label("NES (I4M vs WT; red = up, blue = down)")
    hs = [plt.scatter([], [], s=30 + 90 * v, c="#888888") for v in [1, 2, 3]]
    ax.legend(hs, ["padj 0.1", "padj 0.01", "padj 0.001"], title="dot size (black ring: padj <= 0.05)",
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False, fontsize=8, title_fontsize=8)
    fig.savefig(out / "pathway_themes.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    t.to_csv(rd / "pathway_themes.csv", index=False)
    print(t.pivot(index="label", columns="cell_type", values="NES").round(2).to_string(), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
