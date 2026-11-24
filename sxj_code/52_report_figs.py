"""Figures for the collaborator report deck that no earlier script draws.

  qc_overview.png       nuclei kept at each QC step per mouse + per-nucleus genes / UMI / % mito of the final nuclei
  celltype_counts.png   final high-quality nuclei per major cell type (WT / I4M stacked)
  pw_dot_<ab>.png       GSEA dot plot for cardiomyocytes (cm), endocardium (endo), fibroblasts (fib): per database
                        (GO BP, KEGG, Reactome, WikiPathways, Hallmark) the top 3 up + top 3 down sets by P;
                        x = NES, dot size = -log10 padj, grey = padj > 0.05 (results of 28_pathway_methods, not re-tested)
Usage: python 52_report_figs.py <work_dir> <res_dir>
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
_spec = importlib.util.spec_from_file_location("focus", HERE / "48_hf_focus.py")
focus = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(focus)
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
WT, MUT = list("ABCD"), list("EFGH")
RED, BLUE, GREY = "#B2182B", "#2166AC", "#BDBDBD"
NAME = {"Lymphatic_EC": "Lymphatic EC", "Neural_Schwann": "Neural / Schwann", "SMC": "Smooth muscle"}


def qc_overview(work_dir, res_dir, obs, out):
    steps = pd.read_csv(Path(res_dir) / "tables" / "cells_per_step.csv").set_index("sample")
    ok = obs[obs.flag == "ok"]
    after_cb = obs.groupby("sample", observed=True).size()
    final = ok.groupby("sample", observed=True).size()
    S = list("ABCDEFGH")
    bars = [("Cell Ranger", steps.loc[S, "cellranger_cells"], "#D9D9D9"),
            ("Basic QC + Scrublet", steps.loc[S, "qc_pass"], "#9FB7CF"),
            ("CellBender + doublet clusters", after_cb.loc[S], "#5B84AE"),
            ("Final (low-quality clusters removed)", final.loc[S], "#24466B")]
    fig = plt.figure(figsize=(13.5, 3.7))
    ax = fig.add_axes([0.045, 0.2, 0.3, 0.68])
    w = 0.2
    x = np.arange(8)
    for i, (lab, v, c) in enumerate(bars):
        ax.bar(x + (i - 1.5) * w, np.asarray(v) / 1000, w, color=c, label=lab)
    ax.set_xticks(x, [f"{s}\n{'WT' if s in WT else 'I4M'}" for s in S], fontsize=9)
    for t, s in zip(ax.get_xticklabels(), S):
        t.set_color(GCOL["WT" if s in WT else "I4M"])
    ax.set_ylabel("Nuclei (×1,000)", fontsize=10)
    tot = [int(steps.loc["total", "cellranger_cells"]), int(steps.loc["total", "qc_pass"]), len(obs), len(ok)]
    ax.set_title("Nuclei kept per mouse: " + " → ".join(f"{t:,}" for t in tot), fontsize=10, loc="left")
    ax.legend(fontsize=7.5, frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.0))
    ax.set_ylim(0, 72)
    for j, (col, lab, tr) in enumerate([("n_genes_raw", "Genes per nucleus", np.log10),
                                         ("total_raw", "UMI per nucleus", np.log10),
                                         ("pct_mt", "% mitochondrial reads", None)]):
        ax = fig.add_axes([0.41 + j * 0.2, 0.2, 0.17, 0.68])
        data = [ok.loc[ok["sample"] == s, col].to_numpy(float) for s in S]
        if tr is not None:
            data = [tr(d) for d in data]
        else:
            data = [np.clip(d, 0, 3) for d in data]
        parts = ax.violinplot(data, showextrema=False, widths=0.85)
        for b, s in zip(parts["bodies"], S):
            b.set_facecolor(GCOL["WT" if s in WT else "I4M"]); b.set_alpha(0.75); b.set_edgecolor("none")
        med = [np.median(d) for d in data]
        ax.scatter(range(1, 9), med, s=10, color="black", zorder=3)
        ax.set_xticks(range(1, 9), S, fontsize=9)
        for t, s in zip(ax.get_xticklabels(), S):
            t.set_color(GCOL["WT" if s in WT else "I4M"])
        if tr is not None:
            ticks = [200, 500, 1000, 3000, 10000] if col == "n_genes_raw" else [300, 1000, 3000, 10000, 30000]
            ax.set_yticks(np.log10(ticks), [f"{t:,}" for t in ticks], fontsize=8.5)
            mm = [f"{10 ** m:,.0f}" for m in med]
        else:
            ax.set_ylim(0, 3); ax.tick_params(labelsize=8.5)
            mm = [f"{m:.2f}" for m in med]
        ax.set_title(f"{lab}\n(median {min(mm, key=lambda v: float(v.replace(',', '')))}–"
                     f"{max(mm, key=lambda v: float(v.replace(',', '')))})", fontsize=10)
    fig.text(0.41, 0.03, f"Final {len(ok):,} nuclei; black dot = median of each mouse; brown = WT, green = I4M",
             fontsize=8.5, color="#555555")
    fig.savefig(out / "qc_overview.png", dpi=220)
    plt.close(fig)


def celltype_counts(obs, out):
    ok = obs[obs.flag == "ok"].copy()
    ok["geno"] = np.where(ok["sample"].isin(WT), "WT", "I4M")
    t = pd.crosstab(ok.cell_type_final, ok.geno)[["WT", "I4M"]]
    t = t.loc[t.sum(1).sort_values().index]
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    y = np.arange(len(t))
    ax.barh(y, t.WT / 1000, color=GCOL["WT"], label="WT (4 mice)")
    ax.barh(y, t.I4M / 1000, left=t.WT / 1000, color=GCOL["I4M"], label="I4M (4 mice)")
    tot = t.sum(1)
    for yi, (n, w_, m_) in enumerate(zip(tot, t.WT, t.I4M)):
        ax.text(n / 1000 + 0.6, yi, f"{n:,}  ({100 * n / tot.sum():.1f}%)", va="center", fontsize=8.5)
    ax.set_yticks(y, [NAME.get(c, c) for c in t.index], fontsize=9.5)
    ax.set_xlabel("Nuclei (×1,000)", fontsize=10)
    ax.set_xlim(0, tot.max() / 1000 * 1.35)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title(f"{tot.sum():,} high-quality nuclei, 12 major cell types", fontsize=10.5, loc="left")
    fig.tight_layout()
    fig.savefig(out / "celltype_counts.png", dpi=220)
    plt.close(fig)
    t.assign(total=tot).iloc[::-1].to_csv(out / "celltype_counts.csv")


def pw_dot(wd, out, ct, ab, k=3):
    dbs = list(focus.DBL)
    fig = plt.figure(figsize=(12.3, 4.9))
    for c, dl in enumerate([dbs[:3], dbs[3:]]):
        rows = []
        for db in dl:
            rows.append(("head", db, None, None, None))
            g = pd.read_csv(wd / f"pw_GSEA_{db}_{ct}.csv.gz")
            for d in (g[g.NES > 0].sort_values("p").head(k), g[g.NES < 0].sort_values("p").head(k)):
                rows += [("row", focus.short(db, r.set), r.NES, r.padj, r.size) for r in d.itertuples()]
        ax = fig.add_axes([0.25 + c * 0.5, 0.17, 0.2, 0.78])
        y = 0
        for kind, lab, v, pa, _ in rows:
            if kind == "head":
                ax.text(-0.02, y, {"GO_BP": "GO BP"}.get(lab, lab), transform=ax.get_yaxis_transform(), ha="right",
                        va="center", fontsize=10, weight="bold", color="#333333")
                y += 1
                continue
            sig = pa <= 0.05
            ax.plot([0, v], [y, y], color="#DDDDDD", lw=1, zorder=1)
            ax.scatter(v, y, s=25 + 28 * min(-np.log10(max(pa, 1e-12)), 8),
                       color=(RED if v > 0 else BLUE) if sig else GREY, zorder=2, edgecolor="white", lw=0.5)
            ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5,
                    color="#222222" if sig else "#8A8A8A")
            y += 1
        ax.set_ylim(y - 0.5, -0.6)
        ax.set_xlim(-3.6, 3.6)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_yticks([])
        ax.set_xlabel("NES  (− lower in I4M, + higher in I4M)", fontsize=9)
        ax.tick_params(labelsize=8.5)
    for i, p in enumerate([0.05, 1e-3, 1e-6]):
        fig.gca().scatter([], [], s=25 + 28 * min(-np.log10(p), 8), color="#777777", label=f"padj {p:g}")
    fig.legend(*fig.gca().get_legend_handles_labels(), loc="lower right", fontsize=8, frameon=False, ncol=3,
               bbox_to_anchor=(0.99, 0.0))
    fig.savefig(out / f"pw_dot_{ab}.png", dpi=220)
    plt.close(fig)


def main(work_dir, res_dir):
    out = Path(res_dir) / "genotype"
    wd = Path(work_dir) / "genotype"
    for ct, ab in [("Cardiomyocyte", "cm"), ("Endocardial", "endo"), ("Fibroblast", "fib")]:
        pw_dot(wd, out, ct, ab)
    obs = ad.read_h5ad(Path(work_dir) / "merged" / "final.h5ad", backed="r").obs[
        ["sample", "n_genes_raw", "total_raw", "pct_mt", "flag", "cell_type_final"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    celltype_counts(obs, out)
    qc_overview(work_dir, res_dir, obs, out)
    print("done ->", out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
