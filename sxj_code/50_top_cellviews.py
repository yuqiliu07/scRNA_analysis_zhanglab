"""Single-nucleus views of statistically chosen genes.

  genotype/styled/cellheatmap_top_<type>.png  single-nucleus heatmap (same look as 26_cell_heatmaps.py) of the top
                                              robust DEGs of that cell type (15 up + 15 down, ranked by P, as in
                                              49_top_degs.py); genes detected in < 2% of nuclei are not drawn
  genotype/styled/violin_early_remodel.png    cardiomyocyte nuclei, per mouse (A-D WT, E-H I4M), violins of the
                                              early-change genes of the HF section
Usage: python 50_top_cellviews.py <work_dir> <results_dir>
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
_spec = importlib.util.spec_from_file_location("cellhm", HERE / "26_cell_heatmaps.py")
cellhm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cellhm)
_spec2 = importlib.util.spec_from_file_location("focus", HERE / "48_hf_focus.py")
focus = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(focus)
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = {"Cardiomyocyte": "cardiomyocyte", "Endocardial": "endocardial", "Fibroblast": "fibroblast",
       "Endothelial": "endothelial"}
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
# early-change areas -> candidate genes (robust CM DEGs); per area the first 2 detected in >= 25% of CM nuclei are drawn
VAREAS = [["Kcnd2", "Kcnv2", "Hcn1"], ["Ttn", "Tnni3", "Myl4", "Mybpc2"], ["Cpt1b", "Acadvl", "Slc27a1"],
          ["Bckdhb", "Abat", "Dbt"], ["Hspa5", "Hsp90aa1", "Hspa8", "Calr"], ["Mgst1", "Nqo1", "Nnt"],
          ["Mef2c", "Fgf1", "Fgf9"]]
MIN_VDET = 0.25
N = 15


def top_genes(wd, ct):
    de = pd.read_csv(wd / f"de_{ct}.csv.gz")
    rb = de[de.robust].sort_values(["PValue", "FDR"])
    return {"top up in I4M": list(rb[rb.logFC > 0].gene.head(N)),
            "top down in I4M": list(rb[rb.logFC < 0].gene.head(N))}


def violin(a, meta, wd, out, ct="Cardiomyocyte", areas=None, name="violin_early_remodel.png",
           keylab="key early-change genes", mindet=MIN_VDET):
    areas = areas or VAREAS
    b = a[(a.obs["cell_type_final"].astype(str) == ct).values]
    de = pd.read_csv(wd / f"de_{ct}.csv.gz").set_index("gene")
    order = list(meta.index[meta.genotype == "WT"]) + list(meta.index[meta.genotype == "I4M"])
    rng = np.random.default_rng(0)
    genes = []
    for cand in areas:
        cand = [g for g in cand if g in b.var_names]
        det = np.asarray((b[:, cand].X > 0).mean(0)).ravel()
        genes += [g for g, d_ in zip(cand, det) if d_ >= mindet][:2]
    key = genes[:8]
    # row 1: top robust DEGs (ranked by P) detected in >= 25% of nuclei, 4 up + 4 down
    rb = de[de.robust].sort_values(["PValue", "FDR"])
    top = []
    for sub in (rb[rb.logFC > 0], rb[rb.logFC < 0]):
        cand = [g for g in sub.index if g in b.var_names]
        det = np.asarray((b[:, cand].X > 0).mean(0)).ravel()
        top += [g for g, d_ in zip(cand, det) if d_ >= mindet][:4]
    print("violin top:", top, "key:", key)
    ncol = max(len(top), len(key))
    fig, axes = plt.subplots(2, ncol, figsize=(2.6 * ncol, 6.6))
    smp = b.obs["sample"].astype(str).values
    for ax in list(axes[0][len(top):]) + list(axes[1][len(key):]):
        ax.axis("off")
    rank = {g: i + 1 for i, g in enumerate(rb.index)}
    fig.text(0.005, 0.73, "Top DEGs", rotation=90, fontsize=13, weight="bold", va="center")
    fig.text(0.005, 0.30, "Key genes", rotation=90, fontsize=13, weight="bold", va="center", color="#7A4A00")
    for ax, g in list(zip(axes[0], top)) + list(zip(axes[1], key)):
        X = b[:, g].X
        v = np.asarray(X.toarray() if hasattr(X, "toarray") else X).ravel()
        for j, s in enumerate(order):
            vv = v[smp == s]
            col = GCOL[meta.loc[s, "genotype"]]
            if vv.max() > 0:
                parts = ax.violinplot(vv, positions=[j], widths=0.85, showextrema=False, bw_method=0.3)
                for pc in parts["bodies"]:
                    pc.set_facecolor(col); pc.set_edgecolor("#333333"); pc.set_alpha(0.85); pc.set_lw(0.5)
            k = rng.choice(len(vv), min(300, len(vv)), replace=False)
            ax.scatter(j + rng.uniform(-0.22, 0.22, len(k)), vv[k], s=0.8, c="black", lw=0, rasterized=True)
            ax.scatter(j, vv.mean(), s=28, c="white", edgecolor="black", lw=0.9, zorder=5)
        r = de.loc[g]
        ax.set_title(f"{g}" + (f"  (#{rank[g]})" if g in top else ""), fontsize=12, weight="bold", style="italic")
        ax.text(0.03, 0.97, f"{r.logFC:+.2f}, FDR {r.FDR:.2g}", transform=ax.transAxes, va="top", fontsize=8.5)
        ax.set_xticks(range(8), order, fontsize=8.5)
        ax.axvline(3.5, color="#999999", lw=0.8, ls="--")
        ax.tick_params(axis="y", labelsize=8)
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=GCOL[k]) for k in GCOL] +
               [plt.scatter([], [], s=28, c="white", edgecolor="black")], ["WT (A-D)", "I4M (E-H)", "mean"],
               loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle(f"{ct} nuclei per mouse (log-normalised expression). Row 1: top robust DEGs (rank by P in brackets); row 2: {keylab}. Only genes detected in >= {int(mindet * 100)}% of nuclei; log2FC / FDR from pseudobulk edgeR",
                 fontsize=12)
    fig.tight_layout(rect=(0.015, 0.05, 1, 0.95))
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def hf_cellheatmaps(a, meta, wd, out):
    """Single-nucleus heatmaps of the HF-related genes of the HF section: cardiomyocytes = the first 3 (lead) genes
    genes of each early-change area (48_hf_focus.GENES), endocardium / fibroblasts = 48_hf_focus.KEY groups."""
    de = pd.read_csv(wd / "de_Cardiomyocyte.csv.gz").set_index("gene")
    cm = {area: gl[:3] for area, gl in focus.GENES.items()}  # first 3 of each area (curated order: lead genes first)
    cellhm.draw(a, meta, "Cardiomyocyte", cm, out, "cellheatmap_hf_cm.png")
    cellhm.draw(a, meta, "Endocardial", focus.KEY["Endocardial"], out, "cellheatmap_hf_endo.png")
    cellhm.draw(a, meta, "Fibroblast", focus.KEY["Fibroblast"], out, "cellheatmap_hf_fib.png")


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    wd = work_dir / "genotype"
    out = res_dir / "genotype" / "styled"
    meta = pd.read_csv(HERE / "sample_sheet.csv").set_index("sample")
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    cts = ["Cardiomyocyte", "Endocardial", "Fibroblast"] if ("--violin-only" in sys.argv or "--heatmap-only" in sys.argv) else list(CTS)
    keep = ((a.obs["flag"] == "ok") & a.obs["cell_type_final"].astype(str).isin(cts)).values
    a = a[keep].to_memory()
    hf_cellheatmaps(a, meta, wd, out)
    if "--heatmap-only" in sys.argv:
        return
    violin(a, meta, wd, out)
    for ct, ab in [("Endocardial", "endo"), ("Fibroblast", "fib")]:
        violin(a, meta, wd, out, ct, [list(v) for v in focus.KEY[ct].values()], f"violin_hf_{ab}.png",
               "HF-related genes", mindet=0.15)
    if "--violin-only" not in sys.argv:
        for ct, nm in CTS.items():
            cellhm.draw(a, meta, ct, top_genes(wd, ct), out, f"cellheatmap_top_{nm}.png")


if __name__ == "__main__":
    main(*sys.argv[1:3])
