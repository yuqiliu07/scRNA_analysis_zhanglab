"""Step E (workflow sec. 2.8): QC report tables and figures from steps A-D.

  tables/cells_per_step.csv     per sample: barcodes -> Cell Ranger / EmptyDrops calls -> each filter -> QC pass
  tables/soupx_rho.csv          SoupX contamination estimates (03_soupx.R)
  tables/cluster_composition.csv  nuclei per sample x cluster (fraction of the sample)
  figures/qc_distributions.png  per-sample UMI, genes, % mito with that sample's thresholds
  figures/umap_uncorrected_vs_harmony.png  by sample, before / after Harmony
  figures/umap_clusters_qc.png  clusters, doublet score, % mito, UMI on the Harmony UMAP
  figures/composition.png       cluster composition per sample (sex in the label)

Usage: python 05_qc_report.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SAMPLES = list("ABCDEFGH")


def cells_per_step(work_dir):
    rows = []
    for s in SAMPLES:
        ed = pd.read_csv(work_dir / "emptydrops" / f"{s}_summary.csv").iloc[0]
        t = pd.read_csv(work_dir / "per_sample" / f"{s}_thresholds.csv").iloc[0]
        o = ad.read_h5ad(work_dir / "per_sample" / f"{s}.h5ad", backed="r").obs
        rows.append(dict(sample=s, barcodes_nonzero=int(ed["barcodes_nonzero"]),
                         cellranger_cells=int(ed["cellranger_cells"]),
                         emptydrops_cells=int(ed["emptydrops_cells"]),
                         fail_min_genes=int((~o["pass_min_genes"]).sum()),
                         fail_umi_mad=int((~o["pass_umi"]).sum()),
                         fail_genes_mad=int((~o["pass_genes"]).sum()),
                         fail_mt=int((~o["pass_mt"]).sum()),
                         pass_basic=int(t["n_pass_basic"]),
                         doublets_scrublet_auto=int(t["n_doublet_auto"]),
                         doublets_removed=int(t["n_doublet"]),
                         qc_pass=int(t["n_qc_pass"]),
                         qc_pass_frac_of_cellranger=round(t["n_qc_pass"] / ed["cellranger_cells"], 3)))
    df = pd.DataFrame(rows)
    tot = df.drop(columns=["sample", "qc_pass_frac_of_cellranger"]).sum()
    tot["sample"] = "total"
    tot["qc_pass_frac_of_cellranger"] = round(tot["qc_pass"] / tot["cellranger_cells"], 3)
    return pd.concat([df, tot.to_frame().T], ignore_index=True)


def qc_distributions(work_dir, out):
    fig, axes = plt.subplots(3, 1, figsize=(11, 10))
    data = {k: [] for k in ["total", "n_genes", "pct_mt"]}
    thr = []
    for s in SAMPLES:
        o = ad.read_h5ad(work_dir / "per_sample" / f"{s}.h5ad", backed="r").obs
        o = o[o["pass_min_genes"]]
        for k in data:
            data[k].append(o[k].values)
        thr.append(pd.read_csv(work_dir / "per_sample" / f"{s}_thresholds.csv").iloc[0])
    specs = [("total", "UMI per nucleus", True, "umi_lo", "umi_hi"),
             ("n_genes", "genes per nucleus", True, "genes_lo", "genes_hi"),
             ("pct_mt", "% mitochondrial", False, None, "mt_cut")]
    for ax, (k, lab, logy, lo, hi) in zip(axes, specs):
        vals = [np.log10(v + 1) if logy else v for v in data[k]]
        ax.violinplot(vals, showextrema=False, widths=0.85)
        for i, t in enumerate(thr, start=1):
            for col in [lo, hi]:
                if col:
                    y = np.log10(t[col] + 1) if logy else t[col]
                    ax.hlines(y, i - 0.4, i + 0.4, color="#c0392b", lw=1.5)
        ax.set_xticks(range(1, len(SAMPLES) + 1), SAMPLES)
        ax.set_ylabel(("log10 " if logy else "") + lab)
        if k == "pct_mt":
            ax.set_ylim(0, 6)
    axes[0].set_title("Candidate nuclei (EmptyDrops, >=200 genes); red = per-sample MAD thresholds")
    fig.tight_layout()
    fig.savefig(out / "qc_distributions.png", dpi=150)
    plt.close(fig)


def scatter(ax, xy, c, title, cmap=None, cats=None, s=0.3):
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(xy))
    if cats is not None:
        # one call in shuffled order, so no category is drawn on top of the others
        cmap_ = plt.get_cmap("tab20" if len(cats) > 10 else "tab10")
        code = {k: i for i, k in enumerate(cats)}
        cols = np.array([cmap_(code[k] % cmap_.N) for k in c[idx]])
        ax.scatter(xy[idx, 0], xy[idx, 1], s=s, c=cols, rasterized=True)
        for k in cats:
            ax.scatter([], [], s=20, color=cmap_(code[k] % cmap_.N), label=str(k))
        ax.legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5),
                  ncol=1 if len(cats) <= 10 else 2)
    else:
        sc_ = ax.scatter(xy[idx, 0], xy[idx, 1], s=s, c=c[idx], cmap=cmap, rasterized=True)
        plt.colorbar(sc_, ax=ax, fraction=0.04)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]), ax.set_yticks([])


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    tdir, fdir = res_dir / "tables", res_dir / "figures"
    tdir.mkdir(parents=True, exist_ok=True)
    fdir.mkdir(parents=True, exist_ok=True)

    cps = cells_per_step(work_dir)
    cps.to_csv(tdir / "cells_per_step.csv", index=False)
    print(cps.to_string(index=False), flush=True)
    sx = work_dir / "soupx"
    if sx.exists():
        pd.concat([pd.read_csv(f) for f in sorted(sx.glob("*_soupx.csv"))]).to_csv(tdir / "soupx_rho.csv",
                                                                                   index=False)
    qc_distributions(work_dir, fdir)

    m = ad.read_h5ad(work_dir / "merged" / "merged.h5ad", backed="r")
    o = m.obs.copy()
    um_raw, um_h = m.obsm["X_umap_uncorrected"][:], m.obsm["X_umap"][:]
    m.file.close()
    sex = pd.read_csv(res_dir / "tables" / "sex_check.csv", index_col=0)["inferred_sex"]
    o["cluster"] = o["cluster"].astype(str)
    clusters = sorted(o["cluster"].unique(), key=int)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    scatter(axes[0], um_raw, o["sample"].values, "No integration (by sample)", cats=SAMPLES)
    scatter(axes[1], um_h, o["sample"].values, "Harmony on sample (by sample)", cats=SAMPLES)
    fig.tight_layout()
    fig.savefig(fdir / "umap_uncorrected_vs_harmony.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    scatter(axes[0, 0], um_h, o["cluster"].values, "Leiden clusters (Harmony)", cats=clusters)
    scatter(axes[0, 1], um_h, o["doublet_score"].values, "Scrublet doublet score", cmap="viridis")
    scatter(axes[1, 0], um_h, np.clip(o["pct_mt"].values, 0, 3), "% mitochondrial (clipped at 3)",
            cmap="viridis")
    scatter(axes[1, 1], um_h, np.log10(o["total"].values), "log10 UMI", cmap="viridis")
    fig.tight_layout()
    fig.savefig(fdir / "umap_clusters_qc.png", dpi=150)
    plt.close(fig)

    comp = pd.crosstab(o["sample"], o["cluster"]).reindex(columns=clusters)
    comp.to_csv(tdir / "cluster_composition.csv")
    frac = comp.div(comp.sum(1), axis=0)
    fig, ax = plt.subplots(figsize=(11, 5))
    bottom = np.zeros(len(frac))
    cmap_ = plt.get_cmap("tab20")
    for i, c in enumerate(clusters):
        ax.bar([f"{s} ({sex.get(s, '?')})" for s in frac.index], frac[c].values, bottom=bottom,
               color=cmap_(i % 20), label=c)
        bottom += frac[c].values
    ax.set_ylabel("fraction of QC-pass nuclei")
    ax.legend(title="cluster", fontsize=7, ncol=2, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title("Cluster composition per sample (inferred sex in brackets)")
    fig.tight_layout()
    fig.savefig(fdir / "composition.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:3])
