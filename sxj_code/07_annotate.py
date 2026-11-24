"""Step G: per-cluster review and cell-type annotation of the CellBender re-clustering (06_cellbender_recluster.py).

Each step-06 cluster is classified by rules (evidence columns are written to annotation_clusters.csv):

  1. low quality     median genes per nucleus (CellBender counts) < 200
  2. doublets        > 50% of nuclei are Scrublet-auto doublets AND > 50% carry a second lineage (>= 2 genes of a
                     LINEAGE panel) outside the lineage family of the cluster's top marker panel
  3. epicardial      Muc16 detected in > 30% of nuclei
  4. otherwise       cell type = top marker panel (04_merge_cluster.PANELS, mean log-expression);
                     subtypes: Endothelial -> EC_arterial if Vwf detected in > 15% of nuclei, else EC_capillary;
                     Cardiomyocyte -> CM_conduction if Cntn2 detected in > 20%, else vCM
  The rules reproduce the manual review of the first run (heterotypic doublet clusters with 70-84% two-lineage
  nuclei and 62-74% Scrublet-auto calls; one low-quality cluster with 99 median genes).

The kept nuclei are re-processed (HVG, PCA, Harmony, UMAP, Leiden) to check that no new mixed clusters appear;
the cell-type labels come from the step-06 clusters. Second pass: a new cluster is removed as doublets when > 50%
of its nuclei are Scrublet-auto doublets and > 50% carry a second lineage outside the cluster's own lineage family
(first run: a pericyte cluster of 1,956 nuclei with endothelial markers in 42-69% of nuclei, Scrublet-auto 64%).
The embedding is not recomputed after the second pass.

Usage: python 07_annotate.py <work_dir> <results_dir>
Output: <work_dir>/merged/annotated.h5ad
        <results_dir>/tables/annotation_clusters.csv, celltype_counts.csv, celltype_fractions.csv,
            annot_recluster_review.csv, annot_marker_detection.csv
        <results_dir>/figures/annot_umap.png, annot_markers.png, annot_composition.png
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
import scanpy as sc

_here = Path(__file__).parent


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, _here / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s04 = _load("step04", "04_merge_cluster.py")
s05 = _load("step05", "05_qc_report.py")
SAMPLES, SEED, mixing = s04.SAMPLES, s04.SEED, s04.mixing

MARKERS = {
    "Cardiomyocyte": ["Myh6", "Tnnt2", "Ryr2", "Rbm20", "Myl2", "Mhrt"],
    "CM_conduction": ["Cntn2", "Cacna1g", "Hcn4"],
    "Fibroblast": ["Dcn", "Abca8a", "Pdgfra", "Gsn"],
    "Endothelial": ["Pecam1", "Flt1", "Ptprb", "Kdr"],
    "EC_arterial": ["Vwf", "Eln", "Mecom"],
    "Endocardial": ["Npr3", "Cdh11"],
    "Lymphatic_EC": ["Prox1", "Reln", "Mmrn1"],
    "Pericyte": ["Rgs5", "Kcnj8", "Pdgfrb"],
    "SMC": ["Myh11", "Acta2"],
    "Macrophage": ["F13a1", "Mrc1", "Cd163", "Adgre1"],
    "Epicardial": ["Msln", "Upk3b", "Wt1", "Muc16"],
}
LINEAGE = {  # >= 2 genes detected = lineage present in a nucleus
    "CM": ["Myh6", "Tnnt2", "Ryr2", "Ttn", "Rbm20"],
    "Fib": ["Dcn", "Col1a1", "Pdgfra", "Abca8a", "Gsn"],
    "EC": ["Pecam1", "Cdh5", "Flt1", "Kdr", "Ptprb"],
    "Peri": ["Rgs5", "Kcnj8", "Abcc9", "Pdgfrb"],
    "Mac": ["C1qa", "F13a1", "Mrc1", "Cd163", "Adgre1", "Csf1r"],
    "SMC": ["Myh11", "Acta2", "Cnn1"],
}
FAMILY = {"Cardiomyocyte": {"CM"}, "Fibroblast": {"Fib"}, "Endothelial": {"EC"}, "Endocardial": {"EC"},
          "Lymphatic_EC": {"EC"}, "Pericyte": {"Peri", "SMC"}, "SMC": {"SMC", "Peri"}, "Macrophage": {"Mac"},
          "Epicardial": set()}
ORDER = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
         "Macrophage", "Epicardial"]


PANEL_TO_FAMILY = {"Cardiomyocyte": {"CM"}, "Fibroblast": {"Fib"}, "Endothelial": {"EC"}, "Endocardial": {"EC"},
                   "Lymphatic_EC": {"EC"}, "Pericyte": {"Peri", "SMC"}, "SMC": {"SMC", "Peri"},
                   "Macrophage": {"Mac"}}


def lineage_presence(m):
    out = {}
    for name, genes in LINEAGE.items():
        g = [x for x in genes if x in m.var_names]
        out[name] = np.asarray((m[:, g].layers["counts"] > 0).sum(1)).ravel() >= 2
    return pd.DataFrame(out, index=m.obs_names)


def classify_clusters(m, key="cluster_cb"):
    """Rule-based decision per cluster (see module docstring)."""
    cl = m.obs[key].astype(str)
    scores = {}
    for name, genes in s04.PANELS.items():
        g = [x for x in genes if x in m.var_names]
        if g:
            scores[name] = np.asarray(m[:, g].X.mean(1)).ravel()
    ps = pd.DataFrame(scores, index=m.obs_names).groupby(cl.values).mean()
    top = ps.idxmax(1)
    lin = lineage_presence(m)
    det = {}
    for g in ["Vwf", "Cntn2", "Muc16"]:
        det[g] = pd.Series(np.asarray((m[:, g].layers["counts"] > 0).todense()).ravel(), index=m.obs_names)
    rows = {}
    for c, idx in m.obs.groupby(cl.values, observed=True).groups.items():
        t = top[c]
        fam = PANEL_TO_FAMILY.get(t, set())
        second = lin.loc[idx, [l for l in LINEAGE if l not in fam]].any(axis=1).mean()
        lin_frac = lin.loc[idx].mean()
        r = {"top_panel": t, "n": len(idx),
             "median_genes_cb": float(m.obs.loc[idx, "n_genes_cb"].median()),
             "frac_scrublet_auto": float(m.obs.loc[idx, "predicted_doublet_auto"].astype(bool).mean()),
             "frac_second_lineage": float(second),
             "lineages_present": "+".join(lin_frac.index[lin_frac > 0.3]),
             "Vwf": float(det["Vwf"][idx].mean()), "Cntn2": float(det["Cntn2"][idx].mean()),
             "Muc16": float(det["Muc16"][idx].mean())}
        if r["median_genes_cb"] < 200:
            r.update(cell_type="Low_quality", subtype="", decision="remove")
        elif r["frac_scrublet_auto"] > 0.5 and r["frac_second_lineage"] > 0.5:
            r.update(cell_type="Doublet", subtype=r["lineages_present"], decision="remove")
        elif r["Muc16"] > 0.3:
            r.update(cell_type="Epicardial", subtype="Epicardial", decision="keep")
        else:
            sub = t
            if t == "Endothelial":
                sub = "EC_arterial" if r["Vwf"] > 0.15 else "EC_capillary"
            elif t == "Cardiomyocyte":
                sub = "CM_conduction" if r["Cntn2"] > 0.2 else "vCM"
            r.update(cell_type=t, subtype=sub, decision="keep")
        rows[c] = r
    out = pd.DataFrame(rows).T
    out.index.name = key
    return out.loc[sorted(out.index, key=int)]


def detection_table(m, groupby):
    genes = [g for v in MARKERS.values() for g in v if g in m.var_names]
    x = m[:, genes].layers["counts"]
    df = pd.DataFrame((x > 0).toarray(), columns=genes, index=m.obs_names)
    return df.groupby(m.obs[groupby].values).mean()


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    tdir, fdir = res_dir / "tables", res_dir / "figures"
    m = ad.read_h5ad(work_dir / "merged" / "merged_cellbender.h5ad")
    m.obs["cluster_cb"] = m.obs["cluster"].astype(str)
    a = classify_clusters(m)
    print(a.to_string(), flush=True)
    by_s = pd.crosstab(m.obs["cluster_cb"], m.obs["sample"]).reindex(a.index)
    a = a.join(by_s)
    a.to_csv(tdir / "annotation_clusters.csv")
    for k in ["cell_type", "subtype", "decision"]:
        m.obs[k] = m.obs["cluster_cb"].map(a[k]).astype(str)

    # nuclei removed per sample and reason
    rem = pd.crosstab(m.obs["sample"], m.obs["cell_type"].where(m.obs["decision"] == "remove", "kept"))
    print(rem.to_string(), flush=True)

    k = m[m.obs["decision"] == "keep"].copy()
    del m
    k.obs = k.obs.drop(columns=["cluster", "decision"])
    k.X = k.layers["counts"].copy()
    sc.pp.normalize_total(k, target_sum=1e4)
    sc.pp.log1p(k)
    sc.pp.highly_variable_genes(k, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample")
    sc.pp.pca(k, n_comps=50, mask_var="highly_variable", random_state=SEED)
    sc.external.pp.harmony_integrate(k, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     random_state=SEED, max_iter_harmony=20)
    mix = mixing(k, rep="X_pca_harmony")
    sc.pp.neighbors(k, n_neighbors=20, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(k, random_state=SEED)
    sc.tl.leiden(k, resolution=1.0, key_added="cluster", random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)

    # second pass: doublet clusters among the new clusters
    lin = {}
    for name, genes in LINEAGE.items():
        g = [x for x in genes if x in k.var_names]
        lin[name] = np.asarray((k[:, g].layers["counts"] > 0).sum(1)).ravel() >= 2
    lin = pd.DataFrame(lin, index=k.obs_names)
    own = k.obs["cell_type"].map(lambda t: FAMILY.get(t, set()))
    k.obs["second_lineage"] = [any(lin.at[i, l] for l in LINEAGE if l not in own[i]) for i in k.obs_names]
    g2 = k.obs.groupby("cluster", observed=True)
    p2 = pd.DataFrame({"n": g2.size(), "frac_scrublet_auto": g2["predicted_doublet_auto"].mean(),
                       "frac_second_lineage": g2["second_lineage"].mean()})
    p2["remove"] = (p2["frac_scrublet_auto"] > 0.5) & (p2["frac_second_lineage"] > 0.5)
    p2.to_csv(tdir / "annot_second_pass.csv")
    print(p2.round(3).to_string(), flush=True)
    drop = k.obs["cluster"].isin(p2.index[p2["remove"]])
    print("second pass removes", int(drop.sum()), "nuclei", flush=True)
    k = k[~drop.values].copy()

    # re-clustering check: every new cluster should be dominated by one cell type
    o = k.obs
    xt = pd.crosstab(o["cluster"], o["cell_type"])
    rv = pd.DataFrame({"n": xt.sum(1), "main_cell_type": xt.idxmax(1),
                       "main_frac": (xt.max(1) / xt.sum(1)).round(3),
                       "median_umi_cb": o.groupby("cluster", observed=True)["total_cb"].median(),
                       "mean_doublet_score": o.groupby("cluster", observed=True)["doublet_score"].mean().round(3),
                       "top_sample_frac": o.groupby("cluster", observed=True)["sample"].agg(
                           lambda x: x.value_counts(normalize=True).iloc[0]).round(3)})
    rv = rv.join(xt.add_prefix("n_"))
    rv.to_csv(tdir / "annot_recluster_review.csv")
    print(rv.to_string(), "\n", mix.round(2).to_string(), flush=True)

    det = detection_table(k, "subtype")
    det.to_csv(tdir / "annot_marker_detection.csv")

    order = ORDER + sorted(set(o["cell_type"]) - set(ORDER))
    cnt = pd.crosstab(o["sample"], o["cell_type"]).reindex(columns=order, fill_value=0)
    cnt.loc["total"] = cnt.sum()
    cnt.to_csv(tdir / "celltype_counts.csv")
    frac = cnt.div(cnt.sum(1), axis=0).round(4)
    frac.to_csv(tdir / "celltype_fractions.csv")
    print(cnt.to_string(), "\n", frac.round(3).to_string(), flush=True)

    # figures
    um = k.obsm["X_umap"]
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    s05.scatter(axes[0], um, o["cell_type"].values, "Cell type", cats=order)
    s05.scatter(axes[1], um, o["subtype"].values, "Subtype", cats=sorted(o["subtype"].unique()))
    s05.scatter(axes[2], um, o["sample"].values, "Sample", cats=SAMPLES)
    fig.tight_layout()
    fig.savefig(fdir / "annot_umap.png", dpi=150)
    plt.close(fig)

    sub_order = [s for s in ["vCM", "CM_conduction", "Fibroblast", "EC_capillary", "EC_arterial", "Endocardial",
                             "Lymphatic_EC", "Pericyte", "SMC", "Macrophage", "Epicardial"] if s in det.index]
    d = det.loc[sub_order]
    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(d.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(d.shape[1]), d.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(d.shape[0]), d.index)
    plt.colorbar(im, ax=ax, fraction=0.02, label="fraction of nuclei detected (CellBender counts)")
    ax.set_title("Marker detection per annotated subtype")
    fig.tight_layout()
    fig.savefig(fdir / "annot_markers.png", dpi=150)
    plt.close(fig)

    f = frac.drop(index="total")
    sex = pd.read_csv(tdir / "sex_check.csv", index_col=0)["inferred_sex"]
    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(f))
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(order):
        ax.bar([f"{s} ({sex.get(s, '?')})" for s in f.index], f[c].values, bottom=bottom, color=cmap(i), label=c)
        bottom += f[c].values
    ax.set_ylabel("fraction of nuclei")
    ax.set_title("Cell-type composition per sample (inferred sex in brackets)")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(fdir / "annot_composition.png", dpi=150)
    plt.close(fig)

    k.write_h5ad(work_dir / "merged" / "annotated.h5ad", compression="gzip")


if __name__ == "__main__":
    main(*sys.argv[1:3])
