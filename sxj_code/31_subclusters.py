"""Subclusters of cardiomyocytes, endocardial cells, fibroblasts and endothelial cells (WT vs TNNC1-I4M).

Per cell type (nuclei with flag == ok):
  - log-normalised CellBender counts (10,000 per nucleus);
  - 2,000 highly variable genes (seurat_v3, per mouse), excluding mitochondrial / ribosomal genes, the ambient /
    low-quality signature and ambient genes (>= 5x higher mean expression in another cell type), so that the
    subclusters follow biology rather than library quality or contamination;
  - log(total counts) and the ambient / low-quality fraction regressed out of the HVG matrix, scaling;
  - PCA (30 PCs), Harmony by mouse, 15-NN graph, Leiden, UMAP.
    Note: Harmony by mouse also removes some genotype signal, because genotype = batch here; a state present only in
    I4M mice could be partly merged. Composition per mouse is therefore tested afterwards, and per-cluster quality
    metrics are kept so that quality-driven clusters can be recognised.
  - markers: Wilcoxon, each cluster vs the rest of that cell type.

Outputs (local): <work_dir>/subclusters/sub_<type>.csv.gz      per-nucleus cluster + UMAP + QC
                 <work_dir>/subclusters/markers_<type>.csv.gz  full marker tables
Usage: python 31_subclusters.py <work_dir> [cell_type ...]
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
RES = {"Cardiomyocyte": 0.8, "Endocardial": 0.5, "Fibroblast": 0.5, "Endothelial": 0.5}
EXTRA_RES = [0.3, 0.5, 0.8, 1.2]  # also stored, for review of the cluster granularity
TECH = ["Gm42418", "Cmss1", "Lars2", "Camk1d", "Gphn", "Gm19951", "Cdk8", "AY036118", "Gm26917", "Malat1"]


def is_tech(g):
    g = pd.Index(g)
    return g.isin(TECH) | g.str.startswith("mt-") | g.str.match(r"^Rp[ls]\d")


def run(a, ct, ambient, outd):
    b = a[(a.obs["cell_type_final"].astype(str) == ct).values].copy()
    b.X = b.layers["counts"].copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    try:
        sc.pp.highly_variable_genes(b, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample")
    except ValueError:  # loess fails for small per-mouse batches; use a wider span
        sc.pp.highly_variable_genes(b, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample",
                                    span=0.5)
    bad = is_tech(b.var_names) | b.var_names.isin(ambient[ct])
    hv = b.var["highly_variable_rank"].where(~bad).dropna().sort_values().index[:2000]
    b.var["highly_variable"] = b.var_names.isin(hv)
    # remove library-size / ambient-load variation before PCA (otherwise the first split is high- vs low-RNA nuclei)
    X = b[:, hv].X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float32)
    tot = b.obs["total_cb"].values.astype(float)
    amb = np.asarray(b.layers["counts"][:, b.var_names.isin(TECH)].sum(1)).ravel() / tot
    C = np.column_stack([np.ones(len(tot)), np.log1p(tot), amb])
    X -= (C @ np.linalg.lstsq(C, X, rcond=None)[0]).astype(np.float32)
    X = ((X - X.mean(0)) / X.std(0).clip(1e-6)).clip(-10, 10)
    h = ad.AnnData(X)
    sc.pp.pca(h, n_comps=30, random_state=SEED)
    b.obsm["X_pca"] = h.obsm["X_pca"]
    del X, h
    sc.external.pp.harmony_integrate(b, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     random_state=SEED, max_iter_harmony=20)
    sc.pp.neighbors(b, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.leiden(b, resolution=RES[ct], key_added="sub", random_state=SEED, flavor="igraph", n_iterations=2)
    for r in EXTRA_RES:
        sc.tl.leiden(b, resolution=r, key_added=f"r{r}", random_state=SEED, flavor="igraph", n_iterations=2)
    sc.tl.umap(b, random_state=SEED)
    sc.tl.rank_genes_groups(b, "sub", method="wilcoxon", pts=True)
    m = sc.get.rank_genes_groups_df(b, None)
    m = m[m.logfoldchanges > 0]
    m.to_csv(outd / f"markers_{ct}.csv.gz", index=False)
    techsig = b[:, [g for g in TECH if g in b.var_names]].X.sum(1)
    o = pd.DataFrame({"sample": b.obs["sample"].astype(str).values, "sub": b.obs["sub"].astype(str).values,
                      "umap1": b.obsm["X_umap"][:, 0], "umap2": b.obsm["X_umap"][:, 1],
                      "n_counts": b.obs["total_cb"].values, "n_genes": b.obs["n_genes_cb"].values,
                      "pct_mt": b.obs["pct_mt"].values, "tech_score": np.asarray(techsig).ravel()},
                     index=b.obs_names)
    for r in EXTRA_RES:
        o[f"r{r}"] = b.obs[f"r{r}"].astype(str).values
    o.to_csv(outd / f"sub_{ct}.csv.gz")
    q = o.groupby("sub").agg(n=("sample", "size"), med_genes=("n_genes", "median"), med_tech=("tech_score", "median"))
    comp = pd.crosstab(o["sub"], o["sample"], normalize="columns").round(3)
    top = m.groupby("group").head(12).groupby("group")["names"].apply(lambda x: ",".join(x))
    print(f"==== {ct}: {b.n_obs} nuclei, {o['sub'].nunique()} clusters", flush=True)
    print(pd.concat([q, comp, top.rename("top")], axis=1).to_string(), flush=True)
    for r in EXTRA_RES:
        g = o.groupby(f"r{r}").agg(n=("sample", "size"), med_genes=("n_genes", "median"))
        print(f"  res {r}:", "; ".join(f"{k}:{v.n}/{int(v.med_genes)}" for k, v in g.iterrows()), flush=True)


def main(work_dir, *types):
    work_dir = Path(work_dir)
    outd = work_dir / "subclusters"
    outd.mkdir(exist_ok=True)
    types = list(types) or list(RES)
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    a = a[(a.obs["flag"] == "ok").values].to_memory()
    # ambient genes per cell type: >= 5x higher mean normalised expression in another major type
    major = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "SMC",
             "Lymphatic_EC"]
    C = a.layers["counts"]
    lib = np.asarray(C.sum(1)).ravel()
    means = {}
    for t in major:
        idx = np.where((a.obs["cell_type_final"].astype(str) == t).values)[0]
        x = C[idx].multiply(1e4 / lib[idx][:, None]).tocsr()
        means[t] = np.asarray(x.mean(0)).ravel()
    M = pd.DataFrame(means, index=a.var_names)
    ambient = {t: M.index[(M.drop(columns=t).max(1) >= 5 * M[t]) & (M.drop(columns=t).max(1) > 0.05)]
               for t in RES}
    for t in types:
        run(a, t, ambient, outd)


if __name__ == "__main__":
    main(*sys.argv[1:])
