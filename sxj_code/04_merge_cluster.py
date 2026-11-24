"""Step D (workflow sec. 2.5-2.7): merge QC-pass nuclei of all samples, gene filter, normalise, cluster,
check batch effects without integration, then Harmony; marker-panel scores per cluster for the per-cluster
review; sex check from Xist / chrY genes.

Raw integer counts stay in layers["counts"] (and X before normalisation) for downstream DE.

Usage: python 04_merge_cluster.py <work_dir> <results_dir>
Output: <work_dir>/merged/merged.h5ad
        <results_dir>/tables/cluster_review.csv, sex_check.csv, sample_mixing.csv, top_markers.csv
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

SAMPLES = list("ABCDEFGH")
SEED = 0
MIN_CELLS_PER_GENE = 10
PANELS = {
    "Cardiomyocyte": ["Myh6", "Tnnt2", "Ryr2", "Ttn", "Ankrd1", "Mhrt", "Rbm20", "Myl2", "Myl4", "Nppa"],
    "Fibroblast": ["Dcn", "Col1a1", "Col1a2", "Pdgfra", "Gsn", "Lama2", "Abca8a", "Bicc1"],
    "Endothelial": ["Pecam1", "Cdh5", "Flt1", "Egfl7", "Kdr", "Ptprb", "Emcn"],
    "Endocardial": ["Npr3", "Cdh11", "Vwf"],
    "Lymphatic_EC": ["Prox1", "Lyve1", "Mmrn1", "Reln"],
    "Pericyte": ["Rgs5", "Kcnj8", "Abcc9", "Pdgfrb", "Vtn"],
    "SMC": ["Myh11", "Acta2", "Cnn1", "Tagln", "Lmod1"],
    "Macrophage": ["C1qa", "C1qb", "F13a1", "Mrc1", "Cd163", "Adgre1", "Csf1r"],
    "Neutrophil": ["S100a8", "S100a9", "Retnlg"],
    "T_NK": ["Cd3e", "Cd3g", "Skap1", "Nkg7", "Il7r"],
    "B": ["Cd79a", "Ms4a1", "Pax5", "Bank1"],
    "Epicardial": ["Msln", "Upk3b", "Wt1", "Tcf21"],
    "Adipocyte": ["Adipoq", "Plin1", "Lep", "Cidec"],
    "Neural_Schwann": ["Plp1", "Nrxn1", "Scn7a", "Kcna1"],
    "Hepatocyte": ["Alb", "Apoa1", "Cyp2e1"],
    "Kidney_epi": ["Slc34a1", "Lrp2", "Umod", "Nphs1"],
    "Lung_epi": ["Sftpc", "Ager", "Scgb1a1"],
}
SEX_Y = ["Ddx3y", "Eif2s3y", "Kdm5d", "Uty"]


def load(work_dir):
    parts = []
    for s in SAMPLES:
        a = ad.read_h5ad(Path(work_dir) / "per_sample" / f"{s}.h5ad")
        a = a[a.obs["qc_pass"].values].copy()
        a.obs_names = s + "_" + a.obs_names
        parts.append(a)
    m = ad.concat(parts, join="outer", merge="first")
    m.X = sp.csr_matrix(m.X, dtype=np.float32)
    return m


def mixing(m, key="sample", rep="X_pca", k=30):
    """Fraction of each cell's k neighbours from its own sample vs. that sample's share of all cells."""
    from sklearn.neighbors import NearestNeighbors
    rng = np.random.default_rng(SEED)
    idx = rng.choice(m.n_obs, size=min(20000, m.n_obs), replace=False)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(m.obsm[rep])
    _, nb = nn.kneighbors(m.obsm[rep][idx])
    lab = m.obs[key].values
    same = (lab[nb[:, 1:]] == lab[idx][:, None]).mean(1)
    share = m.obs[key].value_counts(normalize=True)
    df = pd.DataFrame({key: lab[idx], "same": same}).groupby(key)["same"].mean().to_frame("own_sample_neighbour_frac")
    df["expected_if_mixed"] = share.reindex(df.index).values
    df["enrichment"] = df["own_sample_neighbour_frac"] / df["expected_if_mixed"]
    return df


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    (work_dir / "merged").mkdir(parents=True, exist_ok=True)
    (res_dir / "tables").mkdir(parents=True, exist_ok=True)

    m = load(work_dir)
    print("merged", m.shape, flush=True)
    n_cells = np.asarray((m.X > 0).sum(0)).ravel()
    m = m[:, n_cells >= MIN_CELLS_PER_GENE].copy()
    print("after gene filter", m.shape, flush=True)
    m.layers["counts"] = m.X.copy()
    sc.pp.normalize_total(m, target_sum=1e4)
    sc.pp.log1p(m)
    sc.pp.highly_variable_genes(m, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample")
    sc.pp.pca(m, n_comps=50, mask_var="highly_variable", random_state=SEED)

    # 1) no integration
    sc.pp.neighbors(m, n_neighbors=20, n_pcs=30, random_state=SEED)
    sc.tl.umap(m, random_state=SEED)
    m.obsm["X_umap_uncorrected"] = m.obsm["X_umap"].copy()
    mix_raw = mixing(m, rep="X_pca")

    # 2) Harmony on sample
    sc.external.pp.harmony_integrate(m, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     random_state=SEED, max_iter_harmony=20)
    mix_h = mixing(m, rep="X_pca_harmony")
    mix = mix_raw.join(mix_h, rsuffix="_harmony")
    mix.to_csv(res_dir / "tables" / "sample_mixing.csv")
    print(mix.round(3).to_string(), flush=True)

    sc.pp.neighbors(m, n_neighbors=20, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(m, random_state=SEED)
    sc.tl.leiden(m, resolution=1.0, key_added="cluster", random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)

    # marker panel scores (mean log-expression of the panel genes present)
    scores = {}
    for name, genes in PANELS.items():
        g = [x for x in genes if x in m.var_names]
        if g:
            scores[name] = np.asarray(m[:, g].X.mean(1)).ravel()
    sdf = pd.DataFrame(scores, index=m.obs_names)
    for c in sdf:
        m.obs[f"score_{c}"] = sdf[c].values

    o = m.obs
    grp = o.groupby("cluster", observed=True)
    rev = pd.DataFrame({
        "n": grp.size(),
        "median_umi": grp["total"].median(),
        "median_genes": grp["n_genes"].median(),
        "median_pct_mt": grp["pct_mt"].median(),
        "median_pct_malat1": grp["pct_malat1"].median(),
        "mean_doublet_score": grp["doublet_score"].mean(),
        "mean_stress": grp["stress_score"].mean(),
        "frac_in_cellranger": grp["in_cellranger"].mean(),
        "top_sample_frac": grp["sample"].agg(lambda x: x.value_counts(normalize=True).iloc[0]),
        "top_sample": grp["sample"].agg(lambda x: x.value_counts().index[0]),
    })
    ps = sdf.groupby(o["cluster"].values).mean()
    rev["panel_1"] = ps.idxmax(1)
    rev["panel_1_score"] = ps.max(1)
    srt = np.sort(ps.values, 1)
    rev["panel_2"] = ps.apply(lambda r: r.sort_values().index[-2], axis=1)
    rev["panel_2_score"] = srt[:, -2]
    rev = rev.join(ps.add_prefix("score_"))

    # top markers per cluster (Wilcoxon on a subsample for speed)
    rng = np.random.default_rng(SEED)
    keep = o.groupby("cluster", observed=True).apply(
        lambda d: d.sample(min(len(d), 800), random_state=SEED)).index.get_level_values(-1)
    sub = m[keep].copy()
    sc.tl.rank_genes_groups(sub, "cluster", method="wilcoxon", n_genes=30, use_raw=False)
    mk = sc.get.rank_genes_groups_df(sub, None)
    mk = mk[mk["logfoldchanges"] > 0.5].groupby("group", observed=True).head(15)
    mk.to_csv(res_dir / "tables" / "top_markers.csv", index=False)
    rev["top_markers"] = mk.groupby("group", observed=True)["names"].agg(lambda x: ",".join(x[:10])).reindex(rev.index)
    rev.to_csv(res_dir / "tables" / "cluster_review.csv")
    print(rev[["n", "median_umi", "median_genes", "median_pct_mt", "mean_doublet_score", "top_sample_frac",
               "panel_1", "panel_1_score", "panel_2", "panel_2_score", "top_markers"]].round(3).to_string(),
          flush=True)

    # sex check: pseudobulk CPM per sample from raw counts
    cnt = m.layers["counts"]
    rows = []
    for s in SAMPLES:
        mask = (o["sample"] == s).values
        tot = cnt[mask].sum()
        r = {"sample": s}
        for g in ["Xist"] + SEX_Y:
            if g in m.var_names:
                j = m.var_names.get_loc(g)
                col = cnt[mask][:, j]
                r[f"{g}_cpm"] = float(col.sum() / tot * 1e6)
                r[f"{g}_frac_nuclei"] = float((col > 0).mean())
        rows.append(r)
    sex = pd.DataFrame(rows).set_index("sample")
    ycols = [f"{g}_cpm" for g in SEX_Y if f"{g}_cpm" in sex]
    sex["chrY_cpm_sum"] = sex[ycols].sum(1)
    sex["inferred_sex"] = np.where((sex["Xist_cpm"] > 10 * sex["chrY_cpm_sum"]), "F",
                                   np.where(sex["chrY_cpm_sum"] > 10 * sex["Xist_cpm"], "M", "ambiguous"))
    sex.to_csv(res_dir / "tables" / "sex_check.csv")
    print(sex.round(2).to_string(), flush=True)

    m.write_h5ad(work_dir / "merged" / "merged.h5ad", compression="gzip")


if __name__ == "__main__":
    main(*sys.argv[1:3])
