"""Step F: re-cluster the QC-pass nuclei on CellBender-corrected counts and compare with the raw-count clusters.

The nuclei are exactly those of 04_merge_cluster.py (QC-pass from 02_per_sample_qc.py); only their counts are
replaced by CellBender's denoised counts (<work_dir>/cellbender/<S>/<S>_cellbender.h5). Two groups are dropped and
counted: nuclei CellBender did not analyse (outside its droplet window, ~0.5%), and nuclei CellBender calls
non-cells (cell probability < 0.5, ~1.4%). CellBender sets all counts of its non-cells to zero; a first run that
kept them produced ~150 tiny all-zero clusters, mostly from raw-count cluster 11 (low-UMI nuclei near the ambient
plateau). The pipeline is the same as step 04 (gene
filter, log-normalise, 3000 HVGs by sample, 50 PCs, Harmony on sample, Leiden 1.0), so differences come from the
counts, not from settings.

Usage: python 06_cellbender_recluster.py <work_dir> <results_dir>
Output: <work_dir>/merged/merged_cellbender.h5ad
        <results_dir>/tables/cb_cluster_review.csv, cb_top_markers.csv, cb_sample_mixing.csv,
            cb_vs_raw_clusters.csv (cross-table), cb_ambient_check.csv, cb_summary.csv
        <results_dir>/figures/cb_umap.png, cb_vs_raw_ambient.png
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

_spec = importlib.util.spec_from_file_location("step04", Path(__file__).with_name("04_merge_cluster.py"))
step04 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step04)
SAMPLES, SEED, PANELS, MIN_CELLS_PER_GENE = step04.SAMPLES, step04.SEED, step04.PANELS, step04.MIN_CELLS_PER_GENE
mixing = step04.mixing
CM_GENES = ["Ttn", "Myh6", "Ryr2", "Tnnt2", "Mhrt", "Myl2", "Nppa"]


def load(work_dir):
    """CellBender counts for the QC-pass nuclei of each sample, with the step-02/04 annotations."""
    raw = ad.read_h5ad(work_dir / "merged" / "merged.h5ad", backed="r")
    obs_raw = raw.obs.copy()
    raw.file.close()
    parts, missing, noncell = [], {}, {}
    for s in SAMPLES:
        f = work_dir / "cellbender" / s / f"{s}_cellbender.h5"
        a = sc.read_10x_h5(f)
        a.var_names_make_unique()
        with h5py.File(f, "r") as h:
            bc = h["matrix/barcodes"][:].astype(str)
            idx = h["droplet_latents/barcode_indices_for_latents"][:]
            prob = pd.Series(h["droplet_latents/cell_probability"][:], index=bc[idx])
        want = obs_raw.index[obs_raw["sample"] == s].str.slice(2)  # "A_<barcode>" -> "<barcode>"
        have = want[want.isin(prob.index) & want.isin(a.obs_names)]
        missing[s] = len(want) - len(have)
        p = prob.reindex(have).values
        noncell[s] = int((p < 0.5).sum())
        a = a[have[p >= 0.5]].copy()
        a.obs["cb_cell_prob"] = prob.reindex(a.obs_names).values
        a.obs_names = s + "_" + a.obs_names
        parts.append(a)
        print(s, a.shape, "not in CellBender window:", missing[s], "CellBender non-cells:", noncell[s], flush=True)
    m = ad.concat(parts, join="outer", merge="first")
    m.X = sp.csr_matrix(m.X, dtype=np.float32)
    keep = ["sample", "total", "n_genes", "pct_mt", "pct_malat1", "doublet_score", "predicted_doublet_auto",
            "stress_score", "in_cellranger", "cluster"]
    o = obs_raw.loc[m.obs_names, [c for c in keep if c in obs_raw]].rename(
        columns={"cluster": "cluster_raw", "total": "total_raw", "n_genes": "n_genes_raw"})
    m.obs = m.obs[["cb_cell_prob"]].join(o)
    return m, missing, noncell


def ambient_check(work_dir, m, label):
    """CPM of cardiomyocyte genes per cluster label, raw vs CellBender counts, on the same nuclei."""
    genes = [g for g in CM_GENES if g in m.var_names]
    rows = []
    raw_parts = []
    for s in SAMPLES:
        r = ad.read_h5ad(work_dir / "per_sample" / f"{s}.h5ad", backed="r")
        names = m.obs_names[m.obs["sample"] == s].str.slice(2)
        x = r[names.values].to_memory()
        tot = np.asarray(x.X.sum(1)).ravel()
        cm = np.asarray(x[:, genes].X.sum(1)).ravel()
        r.file.close()
        raw_parts.append(pd.DataFrame({"tot": tot, "cm": cm}, index=s + "_" + names))
    rawdf = pd.concat(raw_parts).loc[m.obs_names]
    cnt = m.layers["counts"]
    cbdf = pd.DataFrame({"tot": np.asarray(cnt.sum(1)).ravel(),
                         "cm": np.asarray(cnt[:, [m.var_names.get_loc(g) for g in genes]].sum(1)).ravel()},
                        index=m.obs_names)
    for lab, idx in m.obs.groupby(label, observed=True).groups.items():
        rr, cc = rawdf.loc[idx], cbdf.loc[idx]
        rows.append({label: lab, "n": len(idx),
                     "cm_genes_cpm_raw": rr["cm"].sum() / rr["tot"].sum() * 1e6,
                     "cm_genes_cpm_cellbender": cc["cm"].sum() / cc["tot"].sum() * 1e6,
                     "umi_removed_frac": 1 - cc["tot"].sum() / rr["tot"].sum()})
    df = pd.DataFrame(rows).set_index(label)
    df["cm_genes_fold_reduction"] = df["cm_genes_cpm_raw"] / df["cm_genes_cpm_cellbender"].clip(lower=1e-9)
    return df


def review(m, res_dir):
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
        "median_umi_cb": grp["total_cb"].median(),
        "median_umi_raw": grp["total_raw"].median(),
        "median_genes_cb": grp["n_genes_cb"].median(),
        "median_pct_mt": grp["pct_mt"].median(),
        "mean_doublet_score": grp["doublet_score"].mean(),
        "frac_doublet_auto": grp["predicted_doublet_auto"].mean(),
        "median_cb_cell_prob": grp["cb_cell_prob"].median(),
        "top_sample_frac": grp["sample"].agg(lambda x: x.value_counts(normalize=True).iloc[0]),
        "top_sample": grp["sample"].agg(lambda x: x.value_counts().index[0]),
        "main_raw_cluster": grp["cluster_raw"].agg(lambda x: x.value_counts().index[0]),
        "main_raw_cluster_frac": grp["cluster_raw"].agg(lambda x: x.value_counts(normalize=True).iloc[0]),
    })
    ps = sdf.groupby(o["cluster"].values).mean()
    srt = np.sort(ps.values, 1)
    rev["panel_1"] = ps.idxmax(1)
    rev["panel_1_score"] = srt[:, -1]
    rev["panel_2"] = ps.apply(lambda r: r.sort_values().index[-2], axis=1)
    rev["panel_2_score"] = srt[:, -2]
    rev["panel_margin"] = rev["panel_1_score"] / np.maximum(rev["panel_2_score"], 1e-6)
    rev = rev.join(ps.add_prefix("score_"))

    keep = o.groupby("cluster", observed=True).apply(
        lambda d: d.sample(min(len(d), 800), random_state=SEED)).index.get_level_values(-1)
    sub = m[keep].copy()
    sc.tl.rank_genes_groups(sub, "cluster", method="wilcoxon", n_genes=30, use_raw=False)
    mk = sc.get.rank_genes_groups_df(sub, None)
    mk = mk[mk["logfoldchanges"] > 0.5].groupby("group", observed=True).head(15)
    mk.to_csv(res_dir / "tables" / "cb_top_markers.csv", index=False)
    rev["top_markers"] = mk.groupby("group", observed=True)["names"].agg(
        lambda x: ",".join(x[:10])).reindex(rev.index)
    return rev


def figures(m, fdir):
    import importlib.util as iu
    spec = iu.spec_from_file_location("step05", Path(__file__).with_name("05_qc_report.py"))
    s05 = iu.module_from_spec(spec)
    spec.loader.exec_module(s05)
    o = m.obs
    um = m.obsm["X_umap"]
    clusters = sorted(o["cluster"].astype(str).unique(), key=int)
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    s05.scatter(axes[0, 0], um, o["cluster"].astype(str).values, "Leiden clusters (CellBender counts, Harmony)",
                cats=clusters)
    s05.scatter(axes[0, 1], um, o["panel_1"].astype(str).values, "Top marker panel per cluster",
                cats=sorted(o["panel_1"].astype(str).unique()))
    s05.scatter(axes[1, 0], um, o["sample"].astype(str).values, "Sample", cats=SAMPLES)
    s05.scatter(axes[1, 1], um, np.log10(o["total_cb"].values + 1), "log10 UMI (CellBender)", cmap="viridis")
    fig.tight_layout()
    fig.savefig(fdir / "cb_umap.png", dpi=150)
    plt.close(fig)


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    tdir, fdir = res_dir / "tables", res_dir / "figures"
    tdir.mkdir(parents=True, exist_ok=True)
    fdir.mkdir(parents=True, exist_ok=True)

    m, missing, noncell = load(work_dir)
    print("merged", m.shape, flush=True)
    m.obs["total_cb"] = np.asarray(m.X.sum(1)).ravel()
    m.obs["n_genes_cb"] = np.asarray((m.X > 0).sum(1)).ravel()
    n_cells = np.asarray((m.X > 0).sum(0)).ravel()
    m = m[:, n_cells >= MIN_CELLS_PER_GENE].copy()
    print("after gene filter", m.shape, flush=True)
    m.layers["counts"] = m.X.copy()
    sc.pp.normalize_total(m, target_sum=1e4)
    sc.pp.log1p(m)
    sc.pp.highly_variable_genes(m, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample")
    sc.pp.pca(m, n_comps=50, mask_var="highly_variable", random_state=SEED)
    mix_raw = mixing(m, rep="X_pca")
    sc.external.pp.harmony_integrate(m, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     random_state=SEED, max_iter_harmony=20)
    mix = mix_raw.join(mixing(m, rep="X_pca_harmony"), rsuffix="_harmony")
    mix.to_csv(tdir / "cb_sample_mixing.csv")
    print(mix.round(3).to_string(), flush=True)
    sc.pp.neighbors(m, n_neighbors=20, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(m, random_state=SEED)
    sc.tl.leiden(m, resolution=1.0, key_added="cluster", random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)

    rev = review(m, res_dir)
    m.obs["panel_1"] = m.obs["cluster"].map(rev["panel_1"]).astype(str)
    rev.to_csv(tdir / "cb_cluster_review.csv")
    print(rev[["n", "median_umi_cb", "median_umi_raw", "median_pct_mt", "frac_doublet_auto",
               "median_cb_cell_prob", "top_sample_frac", "main_raw_cluster", "main_raw_cluster_frac",
               "panel_1", "panel_1_score", "panel_2", "panel_margin", "top_markers"]].round(3).to_string(),
          flush=True)

    xt = pd.crosstab(m.obs["cluster"], m.obs["cluster_raw"])
    xt.to_csv(tdir / "cb_vs_raw_clusters.csv")
    ari = adjusted_rand_score(m.obs["cluster"], m.obs["cluster_raw"])

    amb = ambient_check(work_dir, m, "panel_1")
    amb.to_csv(tdir / "cb_ambient_check.csv")
    print(amb.round(2).to_string(), flush=True)

    summ = pd.DataFrame({
        "n_qc_pass": m.obs["sample"].value_counts().reindex(SAMPLES) + pd.Series(missing) + pd.Series(noncell),
        "n_not_in_cellbender_window": pd.Series(missing),
        "n_cellbender_noncell": pd.Series(noncell),
        "n_used": m.obs["sample"].value_counts().reindex(SAMPLES),
        "median_umi_raw": m.obs.groupby("sample", observed=True)["total_raw"].median(),
        "median_umi_cellbender": m.obs.groupby("sample", observed=True)["total_cb"].median(),
        "frac_cb_cell_prob_gt_0.9": m.obs.groupby("sample", observed=True)["cb_cell_prob"].apply(
            lambda x: (x > 0.9).mean()),
    })
    summ.loc["all", "n_used"] = m.n_obs
    summ.loc["all", "ari_vs_raw_clusters"] = ari
    summ.to_csv(tdir / "cb_summary.csv")
    print(summ.round(3).to_string(), "\nARI vs raw-count clusters:", round(ari, 3), flush=True)

    figures(m, fdir)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    a = amb.sort_values("cm_genes_cpm_raw")
    y = np.arange(len(a))
    ax.barh(y - 0.2, a["cm_genes_cpm_raw"], 0.4, label="raw counts")
    ax.barh(y + 0.2, a["cm_genes_cpm_cellbender"], 0.4, label="CellBender")
    ax.set_yticks(y, a.index)
    ax.set_xscale("log")
    ax.set_xlabel("cardiomyocyte genes (" + ",".join(CM_GENES) + "), CPM")
    ax.set_title("Cardiomyocyte transcripts in each cell population, same nuclei")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fdir / "cb_vs_raw_ambient.png", dpi=150)
    plt.close(fig)

    m.write_h5ad(work_dir / "merged" / "merged_cellbender.h5ad", compression="gzip")


if __name__ == "__main__":
    main(*sys.argv[1:3])
