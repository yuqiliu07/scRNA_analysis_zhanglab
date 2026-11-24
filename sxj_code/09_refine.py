"""Step I: two refinements of the step-07 annotation.

0) Lymphocytes are identified per nucleus (see LYMPH_SPECIFIC), including those in removed clusters.
1) Per-nucleus rescue inside the removed doublet clusters. Step 07 removed whole clusters (06 clusters 7, 10, 14,
   15, 17, 19 and one second-pass cluster), which also removes the nuclei in them that show only one lineage. A
   nucleus there is a rescue candidate when it carries exactly one lineage (>= 2 genes of one LINEAGE panel, none of
   the others) and is not a Scrublet-auto doublet. All kept + candidate nuclei are re-embedded (HVG, PCA, Harmony);
   each candidate gets the majority cell type of its 15 nearest kept neighbours, and is kept only when that type
   belongs to its own lineage and >= 70% of the neighbours agree.
2) Subclustering of each compartment (cardiomyocytes, fibroblasts, endothelial/endocardial/lymphatic/epicardial,
   mural, immune) to find populations the global clustering merged (the provider resolved lymphocytes, dendritic
   cells, glia and proliferating cells). Each subcluster gets markers and the fraction of nuclei carrying
   signatures of minor populations; naming is done in 10_name_subclusters.py after review.

Usage: python 09_refine.py <work_dir> <results_dir> [provider_dir]
Output: <work_dir>/merged/refined.h5ad
        <results_dir>/tables/refine_rescue.csv, refine_subclusters.csv, refine_submarkers.csv
        <results_dir>/figures/refine_umap.png, refine_sub_<compartment>.png
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
from sklearn.neighbors import NearestNeighbors

_here = Path(__file__).parent


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, _here / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s04 = _load("step04", "04_merge_cluster.py")
s05 = _load("step05", "05_qc_report.py")
s07 = _load("step07", "07_annotate.py")
SAMPLES, SEED = s04.SAMPLES, s04.SEED
LINEAGE = s07.LINEAGE
LIN_TYPES = {"CM": {"Cardiomyocyte"}, "Fib": {"Fibroblast"}, "EC": {"Endothelial", "Endocardial", "Lymphatic_EC"},
             "Peri": {"Pericyte"}, "Mac": {"Macrophage"}, "SMC": {"SMC"}}
# Lymphocytes: small nuclei with little RNA; most fall into the 06 low-quality cluster or scatter across other
# clusters, so they are identified per nucleus: >= 2 lymphocyte-specific genes (each detected in < 0.5% of
# nuclei overall and of provider macrophage cluster 5, 2-37% of provider lymphocyte cluster 18) plus >= 1
# haematopoietic gene, no non-CM lineage panel (CM transcripts are mostly ambient), not a Scrublet-auto doublet.
LYMPH_SPECIFIC = ["Pax5", "Skap1", "Il7r", "Cd2", "Ms4a1", "Cd19", "Cd79a", "Cd79b", "Cd3e", "Cd3g", "Lck", "Themis",
                  "Ikzf3", "Ms4a4b"]
HAEM = ["Ptprc", "Arhgap15", "Dock2", "Ikzf1"]
COMPARTMENTS = {
    "CM": (["Cardiomyocyte"], 0.4),
    "Fib": (["Fibroblast"], 0.4),
    "EC": (["Endothelial", "Endocardial", "Lymphatic_EC", "Epicardial"], 0.4),
    "Mural": (["Pericyte", "SMC"], 0.4),
    "Immune": (["Macrophage", "Lymphocyte"], 0.8),
}
SIGNATURES = {  # minor populations; >= 2 genes detected = signature present
    "Lymphoid": ["Ptprc", "Bach2", "Aff3", "Arhgap15", "Cd79a", "Ms4a1", "Pax5", "Cd3e", "Skap1", "Il7r", "Nkg7",
                 "Ripor2", "Themis", "Ets1"],
    "B": ["Cd79a", "Ms4a1", "Pax5", "Bank1", "Ebf1", "Cd19"],
    "T_NK": ["Cd3e", "Cd3g", "Skap1", "Themis", "Il7r", "Nkg7", "Klrd1", "Cd247"],
    "DC": ["Wdfy4", "Xcr1", "Flt3", "Clec9a", "Cd209a", "Ccr7", "Cadm1"],
    "Mono": ["Ly6c2", "Ccr2", "Plac8", "Chil3"],
    "Neutrophil": ["S100a8", "S100a9", "Retnlg", "Csf3r"],
    "Mast": ["Cpa3", "Kit", "Cma1", "Mcpt4"],
    "Glia": ["Cdh19", "Plp1", "Sox10", "Nkain2", "Scn7a", "Gfra3", "Mpz"],
    "Prolif": ["Mki67", "Top2a", "Diaph3", "Cenpf", "Hmmr"],
    "Adipo": ["Adipoq", "Plin1", "Cidec", "Lep"],
}


def lineage_matrix(a):
    out = {}
    for name, genes in LINEAGE.items():
        g = [x for x in genes if x in a.var_names]
        out[name] = np.asarray((a[:, g].layers["counts"] > 0).sum(1)).ravel() >= 2
    return pd.DataFrame(out, index=a.obs_names)


def signature_matrix(a):
    out = {}
    for name, genes in SIGNATURES.items():
        g = [x for x in genes if x in a.var_names]
        out[name] = np.asarray((a[:, g].layers["counts"] > 0).sum(1)).ravel() >= 2
    return pd.DataFrame(out, index=a.obs_names)


def embed(a, n_hvg=3000, n_pcs=50, res=None, key="cluster"):
    a.X = a.layers["counts"].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts", batch_key="sample")
    except ValueError as e:  # loess fit can fail when a sample has few nuclei in a small compartment
        print("per-sample seurat_v3 HVG failed (%s); using pooled HVG" % e, flush=True)
        try:
            sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts", span=0.5)
        except ValueError:
            sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor="seurat", batch_key="sample")
    sc.pp.pca(a, n_comps=n_pcs, mask_var="highly_variable", random_state=SEED)
    sc.external.pp.harmony_integrate(a, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     random_state=SEED, max_iter_harmony=20)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    if res is not None:
        sc.tl.leiden(a, resolution=res, key_added=key, random_state=SEED, flavor="igraph", n_iterations=2,
                     directed=False)


def main(work_dir, res_dir, provider_dir=None):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    tdir, fdir = res_dir / "tables", res_dir / "figures"

    # ---------- 1) rescue ----------
    m = ad.read_h5ad(work_dir / "merged" / "merged_cellbender.h5ad")
    m.obs["cluster_cb"] = m.obs["cluster"].astype(str)
    ann = pd.read_csv(tdir / "annotation_clusters.csv", dtype={"cluster_cb": str}).set_index("cluster_cb")
    kept = ad.read_h5ad(work_dir / "merged" / "annotated.h5ad", backed="r")
    kobs = kept.obs[["cell_type", "subtype"]].copy()
    kept.file.close()
    m.obs["cell_type"] = kobs["cell_type"].reindex(m.obs_names).astype(object)
    m.obs["subtype"] = kobs["subtype"].reindex(m.obs_names).astype(object)
    ct06 = m.obs["cluster_cb"].map(ann["cell_type"])
    pool = m.obs_names[(ct06 == "Doublet").values | ((ct06 != "Doublet") & (ct06 != "Low_quality")
                                                    & m.obs["cell_type"].isna()).values]
    lin = lineage_matrix(m[pool])
    n_lin = lin.sum(1)
    single = lin.idxmax(1).where(n_lin == 1)
    cand = pool[(n_lin == 1).values & ~m.obs.loc[pool, "predicted_doublet_auto"].astype(bool).values]
    print("doublet-cluster nuclei", len(pool), "one lineage", int((n_lin == 1).sum()), "candidates", len(cand),
          flush=True)

    # lymphocytes, per nucleus, from all CellBender nuclei (kept, low-quality or doublet clusters)
    def n_det(genes):
        g = [x for x in genes if x in m.var_names]
        return np.asarray((m[:, g].layers["counts"] > 0).sum(1)).ravel()
    lin_all = lineage_matrix(m).drop(columns=["CM"])
    lymph = ((n_det(LYMPH_SPECIFIC) >= 2) & (n_det(HAEM) >= 1) & ~lin_all.any(axis=1).values
             & ~m.obs["predicted_doublet_auto"].astype(bool).values)
    ly = pd.DataFrame({"sample": m.obs["sample"].values, "cluster_cb": m.obs["cluster_cb"].values,
                       "previous": m.obs["cell_type"].fillna("removed_" + ct06.astype(str)).values},
                      index=m.obs_names)[lymph]
    ly.to_csv(tdir / "refine_lymphocytes.csv.gz")
    print("lymphocytes per nucleus", int(lymph.sum()), flush=True)
    print(ly["previous"].value_counts().to_string(), flush=True)
    m.obs.loc[lymph, "cell_type"] = "Lymphocyte"
    m.obs.loc[lymph, "subtype"] = "Lymphocyte"
    cand = cand[~cand.isin(ly.index)]

    keep_mask = m.obs["cell_type"].notna().values | m.obs_names.isin(cand)
    a = m[keep_mask].copy()
    del m
    a.obs = a.obs.drop(columns=["cluster", "panel_1"], errors="ignore")
    embed(a)
    ref = ~a.obs_names.isin(cand)
    emb = a.obsm["X_pca_harmony"]
    nn = NearestNeighbors(n_neighbors=15).fit(emb[ref])
    _, nb = nn.kneighbors(emb[~ref])
    ref_types = a.obs["cell_type"].values[ref]
    ref_sub = a.obs["subtype"].values[ref]
    qidx = a.obs_names[~ref]
    rows = []
    for i, q in enumerate(qidx):
        tt = pd.Series(ref_types[nb[i]]).value_counts()
        st = pd.Series(ref_sub[nb[i]]).value_counts()
        lin_q = single[q]
        ok = (tt.index[0] in LIN_TYPES[lin_q]) and (tt.iloc[0] / 15 >= 0.7)
        rows.append({"barcode": q, "sample": a.obs.at[q, "sample"], "cluster_cb": a.obs.at[q, "cluster_cb"],
                     "lineage": lin_q, "knn_type": tt.index[0], "knn_frac": tt.iloc[0] / 15,
                     "knn_subtype": st.index[0], "rescued": ok})
    rs = pd.DataFrame(rows).set_index("barcode")
    rs.to_csv(tdir / "refine_rescue_cells.csv.gz")
    summ = rs.groupby(["cluster_cb", "lineage"]).agg(candidates=("rescued", "size"), rescued=("rescued", "sum"))
    summ.to_csv(tdir / "refine_rescue.csv")
    print(summ.to_string(), "\nrescued total", int(rs["rescued"].sum()), flush=True)

    a.obs.loc[rs.index[rs["rescued"]], "cell_type"] = rs.loc[rs["rescued"], "knn_type"].values
    a.obs.loc[rs.index[rs["rescued"]], "subtype"] = rs.loc[rs["rescued"], "knn_subtype"].values
    a.obs["rescued"] = a.obs_names.isin(rs.index[rs["rescued"]]) | a.obs_names.isin(ly.index[ly["previous"].str.startswith("removed_")])
    a = a[a.obs["cell_type"].notna().values].copy()
    a.obs["cell_type"] = a.obs["cell_type"].astype(str)
    a.obs["subtype"] = a.obs["subtype"].astype(str)
    print("nuclei after rescue", a.n_obs, flush=True)
    # re-embed the final set (UMAP only for display)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(a, random_state=SEED)

    prov = None
    if provider_dir:
        c = pd.read_csv(Path(provider_dir) / "3.Identification_subpopulation/Clustering/Combination1.clusters_info.xls",
                        sep="\t")
        c.index = c["orig.ident"] + "_" + c["barcode"].str.replace(r"_\d+$", "", regex=True)
        prov = c["seurat_clusters"]
        a.obs["provider_cluster"] = prov.reindex(a.obs_names).astype("Int64").astype(str)

    # ---------- 2) subclustering ----------
    a.obs["subcluster"] = ""
    sig_all = signature_matrix(a)
    rows, mk_all = [], []
    for comp, (types, res) in COMPARTMENTS.items():
        s = a[a.obs["cell_type"].isin(types).values].copy()
        print(comp, s.shape, flush=True)
        embed(s, n_hvg=2000, n_pcs=30, res=res, key="sub")
        s.obs["sub"] = comp + "_" + s.obs["sub"].astype(str)
        a.obs.loc[s.obs_names, "subcluster"] = s.obs["sub"].values
        keep = s.obs.groupby("sub", observed=True).apply(
            lambda d: d.sample(min(len(d), 1000), random_state=SEED)).index.get_level_values(-1)
        t = s[keep].copy()
        sc.tl.rank_genes_groups(t, "sub", method="wilcoxon", n_genes=40, use_raw=False)
        mk = sc.get.rank_genes_groups_df(t, None)
        mk = mk[mk["logfoldchanges"] > 0.5].groupby("group", observed=True).head(20)
        mk_all.append(mk)
        g = s.obs.groupby("sub", observed=True)
        sig = sig_all.loc[s.obs_names].groupby(s.obs["sub"].values).mean()
        info = pd.DataFrame({
            "compartment": comp, "n": g.size(), "frac_of_compartment": g.size() / s.n_obs,
            "main_type": g["cell_type"].agg(lambda x: x.value_counts().index[0]),
            "main_subtype": g["subtype"].agg(lambda x: x.value_counts().index[0]),
            "median_umi_cb": g["total_cb"].median(), "median_pct_mt": g["pct_mt"].median(),
            "frac_scrublet_auto": g["predicted_doublet_auto"].mean(), "frac_rescued": g["rescued"].mean(),
            "top_sample_frac": g["sample"].agg(lambda x: x.value_counts(normalize=True).iloc[0]),
            "top_markers": mk.groupby("group", observed=True)["names"].agg(lambda x: ",".join(x[:12])),
        })
        if prov is not None:
            info["provider_top"] = g["provider_cluster"].agg(
                lambda x: ";".join(f"{k}:{v:.2f}" for k, v in x[x != "<NA>"].value_counts(normalize=True).head(3).items()))
        info = info.join(sig.add_prefix("sig_"))
        rows.append(info)
        # compartment figure
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        subs = sorted(s.obs["sub"].unique(), key=lambda x: int(x.split("_")[-1]))
        s05.scatter(axes[0], s.obsm["X_umap"], s.obs["sub"].values, f"{comp} subclusters", cats=subs)
        um = s.obsm["X_umap"]
        for k in subs:
            mm = (s.obs["sub"] == k).values
            axes[0].text(*np.median(um[mm], 0), k.split("_")[-1], fontsize=10, weight="bold", ha="center")
        s05.scatter(axes[1], um, s.obs["sample"].values, "Sample", cats=SAMPLES)
        fig.tight_layout()
        fig.savefig(fdir / f"refine_sub_{comp}.png", dpi=130)
        plt.close(fig)
        del s, t

    sub = pd.concat(rows)
    sub.index.name = "subcluster"
    sub.to_csv(tdir / "refine_subclusters.csv")
    pd.concat(mk_all).to_csv(tdir / "refine_submarkers.csv", index=False)
    pd.set_option("display.width", 400)
    pd.set_option("display.max_colwidth", 120)
    print(sub.round(3).to_string(), flush=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    types = s07.ORDER + sorted(set(a.obs["cell_type"]) - set(s07.ORDER))
    s05.scatter(ax, a.obsm["X_umap"], a.obs["cell_type"].values, "Cell type after rescue", cats=types)
    fig.tight_layout()
    fig.savefig(fdir / "refine_umap.png", dpi=150)
    plt.close(fig)
    other = (a.obs["subcluster"] == "").values  # types outside the compartments (e.g. Neural_Schwann)
    a.obs.loc[other, "subcluster"] = a.obs.loc[other, "cell_type"].astype(str).values
    a.obs["subcluster"] = a.obs["subcluster"].astype(str)
    a.write_h5ad(work_dir / "merged" / "refined.h5ad", compression="gzip")


if __name__ == "__main__":
    main(*sys.argv[1:4])
