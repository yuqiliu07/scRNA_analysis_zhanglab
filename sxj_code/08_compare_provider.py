"""Step H: compare this pipeline's final annotated nuclei with the sequencing provider's delivered results.

Provider files (local only, not in the repo): <provider_dir>/1.Data_process/<S>/<S>_BarcodeMatrix_{meta,filter,
filter_single_double}.xls (per-sample Seurat QC and DoubletFinder), 3.Identification_subpopulation/{Clustering,UMAP}/
Combination1.* (merged clustering, 26 clusters), 4.Detection_Markergene/Combination1/Combination1.cluster.diffexp.
significant.xls (cluster markers).

Usage: python 08_compare_provider.py <provider_dir> <work_dir> <out_dir>
Output (<out_dir>): tables (cmp_*.csv), figures (cmp_*.png), cmp_summary.json
"""
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

SAMPLES = list("ABCDEFGH")
SEED = 0
TYPES = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
         "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"]
TYPE_ZH = {"Cardiomyocyte": "心肌", "Fibroblast": "成纤维", "Endothelial": "内皮", "Endocardial": "心内膜",
           "Lymphatic_EC": "淋巴管内皮", "Pericyte": "周细胞", "SMC": "平滑肌", "Macrophage": "巨噬",
           "Epicardial": "心外膜", "Neural_Schwann": "神经 / 施万", "Lymphocyte": "淋巴细胞", "Dendritic": "树突状细胞"}
COL = dict(zip(TYPES, plt.get_cmap("tab20").colors[::2] + plt.get_cmap("tab20").colors[1::2]))
MARKERS = ["Myh6", "Tnnt2", "Ryr2", "Ttn", "Dcn", "Abca8a", "Pdgfra", "Pecam1", "Flt1", "Ptprb", "Npr3", "Cdh11",
           "Prox1", "Reln", "Rgs5", "Kcnj8", "Pdgfrb", "Myh11", "Acta2", "F13a1", "Mrc1", "Cd163", "Msln", "Upk3b"]
CM_GENES = ["Myh6", "Tnnt2", "Ryr2", "Ttn"]
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def provider_cells(pdir):
    rows = []
    for s in SAMPLES:
        base = pdir / "1.Data_process" / s
        meta = pd.read_csv(base / f"{s}_BarcodeMatrix_meta.xls", sep="\t").set_index("Barcode")
        filt = pd.read_csv(base / f"{s}_BarcodeMatrix_filter.xls", sep="\t").set_index("Barcode")
        dd = pd.read_csv(base / f"{s}_BarcodeMatrix_filter_single_double.xls", sep="\t").set_index("Barcode")
        dc = [c for c in dd.columns if c.startswith("DF.class")][-1]
        d = pd.DataFrame({"sample": s, "pct_mt_provider": meta["percent.MT"]})
        d["provider_filter_pass"] = d.index.isin(filt.index)
        d["provider_doublet"] = dd[dc].reindex(d.index).eq("Doublet")
        d["provider_final"] = d.index.isin(dd.index[dd[dc] == "Singlet"])
        d.index = s + "_" + d.index
        rows.append(d)
    p = pd.concat(rows)
    cl = pd.read_csv(pdir / "3.Identification_subpopulation/Clustering/Combination1.clusters_info.xls", sep="\t")
    um = pd.read_csv(pdir / "3.Identification_subpopulation/UMAP/Combination1.cluster.UMAP.xls", sep="\t")
    key = (cl["orig.ident"] + "_" + cl["barcode"].str.replace(r"_\d+$", "", regex=True)).values
    um = um.set_index("barcode").reindex(cl["barcode"].values)
    cl.index = key
    um.index = key
    p["provider_cluster"] = cl["seurat_clusters"].reindex(p.index)
    p[["provider_umap1", "provider_umap2"]] = um[["UMAP_1", "UMAP_2"]].reindex(p.index).values
    return p


def our_cells(work_dir):
    per = []
    for s in SAMPLES:
        o = ad.read_h5ad(work_dir / "per_sample" / f"{s}.h5ad", backed="r").obs
        o = o[["pass_min_genes", "pass_umi", "pass_genes", "pass_mt", "predicted_doublet", "qc_pass", "total",
               "n_genes", "pct_mt"]].copy()
        o.index = s + "_" + o.index
        per.append(o)
    per = pd.concat(per)
    cb = ad.read_h5ad(work_dir / "merged" / "merged_cellbender.h5ad", backed="r")
    cbo = cb.obs[["cluster"]].astype(str)
    cb.file.close()
    ann = pd.read_csv(OUT_TABLES_SRC / "annotation_clusters.csv", dtype={"cluster_cb": str}).set_index("cluster_cb")
    a = read_final(work_dir, backed="r")
    fin = a.obs[["sample", "cell_type", "subtype"]].copy()
    fin[["umap1", "umap2"]] = a.obsm["X_umap"][:]
    a.file.close()
    return per, cbo, ann, fin


def reason_ours(idx, per, cbo, ann, fin):
    r = pd.Series("kept", index=idx)
    inper = idx.isin(per.index)
    r[~inper] = "EmptyDrops 未判为细胞"
    p = per.reindex(idx)
    m = inper & (r == "kept")
    r[m & ~p["pass_min_genes"].fillna(True).astype(bool)] = "基因数 < 200"
    m = r == "kept"
    r[m & inper & ~p["pass_mt"].fillna(True).astype(bool)] = "线粒体比例超阈值"
    m = r == "kept"
    r[m & inper & ~(p["pass_umi"].fillna(True).astype(bool) & p["pass_genes"].fillna(True).astype(bool))] = "UMI/基因 MAD"
    m = r == "kept"
    r[m & inper & p["predicted_doublet"].fillna(False).astype(bool)] = "Scrublet 双细胞"
    m = (r == "kept") & inper
    notcb = m & ~idx.isin(cbo.index)
    r[notcb] = "CellBender 非细胞/窗口外"
    m = (r == "kept") & idx.isin(cbo.index)
    cbt = cbo["cluster"].reindex(idx).map(ann["cell_type"])
    r[m & (cbt == "Doublet")] = "双细胞簇（07 第一轮）"
    r[m & (cbt == "Low_quality")] = "低质量簇（07）"
    m = (r == "kept") & idx.isin(cbo.index) & ~idx.isin(fin.index)
    r[m] = "双细胞（07 第二轮 / 10 亚群）"
    return r


def detection(raw_dir, cells, genes, group):
    """fraction of nuclei with >=1 raw count, per group; cells index 'S_barcode'."""
    out = []
    for s in SAMPLES:
        sub = cells[cells.index.str.startswith(s + "_")]
        if not len(sub):
            continue
        a = sc.read_10x_h5(raw_dir / s / "raw_feature_bc_matrix.h5")
        a.var_names_make_unique()
        bc = sub.index.str.slice(2)
        x = a[bc.values, genes].X
        df = pd.DataFrame((x > 0).toarray(), columns=genes, index=sub.index)
        df["group"] = sub[group].values
        out.append(df)
    df = pd.concat(out)
    return df.groupby("group").mean()


def final_h5ad(work_dir):
    """Final nucleus set: 10_name_subclusters.py, else 09_refine.py (rescue), else 07 annotation."""
    for f in ("final.h5ad", "refined.h5ad"):
        if (work_dir / "merged" / f).exists():
            return work_dir / "merged" / f
    return work_dir / "merged" / "annotated.h5ad"


def read_final(work_dir, **kw):
    a = ad.read_h5ad(final_h5ad(work_dir), **kw)
    if "cell_type_final" in a.obs:
        a.obs["cell_type"] = a.obs["cell_type_final"].astype(str)
    return a


def main(pdir, work_dir, out):
    global OUT_TABLES_SRC
    pdir, work_dir, out = Path(pdir), Path(work_dir), Path(out)
    OUT_TABLES_SRC = Path(__file__).resolve().parents[2] / "results" / "mouse_qc" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    S = {}

    prov = provider_cells(pdir)
    per, cbo, ann, fin = our_cells(work_dir)
    pf = prov.index[prov["provider_final"]]
    of = fin.index
    both = pf.intersection(of)
    S["n_provider_final"], S["n_ours_final"], S["n_both"] = len(pf), len(of), len(both)
    S["n_provider_only"], S["n_ours_only"] = len(pf.difference(of)), len(of.difference(pf))

    # 1) per-sample cell sets
    ps = pd.DataFrame({"provider": prov[prov.provider_final].groupby("sample").size(),
                       "ours": fin.groupby("sample").size(),
                       "both": pd.Series([s[0] for s in both]).value_counts()}).reindex(SAMPLES)
    ps["provider_only"] = ps["provider"] - ps["both"]
    ps["ours_only"] = ps["ours"] - ps["both"]
    ps["cellranger"] = prov.groupby("sample").size()
    ps.to_csv(out / "cmp_cells_per_sample.csv")

    # 2) why provider-only nuclei were removed here; why our-only nuclei were removed by the provider
    po = pf.difference(of)
    r1 = reason_ours(po, per, cbo, ann, fin).value_counts()
    r1.to_csv(out / "cmp_provider_only_reasons.csv")
    oo = of.difference(pf)
    r2 = pd.Series("其他", index=oo)
    r2[~oo.isin(prov.index)] = "Cell Ranger 未判为细胞"
    inprov = oo.isin(prov.index)
    pv = prov.reindex(oo)
    r2[inprov & ~pv["provider_filter_pass"].fillna(False).astype(bool)] = "公司 Seurat 过滤"
    r2[inprov & pv["provider_doublet"].fillna(False).astype(bool)] = "DoubletFinder 双细胞"
    r2 = r2.value_counts()
    r2.to_csv(out / "cmp_ours_only_reasons.csv")
    pm = per.reindex(po)
    S["provider_only_median_pct_mt"] = float(pm["pct_mt"].median())
    S["both_median_pct_mt"] = float(per.reindex(both)["pct_mt"].median())
    S["provider_only_median_umi"] = float(pm["total"].median())
    S["both_median_umi"] = float(per.reindex(both)["total"].median())

    # 3) DoubletFinder calls by our final cell type
    f2 = fin.join(prov[["provider_doublet"]])
    dfr = f2.groupby("cell_type")["provider_doublet"].mean().reindex(TYPES)
    dfr.to_csv(out / "cmp_doubletfinder_by_type.csv")

    # 4) cluster correspondence on shared nuclei
    j = fin.loc[both].join(prov[["provider_cluster"]])
    S["ari_shared"] = float(adjusted_rand_score(j["provider_cluster"], j["subtype"]))
    xt = pd.crosstab(j["provider_cluster"], j["cell_type"]).reindex(columns=TYPES, fill_value=0)
    fr = xt.div(xt.sum(1), axis=0)
    corr = pd.DataFrame({"n_shared": xt.sum(1), "our_type": fr.idxmax(1), "frac": fr.max(1).round(3),
                         "second": fr.apply(lambda r: r.sort_values().index[-2], axis=1),
                         "second_frac": fr.apply(lambda r: r.sort_values().iloc[-2], axis=1).round(3)})
    corr["n_provider"] = prov[prov.provider_final].groupby("provider_cluster").size()
    corr["mixed"] = corr["frac"] < 0.7
    corr.to_csv(out / "cmp_cluster_correspondence.csv")
    xt.to_csv(out / "cmp_cluster_crosstab.csv")

    # 5) composition: provider clusters mapped to our types by majority on shared nuclei
    pfin = prov[prov.provider_final].copy()
    pfin["type"] = pfin["provider_cluster"].map(corr["our_type"])
    pfin.loc[pfin["provider_cluster"].map(corr["mixed"]).fillna(False).astype(bool), "type"] = "Mixed"
    comp_p = pd.crosstab(pfin["sample"], pfin["type"], normalize="index")
    comp_o = pd.crosstab(fin["sample"], fin["cell_type"], normalize="index")
    comp = pd.concat({"provider": comp_p, "ours": comp_o}, axis=1).fillna(0)
    comp.to_csv(out / "cmp_composition.csv")
    tot = pd.DataFrame({"provider_n": pfin["type"].value_counts(), "ours_n": fin["cell_type"].value_counts()})
    tot["provider_frac"] = tot["provider_n"] / tot["provider_n"].sum()
    tot["ours_frac"] = tot["ours_n"] / tot["ours_n"].sum()
    tot.to_csv(out / "cmp_composition_total.csv")

    # 6) marker detection: provider clusters on raw counts vs our types on CellBender counts
    raw_dir = pdir / "1.Data_process"
    det_p = detection(raw_dir, pfin.assign(grp=pfin["provider_cluster"].astype(int)), MARKERS, "grp")
    det_p.to_csv(out / "cmp_marker_detection_provider.csv")
    a = read_final(work_dir)
    x = a[:, MARKERS].layers["counts"]
    det_o = pd.DataFrame((x > 0).toarray(), columns=MARKERS, index=a.obs_names).groupby(
        a.obs["cell_type"].values).mean().reindex(TYPES)
    det_o.to_csv(out / "cmp_marker_detection_ours.csv")
    # same nuclei, our types, raw counts (to separate counts from cell set)
    det_or = detection(raw_dir, fin, MARKERS, "cell_type").reindex(TYPES)
    det_or.to_csv(out / "cmp_marker_detection_ours_rawcounts.csv")

    # 7) marker sharpness: provider top-10 markers per cluster vs ours per type (one-vs-rest Wilcoxon)
    mk = pd.read_csv(pdir / "4.Detection_Markergene/Combination1/Combination1.cluster.diffexp.significant.xls",
                     sep="\t")
    mk = mk[mk["avg_log2FC"] > 0].sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
    top_p = mk.groupby("cluster").head(10)
    sp = top_p.groupby("cluster").agg(top_markers=("gene", lambda g: ", ".join(g)),
                                      median_log2FC=("avg_log2FC", "median"), median_pct1=("pct.1", "median"),
                                      median_pct2=("pct.2", "median"))
    sp = sp.join(corr[["our_type", "frac", "mixed", "n_provider"]])
    sp.to_csv(out / "cmp_markers_provider.csv")
    rng = np.random.default_rng(SEED)
    keep = np.concatenate([rng.choice(np.where(a.obs["cell_type"] == t)[0],
                                      min(2000, int((a.obs["cell_type"] == t).sum())), replace=False)
                           for t in TYPES])
    sub = a[np.sort(keep)].copy()
    del a
    sc.tl.rank_genes_groups(sub, "cell_type", method="wilcoxon", n_genes=50, use_raw=False)
    rows = []
    cnt = sub.layers["counts"]
    for t in TYPES:
        r = sc.get.rank_genes_groups_df(sub, t)
        r = r[r["logfoldchanges"] > 0].head(10)
        ingrp = (sub.obs["cell_type"] == t).values
        gi = [sub.var_names.get_loc(g) for g in r["names"]]
        det = (cnt[:, gi] > 0).toarray()
        rows.append({"cell_type": t, "top_markers": ", ".join(r["names"]),
                     "median_log2FC": float(r["logfoldchanges"].median()),
                     "median_pct1": float(np.median(det[ingrp].mean(0))),
                     "median_pct2": float(np.median(det[~ingrp].mean(0)))})
    so = pd.DataFrame(rows).set_index("cell_type")
    so.to_csv(out / "cmp_markers_ours.csv")

    # 8) CM transcripts in non-CM populations: provider (raw counts, provider clusters) vs ours (CellBender)
    amb = pd.read_csv(OUT_TABLES_SRC / "cb_ambient_check.csv", index_col=0)
    amb.to_csv(out / "cmp_ambient_cm_genes.csv")

    # ---------------- figures ----------------
    fig, ax = plt.subplots(figsize=(9, 4.2))
    xs = np.arange(len(SAMPLES))
    ax.bar(xs, ps["both"], color="#7f7f7f", label="两边共有")
    ax.bar(xs, ps["provider_only"], bottom=ps["both"], color="#1f77b4", label="仅公司保留")
    ax.bar(xs, -ps["ours_only"], color="#ff7f0e", label="仅本流程保留（向下）")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs, SAMPLES)
    ax.set_ylabel("细胞核数")
    ax.set_title("最终细胞集合：公司 vs 本流程")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    fig.savefig(out / "cmp_fig1_cellsets.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    r1.sort_values().plot.barh(ax=axes[0], color="#1f77b4")
    axes[0].set_title(f"仅公司保留的 {len(po):,} 个核：本流程去除原因", fontsize=10)
    r2.sort_values().plot.barh(ax=axes[1], color="#ff7f0e")
    axes[1].set_title(f"仅本流程保留的 {len(oo):,} 个核：公司去除原因", fontsize=10)
    for axx in axes:
        axx.set_xlabel("细胞核数")
    fig.tight_layout()
    fig.savefig(out / "cmp_fig2_reasons.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 3.8))
    v = dfr.dropna()
    ax.bar([TYPE_ZH[t] for t in v.index], v.values * 100, color=[COL[t] for t in v.index])
    for i, val in enumerate(v.values):
        ax.text(i, val * 100 + 1, f"{val*100:.0f}%", ha="center", fontsize=8)
    ax.set_ylabel("被 DoubletFinder 判为双细胞 (%)")
    ax.set_title("本流程最终细胞中，被公司判为双细胞的比例（按细胞类型）")
    fig.tight_layout()
    fig.savefig(out / "cmp_fig3_doubletfinder_by_type.png", dpi=200)
    plt.close(fig)

    shared_p = prov.loc[both]
    rng = np.random.default_rng(SEED)
    idx_p = rng.permutation(len(pfin))
    idx_o = rng.permutation(len(fin))
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    pcl = pfin["provider_cluster"].astype(int)
    cm20 = plt.get_cmap("tab20")
    cols = np.array([cm20(c % 20) for c in pcl])
    axes[0, 0].scatter(pfin["provider_umap1"].values[idx_p], pfin["provider_umap2"].values[idx_p], s=0.2,
                       c=cols[idx_p], rasterized=True)
    for c in sorted(pcl.unique()):
        mm = pcl == c
        axes[0, 0].text(pfin.loc[mm, "provider_umap1"].median(), pfin.loc[mm, "provider_umap2"].median(), str(c),
                        fontsize=8, ha="center", weight="bold")
    axes[0, 0].set_title("公司 UMAP：公司 26 个簇")
    sj = shared_p.join(fin[["cell_type"]])
    ii = rng.permutation(len(sj))
    axes[0, 1].scatter(sj["provider_umap1"].values[ii], sj["provider_umap2"].values[ii], s=0.2,
                       c=[COL[t] for t in sj["cell_type"].values[ii]], rasterized=True)
    axes[0, 1].set_title("公司 UMAP：按本流程细胞类型着色（共有核）")
    axes[1, 0].scatter(fin["umap1"].values[idx_o], fin["umap2"].values[idx_o], s=0.2,
                       c=[COL[t] for t in fin["cell_type"].values[idx_o]], rasterized=True)
    axes[1, 0].set_title("本流程 UMAP：细胞类型")
    fj = fin.loc[both].join(prov[["provider_cluster"]])
    ii = rng.permutation(len(fj))
    axes[1, 1].scatter(fj["umap1"].values[ii], fj["umap2"].values[ii], s=0.2,
                       c=[cm20(int(c) % 20) for c in fj["provider_cluster"].values[ii]], rasterized=True)
    axes[1, 1].set_title("本流程 UMAP：按公司簇着色（共有核）")
    for axx in axes.ravel():
        axx.set_xticks([]), axx.set_yticks([])
    handles = [plt.Line2D([], [], marker="o", ls="", color=COL[t], label=TYPE_ZH[t]) for t in TYPES]
    fig.legend(handles=handles, loc="lower center", ncol=9, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out / "cmp_fig4_umaps.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(fr.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(TYPES)), [TYPE_ZH[t] for t in TYPES], rotation=45, ha="right")
    ax.set_yticks(range(len(fr)), [f"簇 {c} (n={corr.loc[c,'n_provider']:,})" for c in fr.index], fontsize=8)
    for i in range(fr.shape[0]):
        for k2 in range(fr.shape[1]):
            if fr.values[i, k2] >= 0.1:
                ax.text(k2, i, f"{fr.values[i,k2]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if fr.values[i, k2] > 0.6 else "black")
    plt.colorbar(im, ax=ax, fraction=0.04, label="公司簇中落在各本流程类型的比例（共有核）")
    ax.set_title("公司簇 → 本流程细胞类型")
    fig.tight_layout()
    fig.savefig(out / "cmp_fig5_correspondence.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for axx, key, title in [(axes[0], "provider", "公司（簇按多数映射；混合簇单列）"), (axes[1], "ours", "本流程")]:
        c = comp[key]
        bottom = np.zeros(len(c))
        for t in TYPES + ["Mixed"]:
            if t in c:
                axx.bar(c.index, c[t].values, bottom=bottom, color=COL.get(t, "#cccccc"),
                        label=TYPE_ZH.get(t, "混合簇"))
                bottom += c[t].values
        axx.set_title(title)
    axes[0].set_ylabel("细胞比例")
    axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(out / "cmp_fig6_composition.png", dpi=200)
    plt.close(fig)

    order_p = corr.sort_values(["our_type", "n_provider"], ascending=[True, False]).index
    fig, axes = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={"width_ratios": [1, 0.45]})
    dp = det_p.loc[order_p]
    axes[0].imshow(dp.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    axes[0].set_xticks(range(len(MARKERS)), MARKERS, rotation=90, fontsize=8)
    axes[0].set_yticks(range(len(dp)), [f"簇 {c} → {TYPE_ZH.get(corr.loc[c,'our_type'],'')}"
                                        + ("（混合）" if corr.loc[c, "mixed"] else "") for c in dp.index], fontsize=8)
    axes[0].set_title("公司：各簇 marker 检出比例（原始计数）")
    im = axes[1].imshow(det_o.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    axes[1].set_xticks(range(len(MARKERS)), MARKERS, rotation=90, fontsize=8)
    axes[1].set_yticks(range(len(TYPES)), [TYPE_ZH[t] for t in TYPES], fontsize=9)
    axes[1].set_title("本流程：各类型 marker 检出比例（CellBender）")
    plt.colorbar(im, ax=axes[1], fraction=0.05, label="检出比例")
    fig.tight_layout()
    fig.savefig(out / "cmp_fig7_marker_detection.png", dpi=180)
    plt.close(fig)

    # CM genes detected in non-CM types: provider (raw counts, mapped types) vs ours (CellBender)
    nonp = det_p.join(corr[["our_type", "n_provider", "mixed"]])
    nonp = nonp[(nonp["our_type"] != "Cardiomyocyte") & ~nonp["mixed"].fillna(True).astype(bool)]
    wp = nonp.groupby("our_type").apply(lambda d: (d[CM_GENES].mul(d["n_provider"], axis=0).sum()
                                                   / d["n_provider"].sum()))
    cmp_cm = pd.concat({"provider_raw": wp, "ours_raw_same_types": det_or.loc[wp.index, CM_GENES],
                        "ours_cellbender": det_o.loc[wp.index, CM_GENES]}, axis=1)
    cmp_cm.to_csv(out / "cmp_cm_gene_detection_nonCM.csv")
    fig, ax = plt.subplots(figsize=(9, 4))
    tt = list(wp.index)
    w = 0.27
    xs = np.arange(len(tt))
    ax.bar(xs - w, wp["Ttn"].values * 100, w, label="公司（原始计数）", color="#1f77b4")
    ax.bar(xs, det_or.loc[tt, "Ttn"].values * 100, w, label="本流程细胞 + 原始计数", color="#aec7e8")
    ax.bar(xs + w, det_o.loc[tt, "Ttn"].values * 100, w, label="本流程（CellBender）", color="#ff7f0e")
    ax.set_xticks(xs, [TYPE_ZH[t] for t in tt])
    ax.set_ylabel("检出 Ttn 的细胞核 (%)")
    ax.set_title("非心肌细胞中心肌基因 Ttn 的检出比例")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "cmp_fig8_ttn_nonCM.png", dpi=200)
    plt.close(fig)

    S["doubletfinder_rate_by_type"] = dfr.astype(float).round(4).to_dict()
    S["provider_only_reasons"] = r1.to_dict()
    S["ours_only_reasons"] = r2.to_dict()
    (out / "cmp_summary.json").write_text(json.dumps(S, ensure_ascii=False, indent=1, default=float),
                                          encoding="utf8")
    print(json.dumps(S, ensure_ascii=False, indent=1, default=float))
    print(ps.to_string(), "\n", corr.to_string(), "\n", tot.round(3).to_string(), "\n",
          sp.round(3).to_string(), "\n", so.round(3).to_string(), "\n", cmp_cm.round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
