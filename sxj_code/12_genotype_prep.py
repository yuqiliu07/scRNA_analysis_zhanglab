"""WT vs TNNC1-I4M (+/-): figures that need no DE test, and pseudobulk matrices for 13_pseudobulk_de.R.

Nuclei: final.h5ad (10_name_subclusters.py) with flag == "ok" (low-quality / ambient-high subclusters excluded).
Groups: sample_sheet.csv (4 WT: A-D; 4 I4M: E-H; 2 male + 2 female per group).

Outputs (results_dir/genotype/):
  fig_umap_by_genotype.png       UMAP split by genotype, coloured by cell type
  fig_composition.png            per-sample proportion of each cell type, WT vs I4M
  composition_test.csv           proportion test per type: asin(sqrt(p)) ~ sex + genotype (OLS, BH), as in propeller
  fig_marker_dotplot.png         top marker genes of each cell type
  fig_subcluster_dotplot_<X>.png top markers of the subclusters of CM / EC / Fib / Mural / Immune
  cell_counts.csv                nuclei per sample and type
Pseudobulk (work_dir/genotype/pb_<type>.csv.gz): summed CellBender counts, genes x samples, for every type with
>= MIN_NUC nuclei in every sample.

Usage: python 12_genotype_prep.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

plt.rcParams["font.family"] = "Arial"
ORDER = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
         "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"]
LABEL = {"Cardiomyocyte": "Cardiomyocyte", "Fibroblast": "Fibroblast", "Endothelial": "Endothelial",
         "Endocardial": "Endocardial", "Lymphatic_EC": "Lymphatic EC", "Pericyte": "Pericyte", "SMC": "SMC",
         "Macrophage": "Macrophage", "Lymphocyte": "Lymphocyte", "Dendritic": "Dendritic (cDC1)",
         "Epicardial": "Epicardial", "Neural_Schwann": "Neural / Schwann"}
COLORS = dict(zip(ORDER, [plt.get_cmap("tab20")(i) for i in list(range(0, 20, 2)) + [1, 3]]))
GCOL = {"WT": "#7F7F7F", "I4M": "#C0392B"}
MIN_NUC = 20
AMBIENT = {"Cmss1", "Lars2", "Camk1d", "Gphn", "Cdk8", "AY036118", "Malat1"}  # ambient / low-quality signature
SEED = 0
# canonical markers shown for every type (plus the top data-driven markers)
CANON = {"Cardiomyocyte": ["Ttn", "Myh6", "Tnnt2"], "Fibroblast": ["Dcn", "Pdgfra", "Abca8a"],
         "Endothelial": ["Pecam1", "Flt1", "Ptprb"], "Endocardial": ["Npr3", "Tmem108", "Cdh11"],
         "Lymphatic_EC": ["Prox1", "Mmrn1", "Reln"], "Pericyte": ["Pdgfrb", "Abcc9", "Kcnj8"],
         "SMC": ["Myh11", "Acta2", "Kcnab1"], "Macrophage": ["F13a1", "Mrc1", "Cd163"],
         "Lymphocyte": ["Ptprc", "Skap1", "Pax5"], "Dendritic": ["Clec9a", "Flt3", "Wdfy4"],
         "Epicardial": ["Msln", "Upk3b", "Muc16"], "Neural_Schwann": ["Plp1", "Cdh19", "Scn7a"]}


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    out = res_dir / "genotype"
    out.mkdir(parents=True, exist_ok=True)
    pbdir = work_dir / "genotype"
    pbdir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")

    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad")
    a = a[(a.obs["flag"] == "ok").values].copy()
    a.obs["cell_type"] = a.obs["cell_type_final"].astype(str)
    a.obs["genotype"] = a.obs["sample"].map(meta["genotype"]).astype(str)
    a.obs["sex"] = a.obs["sample"].map(meta["sex"]).astype(str)
    print(a.n_obs, "nuclei", flush=True)

    # ---- counts and composition ----
    cnt = pd.crosstab(a.obs["sample"], a.obs["cell_type"]).reindex(columns=ORDER, fill_value=0)
    cnt.to_csv(out / "cell_counts.csv")
    prop = cnt.div(cnt.sum(1), axis=0)
    rows = []
    for t in ORDER:
        d = meta.assign(p=prop[t]).assign(y=lambda x: np.arcsin(np.sqrt(x["p"])))
        f = smf.ols("y ~ C(sex) + C(genotype, Treatment('WT'))", data=d).fit()
        k = [c for c in f.params.index if "genotype" in c][0]
        rows.append({"cell_type": t, "mean_WT": d.loc[d.genotype == "WT", "p"].mean(),
                     "mean_I4M": d.loc[d.genotype == "I4M", "p"].mean(), "coef_asin": f.params[k],
                     "p_value": f.pvalues[k]})
    ct = pd.DataFrame(rows)
    ct["ratio_I4M_WT"] = ct["mean_I4M"] / ct["mean_WT"]
    ct["fdr"] = multipletests(ct["p_value"], method="fdr_bh")[1]
    ct.to_csv(out / "composition_test.csv", index=False)
    print(ct.round(4).to_string(), flush=True)

    fig, ax = plt.subplots(figsize=(12, 4.6))
    w = 0.36
    for i, t in enumerate(ORDER):
        for j, g in enumerate(["WT", "I4M"]):
            v = prop.loc[meta.index[meta.genotype == g], t].values
            x = i + (j - 0.5) * w * 1.1
            ax.boxplot(v, positions=[x], widths=w, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=GCOL[g], alpha=0.35, color=GCOL[g]),
                       medianprops=dict(color=GCOL[g]), whiskerprops=dict(color=GCOL[g]),
                       capprops=dict(color=GCOL[g]))
            sx = meta.loc[meta.genotype == g, "sex"].values
            for vv, s in zip(v, sx):
                ax.scatter(x + np.random.default_rng(i * 10 + j).uniform(-0.07, 0.07), vv, s=16,
                           color=GCOL[g], marker="o" if s == "M" else "^", zorder=3, edgecolor="k", lw=0.3)
        r = ct.set_index("cell_type").loc[t]
        if r["fdr"] < 0.1:
            top = prop[t].max()
            ax.text(i, top * 1.04 + 0.005, "*" if r["fdr"] < 0.05 else "†", ha="center", fontsize=13,
                    color="#C0392B")
    for i in range(len(ORDER) - 1):
        ax.axvline(i + 0.5, color="#e0e0e0", lw=0.8)
    ax.set_xticks(range(len(ORDER)), [LABEL[t] for t in ORDER], rotation=40, ha="right")
    ax.set_ylabel("Proportion of nuclei")
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    h = [plt.Rectangle((0, 0), 1, 1, color=GCOL[g], alpha=0.5) for g in ["WT", "I4M"]]
    h += [plt.Line2D([], [], marker="o", ls="", color="k", label="male"),
          plt.Line2D([], [], marker="^", ls="", color="k", label="female")]
    ax.legend(h, ["WT (n=4)", "TNNC1-I4M (n=4)", "male", "female"], frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "fig_composition.png", dpi=200)
    # zoom on the smaller types
    ax.set_ylim(0, 0.1)
    ax.set_title("Zoom: proportion 0-0.10")
    fig.savefig(out / "fig_composition_zoom.png", dpi=200)
    plt.close(fig)

    # ---- UMAP by genotype ----
    um = a.obsm["X_umap"]
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.4))
    for ax, g in zip(axes, ["WT", "I4M"]):
        m = np.where((a.obs["genotype"] == g).values)[0]
        ax.scatter(um[:, 0], um[:, 1], s=0.2, c="#eeeeee", rasterized=True)
        m = rng.permutation(m)
        ax.scatter(um[m, 0], um[m, 1], s=0.3, c=[COLORS[t] for t in a.obs["cell_type"].values[m]],
                   rasterized=True)
        ax.set_title(f"{'WT' if g == 'WT' else 'TNNC1-I4M (+/-)'}   n = {len(m):,} nuclei", fontsize=13)
        ax.set_xticks([]), ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for i, t in enumerate(ORDER):
        mm = (a.obs["cell_type"] == t).values
        xy = np.median(um[mm], 0)
        for ax in axes:
            ax.text(*xy, str(i + 1), fontsize=10, weight="bold", ha="center", va="center",
                    bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec=COLORS[t], alpha=0.85))
    h = [plt.Line2D([], [], marker="o", ls="", color=COLORS[t], label=f"{i + 1} {LABEL[t]}")
         for i, t in enumerate(ORDER)]
    fig.legend(handles=h, loc="center right", frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.84, 1))
    fig.savefig(out / "fig_umap_by_genotype.png", dpi=180)
    plt.close(fig)

    # ---- marker dot plot (cell types) ----
    mk = pd.read_csv(res_dir / "provider_compare" / "cmp_markers_ours.csv").set_index("cell_type")
    genes = {}
    for t in ORDER:
        top = [g for g in mk.loc[t, "top_markers"].split(", ") if g in a.var_names][:3]
        gl = []
        for g in CANON[t] + top:
            if g in a.var_names and g not in gl and all(g not in v for v in genes.values()):
                gl.append(g)
        genes[LABEL[t]] = gl[:5]
    a.obs["ct_label"] = pd.Categorical(a.obs["cell_type"].map(LABEL), categories=[LABEL[t] for t in ORDER])
    dp = sc.pl.dotplot(a, genes, groupby="ct_label", layer=None, standard_scale="var", show=False,
                       return_fig=True, figsize=(17, 4.8), cmap="Reds")
    dp.savefig(out / "fig_marker_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close("all")

    # ---- subcluster marker dot plots (5 compartments) ----
    smk = pd.read_csv(res_dir / "tables" / "refine_submarkers.csv")
    for comp in ["CM", "EC", "Fib", "Mural", "Immune"]:
        b = a[a.obs["subcluster"].astype(str).str.startswith(comp + "_").values].copy()
        b.obs["st"] = b.obs["state"].astype(str)
        states = b.obs.groupby("st", observed=True).size()
        states = states[states >= 50].index.tolist()
        b = b[b.obs["st"].isin(states).values].copy()
        gd, used = {}, set()
        for st in states:
            subs = b.obs.loc[b.obs["st"] == st, "subcluster"].astype(str).value_counts().index[0]
            gl = [g for g in smk.loc[smk["group"] == subs, "names"] if g in b.var_names and g not in used
                  and not g.startswith(("mt-", "Gm", "Rpl", "Rps")) and "Rik" not in g and g not in AMBIENT][:5]
            used.update(gl)
            gd[st] = gl
        dp = sc.pl.dotplot(b, gd, groupby="st", standard_scale="var", show=False, return_fig=True,
                           figsize=(max(8, 0.42 * sum(len(v) for v in gd.values()) + 3), 0.45 * len(gd) + 2.2),
                           cmap="Reds")
        dp.savefig(out / f"fig_subcluster_dotplot_{comp}.png", dpi=200, bbox_inches="tight")
        plt.close("all")

    # ---- pseudobulk ----
    X = a.layers["counts"]
    keep = [t for t in ORDER if (cnt[t] >= MIN_NUC).all()]
    print("pseudobulk types:", keep, flush=True)
    for t in keep:
        cols = {}
        for s in meta.index:
            m = np.where(((a.obs["cell_type"] == t) & (a.obs["sample"] == s)).values)[0]
            cols[s] = np.asarray(X[m].sum(0)).ravel()
        pd.DataFrame(cols, index=a.var_names).astype(int).to_csv(pbdir / f"pb_{t}.csv.gz")
    meta.assign(n_nuclei_total=cnt.sum(1)).to_csv(pbdir / "meta.csv")
    pd.Series(a.var["gene_ids"].values, index=a.var_names, name="ensembl").to_csv(pbdir / "genes.csv")


if __name__ == "__main__":
    main(*sys.argv[1:3])
