"""Fibroblasts, WT vs TNNC1-I4M: state scores (project definitions) and an extracellular-matrix heatmap.

1) Fibroblast states: the reference states of the project's fibroblast-state method (src/v3/fib_states_v2.py, REF),
   as mouse orthologs. Every fibroblast nucleus (final.h5ad, flag == "ok") is scored per state with
   scanpy.tl.score_genes (mean expression of the state genes minus a matched random control set). Per mouse:
   mean score per state, and the share of fibroblasts with a high "activated" score (above the 95th percentile of
   all fibroblasts). WT vs I4M: OLS ~ sex + genotype, BH over states.
2) Matrisome heatmap: fibroblast pseudobulk, genes grouped by the core-matrisome / matrisome-associated
   categories (collagens, ECM glycoproteins, proteoglycans, ECM regulators, secreted pro-fibrotic factors),
   drawn with 15_genotype_figures_styled.heatmap.

Outputs (results_dir/genotype/): fib_state_scores_per_mouse.csv, fib_state_test.csv,
  styled/fib_state_scores.png, styled/heatmap_fibroblast_matrisome.png

Usage: python 22_fibroblast_states.py <work_dir> <results_dir>
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
# src/v3/fib_states_v2.py REF, mouse orthologs (WISP2 = Ccn5; MXRA5 has no mouse ortholog)
STATES = {
    "Universal Pi16+": ["Pi16", "Mfap5", "Cd34", "Dpp4", "Sema3c", "Sfrp2", "Ccn5", "Pcolce2"],
    "Universal Col15a1+": ["Col15a1", "Col4a1", "Col4a2", "Lamb1", "Npnt", "Itga8", "Gpx3", "Scn7a"],
    "Activated Cthrc1+/Postn+": ["Cthrc1", "Postn", "Lrrc15", "Col11a1", "Comp", "Fap", "Tnc", "Col1a1", "Col3a1",
                                 "Adam12", "Aspn"],
    "Inflammatory": ["Ccl19", "Cxcl9", "Cxcl10", "Ccl2", "Il6", "Mmp3", "Cxcl1", "Sod2"],
    "Contractile / myofibroblast": ["Acta2", "Myh11", "Tagln", "Des", "Cnn1", "Myl9"],
}
MATRISOME = {
    "collagens": ["Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col5a1", "Col5a2", "Col6a1", "Col6a2", "Col6a3", "Col8a1",
                  "Col14a1", "Col15a1"],
    "ECM glycoproteins": ["Fn1", "Postn", "Tnc", "Fbn1", "Eln", "Lama2", "Lamb1", "Nid1", "Thbs1", "Thbs2", "Sparc",
                          "Comp", "Fbln1", "Fbln2", "Mfap5"],
    "proteoglycans": ["Dcn", "Lum", "Bgn", "Vcan", "Fmod", "Aspn", "Hspg2", "Prelp", "Ogn"],
    "ECM regulators": ["Lox", "Loxl1", "Loxl2", "Mmp2", "Mmp14", "Timp1", "Timp2", "Timp3", "Serpine1", "Serpinh1",
                       "Ctsk", "Adamts2", "Pcolce"],
    "pro-fibrotic secreted factors": ["Tgfb1", "Tgfb2", "Tgfb3", "Ccn2", "Pdgfa", "Ltbp2", "Fgf2"],
}
SEED = 0


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    out = res_dir / "genotype" / "styled"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    keep = ((a.obs["flag"] == "ok") & (a.obs["cell_type_final"].astype(str) == "Fibroblast")).values
    f = a[keep].to_memory()
    a.file.close()
    print(f.n_obs, "fibroblast nuclei", flush=True)

    used = {}
    for k, gl in STATES.items():
        g = [x for x in gl if x in f.var_names]
        used[k] = g
        sc.tl.score_genes(f, g, score_name=k, random_state=SEED, use_raw=False)
    thr = np.percentile(f.obs["Activated Cthrc1+/Postn+"], 95)
    f.obs["high_activated"] = f.obs["Activated Cthrc1+/Postn+"] > thr
    per = f.obs.groupby("sample", observed=True)[list(STATES)].mean()
    per["% high-activated fibroblasts"] = f.obs.groupby("sample", observed=True)["high_activated"].mean() * 100
    per = per.reindex(meta.index)
    per.to_csv(res_dir / "genotype" / "fib_state_scores_per_mouse.csv")
    rows = []
    for c in per.columns:
        d = meta.assign(y=per[c].values)
        fit = smf.ols("y ~ C(sex) + C(genotype, Treatment('WT'))", data=d).fit()
        k = [x for x in fit.params.index if "genotype" in x][0]
        rows.append({"measure": c, "WT_mean": d.loc[d.genotype == "WT", "y"].mean(),
                     "I4M_mean": d.loc[d.genotype == "I4M", "y"].mean(), "p_value": fit.pvalues[k],
                     "genes_used": ", ".join(used.get(c, []))})
    tst = pd.DataFrame(rows)
    tst["fdr"] = multipletests(tst["p_value"], method="fdr_bh")[1]
    tst.to_csv(res_dir / "genotype" / "fib_state_test.csv", index=False)
    print(tst.drop(columns="genes_used").round(4).to_string(index=False), flush=True)

    cols = list(per.columns)
    fig, axes = plt.subplots(1, len(cols), figsize=(2.75 * len(cols), 3.9))
    for ax, c in zip(axes, cols):
        for j, g in enumerate(["WT", "I4M"]):
            mm = meta.index[meta.genotype == g]
            v = per.loc[mm, c].values
            ax.boxplot(v, positions=[j], widths=0.55, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=GCOL[g], color="#333333", lw=0.8), medianprops=dict(color="#333333"),
                       whiskerprops=dict(color="#333333", lw=0.8), capprops=dict(color="#333333", lw=0.8))
            for x_, v_, s_ in zip(np.linspace(-0.12, 0.12, len(v)), v, meta.loc[mm, "sex"]):
                ax.scatter(j + x_, v_, s=22, color="black", marker="o" if s_ == "M" else "^", zorder=3)
        r = tst.set_index("measure").loc[c]
        lab = "ns" if r["fdr"] >= 0.05 else f"FDR = {r['fdr']:.2g}"
        ax.set_title(c.replace(" / ", " /\n").replace("% high-activated fibroblasts", "% fibroblasts with\nhigh activated score")
                     + f"\n({lab})", fontsize=10)
        ax.set_xticks([0, 1], ["WT", "I4M"])
    axes[0].set_ylabel("State score (per mouse)")
    axes[-1].set_ylabel("% of fibroblasts")
    h = [Patch(facecolor=GCOL[g], edgecolor="#333333", label=g) for g in ["WT", "I4M"]]
    h += [Line2D([], [], marker="o", ls="", color="black", label="male"),
          Line2D([], [], marker="^", ls="", color="black", label="female")]
    fig.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fib_state_scores.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    s15 = _load("step15", "15_genotype_figures_styled.py")
    pbf = pd.read_csv(work_dir / "genotype" / "pb_Fibroblast.csv.gz", index_col=0)[meta.index]
    de = pd.read_csv(work_dir / "genotype" / "de_Fibroblast.csv.gz")
    s15.heatmap(pbf, de, meta, MATRISOME, out, "heatmap_fibroblast_matrisome.png",
                "Fibroblasts: extracellular matrix (matrisome) genes, per mouse")
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
