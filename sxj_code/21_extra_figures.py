"""Extra figures for the WT vs I4M deck (style as 15_genotype_figures_styled.py):

  styled/confound_metrics.png   per-mouse technical metrics (from 14_genotype_figures.py: sample_technical_metrics.csv)
                                in the deck's genotype colours, with plain axis titles
  styled/ec_subtype_composition.png  share of each endothelial / endocardial subtype within the endothelial
                                compartment per mouse, WT vs I4M (asin-sqrt ~ sex + genotype, BH)
  ec_subtype_composition_test.csv

Usage: python 21_extra_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
# capillary high/low-UMI and endocardial 1/2 subclusters differ mainly by RNA quality, so they are merged
EC_STATES = {"EC_arterial": "Arterial EC", "EC_capillary_high_UMI": "Capillary EC", "EC_capillary_low_UMI": "Capillary EC",
             "Endocardial_1": "Endocardial", "Endocardial_2": "Endocardial", "Lymphatic_EC": "Lymphatic EC"}
EC_ORDER = ["Capillary EC", "Arterial EC", "Endocardial", "Lymphatic EC"]


def box_by_genotype(ax, prop, meta, cols, fdr=None):
    for i, t in enumerate(cols):
        for j, g in enumerate(["WT", "I4M"]):
            v = prop.loc[meta.index[meta.genotype == g], t].values
            x = i + (j - 0.5) * 0.38
            ax.boxplot(v, positions=[x], widths=0.34, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=GCOL[g], color="#333333", lw=0.8), medianprops=dict(color="#333333"),
                       whiskerprops=dict(color="#333333", lw=0.8), capprops=dict(color="#333333", lw=0.8))
            ax.scatter(np.full(len(v), x) + np.linspace(-0.06, 0.06, len(v)), v, s=12, color="black", zorder=3)
        if fdr is not None and fdr[t] < 0.05:
            ax.text(i, prop[t].max() * 1.05, "*", color="red", ha="center", fontsize=14)
        if i:
            ax.axvline(i - 0.5, color="#d9d9d9", lw=0.8)
    ax.set_xticks(range(len(cols)), cols, rotation=30, ha="right")
    ax.set_xlim(-0.6, len(cols) - 0.4)


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    out = res_dir / "genotype" / "styled"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")

    # ---- technical metrics ----
    m = pd.read_csv(res_dir / "genotype" / "sample_technical_metrics.csv", index_col=0)
    titles = ["Reads in introns (%)", "Reads between genes (%)", "Background (ambient) RNA\n(CPM, all nuclei)",
              "Atrial genes Myl7 + Sln\n(CPM, cardiomyocytes)"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, c, ttl in zip(axes, m.columns, titles):
        ax.bar(m.index, m[c], color=[GCOL[meta.loc[s, "genotype"]] for s in m.index], edgecolor="#333333", lw=0.4)
        ax.set_title(ttl, fontsize=11)
        ax.set_xlabel("Mouse")
    fig.legend([Patch(facecolor=GCOL[g], edgecolor="#333333") for g in ["WT", "I4M"]], ["WT (A-D)", "I4M (E-H)"],
               loc="upper right", ncol=2, frameon=False, bbox_to_anchor=(1.0, 1.08))
    fig.tight_layout()
    fig.savefig(out / "confound_metrics.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- endothelial subtype composition ----
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    obs = a.obs[["sample", "state", "flag"]].copy()
    a.file.close()
    obs = obs[(obs["flag"] == "ok").values]
    obs["sub"] = obs["state"].astype(str).map(EC_STATES)
    obs = obs.dropna(subset=["sub"])
    cnt = pd.crosstab(obs["sample"], obs["sub"]).reindex(columns=EC_ORDER, fill_value=0)
    prop = cnt.div(cnt.sum(1), axis=0)
    rows = []
    for t in EC_ORDER:
        d = meta.assign(y=np.arcsin(np.sqrt(prop[t])))
        f = smf.ols("y ~ C(sex) + C(genotype, Treatment('WT'))", data=d).fit()
        k = [c for c in f.params.index if "genotype" in c][0]
        rows.append({"subtype": t, "WT_mean": prop.loc[meta.genotype == "WT", t].mean(),
                     "I4M_mean": prop.loc[meta.genotype == "I4M", t].mean(), "p_value": f.pvalues[k],
                     "n_nuclei": int(cnt[t].sum())})
    ct = pd.DataFrame(rows)
    ct["fdr"] = multipletests(ct["p_value"], method="fdr_bh")[1]
    ct.to_csv(res_dir / "genotype" / "ec_subtype_composition_test.csv", index=False)
    print(ct.round(4).to_string(index=False), flush=True)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    box_by_genotype(ax, prop * 100, meta, EC_ORDER, fdr=ct.set_index("subtype")["fdr"])
    ax.set_ylabel("% of endothelial-lineage nuclei")
    ax.legend([Patch(facecolor=GCOL[g], edgecolor="#333333") for g in ["WT", "I4M"]], ["WT", "I4M"], frameon=True)
    fig.savefig(out / "ec_subtype_composition.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:3])
