"""WT vs TNNC1-I4M by sex: composition interaction test and figures (style as 15_genotype_figures_styled.py).

Composition: per-sample proportion, asin(sqrt(p)) ~ sex * genotype (OLS; 2 mice per group, 4 residual df);
reports the genotype effect within each sex and the sex x genotype interaction.
DE: 17_sex_genotype_de.R outputs.

Outputs (results_dir/genotype/):
  sex_composition_test.csv
  styled/sex_composition_interaction.png   interaction plots for the main cell types
  styled/sex_logfc_scatter.png             genotype logFC in males vs females (CM, endocardial, fibroblast)
  styled/sex_deg_counts.png                DEGs within each sex and sex x genotype interaction genes

Usage: python 18_sex_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
NAME = {"Cardiomyocyte": "Cardiomyocyte", "Fibroblast": "Fibroblast", "Endothelial": "Endothelial",
        "Endocardial": "Endocardial", "Lymphatic_EC": "Lymphatic endothelial", "Pericyte": "Pericyte", "SMC": "VSMC",
        "Macrophage": "Macrophage", "Lymphocyte": "Lymphocyte", "Neural_Schwann": "Neural / Schwann",
        "Dendritic": "Dendritic (cDC1)", "Epicardial": "Epicardial"}
SCOL = {"M": "#2B6CB0", "F": "#D53F8C"}
SLAB = {"M": "Male", "F": "Female"}


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    res_dir = Path(res_dir)
    out = res_dir / "genotype" / "styled"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    cnt = pd.read_csv(res_dir / "genotype" / "cell_counts.csv", index_col=0)
    prop = cnt.div(cnt.sum(1), axis=0)

    rows = []
    for t in prop.columns:
        d = meta.assign(p=prop[t], y=np.arcsin(np.sqrt(prop[t])))
        d["g"] = (d["genotype"] == "I4M").astype(int)
        d["f"] = (d["sex"] == "F").astype(int)
        fit = smf.ols("y ~ f * g", data=d).fit()
        r = {"cell_type": t, "p_interaction": fit.pvalues["f:g"],
             # genotype effect within each sex, from the same model (pooled residual variance, 4 df)
             "p_male": fit.pvalues["g"], "p_female": float(fit.t_test("g + f:g = 0").pvalue)}
        for sx in ["M", "F"]:
            m = d[d.sex == sx]
            r[f"WT_{sx}"] = m.loc[m.genotype == "WT", "p"].mean()
            r[f"I4M_{sx}"] = m.loc[m.genotype == "I4M", "p"].mean()
        rows.append(r)
    ct = pd.DataFrame(rows)
    ct["fdr_interaction"] = multipletests(ct["p_interaction"], method="fdr_bh")[1]
    ct["fdr_male"] = multipletests(ct["p_male"], method="fdr_bh")[1]
    ct["fdr_female"] = multipletests(ct["p_female"], method="fdr_bh")[1]
    ct["ratio_male"] = ct["I4M_M"] / ct["WT_M"]
    ct["ratio_female"] = ct["I4M_F"] / ct["WT_F"]

    # ---- sex-stratified composition: all cell types, WT vs I4M within males and within females ----
    gcol = {"WT": "#B8682E", "I4M": "#3F8F6E"}
    order = [c for c in NAME if c in prop.columns]
    fig, axes = plt.subplots(2, 1, figsize=(12, 7.6), sharex=True)
    for ax, sx in zip(axes, ["M", "F"]):
        for i, t in enumerate(order):
            for j, g in enumerate(["WT", "I4M"]):
                v = prop.loc[meta.index[(meta.genotype == g) & (meta.sex == sx)], t].values * 100
                x = i + (j - 0.5) * 0.38
                ax.bar(x, v.mean(), 0.34, color=gcol[g], edgecolor="#333333", lw=0.6)
                ax.scatter([x - 0.06, x + 0.06][:len(v)], v, s=14, color="black", zorder=3)
            rr = ct.loc[ct.cell_type == t].iloc[0]
            f_ = rr["fdr_male" if sx == "M" else "fdr_female"]
            top = prop.loc[meta.index[meta.sex == sx], t].max() * 100
            ax.text(i, top * 1.08 + 0.3, "*" if f_ <= 0.05 else "", ha="center", color="red", fontsize=13)
            if i:
                ax.axvline(i - 0.5, color="#e0e0e0", lw=0.7)
        ax.set_yscale("symlog", linthresh=0.1)
        ticks = [0, 0.1, 0.3, 1, 3, 10, 30]
        ax.set_yticks(ticks, [f"{v:g}%" for v in ticks])
        ax.set_ylim(0, 80)
        ax.set_ylabel("% of nuclei (log scale)")
        ax.set_title(f"{SLAB[sx]}s: WT (n = 2) vs I4M (n = 2)", fontsize=12)
    axes[-1].set_xticks(range(len(order)), [NAME[t] for t in order], rotation=40, ha="right")
    axes[0].legend([plt.Rectangle((0, 0), 1, 1, color=gcol[g]) for g in ["WT", "I4M"]], ["WT", "I4M"],
                   frameon=False, loc="upper right")
    fig.text(0.01, -0.02, "* FDR <= 0.05 for the genotype effect within that sex (asin-sqrt OLS ~ sex * genotype, "
             "BH over cell types); y axis: symmetric log scale", fontsize=8, color="#555555")
    fig.tight_layout()
    fig.savefig(out / "sex_composition_by_sex.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- same figure types as the pooled composition slide (box + points; stacked bar), split by sex ----
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.4), sharex=True)
    for ax, sx in zip(axes, ["M", "F"]):
        for i, t in enumerate(order):
            for j, g in enumerate(["WT", "I4M"]):
                v = prop.loc[meta.index[(meta.genotype == g) & (meta.sex == sx)], t].values
                x = i + (j - 0.5) * 0.38
                ax.boxplot(v, positions=[x], widths=0.34, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=gcol[g], color="#333333", lw=0.8), medianprops=dict(color="#333333"),
                           whiskerprops=dict(color="#333333", lw=0.8), capprops=dict(color="#333333", lw=0.8))
                ax.scatter(np.full(len(v), x) + np.linspace(-0.05, 0.05, len(v)), v, s=12, color="black", zorder=3)
            rr = ct.loc[ct.cell_type == t].iloc[0]
            if rr["fdr_male" if sx == "M" else "fdr_female"] <= 0.05:
                top = prop.loc[meta.index[meta.sex == sx], t].max() + 0.015
                ax.plot([i - 0.19, i + 0.19], [top, top], color="red", lw=1)
                ax.text(i, top, "*", color="red", ha="center", fontsize=14)
            if i:
                ax.axvline(i - 0.5, color="#d9d9d9", lw=0.8)
        ax.set_ylabel("Proportion")
        ax.set_title(f"{SLAB[sx]}s", fontsize=12)
        ax.set_xlim(-0.6, len(order) - 0.4)
    axes[-1].set_xticks(range(len(order)), [NAME[t] for t in order], rotation=40, ha="right")
    axes[0].legend([plt.Rectangle((0, 0), 1, 1, facecolor=gcol[g], edgecolor="#333333") for g in ["WT", "I4M"]],
                   ["WT", "I4M"], loc="upper right", frameon=True)
    fig.tight_layout()
    fig.savefig(out / "sex_composition_box.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    hue = ["#F8766D", "#DE8C00", "#B79F00", "#7CAE00", "#00BA38", "#00C08B", "#00BFC4", "#00B4F0", "#619CFF",
           "#C77CFF", "#F564E3", "#FF64B0"]
    tco = dict(zip(["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
                    "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"], hue))
    grp = [("WT", "M"), ("I4M", "M"), ("WT", "F"), ("I4M", "F")]
    pc = pd.DataFrame({f"{g}\n{SLAB[s]}": cnt.loc[meta.index[(meta.genotype == g) & (meta.sex == s)]].sum()
                       for g, s in grp}).T
    pc = pc.div(pc.sum(1), axis=0)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    bottom = np.zeros(len(pc))
    for t in list(tco)[::-1]:
        if t in pc.columns:
            ax.bar(range(len(pc)), pc[t].values, bottom=bottom, color=tco[t], width=0.75, label=NAME[t])
            bottom += pc[t].values
    ax.set_xticks(range(len(pc)), pc.index)
    ax.set_ylabel("Cell proportion")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"])
    hd, lb = ax.get_legend_handles_labels()
    ax.legend(hd[::-1], lb[::-1], title="Cell type", loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    fig.savefig(out / "sex_composition_stacked.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    ct.to_csv(res_dir / "genotype" / "sex_composition_test.csv", index=False)
    print(ct.round(4).to_string(), flush=True)
    ct = ct.set_index("cell_type")

    # ---- interaction plots ----
    types = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Pericyte", "Macrophage"]
    fig, axes = plt.subplots(1, len(types), figsize=(3.0 * len(types), 3.6))
    for ax, t in zip(axes, types):
        for sx in ["M", "F"]:
            means = []
            for j, g in enumerate(["WT", "I4M"]):
                v = prop.loc[meta.index[(meta.genotype == g) & (meta.sex == sx)], t].values * 100
                ax.scatter(np.full(len(v), j) + (-0.06 if sx == "M" else 0.06), v, color=SCOL[sx], s=28, zorder=3,
                           edgecolor="white", lw=0.5)
                means.append(v.mean())
            ax.plot([0, 1], means, color=SCOL[sx], lw=2)
        ax.set_xticks([0, 1], ["WT", "I4M"])
        ax.set_xlim(-0.4, 1.4)
        ax.set_title(NAME[t] + "\ninteraction p = " + f"{ct.loc[t, 'p_interaction']:.2g}", fontsize=11)
    axes[0].set_ylabel("Proportion of nuclei (%)")
    h = [Line2D([], [], color=SCOL[s], marker="o", lw=2, label=SLAB[s]) for s in ["M", "F"]]
    fig.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(out / "sex_composition_interaction.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- logFC male vs female ----
    summ = pd.read_csv(res_dir / "genotype" / "sex_summary.csv").set_index("cell_type")
    st = ["Cardiomyocyte", "Endocardial", "Fibroblast"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    for ax, t in zip(axes, st):
        d = pd.read_csv(wd / f"sex_{t}.csv.gz")
        bg = d[~d.overall_robust & ~d.flagged]
        ax.scatter(bg.logFC_male, bg.logFC_female, s=3, c="#D0D0D0", lw=0, rasterized=True)
        r = d[d.overall_robust]
        ax.scatter(r.logFC_male, r.logFC_female, s=10, c=np.where(r.overall_logFC > 0, "#D7263D", "#1B2A7A"), lw=0)
        ia = d[(d.FDR_interaction <= 0.1) & ~d.flagged]
        ax.scatter(ia.logFC_male, ia.logFC_female, s=40, facecolors="none", edgecolors="#E69F00", lw=1.2)
        for g in ia.itertuples():
            ax.annotate(g.gene, (g.logFC_male, g.logFC_female), fontsize=7, xytext=(3, 3), textcoords="offset points")
        lim = np.nanpercentile(np.abs(np.r_[d.logFC_male, d.logFC_female]), 99.8)
        ax.plot([-lim, lim], [-lim, lim], color="#888888", ls="--", lw=0.8)
        ax.axhline(0, color="#BBBBBB", lw=0.6); ax.axvline(0, color="#BBBBBB", lw=0.6)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        s_ = summ.loc[t]
        ax.set_title(f"{NAME[t]}\n{int(s_.n_overall_robust)} DEGs: same direction {s_.same_direction_frac * 100:.0f}%, "
                     f"r = {s_.r_logFC_overall_degs:.2f}", fontsize=11)
        ax.set_xlabel("log2FoldChange in males (I4M vs WT)")
        ax.set_ylabel("log2FoldChange in females (I4M vs WT)")
    h = [Line2D([], [], marker="o", ls="", color="#D7263D", label="DEG up in I4M (both sexes pooled)"),
         Line2D([], [], marker="o", ls="", color="#1B2A7A", label="DEG down in I4M"),
         Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="#E69F00",
                label="sex x genotype interaction (FDR <= 0.1)"),
         Line2D([], [], marker="o", ls="", color="#D0D0D0", label="other genes")]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    fig.savefig(out / "sex_logfc_scatter.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---- DEG counts within each sex + interaction ----
    order = [t for t in NAME if t in summ.index]
    s = summ.reindex(order)
    x = np.arange(len(s))
    w = 0.27
    fig, ax = plt.subplots(figsize=(10, 4.4))
    ax.bar(x - w, s.male_up + s.male_down, w, color=SCOL["M"], label="Male: I4M vs WT")
    ax.bar(x, s.female_up + s.female_down, w, color=SCOL["F"], label="Female: I4M vs WT")
    ax.bar(x + w, s.interaction, w, color="#E69F00", label="Sex x genotype interaction")
    for i, t in enumerate(order):
        ax.text(i - w, s.loc[t, "male_up"] + s.loc[t, "male_down"] + 1, int(s.loc[t, "male_up"] + s.loc[t, "male_down"]),
                ha="center", fontsize=8)
        ax.text(i, s.loc[t, "female_up"] + s.loc[t, "female_down"] + 1,
                int(s.loc[t, "female_up"] + s.loc[t, "female_down"]), ha="center", fontsize=8)
        ax.text(i + w, s.loc[t, "interaction"] + 1, int(s.loc[t, "interaction"]), ha="center", fontsize=8)
    ax.set_xticks(x, [NAME[t] for t in order], rotation=40, ha="right")
    ax.set_ylabel("# differentially\nexpressed genes (FDR <= 0.1)")
    ax.legend(frameon=True)
    fig.savefig(out / "sex_deg_counts.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
