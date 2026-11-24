"""WT vs TNNC1-I4M (+/-): figures from the pseudobulk DE (13_pseudobulk_de.R).

Outputs (results_dir/genotype/):
  fig_confound.png        per-sample technical metrics by genotype (intronic / intergenic reads, ambient signature,
                          atrial genes) - genotype is confounded with library quality and dissection
  fig_deg_counts.png      number of DEGs per cell type (FDR < 0.05), robust vs flagged as technical / ambient / shared
  fig_volcano.png         volcano plots for six cell types
  fig_go.png              GO BP over-representation of robust DEGs (cell types with >= 15 robust DEGs)
  fig_panel_cm.png        curated cardiomyocyte genes (sarcomere, calcium, stress, metabolism): logFC and FDR
  fig_panel_fib.png       fibroblast-activation genes in fibroblasts: per-sample CPM
  lncrna_check.csv        share of Gm / Rik (lncRNA-like) genes among robust up / down CM DEGs vs background

Usage: python 14_genotype_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "Arial"
GCOL = {"WT": "#7F7F7F", "I4M": "#C0392B"}
UP, DOWN, FLAG = "#C0392B", "#2E6DA4", "#E0B060"
LABEL = {"Cardiomyocyte": "Cardiomyocyte", "Fibroblast": "Fibroblast", "Endothelial": "Endothelial",
         "Endocardial": "Endocardial", "Lymphatic_EC": "Lymphatic EC", "Pericyte": "Pericyte", "SMC": "SMC",
         "Macrophage": "Macrophage", "Lymphocyte": "Lymphocyte", "Neural_Schwann": "Neural / Schwann"}
ORDER = list(LABEL)
# Cell Ranger summary per sample (provider delivery, as listed in the project slides)
CR = pd.DataFrame({"intronic": [23.1, 24.8, 23.1, 20.8, 32.2, 31.3, 29.4, 36.0],
                   "intergenic": [22.2, 20.1, 22.8, 23.9, 17.7, 18.2, 17.6, 14.3]}, index=list("ABCDEFGH"))
AMB = ["Gm42418", "Cmss1", "Lars2", "Camk1d", "Gphn", "Gm19951", "Cdk8", "AY036118"]


def cpm(x):
    return x / x.sum() * 1e6


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype"
    meta = pd.read_csv(wd / "meta.csv", index_col=0)
    gt = meta["genotype"]
    de = {t: pd.read_csv(wd / f"de_{t}.csv.gz") for t in ORDER if (wd / f"de_{t}.csv.gz").exists()}
    summ = pd.read_csv(out / "de_summary.csv").set_index("cell_type")

    # ---- confound panel ----
    pbc = pd.read_csv(wd / "pb_Cardiomyocyte.csv.gz", index_col=0)
    pbf = pd.read_csv(wd / "pb_Fibroblast.csv.gz", index_col=0)
    allc = sum(pd.read_csv(wd / f"pb_{t}.csv.gz", index_col=0) for t in de)
    metrics = pd.DataFrame({
        "Intronic reads (%)": CR["intronic"], "Intergenic reads (%)": CR["intergenic"],
        "Ambient / low-quality signature\n(CPM, all nuclei)": cpm(allc).loc[AMB].sum(),
        "Atrial genes Myl7 + Sln\n(CPM, cardiomyocytes)": cpm(pbc).loc[["Myl7", "Sln"]].sum()})
    metrics.to_csv(out / "sample_technical_metrics.csv")
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, c in zip(axes, metrics.columns):
        ax.bar(metrics.index, metrics[c], color=[GCOL[gt[s]] for s in metrics.index])
        ax.set_title(c, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    h = [plt.Rectangle((0, 0), 1, 1, color=GCOL[g]) for g in ["WT", "I4M"]]
    fig.legend(h, ["WT (A-D)", "TNNC1-I4M (E-H)"], loc="upper right", frameon=False, ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(out / "fig_confound.png", dpi=200)
    plt.close(fig)

    # ---- DEG counts ----
    s = summ.reindex([t for t in ORDER if t in summ.index])
    x = np.arange(len(s))
    fig, ax = plt.subplots(figsize=(10, 4.4))
    ax.bar(x, s["up_robust"], color=UP, label="Up in I4M (robust)")
    ax.bar(x, s["up"] - s["up_robust"], bottom=s["up_robust"], color=UP, alpha=0.3, hatch="//",
           edgecolor="white", label="Up, flagged")
    ax.bar(x, -s["down_robust"], color=DOWN, label="Down in I4M (robust)")
    ax.bar(x, -(s["down"] - s["down_robust"]), bottom=-s["down_robust"], color=DOWN, alpha=0.3, hatch="//",
           edgecolor="white", label="Down, flagged")
    for i, t in enumerate(s.index):
        ax.text(i, s.loc[t, "up"] + 8, str(s.loc[t, "up"]), ha="center", fontsize=9)
        ax.text(i, -s.loc[t, "down"] - 8, str(s.loc[t, "down"]), ha="center", va="top", fontsize=9)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x, [LABEL[t] for t in s.index], rotation=35, ha="right")
    ax.set_ylabel("Number of DEGs (FDR < 0.05)")
    lim = max(s["up"].max(), s["down"].max()) * 1.25
    ax.set_ylim(-lim, lim)
    ax.set_yticks(ax.get_yticks(), [str(abs(int(v))) for v in ax.get_yticks()])
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "fig_deg_counts.png", dpi=200)
    plt.close(fig)

    # ---- volcano ----
    vt = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Pericyte", "Macrophage"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, t in zip(axes.ravel(), vt):
        d = de[t].copy()
        d["nl"] = -np.log10(d["PValue"].clip(1e-300))
        sig = d["FDR"] < 0.05
        flag = sig & ~d["robust"]
        ax.scatter(d.loc[~sig, "logFC"], d.loc[~sig, "nl"], s=3, c="#cfcfcf", rasterized=True)
        ax.scatter(d.loc[flag, "logFC"], d.loc[flag, "nl"], s=8, c=FLAG, label="FDR<0.05, flagged")
        r = d[d["robust"]]
        ax.scatter(r["logFC"], r["nl"], s=9, c=np.where(r["logFC"] > 0, UP, DOWN), label="FDR<0.05, robust")
        lab = pd.concat([r[r.logFC > 0].nsmallest(6, "PValue"), r[r.logFC < 0].nsmallest(6, "PValue")])
        lab = lab[~lab["gene"].str.match(r"^(Gm\d|.*Rik$)")]
        for _, g in lab.iterrows():
            ax.annotate(g["gene"], (g["logFC"], g["nl"]), fontsize=7, xytext=(3, 2), textcoords="offset points")
        n_up, n_dn = int((r.logFC > 0).sum()), int((r.logFC < 0).sum())
        ax.set_title(f"{LABEL[t]}   robust up {n_up} / down {n_dn}", fontsize=11)
        ax.set_xlabel("log2 fold change (I4M vs WT)")
        ax.set_ylabel("-log10 P")
        ax.axvline(0, color="#999999", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "fig_volcano.png", dpi=200)
    plt.close(fig)

    # ---- GO ----
    go = pd.read_csv(out / "go_top.csv")
    blocks = [(t, d) for t, d in go[["cell_type", "direction"]].drop_duplicates().itertuples(index=False)
              if (go[(go.cell_type == t) & (go.direction == d)]["FDR"] < 0.25).any()]
    fig, axes = plt.subplots(1, len(blocks), figsize=(6.5 * len(blocks), 5.2), gridspec_kw={"wspace": 1.5})
    axes = np.atleast_1d(axes)
    for ax, (t, d) in zip(axes, blocks):
        g = go[(go.cell_type == t) & (go.direction == d)].head(10).iloc[::-1]
        y = np.arange(len(g))
        sc_ = ax.scatter(g["DE"] / g["N"], y, s=g["DE"] * 22, c=-np.log10(g["FDR"]), cmap="Reds",
                         vmin=0, vmax=3, edgecolor="k", lw=0.3)
        ax.set_yticks(y, [s_[:48] + ("..." if len(s_) > 48 else "") for s_ in g["term"]], fontsize=8)
        ax.set_xlabel("Gene ratio (DEGs / genes in term)")
        ax.set_title(f"{LABEL[t]} - {d} in I4M\n({int(g['n_deg'].iloc[0])} robust DEGs)", fontsize=10)
        ax.margins(x=0.25, y=0.06)
    cb = fig.colorbar(sc_, ax=axes, fraction=0.02)
    cb.set_label("-log10 FDR")
    fig.savefig(out / "fig_go.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- curated CM panel ----
    groups = {"Thin filament / sarcomere": ["Tnnc1", "Tnni3", "Tnnt2", "Tpm1", "Actc1", "Myh6", "Myh7", "Myl2",
                                             "Myl3", "Mybpc3", "Ttn", "Actn2"],
              "Calcium handling": ["Atp2a2", "Pln", "Ryr2", "Cacna1c", "Slc8a1", "Casq2", "Camk2d", "Trdn"],
              "Stress / hypertrophy": ["Nppa", "Nppb", "Acta1", "Ankrd1", "Ccn2", "Rcan1", "Fhl1", "Atf3", "Myh7b"],
              "Metabolism": ["Ppargc1a", "Ppara", "Cpt1b", "Acadm", "Hadha", "Pdk4", "Ucp3", "Slc2a4", "Slc2a1",
                             "Hk2", "Ldha"]}
    d = de["Cardiomyocyte"].set_index("gene")
    rows = [(grp, g) for grp, gl in groups.items() for g in gl if g in d.index]
    fig, ax = plt.subplots(figsize=(13, 3.6))
    xs = np.arange(len(rows))
    lf = np.array([d.loc[g, "logFC"] for _, g in rows])
    fd = np.array([d.loc[g, "FDR"] for _, g in rows])
    ax.bar(xs, lf, color=np.where(lf > 0, UP, DOWN), alpha=np.where(fd < 0.05, 1, 0.35).tolist() and 1)
    for i, (f, l) in enumerate(zip(fd, lf)):
        if f < 0.05:
            ax.text(i, l + (0.05 if l > 0 else -0.05), "*", ha="center", va="bottom" if l > 0 else "top", fontsize=12)
        elif f < 0.1:
            ax.text(i, l + (0.05 if l > 0 else -0.05), "·", ha="center", va="bottom" if l > 0 else "top", fontsize=14)
    ax.set_xticks(xs, [g for _, g in rows], rotation=60, ha="right", fontsize=8, style="italic")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_ylabel("log2FC in cardiomyocytes\n(I4M vs WT)")
    pos = 0
    for grp, gl in groups.items():
        n = sum(g in d.index for g in gl)
        ax.text(pos + n / 2 - 0.5, ax.get_ylim()[1] * 0.98, grp, ha="center", va="top", fontsize=9, weight="bold")
        if pos:
            ax.axvline(pos - 0.5, color="#bbbbbb", lw=0.8)
        pos += n
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(min(lf.min() * 1.3, -0.5), max(lf.max() * 1.6, 0.6))
    fig.tight_layout()
    fig.savefig(out / "fig_panel_cm.png", dpi=200)
    plt.close(fig)

    # ---- fibroblast activation genes: per-sample CPM ----
    fg = ["Postn", "Cthrc1", "Col1a1", "Col3a1", "Fn1", "Acta2", "Lox", "Timp1", "Ccn2", "Comp", "Fap", "Meox1"]
    fg = [g for g in fg if g in pbf.index]
    c = cpm(pbf).loc[fg]
    dfi = de["Fibroblast"].set_index("gene")
    fig, axes = plt.subplots(2, 6, figsize=(15, 5.4))
    for ax, g in zip(axes.ravel(), fg):
        for j, grp in enumerate(["WT", "I4M"]):
            v = c.loc[g, meta.index[gt == grp]].values
            ax.bar(j, v.mean(), color=GCOL[grp], alpha=0.4, width=0.6)
            ax.scatter(np.full(len(v), j) + np.linspace(-0.12, 0.12, len(v)), v, color=GCOL[grp], s=18, zorder=3)
        fdr = dfi.loc[g, "FDR"] if g in dfi.index else np.nan
        ax.set_title(f"{g}\nFDR {fdr:.2f}" if fdr == fdr else f"{g}\n(not tested)", fontsize=9, style="italic")
        ax.set_xticks([0, 1], ["WT", "I4M"])
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_ylabel("CPM (fibroblast pseudobulk)")
    axes[1, 0].set_ylabel("CPM (fibroblast pseudobulk)")
    fig.tight_layout()
    fig.savefig(out / "fig_panel_fib.png", dpi=200)
    plt.close(fig)

    # ---- lncRNA-like gene check (CM) ----
    d = de["Cardiomyocyte"]
    lnc = d["gene"].str.match(r"^(Gm\d|.*Rik$|.*os$)")
    rows = []
    for name, m in [("background (tested)", np.ones(len(d), bool)), ("robust up", d.robust & (d.logFC > 0)),
                    ("robust down", d.robust & (d.logFC < 0))]:
        rows.append({"set": name, "n": int(m.sum()), "frac_lncRNA_like": float(lnc[m].mean())})
    pd.DataFrame(rows).to_csv(out / "lncrna_check.csv", index=False)
    print(pd.DataFrame(rows), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
