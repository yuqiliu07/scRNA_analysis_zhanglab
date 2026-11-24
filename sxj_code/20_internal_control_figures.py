"""Figures for the internal-control analysis (19_internal_control.R).

For the robust DEGs of a target cell type (cardiomyocyte, endocardial):
  styled/ic_scatter.png     genotype logFC in the target vs in the reference cell type (fibroblast, endothelial),
                            testable genes only; points near the diagonal change equally in both cell types
  styled/ic_categories.png  robust DEGs split into: not testable (not expressed by the references),
                            same-direction change in the references (batch-like), target-specific evidence
                            (interaction P < 0.05 against every testable reference, same sign; none reaches FDR < 0.05)
  ic_categories.csv, ic_candidates.csv

Usage: python 20_internal_control_figures.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
TG = ["Cardiomyocyte", "Endocardial"]
RF = ["Fibroblast", "Endothelial"]
CAT_COL = {"target-specific evidence": "#D7263D", "same direction in references": "#8C8C8C",
           "other / inconclusive": "#C9C9C9", "not testable": "#E8C27A"}


def categorise(d):
    r = d[d.robust].copy()
    tested = [r[f"P_int_vs_{f}"].notna() for f in RF]
    anyt = tested[0] | tested[1]
    ok = np.ones(len(r), bool)
    for f, t in zip(RF, tested):
        good = (r[f"P_int_vs_{f}"] < 0.05) & (np.sign(r[f"logFC_int_vs_{f}"]) == np.sign(r["logFC_target"]))
        ok &= (~t) | good
    same = np.zeros(len(r), bool)
    for f in RF:
        lf = r[f"logFC_in_{f}"]
        same |= r[f"P_int_vs_{f}"].notna() & (np.sign(lf) == np.sign(r["logFC_target"])) & (lf.abs() >= 0.5 * r["logFC_target"].abs())
    cat = np.where(~anyt, "not testable", np.where(ok, "target-specific evidence",
                   np.where(same, "same direction in references", "other / inconclusive")))
    r["ic_category"] = cat
    return r


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    res_dir = Path(res_dir)
    out = res_dir / "genotype" / "styled"
    cats, cands, rows = {}, [], []
    for tg in TG:
        r = categorise(pd.read_csv(wd / f"ic_{tg}.csv.gz"))
        cats[tg] = r
        vc = r["ic_category"].value_counts()
        rows.append({"target": tg, "n_robust": len(r), **{k: int(vc.get(k, 0)) for k in CAT_COL}})
        c = r[r.ic_category == "target-specific evidence"].copy()
        c["target"] = tg
        cands.append(c[["target", "gene", "logFC_target", "FDR_target", "logFC_in_Fibroblast", "logFC_in_Endothelial",
                        "P_int_vs_Fibroblast", "P_int_vs_Endothelial", "FDR_int_vs_Fibroblast", "FDR_int_vs_Endothelial"]])
    pd.DataFrame(rows).to_csv(res_dir / "genotype" / "ic_categories.csv", index=False)
    cand = pd.concat(cands).sort_values(["target", "P_int_vs_Fibroblast"])
    cand.to_csv(res_dir / "genotype" / "ic_candidates.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    print(cand.round(3).to_string(index=False), flush=True)

    # scatter target vs reference
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))
    k = 0
    for tg in TG:
        r = cats[tg]
        for f in RF:
            ax = axes[k]; k += 1
            t = r[r[f"P_int_vs_{f}"].notna()]
            col = [CAT_COL[c] if c != "not testable" else "#C9C9C9" for c in t.ic_category]
            ax.scatter(t.logFC_target, t[f"logFC_in_{f}"].fillna(0), s=14, c=col, lw=0)
            lim = np.nanmax(np.abs(np.r_[t.logFC_target, t[f"logFC_in_{f}"].fillna(0)])) * 1.05
            ax.plot([-lim, lim], [-lim, lim], ls="--", color="#666666", lw=0.8)
            ax.axhline(0, color="#BBBBBB", lw=0.6); ax.axvline(0, color="#BBBBBB", lw=0.6)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            rr = np.corrcoef(t.logFC_target, t[f"logFC_in_{f}"].fillna(0))[0, 1]
            same = (np.sign(t.logFC_target) == np.sign(t[f"logFC_in_{f}"].fillna(0))).mean()
            ax.set_title(f"{tg} DEGs vs {f}\n{len(t)} genes, same direction {same * 100:.0f}%, r = {rr:.2f}", fontsize=10.5)
            ax.set_xlabel(f"log2FC in {tg} (I4M vs WT)")
            ax.set_ylabel(f"log2FC in {f} (I4M vs WT)")
            for g in t[t.ic_category == "target-specific evidence"].itertuples():
                ax.annotate(g.gene, (g.logFC_target, getattr(g, f"logFC_in_{f}")), fontsize=7, xytext=(3, 3),
                            textcoords="offset points")
    h = [Line2D([], [], marker="o", ls="", color=CAT_COL[c], label=c) for c in
         ["target-specific evidence", "same direction in references", "other / inconclusive"]]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout()
    fig.savefig(out / "ic_scatter.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # category bars
    df = pd.DataFrame(rows).set_index("target")
    fig, ax = plt.subplots(figsize=(8, 3.2))
    left = np.zeros(len(df))
    for c in ["target-specific evidence", "same direction in references", "other / inconclusive", "not testable"]:
        ax.barh(range(len(df)), df[c], left=left, color=CAT_COL[c], label=c, edgecolor="white")
        for i, v in enumerate(df[c]):
            if v >= 8:
                ax.text(left[i] + v / 2, i, str(v), ha="center", va="center", fontsize=9)
        left += df[c].values
    ax.set_yticks(range(len(df)), [f"{t}\n({n} DEGs)" for t, n in zip(df.index, df.n_robust)])
    ax.invert_yaxis()
    ax.set_xlabel("Robust DEGs (I4M vs WT)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2)
    fig.savefig(out / "ic_categories.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:3])
