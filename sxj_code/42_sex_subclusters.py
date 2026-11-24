"""Subcluster composition by sex: is the genotype effect on subcluster proportions different in males and females?

For each of CM / endocardial / fibroblast / endothelial, per-mouse proportion of every named subcluster within its cell
type (31_subclusters.py, names from 32_subcluster_figures.py; low-quality / doublet-like / "Minor" groups left out).
Test per subcluster (as 18_sex_figures.py for cell types): asin-sqrt(proportion) ~ sex * genotype, OLS on the 8 mice;
reports the WT -> I4M change within each sex and the sex x genotype interaction p (BH within cell type).
With 2 mice per sex and genotype these are descriptive; only large, opposite changes can reach significance.

Outputs: results_dir/genotype/styled/sex_subcluster_interaction.png   interaction plots (male / female lines)
         results_dir/genotype/sex_subcluster_test.csv
Usage: python 42_sex_subclusters.py <work_dir> <results_dir>
"""
import importlib.util
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
TYPES = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
SKIP = ("Low quality", "doublet", "Minor", "Glia-like")
SCOL = {"M": "#2E6DA4", "F": "#D6609A"}
SLAB = {"M": "Male", "F": "Female"}

spec = importlib.util.spec_from_file_location("subfig", Path(__file__).with_name("32_subcluster_figures.py"))
subfig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subfig)


def main(work_dir, res_dir):
    sd = Path(work_dir) / "subclusters"
    rd = Path(res_dir) / "genotype"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    rows, panels = [], []
    for t in TYPES:
        s = pd.read_csv(sd / f"sub_{t}.csv.gz", index_col=0)
        s["sub"] = s["sub"].astype(str)
        s["name"] = [subfig.NAMES[t].get(k, (f"{t[:3]}_{k}", []))[0] for k in s["sub"]]
        prop = pd.crosstab(s["sample"], s["name"], normalize="index").reindex(meta.index)
        keep = [n for n in prop.columns if not any(w in n for w in SKIP)]
        keep = sorted(keep, key=lambda n: -prop[n].mean())
        res = []
        for n in keep:
            d = meta.join(prop[n].rename("f"))
            d["y"] = np.arcsin(np.sqrt(d.f))
            d["g"] = (d.genotype == "I4M").astype(int); d["fm"] = (d.sex == "F").astype(int)
            fit = smf.ols("y ~ fm * g", d).fit()
            r = {"cell_type": t, "subcluster": n, "p_interaction": fit.pvalues["fm:g"]}
            for sx in ["M", "F"]:
                for g in ["WT", "I4M"]:
                    r[f"mean_{g}_{sx}"] = d.f[(d.genotype == g) & (d.sex == sx)].mean()
            res.append(r)
        res = pd.DataFrame(res)
        res["fdr_interaction"] = multipletests(res.p_interaction, method="fdr_bh")[1]
        rows.append(res)
        panels += [(t, n, prop[n]) for n in keep]
    R = pd.concat(rows)
    R.to_csv(rd / "sex_subcluster_test.csv", index=False)

    # ---- one row per cell type, one small interaction plot per subcluster ----
    ncol = max(sum(1 for p in panels if p[0] == t) for t in TYPES)
    fig, axes = plt.subplots(len(TYPES), ncol, figsize=(2.15 * ncol, 2.45 * len(TYPES)), squeeze=False)
    Ri = R.set_index(["cell_type", "subcluster"])
    for i, t in enumerate(TYPES):
        pt = [p for p in panels if p[0] == t]
        for j in range(ncol):
            ax = axes[i, j]
            if j >= len(pt):
                ax.axis("off"); continue
            _, n, v = pt[j]
            for sx in ["M", "F"]:
                means = []
                for k, g in enumerate(["WT", "I4M"]):
                    vv = v[meta.index[(meta.genotype == g) & (meta.sex == sx)]].values * 100
                    ax.scatter(np.full(len(vv), k) + (-0.07 if sx == "M" else 0.07), vv, color=SCOL[sx], s=18,
                               zorder=3, edgecolor="white", lw=0.4)
                    means.append(vv.mean())
                ax.plot([0, 1], means, color=SCOL[sx], lw=1.8)
            ax.set_xticks([0, 1], ["WT", "I4M"], fontsize=8.5)
            ax.tick_params(axis="y", labelsize=8)
            ax.set_xlim(-0.4, 1.4)
            p = Ri.loc[(t, n), "p_interaction"]
            nm = n.replace(" (main)", "").replace(" (EC-like)", "")
            nm = nm if len(nm) <= 22 else nm[:20] + "..."
            ax.set_title(f"{nm}\np = {p:.2g}", fontsize=9, color="#C0392B" if p <= 0.05 else "#1F1F1F")
        axes[i, 0].set_ylabel(f"{t}\n% of nuclei", fontsize=10)
    h = [Line2D([], [], color=SCOL[s], marker="o", lw=2, label=SLAB[s]) for s in ["M", "F"]]
    fig.legend(handles=h, loc="upper right", ncol=2, frameon=False, fontsize=10)
    fig.suptitle("Subcluster proportions by sex (p = sex x genotype interaction; red: p <= 0.05)", fontsize=12, x=0.4)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(rd / "styled" / "sex_subcluster_interaction.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(R[["cell_type", "subcluster", "p_interaction", "fdr_interaction"]].round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
