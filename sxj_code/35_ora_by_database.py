"""ORA results per database (KEGG, Reactome, WikiPathways, MSigDB Hallmark; GO BP for comparison): one figure per database.

Input: work_dir/genotype/pw_ORA_<db>_<type>_<Up|Down>.csv.gz (23_pathways.R; hypergeometric test of robust DEGs,
       FDR <= 0.1 and |log2FC| >= 0.25, against the genes tested in that cell type; BH within each comparison).
Plot: clusterProfiler-style dot plot, columns = comparisons (cell type x direction with >= 15 DEGs), rows = the top 5
      sets of each comparison (by p); dot size = gene ratio (overlap / input DEGs), colour = padj; ring = padj <= 0.05.
Outputs: results_dir/genotype/styled/ora_<db>.png; results_dir/genotype/ora_by_database_top.csv
Usage: python 35_ora_by_database.py <work_dir> <results_dir>
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial"})
DBS = {"KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WikiPathways", "GO_BP": "GO biological process",
       "Hallmark": "MSigDB Hallmark"}
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
SHORT = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC"}
TOP = 5


def clean(s, db):
    s = re.sub(r"^(REACTOME|WP|HALLMARK)_", "", s) if db in ("Reactome", "WikiPathways", "Hallmark") else s
    if db in ("Reactome", "WikiPathways", "Hallmark"):
        s = s.replace("_", " ").lower()
        s = s[0].upper() + s[1:]
    return s if len(s) <= 55 else s[:53] + "..."


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    rd = Path(res_dir) / "genotype"
    alltop = []
    for db, dbl in DBS.items():
        rows, groups = [], []
        for ct in CTS:
            de = pd.read_csv(wd / f"de_{ct}.csv.gz")
            for d in ["Up", "Down"]:
                f = wd / f"pw_ORA_{db}_{ct}_{d}.csv.gz"
                if not f.exists():
                    continue
                n_in = int((de.robust & ((de.logFC > 0) if d == "Up" else (de.logFC < 0))).sum())
                g = f"{SHORT[ct]} {'up' if d == 'Up' else 'down'}\n(n={n_in})"
                groups.append(g)
                t = pd.read_csv(f).sort_values("p")
                for r in t.head(TOP).itertuples():
                    rows.append({"group": g, "cell_type": ct, "direction": d, "set": r.set, "term": clean(r.set, db),
                                 "ratio": r.overlap / n_in, "overlap": r.overlap, "size": r.size, "p": r.p,
                                 "padj": r.padj, "genes": r.genes})
        o = pd.DataFrame(rows)
        o.insert(0, "database", db)
        alltop.append(o)
        terms = list(dict.fromkeys(o.term))
        # also show each top term's value in the other comparisons (so a row can be read across)
        full = []
        for g in groups:
            ct = [c for c in CTS if g.startswith(SHORT[c] + " ")][0]
            d = "Up" if " up" in g else "Down"
            t = pd.read_csv(wd / f"pw_ORA_{db}_{ct}_{d}.csv.gz")
            t["term"] = [clean(s, db) for s in t.set]
            t = t[t.term.isin(terms)]
            n_in = int(g.split("n=")[1].rstrip(")"))
            for r in t.itertuples():
                if r.overlap > 0:
                    full.append({"group": g, "term": r.term, "ratio": r.overlap / n_in, "padj": r.padj})
        F = pd.DataFrame(full)
        fig, ax = plt.subplots(figsize=(1.25 * len(groups) + 5.5, 0.3 * len(terms) + 1.8))
        cm = LinearSegmentedColormap.from_list("ora", ["#D7263D", "#F2A3A8", "#3B5BA9"])
        x = [groups.index(g) for g in F.group]
        y = [terms.index(t) for t in F.term]
        sc_ = ax.scatter(x, y, s=(F.ratio * 1800).clip(6, 400), c=np.log10(F.padj.clip(1e-6)), cmap=cm, vmin=-4, vmax=0,
                         edgecolors=["black" if p <= 0.05 else "#BBBBBB" for p in F.padj],
                         linewidths=[1.1 if p <= 0.05 else 0.4 for p in F.padj])
        ax.set_xticks(range(len(groups)), groups, fontsize=9)
        ax.set_yticks(range(len(terms)), terms, fontsize=8.5)
        ax.set_ylim(len(terms) - 0.5, -0.5)
        ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.grid(color="#EEEEEE", lw=0.6); ax.set_axisbelow(True)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.03, pad=0.02)
        cb.set_ticks([-4, -3, -2, -1, 0], labels=["1e-4", "1e-3", "0.01", "0.1", "1"])
        cb.set_label("p.adjust")
        hs = [plt.scatter([], [], s=v * 1800, c="#999999") for v in [0.02, 0.05, 0.1]]
        hs.append(plt.scatter([], [], s=80, c="white", edgecolors="black", linewidths=1.1))
        ax.legend(hs, ["0.02", "0.05", "0.10", "padj <= 0.05"], title="GeneRatio", loc="upper left",
                  bbox_to_anchor=(1.2, 1), frameon=False, labelspacing=1.2)
        ax.set_title(f"{dbl}: over-representation of DEGs (top {TOP} per comparison)", fontsize=11)
        fig.savefig(rd / "styled" / f"ora_{db}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)
        sig = o[o.padj <= 0.05]
        print(db, "groups:", len(groups), "terms:", len(terms), "| significant among top:",
              sig.groupby("group").size().to_dict(), flush=True)
    pd.concat(alltop).drop(columns="term").to_csv(rd / "ora_by_database_top.csv", index=False)


if __name__ == "__main__":
    main(*sys.argv[1:3])
