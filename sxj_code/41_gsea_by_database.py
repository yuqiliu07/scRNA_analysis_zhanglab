"""GSEA results per database (GO BP, KEGG, Reactome, WikiPathways, MSigDB Hallmark): one dot plot per database.

Input: work_dir/genotype/pw_GSEA_<db>_<type>.csv.gz (23_pathways.R; fgsea on all tested genes ranked by
       sign(logFC) * -log10(P), technical / ambient genes removed; sets of 15-500 genes; BH within each cell type).
Plot: columns = the four main cell types; rows = the 3 most significant shared sets, then the top 4 cell-type-specific
      significant sets (padj <= 0.05, by P) of each cell type;
      colour = NES (red = up in I4M, blue = down), size = -log10(padj), ring = padj <= 0.05.
      Sets significant in the same direction in >= 3 of the 10 tested cell types are labelled grey with "(shared)":
      they behave like a batch / systemic effect rather than a cell-type-specific one (marked "*", grey).
Outputs: results_dir/genotype/styled/gsea_<db>.png; results_dir/genotype/gsea_by_database_top.csv
Usage: python 41_gsea_by_database.py <work_dir> <results_dir>
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial"})
DBS = {"GO_BP": "GO biological process", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WikiPathways",
       "Hallmark": "MSigDB Hallmark"}
ALL = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "Pericyte", "SMC", "Macrophage",
       "Lymphocyte", "Neural_Schwann"]
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
SHORT = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC"}
TOP, NSHARED = 4, 3


def clean(s, db):
    if db in ("Reactome", "WikiPathways", "Hallmark"):
        s = re.sub(r"^(REACTOME|WP|HALLMARK)_", "", s).replace("_", " ").lower()
        s = s[0].upper() + s[1:]
    s = re.sub(r" \(GO:\d+\)$", "", s)
    return s if len(s) <= 38 else s[:36] + "..."


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    rd = Path(res_dir) / "genotype"
    alltop = []
    for db, dbl in DBS.items():
        tabs = {ct: pd.read_csv(wd / f"pw_GSEA_{db}_{ct}.csv.gz").set_index("set") for ct in ALL
                if (wd / f"pw_GSEA_{db}_{ct}.csv.gz").exists()}
        sig = pd.concat([t[t.padj <= 0.05].assign(ct=ct, sgn=np.sign(t[t.padj <= 0.05].NES)) for ct, t in tabs.items()])
        nshare = sig.reset_index().groupby(["set", "sgn"]).ct.nunique()
        is_sh = lambda s: max(nshare.get((s, 1.0), 0), nshare.get((s, -1.0), 0)) >= 3
        sh_p = {}
        for ct in CTS:
            t = tabs[ct]
            for s_, p_ in t[t.padj <= 0.05].p.items():
                if is_sh(s_):
                    sh_p[s_] = min(p_, sh_p.get(s_, 1))
        rows = sorted(sh_p, key=sh_p.get)[:NSHARED]
        for ct in CTS:
            t = tabs[ct]
            t = t[(t.padj <= 0.05) & ~t.index.map(is_sh)].sort_values("p")
            rows += [s_ for s_ in t.index if s_ not in rows][:TOP]
        if not rows:
            print(db, "no significant sets", flush=True)
            continue
        shared = {s_: is_sh(s_) for s_ in rows}
        recs = []
        for s in rows:
            for j, ct in enumerate(CTS):
                if s in tabs[ct].index:
                    r = tabs[ct].loc[s]
                    recs.append({"set": s, "x": j, "NES": r.NES, "padj": r.padj})
                    alltop.append({"database": db, "set": s, "cell_type": ct, "NES": r.NES, "padj": r.padj,
                                   "shared_3plus": shared[s]})
        F = pd.DataFrame(recs)
        F["y"] = F.set.map({s: i for i, s in enumerate(rows)})
        fig, ax = plt.subplots(figsize=(3.2, 0.3 * len(rows) + 1.2))
        sc_ = ax.scatter(F.x, F.y, c=F.NES.clip(-3, 3), cmap="RdBu_r", vmin=-3, vmax=3,
                         s=(-np.log10(F.padj.clip(1e-6)) * 30 + 6).clip(6, 200),
                         edgecolors=["black" if p <= 0.05 else "#BBBBBB" for p in F.padj],
                         linewidths=[1.0 if p <= 0.05 else 0.4 for p in F.padj])
        ax.set_xticks(range(len(CTS)), [SHORT[c] for c in CTS], fontsize=11)
        labs = [clean(s, db) + (" *" if shared[s] else "") for s in rows]
        ax.set_yticks(range(len(rows)), labs, fontsize=10.5)
        for tl, s in zip(ax.get_yticklabels(), rows):
            tl.set_color("#999999" if shared[s] else "#1F1F1F")
        ax.set_ylim(len(rows) - 0.5, -0.5); ax.set_xlim(-0.6, len(CTS) - 0.4)
        ax.grid(color="#EEEEEE", lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        cax = ax.inset_axes([1.07, 0.62, 0.07, 0.36])
        cb = fig.colorbar(sc_, cax=cax)
        cb.set_label("NES", fontsize=10.5)
        hs = [plt.scatter([], [], s=v * 30 + 6, c="#999999") for v in [1.3, 2, 4]]
        hs.append(plt.scatter([], [], s=45, c="white", edgecolors="black", linewidths=1.0))
        ax.legend(hs, ["0.05", "0.01", "1e-4", "sig."], title="padj", loc="upper left", bbox_to_anchor=(1.0, 0.52),
                  frameon=False, fontsize=9.5, title_fontsize=10, labelspacing=0.8, handletextpad=0.2)
        ax.set_title(f"GSEA: {dbl}", fontsize=12)
        fig.savefig(rd / "styled" / f"gsea_{db}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)
        print(db, "rows", len(rows), "shared", sum(shared.values()), flush=True)
    pd.DataFrame(alltop).to_csv(rd / "gsea_by_database_top.csv", index=False)


if __name__ == "__main__":
    main(*sys.argv[1:3])
