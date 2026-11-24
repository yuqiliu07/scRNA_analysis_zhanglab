"""Name the subclusters from 31_subclusters.py, test their composition WT vs I4M and draw the figures.

Names were assigned by reviewing the marker tables (markers_<type>.csv.gz) and per-cluster quality metrics
(median genes, ambient / low-quality signature); clusters whose markers are mainly ambient / low-quality genes or
that are driven by library size are labelled "low quality" and kept in the figures but not interpreted.
Every name is checked against marker genes that must be among the cluster's top 40 markers (stops if a re-run
renumbers the clusters).

Composition: fraction of each subcluster per mouse (of that cell type); arcsine-square-root OLS ~ sex + genotype;
BH over the subclusters of a cell type.

Outputs: results_dir/genotype/subcluster_composition_test.csv, subcluster_markers_top.csv, subcluster_qc.csv
         results_dir/genotype/styled/sub_<type>_umap.png, sub_<type>_dotplot.png, sub_<type>_composition.png
Usage: python 32_subcluster_figures.py <work_dir> <results_dir> [cell_type ...]  (tables written only for a full run)
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
HUE = ["#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7", "#8C564B", "#7F7F7F", "#17BECF",
       "#BCBD22", "#E377C2"]
TYPES = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
# cluster id -> (name, marker genes that must be among its top-40 markers); filled after review
NAMES = {t: {} for t in TYPES}
NAMES["Cardiomyocyte"] = {
    "0": ("vCM Ankrd1+ Xirp2+", ["Ankrd1", "Xirp2", "Nrap"]),
    "1": ("vCM Grm1+ Prune2+", ["Grm1", "Prune2", "Mylk4"]),
    "2": ("vCM Kcnj3+ Corin+", ["Kcnj3", "Corin", "Hcn1"]),
    "3": ("Conduction-like Etv1+", ["Etv1", "Slit3", "Kirrel3"]),
    "4": ("vCM Scn5a+ Fgf14+", ["Scn5a", "Fgf14", "Chrm2"]),
    "5": ("vCM Mb+ OXPHOS-high", ["Mb", "Cox6a2"]),
    "6": ("Minor (19 nuclei)", []),
}
NAMES["Endocardial"] = {
    "0": ("Endo Nrg1+", ["Nrg1", "Chrm3", "Bmp6"]),
    "1": ("Endo Stox2+ Sulf1+", ["Stox2", "Klf12", "Sulf1"]),
    "2": ("Endo Ptprb+ (EC-like)", ["Ptprb", "Actn1", "Calcrl"]),
    "3": ("Glia-like (doublet?)", ["Csmd1", "Cdh19"]),
    "4": ("Proliferating", ["Mki67", "Top2a"]),
    "5": ("Valve endocardium", ["Lsamp", "Postn", "Prox1"]),
}
NAMES["Fibroblast"] = {
    "0": ("Fib Pi16+", ["Pi16", "Gfpt2", "Fbn1"]),
    "1": ("Fib Abca8a+ (main)", ["Abca8a", "Rora", "Gpc6"]),
    "2": ("Fib Cfh+ Tbx20+", ["Cfh", "Tbx20", "Fbln1"]),
    "3": ("Fib Slit3+ Aldh1a1+", ["Slit3", "Svep1", "Aldh1a1"]),
    "4": ("Low quality (CM RNA)", ["Mb", "Ttn", "mt-Co1"]),
    "5": ("Glia-like (doublet?)", ["Csmd1", "Cadm2"]),
}
NAMES["Endothelial"] = {
    "0": ("Arterial EC", ["Smad6", "Prdm16", "Efnb2"]),
    "1": ("IFN-high, low RNA", ["Iigp1", "Gbp4"]),
    "2": ("Capillary EC Cd36+", ["Cd36", "Xdh", "Ccdc85a"]),
    "3": ("Low quality (CM RNA)", ["Mb", "Hbb-bs", "Ttn"]),
    "4": ("Venous EC", ["Vwf", "Vegfc", "Calcrl"]),
    "5": ("Capillary EC Col4a1+ Kit+", ["Kit", "Col4a1", "Nrp2"]),
    "6": ("Proliferating", ["Mki67", "Top2a"]),
}


def names_for(ct, m):
    top = m.groupby("group").head(40).groupby("group")["names"].apply(set)
    out = {}
    for c in sorted(top.index, key=int):
        if c in NAMES[ct]:
            name, req = NAMES[ct][c]
            miss = [g for g in req if g not in top[c]]
            if miss:
                sys.exit(f"{ct} cluster {c} ({name}): expected markers missing {miss}; check the clustering")
            out[c] = name
        else:
            out[c] = f"{ct[:3]}_{c}"
    return out


def comp_test(o, meta, ct):
    p = pd.crosstab(o["sample"], o["name"], normalize="index")
    rows = []
    for n in p.columns:
        d = meta.join(p[n].rename("f"))
        d["y"] = np.arcsin(np.sqrt(d.f))
        f = smf.ols("y ~ sex + genotype", d).fit()
        rows.append({"cell_type": ct, "subcluster": n, "mean_WT": d.f[d.genotype == "WT"].mean(),
                     "mean_I4M": d.f[d.genotype == "I4M"].mean(), "coef": f.params["genotype[T.WT]"] * -1,
                     "p": f.pvalues["genotype[T.WT]"]})
    r = pd.DataFrame(rows)
    r["fdr"] = multipletests(r.p, method="fdr_bh")[1]
    return r, p


def main(work_dir, res_dir, *types):
    wd, rd = Path(work_dir) / "subclusters", Path(res_dir) / "genotype"
    out = rd / "styled"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")[["genotype", "sex"]]
    a = ad.read_h5ad(Path(work_dir) / "merged" / "final.h5ad", backed="r")
    tests, tops, qcs = [], [], []
    for ct in (list(types) or TYPES):
        o = pd.read_csv(wd / f"sub_{ct}.csv.gz", index_col=0, dtype={"sub": str})
        m = pd.read_csv(wd / f"markers_{ct}.csv.gz", dtype={"group": str})
        nm = names_for(ct, m)
        o["name"] = o["sub"].map(nm)
        order = [nm[c] for c in sorted(nm, key=int)]
        col = {n: HUE[i % len(HUE)] for i, n in enumerate(order)}
        o["genotype"] = meta.loc[o["sample"], "genotype"].values
        r, p = comp_test(o, meta, ct)
        tests.append(r)
        q = o.groupby("name").agg(n=("sample", "size"), median_genes=("n_genes", "median"),
                                  median_ambient_signature=("tech_score", "median"), median_pct_mt=("pct_mt", "median"))
        q.insert(0, "cell_type", ct)
        qcs.append(q.reindex(order).reset_index())

        # ---- UMAP: all / WT / I4M ----
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
        rng = np.random.default_rng(0)
        idx = rng.permutation(len(o))
        oo = o.iloc[idx]
        for ax, sel, ttl in [(axes[0], None, f"{ct}: {len(o):,} nuclei"), (axes[1], "WT", "WT"), (axes[2], "I4M", "I4M")]:
            if sel:
                ax.scatter(oo.umap1, oo.umap2, s=0.5, c="#E5E5E5", rasterized=True, lw=0)
                d = oo[oo.genotype == sel]
            else:
                d = oo
            ax.scatter(d.umap1, d.umap2, s=0.8 if len(o) > 20000 else 2, c=d.name.map(col), rasterized=True, lw=0)
            ax.set_title(ttl, fontsize=13, color=GCOL.get(sel, "black"))
            ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
        for n in order:
            if (o.name == n).sum() < 100:  # too small / scattered to label
                continue
            med = o[o.name == n][["umap1", "umap2"]].median()
            axes[0].text(med.umap1, med.umap2, n.split(" (")[0], fontsize=8, weight="bold",
                         ha="center", va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        hs = [plt.scatter([], [], s=40, c=col[n]) for n in order]
        fig.legend(hs, order, loc="center right", frameon=False, fontsize=10, bbox_to_anchor=(1.12, 0.5))
        fig.tight_layout()
        fig.savefig(out / f"sub_{ct}_umap.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        # ---- dot plot of top markers ----
        mk = []
        for c in sorted(nm, key=int):
            g = m[(m.group == c) & (m.pvals_adj < 0.05)].sort_values("scores", ascending=False)
            g = g[~g.names.str.match(r"^(mt-|Rp[ls]\d|Gm\d|ENSMUSG)") & ~g.names.str.endswith("Rik")]
            req = NAMES[ct].get(c, ("", []))[1]
            for x in (list(req) + list(g.names))[:4] if req else g.names.head(4):
                if x not in mk:
                    mk.append(x)
            tops.append({"cell_type": ct, "subcluster": nm[c], "top_markers": ",".join(g.names.head(10))})
        X = a[o.index.values, mk].to_memory()
        Xc = X.layers["counts"]
        lib = np.asarray(o.n_counts.values, dtype=float)
        L = np.log1p(Xc.multiply(1e4 / lib[:, None]).toarray())
        E = pd.DataFrame(L, columns=mk); E["name"] = o.name.values
        D = pd.DataFrame((Xc.toarray() > 0), columns=mk); D["name"] = o.name.values
        avg = E.groupby("name").mean().reindex(order); pct = D.groupby("name").mean().reindex(order)
        z = (avg - avg.mean()) / avg.std().replace(0, 1)
        fig, ax = plt.subplots(figsize=(0.32 * len(mk) + 3, 0.42 * len(order) + 1.6))
        yy, xx = np.meshgrid(range(len(order)), range(len(mk)), indexing="ij")
        sc_ = ax.scatter(xx.ravel(), yy.ravel(), s=(pct.values.ravel() * 140), c=z.values.ravel().clip(-2, 2.5),
                         cmap="RdBu_r", vmin=-2, vmax=2.5, edgecolors="#555555", lw=0.3)
        ax.set_xticks(range(len(mk)), mk, rotation=90, style="italic", fontsize=9)
        ax.set_yticks(range(len(order)), order, fontsize=10)
        ax.set_ylim(len(order) - 0.5, -0.5); ax.set_xlim(-0.6, len(mk) - 0.4)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.025, pad=0.01); cb.set_label("Average expression (scaled)")
        hs = [plt.scatter([], [], s=v * 140, c="#888888") for v in [0.25, 0.5, 0.75]]
        ax.legend(hs, ["25%", "50%", "75%"], title="Percent expressed", loc="upper left", bbox_to_anchor=(1.13, 0.45),
                  frameon=False, fontsize=8, title_fontsize=8)
        fig.savefig(out / f"sub_{ct}_dotplot.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        # ---- composition per mouse ----
        fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1, 2.2]})
        ax = axes[0]
        samples = list(meta.index)
        bottom = np.zeros(len(samples))
        for n in order:
            v = p.reindex(samples)[n].fillna(0).values * 100
            ax.bar(range(len(samples)), v, bottom=bottom, color=col[n], width=0.8, label=n)
            bottom += v
        ax.set_xticks(range(len(samples)), [f"{s}\n{meta.loc[s, 'sex']}" for s in samples])
        for i, s in enumerate(samples):
            ax.get_xticklabels()[i].set_color(GCOL[meta.loc[s, "genotype"]])
        ax.set_ylabel(f"% of {ct.lower()} nuclei"); ax.set_ylim(0, 100)
        ax.set_title("Per mouse (WT A-D, I4M E-H)", fontsize=11)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=8)
        ax = axes[1]
        for i, n in enumerate(order):
            for j, g in enumerate(["WT", "I4M"]):
                v = p.loc[meta.index[meta.genotype == g], n] * 100
                x0 = i + (j - 0.5) * 0.38
                ax.boxplot(v, positions=[x0], widths=0.3, patch_artist=True, showfliers=False,
                           boxprops=dict(facecolor=GCOL[g], alpha=0.35, color=GCOL[g]), medianprops=dict(color=GCOL[g]),
                           whiskerprops=dict(color=GCOL[g]), capprops=dict(color=GCOL[g]))
                for sx, mk_ in [("M", "o"), ("F", "^")]:
                    vv = v[meta.loc[v.index, "sex"] == sx]
                    ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(vv)), vv, s=24, c=GCOL[g], marker=mk_,
                               edgecolors="black", lw=0.4, zorder=3)
            rr = r[r.subcluster == n].iloc[0]
            ymax = p[n].max() * 100
            ax.text(i, ymax + p.values.max() * 100 * 0.04, f"FDR {rr.fdr:.2f}" if rr.fdr >= 0.01 else "FDR<0.01", ha="center", fontsize=8,
                    color="#C0392B" if rr.fdr <= 0.05 else "#555555")
        ax.set_xticks(range(len(order)), order, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("% of nuclei")
        ax.set_ylim(top=p.values.max() * 100 * 1.13)
        ax.set_title("WT vs I4M (asin-sqrt OLS ~ sex + genotype, BH)", fontsize=11)
        hs = [plt.Rectangle((0, 0), 1, 1, color=GCOL[g], alpha=0.6) for g in GCOL] + \
             [plt.scatter([], [], s=24, c="white", marker=k, edgecolors="black") for k in "o^"]
        ax.legend(hs, list(GCOL) + ["male", "female"], frameon=False, loc="upper right", ncol=2, fontsize=9)
        fig.tight_layout()
        fig.savefig(out / f"sub_{ct}_composition.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(ct, r.round(4).to_string(), flush=True)

    if not types:
        pd.concat(tests).to_csv(rd / "subcluster_composition_test.csv", index=False)
        pd.DataFrame(tops).to_csv(rd / "subcluster_markers_top.csv", index=False)
        pd.concat(qcs).to_csv(rd / "subcluster_qc.csv", index=False)


if __name__ == "__main__":
    main(*sys.argv[1:])
