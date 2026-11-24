"""Heart-failure-focused figures: is there an HF signature in I4M hearts at 4 weeks, and which early changes exist?

  genotype/styled/hf_gene_scorecard.png      classic HF genes (CM, fibroblast, macrophage): log2FC I4M vs WT, coloured
                                             by whether the change goes the way it does in failing hearts
  genotype/styled/heatmap_hf_panel.png       per-mouse heatmap of the cardiomyocyte HF genes
  genotype/styled/early_remodel_mice.png     per-mouse CPM of the early-remodelling genes (Kcnd2, Cpt1b, Ttn, Tnni3,
                                             Myh6, Myh7, Fhl1) and the Myh7 / Myh6 ratio
  genotype/styled/heatmap_early_remodel.png  per-mouse heatmap of robust CM DEGs in ion-channel, sarcomere, fatty-acid
                                             and BCAA gene sets
  genotype/styled/hf_pathways_function.png   HF pathways (stress / cardiomyopathy, contraction, Ca2+, electrical) in
                                             5 databases x 2 methods (GSEA in 5 cell types, ORA on robust DEGs)
  genotype/styled/hf_pathways_metabolism.png same for energy metabolism
  genotype/styled/hf_pathways_stress_fibrosis.png  same for protein stress / hypoxia, fibrosis, inflammation
  genotype/styled/early_remodel_gsea.png     GSEA running-score curves in cardiomyocytes for the early-remodelling sets
  genotype/hf_pathways.csv                   the pathway values plotted (NES / padj / ORA overlap / padj)
GSEA / ORA values are read from the main pathway analysis (pw_GSEA_*, pw_ORA_*); nothing is re-tested except the
running-score curves, which are recomputed from the same ranking (sign(log2FC) x -log10 P, technical / ambient genes
removed) to show where the genes sit.
Usage: python 46_heart_failure.py <work_dir> <results_dir>
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

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from hf_sets import GENES, THEMES  # noqa: E402

spec = importlib.util.spec_from_file_location("sty", HERE / "15_genotype_figures_styled.py")
sty = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sty)
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})

DBL = {"GO_BP": "GO", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WP", "Hallmark": "Hallmark"}
GCT = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Macrophage"]
GLAB = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC", "Macrophage": "Mac"}
OCT = ["Cardiomyocyte", "Endocardial", "Fibroblast"]
WT, MUT = list("ABCD"), list("EFGH")
RED, BLUE, GREY = "#B2182B", "#2166AC", "#9E9E9E"


def nice(db, s):
    if db == "GO_BP":
        s = s.rsplit(" (GO:", 1)[0]
    elif db in ("Reactome", "WikiPathways", "Hallmark"):
        s = s.split("_", 1)[1].replace("_", " ").lower()
        for a, b in [("tca", "TCA"), ("ppar", "PPAR"), ("tgf beta", "TGF-beta"), ("tgfbeta", "TGF-beta"), ("upr", "UPR"),
                     ("mir208", "miR-208"), ("micrornas", "microRNAs"), ("tnfa", "TNF-a"), ("nfkb", "NF-kB"),
                     ("il6 jak stat3", "IL6-JAK-STAT3"), ("ecm", "ECM"), ("betaoxidation", "beta-oxidation"),
                     ("interleukin 1", "interleukin-1"), ("pparalpha", "PPAR-alpha")]:
            s = s.replace(a, b)
    s = s[0].upper() + s[1:]
    return s if len(s) <= 52 else s[:50] + "…"


def load_gmt(wd):
    g = {}
    for db in DBL:
        for line in open(wd / f"gmt_{db}.gmt", encoding="utf-8"):
            x = line.rstrip("\n").split("\t")
            g[(db, x[0])] = set(x[2:])
    return g


# ---------------------------------------------------------------- genes
def scorecard(wd, out):
    rows = []
    for grp, L in GENES.items():
        for g, exp, ct in L:
            d = pd.read_csv(wd / f"de_{ct}.csv.gz").set_index("gene")
            if g not in d.index:
                continue
            r = d.loc[g]
            rows.append(dict(group=grp, gene=g, exp=exp, ct=ct, lfc=r.logFC, p=r.PValue, fdr=r.FDR,
                             sig=(r.FDR <= 0.1) and abs(r.logFC) >= 0.25))
    T = pd.DataFrame(rows)
    groups = list(GENES)
    cols = [groups[:4], groups[4:]]
    fig, axes = plt.subplots(1, 2, figsize=(13, 9.2))
    for ax, gl in zip(axes, cols):
        y, ticks, labs = 0, [], []
        for grp in gl:
            t = T[T.group == grp]
            ax.text(-2.95, y - 0.2, grp, fontsize=10.5, weight="bold", va="bottom")
            y += 0.9
            for r in t.itertuples():
                same = np.sign(r.lfc) == r.exp
                col = (RED if same else BLUE) if r.p < 0.05 else ("#F4A6A6" if same else "#A6C4E4")
                ax.barh(y, r.lfc, color=col, height=0.7)
                mark = "**" if r.sig else ("*" if r.p < 0.05 else "")
                ax.text(r.lfc + (0.06 if r.lfc >= 0 else -0.06), y, mark, va="center",
                        ha="left" if r.lfc >= 0 else "right", fontsize=11)
                ax.text(2.95, y, "↑" if r.exp > 0 else "↓", va="center", ha="right", fontsize=11, color="#444444")
                ticks.append(y); labs.append(r.gene)
                y += 1
            y += 0.6
        ax.set_yticks(ticks, labs, fontsize=9.5, style="italic")
        ax.set_ylim(y - 0.4, -1.2)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_xlim(-3, 3)
        ax.set_xlabel("log2FC (I4M vs WT)")
        ax.text(2.95, -1.0, "in HF", ha="right", fontsize=9, color="#444444")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in [RED, "#F4A6A6", BLUE, "#A6C4E4"]]
    fig.legend(h, ["same direction as in HF, P < 0.05", "same direction, n.s.", "opposite to HF, P < 0.05",
                   "opposite, n.s."], loc="lower center", ncol=4, frameon=False, fontsize=10)
    fig.suptitle("Classic heart-failure genes in I4M vs WT (pseudobulk edgeR ~ sex + genotype)\n"
                 "** FDR <= 0.1 and |log2FC| >= 0.25;  * P < 0.05;  arrow = usual direction in failing hearts; "
                 "fibrosis genes read in fibroblasts, inflammation genes in macrophages", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out / "hf_gene_scorecard.png", dpi=220)
    plt.close(fig)
    return T


def per_mouse(wd, out, meta):
    pb = pd.read_csv(wd / "pb_Cardiomyocyte.csv.gz", index_col=0)[meta.index]
    de = pd.read_csv(wd / "de_Cardiomyocyte.csv.gz")
    cm_groups = {k: [g for g, _, ct in v if ct == "Cardiomyocyte"] for k, v in GENES.items()
                 if all(ct == "Cardiomyocyte" for _, _, ct in v)}
    sty.heatmap(pb, de, meta, cm_groups, out, "heatmap_hf_panel.png",
                "Cardiomyocytes: classic heart-failure genes per mouse")
    early = {"ion channels": ["Kcnd2", "Kcnv2", "Kcng2", "Hcn1", "Kcnj2", "Cacna1s", "Cacna2d1", "Slc4a3", "Rnf207",
                              "Irx5"],
             "sarcomere / Z-disc": ["Ttn", "Tnni3", "Tnnt2", "Myl4", "Myl1", "Mybpc2", "Myh7b", "Flnc", "Ankrd23",
                                    "Cryab", "Hspb8", "Dsp"],
             "fatty-acid uptake / oxidation": ["Cpt1b", "Acadvl", "Slc27a1", "Slc27a4", "Plin4", "Echdc2", "Acox1",
                                               "Ehhadh", "Etfa", "Acad11", "Slc22a5", "Plin2"],
             "BCAA catabolism": ["Bckdhb", "Dbt", "Mccc2", "Acadsb", "Abat"]}
    d = de.set_index("gene")
    for k, v in early.items():
        bad = [g for g in v if not d.loc[g, "robust"]]
        assert not bad, (k, bad)
    sty.heatmap(pb, de, meta, early, out, "heatmap_early_remodel.png",
                "Cardiomyocytes: robust DEGs in ion-channel, sarcomere, fatty-acid and BCAA gene sets")
    cpm = pb / pb.sum() * 1e6
    genes = ["Kcnd2", "Cpt1b", "Ttn", "Tnni3", "Myh6", "Myh7", "Fhl1"]
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.2))
    rng = np.random.default_rng(1)
    for ax, g in zip(axes.ravel(), genes + ["ratio"]):
        if g == "ratio":
            v = cpm.loc["Myh7"] / cpm.loc["Myh6"] * 100
            ttl, sub, ylab = "Myh7 / Myh6", "ratio x 100 (not tested)", "Myh7 / Myh6 (%)"
        else:
            v = cpm.loc[g]
            r = d.loc[g]
            ttl, ylab = g, "CPM"
            sub = f"log2FC {r.logFC:+.2f}, FDR {r.FDR:.2g}" + ("  (robust)" if r.robust else "")
        for k_, grp in enumerate([WT, MUT]):
            x = k_ + rng.uniform(-0.12, 0.12, 4)
            cols = ["#4D4D4D" if k_ == 0 else RED] * 4
            mk = ["o" if meta.loc[s, "sex"] == "M" else "^" for s in grp]
            for xi, s, c, m in zip(x, grp, cols, mk):
                ax.scatter(xi, v[s], color=c, marker=m, s=46, zorder=3)
            ax.hlines(v[grp].mean(), k_ - 0.25, k_ + 0.25, color="#222222", lw=2)
        ax.set_xticks([0, 1], ["WT", "I4M"])
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylim(0, v.max() * 1.25)
        ax.set_title(f"{ttl}\n{sub}", fontsize=10, style="italic" if g != "ratio" else "normal")
        ax.set_ylabel(ylab, fontsize=9)
    fig.legend([Line2D([], [], marker="o", ls="", color="#666666"), Line2D([], [], marker="^", ls="", color="#666666")],
               ["male", "female"], loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Early-remodelling genes in cardiomyocytes: each point one mouse (pseudobulk CPM), bar = group mean",
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(out / "early_remodel_mice.png", dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------- pathways
def pathway_table(wd):
    G = {(db, ct): pd.read_csv(wd / f"pw_GSEA_{db}_{ct}.csv.gz").set_index("set") for db in DBL for ct in GCT}
    O = {(db, ct, d): pd.read_csv(wd / f"pw_ORA_{db}_{ct}_{d}.csv.gz").set_index("set")
         for db in DBL for ct in OCT for d in ["Up", "Down"]}
    rows = []
    for th, L in THEMES.items():
        for db, s in L:
            r = dict(theme=th, db=db, set=s, label=f"{nice(db, s)} [{DBL[db]}]")
            ok = False
            for ct in GCT:
                g = G[(db, ct)]
                if s in g.index:
                    r[f"NES_{ct}"], r[f"padj_{ct}"], r[f"p_{ct}"] = g.loc[s, "NES"], g.loc[s, "padj"], g.loc[s, "p"]
                    ok = True
            for ct in OCT:
                for d in ["Up", "Down"]:
                    o = O[(db, ct, d)]
                    if s in o.index:
                        r[f"ORAn_{ct}_{d}"], r[f"ORApadj_{ct}_{d}"] = o.loc[s, "overlap"], o.loc[s, "padj"]
                        r[f"ORAp_{ct}_{d}"] = o.loc[s, "p"]
            if ok:
                rows.append(r)
    return pd.DataFrame(rows)


def pathway_plot(T, themes, out, name, title):
    T = T[T.theme.isin(themes)].reset_index(drop=True)
    n = len(T)
    H = 0.27 * n + 2.2
    fig = plt.figure(figsize=(11.5, H))
    gs = fig.add_gridspec(1, 2, width_ratios=[5, 6], wspace=0.05, left=0.4, right=0.85, top=1 - 1.05 / H,
                          bottom=0.6 / H)
    ax1, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    y = np.arange(n)
    # GSEA
    for j, ct in enumerate(GCT):
        nes, pa, p = T.get(f"NES_{ct}"), T.get(f"padj_{ct}"), T.get(f"p_{ct}")
        for i in range(n):
            if pd.isna(nes[i]):
                continue
            sz = 18 + 110 * min(-np.log10(p[i]), 4) / 4
            sc = ax1.scatter(j, i, s=sz, c=[nes[i]], cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                             edgecolors="black" if pa[i] <= 0.05 else "none", linewidths=1.4, zorder=3)
    ax1.set_xticks(range(len(GCT)), [GLAB[c] + "\n" for c in GCT])
    ax1.set_title("GSEA (all genes ranked)", fontsize=12)
    # ORA
    oc = [(ct, d) for ct in OCT for d in ["Up", "Down"]]
    for j, (ct, d) in enumerate(oc):
        nn, pa, p = T.get(f"ORAn_{ct}_{d}"), T.get(f"ORApadj_{ct}_{d}"), T.get(f"ORAp_{ct}_{d}")
        if nn is None:
            continue
        for i in range(n):
            if pd.isna(nn[i]):
                continue
            col = RED if d == "Up" else BLUE
            alpha = 1 if p[i] < 0.05 else 0.3
            ax2.scatter(j, i, s=12 + 9 * min(nn[i], 16), color=col, alpha=alpha,
                        edgecolors="black" if pa[i] <= 0.05 else "none", linewidths=1.4, zorder=3)
    ax2.set_xticks(range(len(oc)), [GLAB[c] + "\n" + ("↑" if d == "Up" else "↓") for c, d in oc])
    ax2.set_title("ORA (robust DEGs)", fontsize=12)
    for ax, k in [(ax1, len(GCT)), (ax2, len(oc))]:
        ax.set_ylim(n - 0.5, -0.5); ax.set_xlim(-0.6, k - 0.4)
        ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=2, labelsize=11)
        ax.grid(axis="y", color="#EEEEEE", zorder=0)
        for s_ in ax.spines.values():
            s_.set_visible(False)
    ax1.set_yticks(y, T.label, fontsize=10.5)
    ax2.set_yticks(y, [""] * n)
    ax1.tick_params(axis="y", length=0); ax2.tick_params(axis="y", length=0)
    # theme separators
    starts = {th: T.index[T.theme == th].min() for th in T.theme.unique()}
    for th, i0 in starts.items():
        if i0 > 0:
            for ax in (ax1, ax2):
                ax.axhline(i0 - 0.5, color="#888888", lw=0.8)
        i1 = T.index[T.theme == th].max()
        yf = fig.transFigure.inverted().transform(ax1.transData.transform((0, (i0 + i1) / 2)))[1]
        fig.text(0.012, yf, th.replace(" / ", " /\n") if i1 - i0 < 6 else th, rotation=90, ha="center", va="center",
                 fontsize=10.5, weight="bold", color="#444444", linespacing=1.0)
    cax = fig.add_axes([0.885, 0.55, 0.012, 0.25])
    cb = fig.colorbar(sc, cax=cax); cb.set_label("GSEA NES", fontsize=9)
    hs = [plt.scatter([], [], s=18 + 110 * v / 4, color="#777777") for v in [1, 2, 4]]
    fig.legend(hs, ["P 0.1", "0.01", "<= 1e-4"], title="GSEA P", loc="center left", bbox_to_anchor=(0.875, 0.38),
               frameon=False, fontsize=8, title_fontsize=8.5)
    hs = [plt.scatter([], [], s=12 + 9 * v, color="#777777") for v in [2, 8, 16]]
    fig.legend(hs, ["2", "8", ">= 16"], title="ORA genes", loc="center left", bbox_to_anchor=(0.875, 0.22),
               frameon=False, fontsize=8, title_fontsize=8.5)
    fig.text(0.02, 0.15 / H, "Black ring = padj <= 0.05.  ORA: red = up-regulated list, blue = down-regulated list; faint = "
             "P >= 0.05.  Blank = set not tested (size filter) or no overlap.", fontsize=8.5, color="#444444")
    fig.suptitle(title, fontsize=13, x=0.5, y=1 - 0.2 / H)
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def gsea_curves(wd, out, gmt):
    de = pd.read_csv(wd / "de_Cardiomyocyte.csv.gz")
    de = de[~(de.technical | de.ambient)]
    rk = (np.sign(de.logFC) * -np.log10(de.PValue.clip(1e-300))).values
    genes = de.gene.values
    o = np.argsort(-rk); rk, genes = rk[o], genes[o]
    G = {db: pd.read_csv(wd / f"pw_GSEA_{db}_Cardiomyocyte.csv.gz").set_index("set") for db in DBL}
    sets = [("Reactome", "REACTOME_POTASSIUM_CHANNELS", ["Kcnd2", "Kcnv2", "Kcng2", "Kcnj2", "Hcn1"]),
            ("KEGG", "Fatty acid degradation", ["Cpt1b", "Acadvl", "Ehhadh", "Acox1"]),
            ("Reactome", "REACTOME_STRIATED_MUSCLE_CONTRACTION", ["Ttn", "Tnni3", "Tnnt2", "Myl4", "Myh6", "Myh7"]),
            ("KEGG", "Hypertrophic cardiomyopathy", ["Ttn", "Tnni3", "Myh6", "Myh7", "Cacna1s"]),
            ("KEGG", "Valine, leucine and isoleucine degradation", ["Bckdhb", "Dbt", "Mccc2", "Acadsb", "Abat"]),
            ("Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", ["Hspa5", "Calr", "Hsp90b1", "Hyou1"])]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.2))
    N = len(genes)
    for ax, (db, s, lab) in zip(axes.ravel(), sets):
        S = gmt[(db, s)]
        hit = np.isin(genes, list(S))
        w = np.abs(rk) * hit
        ph = np.cumsum(w) / w.sum()
        pm = np.cumsum(~hit) / (~hit).sum()
        es = ph - pm
        col = RED if G[db].loc[s, "NES"] > 0 else BLUE
        ax.plot(np.arange(N), es, color=col, lw=1.6)
        ax.axhline(0, color="#999999", lw=0.6)
        lo = es.min() - 0.12 * (es.max() - es.min() + 1e-9)
        hi = es.max()
        ax.vlines(np.where(hit)[0], lo - 0.08, lo, color="#333333", lw=0.5)
        ax.set_ylim(lo - 0.25, hi + 0.12)
        pos = sorted((int(np.where(genes == g)[0][0]), g) for g in lab if g in S and g in genes)
        gap, xs = 0.13 * N, [p for p, _ in pos]
        for k_ in range(1, len(xs)):            # spread labels left to right with a minimum gap ...
            xs[k_] = max(xs[k_], xs[k_ - 1] + gap)
        lim = N * 0.97
        for k_ in range(len(xs) - 1, -1, -1):   # ... then pull them back inside the axis
            xs[k_] = min(xs[k_], lim - (len(xs) - 1 - k_) * gap)
        xs = [max(x, N * 0.03 + k_ * gap) for k_, x in enumerate(xs)]
        for (i, g), xt in zip(pos, xs):
            ax.annotate(g, xy=(i, lo - 0.08), xytext=(xt, lo - 0.2), fontsize=7.5, style="italic",
                        ha="center", va="top", color="#B2182B" if rk[i] > 0 else "#2166AC",
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
        g_ = G[db].loc[s]
        ax.set_title(f"{nice(db, s)} [{DBL[db]}]\nNES {g_.NES:+.2f}, P {g_.p:.2g}, padj {g_.padj:.2g}", fontsize=10)
        ax.set_xlim(0, N); ax.set_xticks([0, N], ["up in I4M", "down in I4M"], fontsize=8)
        ax.set_ylabel("running enrichment score", fontsize=8.5)
    fig.suptitle("Cardiomyocytes: GSEA running-score curves of early-remodelling and HF-related pathways\n"
                 "genes ranked from most up (left) to most down (right) in I4M; ticks = genes of the pathway; "
                 "red curve = pathway up, blue = down", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out / "early_remodel_gsea.png", dpi=220)
    plt.close(fig)


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    meta = pd.read_csv(HERE / "sample_sheet.csv").set_index("sample")
    S = scorecard(wd, out)
    print(S.round(3).to_string(index=False))
    per_mouse(wd, out, meta)
    T = pathway_table(wd)
    T.to_csv(Path(res_dir) / "genotype" / "hf_pathways.csv", index=False)
    th = list(THEMES)
    pathway_plot(T, th[:4], out, "hf_pathways_function.png",
                 "Heart-failure pathways (1): stress / cardiomyopathy, contraction, Ca2+ handling, electrical")
    pathway_plot(T, th[4:5], out, "hf_pathways_metabolism.png", "Heart-failure pathways (2): energy metabolism")
    pathway_plot(T, th[5:], out, "hf_pathways_stress_fibrosis.png",
                 "Heart-failure pathways (3): protein stress / hypoxia, fibrosis, inflammation")
    gsea_curves(wd, out, load_gmt(wd))
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
