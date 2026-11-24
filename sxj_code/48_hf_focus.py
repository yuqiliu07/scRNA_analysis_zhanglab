"""Focused figures for the HF section: which processes change early in I4M cardiomyocytes, which genes, which pathways.

  genotype/styled/early_remodel_genes.png     robust CM DEGs of the seven early-change areas, log2FC bars
  genotype/styled/early_remodel_pathways.png  CM pathways of those areas in the 5 databases: GSEA NES (bar, coloured
                                              when padj <= 0.05), ORA on robust up / down DEGs, and the robust genes
                                              that drive each set (GSEA leading edge, or robust DEGs in the set when
                                              the set is not significant)
Areas: significant pathways from the main CM GSEA / ORA (padj <= 0.05) that relate to HF biology, plus the electrical,
sarcomere and fatty-acid sets whose key genes change although the set as a whole does not (shown greyed).
All statistics are read from the main analysis (de_Cardiomyocyte, pw_GSEA_*, pw_ORA_*); nothing is re-tested.
Usage: python 48_hf_focus.py <work_dir> <results_dir>
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import blended_transform_factory as blended

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
RED, BLUE, LRED, LBLUE, GREY = "#B2182B", "#2166AC", "#F4A6A6", "#A6C4E4", "#BDBDBD"
CT = "Cardiomyocyte"

# area -> robust CM DEGs shown in the gene figure
GENES = {
    "Electrical: K+ channels": ["Kcnd2", "Kcnv2", "Kcng2", "Kcnj2", "Hcn1"],
    "Sarcomere / cytoskeleton": ["Ttn", "Tnni3", "Tnnt2", "Myl4", "Myl1", "Mybpc2", "Flnc", "Cryab"],
    "Fatty-acid uptake / oxidation": ["Cpt1b", "Acadvl", "Slc27a1", "Acox1", "Ehhadh"],
    "BCAA catabolism": ["Bckdhb", "Dbt", "Mccc2", "Acadsb", "Abat"],
    "Protein folding / ER stress": ["Hspa5", "Hsp90aa1", "Hsp90b1", "Hspa8", "Calr", "Pdia3", "Hyou1", "Hspd1"],
    "Antioxidant / detoxification": ["Nqo1", "Mgst1", "Nnt", "Sod1", "Gstk1"],
    "Cardiac growth signals": ["Mef2c", "Fgf1", "Fgf9", "Adra1b"],
}
# area -> [(db, set)]
SETS = {
    "Electrical: K+ channels": [("Reactome", "REACTOME_POTASSIUM_CHANNELS")],
    "Sarcomere / cytoskeleton": [("GO_BP", "supramolecular fiber organization (GO:0097435)"),
                                 ("Reactome", "REACTOME_STRIATED_MUSCLE_CONTRACTION"),
                                 ("WikiPathways", "WP_STRIATED_MUSCLE_CONTRACTION")],
    "Fatty-acid uptake / oxidation": [("KEGG", "Fatty acid degradation"),
                                      ("WikiPathways", "WP_FATTY_ACID_BETAOXIDATION")],
    "BCAA catabolism": [("GO_BP", "branched-chain amino acid catabolic process (GO:0009083)"),
                        ("KEGG", "Valine, leucine and isoleucine degradation"),
                        ("Reactome", "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM")],
    "Protein folding / ER stress": [("GO_BP", "protein folding (GO:0006457)"),
                                    ("GO_BP", "response to endoplasmic reticulum stress (GO:0034976)"),
                                    ("KEGG", "Protein processing in endoplasmic reticulum"),
                                    ("KEGG", "Integrated stress response (ISR) signaling pathway"),
                                    ("Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")],
    "Antioxidant / detoxification": [("GO_BP", "cellular oxidant detoxification (GO:0098869)"),
                                     ("WikiPathways", "WP_OXIDATIVE_STRESS_RESPONSE")],
    "Cardiac growth signals": [("GO_BP", "cardiac muscle tissue growth (GO:0055017)"),
                               ("GO_BP", "cellular response to fibroblast growth factor stimulus (GO:0044344)")],
}
SHORT = {"Electrical: K+ channels": "Electrical", "Sarcomere / cytoskeleton": "Sarcomere /\ncytoskeleton",
         "Fatty-acid uptake / oxidation": "Fatty acids", "BCAA catabolism": "BCAA\ncatabolism",
         "Protein folding / ER stress": "Protein folding /\nER stress", "Antioxidant / detoxification": "Antioxidant",
         "Cardiac growth signals": "Cardiac growth\nsignals"}
# endocardium / fibroblasts: robust DEGs that belong to the HF gene sets (hf_sets THEMES / classic genes, see
# 47_hf_deg_classes.py) or to HF-relevant endothelial / matrix functions, grouped by function
KEY = {
    "Endocardial": {
        "Signals to the myocardium": ["Nrg1", "Pdgfc", "Bmp6", "Gdf6", "Nrp1", "Il6st"],
        "Flow / nitric oxide": ["Nos3", "Klf4", "Piezo1", "Hes1", "Kdr", "Bcl6b"],
        "Ca2+ / hypertrophy signalling": ["Nfatc1", "Itpr3", "Plcd1", "Plcb1", "Gata6"],
        "Matrix": ["Fn1", "Vim", "Reck", "Dpp4", "Cyp1b1"],
        "Stress": ["Hsp90b1", "Mb"],
    },
    "Fibroblast": {
        "Collagen / matrix down": ["Col6a1", "Col5a1", "Eln", "Adamts14", "Adamts15", "Sh3pxd2b"],
        "Matrix / anti-fibrotic up": ["Dcn", "Timp3", "Fbln1", "Pcolce2", "Svep1", "Ccn5"],
        "Growth-factor signalling": ["Igf1r", "Akt3", "Irs1", "Vegfc", "Runx2"],
        "Protein stress": ["Hyou1", "Pdia6", "Hspa8", "Prdx4", "Selenbp1"],
        "Mitochondria / metabolism": ["Ndufa3", "Taldo1"],
    },
}
DBL = {"GO_BP": "GO", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WP", "Hallmark": "Hallmark"}


def nice(db, s):
    if db == "GO_BP":
        s = s.rsplit(" (GO:", 1)[0]
    elif db in ("Reactome", "Hallmark", "WikiPathways"):
        s = s.split("_", 1)[1].replace("_", " ").lower()
    s = s.replace("bcaa", "BCAA").replace("cellular response to fibroblast growth factor stimulus",
                                          "response to FGF stimulus")
    s = s.replace("Integrated stress response (ISR) signaling pathway", "Integrated stress response")
    return (s[0].upper() + s[1:]) + f"  [{DBL[db]}]"


def gene_fig(de, out, groups=None, name="early_remodel_genes.png",
             title="Cardiomyocytes: robust DEGs of the early-change areas"):
    groups = groups or GENES
    rows = []
    for area, gl in groups.items():
        for g in gl:
            r = de.loc[g]
            assert r.robust, g
            rows.append((area, g, r.logFC, r.FDR))
    T = pd.DataFrame(rows, columns=["area", "gene", "lfc", "fdr"])
    areas = list(groups)
    half = (len(T) + 1) // 2
    cum, cut = 0, len(areas)
    for i_, a_ in enumerate(areas):
        cum += (T.area == a_).sum()
        if cum >= half:
            cut = i_ + 1
            break
    cols = [areas[:cut], areas[cut:]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.6))
    ymax = 0
    lim = max(2.0, T.lfc.abs().max() * 1.3)
    for ax, al in zip(axes, cols):
        y, ticks, labs = 0, [], []
        for a in al:
            t = T[T.area == a]
            ax.text(-lim * 0.97, y - 0.15, a, fontsize=11.5, weight="bold", va="bottom")
            y += 0.85
            for r in t.itertuples():
                ax.barh(y, r.lfc, color=RED if r.lfc > 0 else BLUE, height=0.7)
                ax.text(r.lfc + (0.04 if r.lfc > 0 else -0.04), y, f"{r.lfc:+.2f}", va="center",
                        ha="left" if r.lfc > 0 else "right", fontsize=9, color="#444444")
                ticks.append(y); labs.append(r.gene)
                y += 1
            y += 1.0
        ax.set_yticks(ticks, labs, fontsize=10.5, style="italic")
        ymax = max(ymax, y)
        ax.set_xlim(-lim, lim)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_xlabel("log2FC (I4M vs WT)", fontsize=11)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    for ax in axes:
        ax.set_ylim(ymax - 0.8, -1.0)
    fig.suptitle(f"{title} (all FDR <= 0.1, |log2FC| >= 0.25, not flagged)\n"
                 "blue = lower in I4M, red = higher", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def path_fig(wd, de, out):
    rb = set(de.index[de.robust])
    rows = []
    for area, L in SETS.items():
        for db, s in L:
            g = pd.read_csv(wd / f"pw_GSEA_{db}_{CT}.csv.gz").set_index("set")
            gmt = {}
            for line in open(wd / f"gmt_{db}.gmt", encoding="utf-8"):
                x = line.rstrip("\n").split("\t")
                if x[0] == s:
                    gmt = set(x[2:])
            r = g.loc[s]
            ora = {}
            for d in ["Up", "Down"]:
                o = pd.read_csv(wd / f"pw_ORA_{db}_{CT}_{d}.csv.gz").set_index("set")
                ora[d] = (o.loc[s, "padj"], o.loc[s, "overlap"]) if s in o.index else (np.nan, 0)
            if r.padj <= 0.05:
                le = [x.strip() for x in str(r.leading_edge).split(",")]
                genes = [x for x in le if x in rb]
            else:
                sub = de[de.index.isin(gmt) & de.robust].sort_values("FDR")
                genes = list(sub.index)
            rows.append(dict(area=area, label=nice(db, s), NES=r.NES, padj=r.padj, p=r.p, genes=genes,
                             ora_up=ora["Up"][0], ora_dn=ora["Down"][0], n_up=ora["Up"][1], n_dn=ora["Down"][1]))
    T = pd.DataFrame(rows)
    n = len(T)
    W, H = 12.3, 0.26 * n + 1.25
    fig = plt.figure(figsize=(W, H))
    bot, top = 0.55 / H, 1 - 0.7 / H
    ax = fig.add_axes([4.75 / W, bot, 2.0 / W, top - bot])
    ax2 = fig.add_axes([6.85 / W, bot, 0.75 / W, top - bot])
    y, ys, labs = 0, [], []
    for k, area in enumerate(SETS):
        t = T[T.area == area]
        y0 = y
        for r in t.itertuples():
            sig = r.padj <= 0.05
            ax.barh(y, r.NES, color=(RED if r.NES > 0 else BLUE) if sig else GREY, height=0.66)
            lab = (f"{r.padj:.1e}" if r.padj < 0.001 else f"{r.padj:.3f}") if sig else "n.s."
            ax.text(r.NES + (0.07 if r.NES > 0 else -0.07), y, lab, va="center",
                    ha="left" if r.NES > 0 else "right", fontsize=8.5, color="#333333")
            for j, (pa, col) in enumerate([(r.ora_up, RED), (r.ora_dn, BLUE)]):
                if pd.notna(pa):
                    ax2.scatter(j, y, s=70 if pa <= 0.05 else 28, color=col if pa <= 0.05 else "white",
                                edgecolor=col, lw=1.1)
            gl = r.genes[:6]
            txt = ", ".join(gl) + (" …" if len(r.genes) > 6 else "")
            ax2.text(1.6, y, txt, va="center", ha="left", fontsize=10, style="italic",
                     color=("#222222" if sig else "#8A8A8A"), clip_on=False)
            ys.append(y); labs.append(r.label)
            y += 1
        if k % 2 == 0:
            ax.axhspan(y0 - 0.5, y - 0.5, xmin=-2.3, xmax=1.0, color="#F3F3F3", zorder=0, clip_on=False)
        ax.text(0.012, (y0 + y - 1) / 2, SHORT[area], fontsize=10.5, weight="bold", va="center", ha="left",
                color="#333333", clip_on=False, transform=blended(fig.transFigure, ax.transData))
        y += 0.0
    for a_ in (ax, ax2):
        a_.set_ylim(y - 0.5, -0.5)
    ax.set_yticks(ys, labs, fontsize=10)
    ax.set_xlim(-3.4, 3.7)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_xlabel("GSEA NES (label = padj)", fontsize=9.5)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax2.set_xlim(-0.6, 1.5)
    ax2.set_xticks([0, 1], ["up", "down"], fontsize=9.5)
    ax2.set_xlabel("ORA", fontsize=9.5)
    ax2.set_yticks([])
    for sp in ["left", "bottom"]:
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(axis="x", length=0)
    ax2.text(1.6, -1.0, "robust DEGs driving the set (grey = set n.s.)", fontsize=10, weight="bold",
             va="center", clip_on=False)
    fig.text(0.01, 1 - 0.28 / H, "Cardiomyocytes: pathways of the early-change areas.  Bar coloured = GSEA padj <= 0.05 "
             "(red up, blue down), grey = n.s.;  ORA dot filled = padj <= 0.05, open = n.s.", fontsize=10.5)
    fig.savefig(out / "early_remodel_pathways.png", dpi=220)
    plt.close(fig)
    return T


def top_fig(de, out, n=15, groups=None, name="hf_top_cm.png"):
    """Top robust DEGs (ranked by P), key genes marked: the unbiased starting point of the HF slides."""
    key = {g: a for a, gl in (groups or GENES).items() for g in gl}
    rb = de[de.robust].sort_values(["PValue", "FDR"])
    up, dn = rb[rb.logFC > 0].head(n), rb[rb.logFC < 0].head(n)
    fig, axes = plt.subplots(2, 1, figsize=(4.4, 7.4))
    for ax, d, col, lab in [(axes[0], up, RED, f"Top {n} higher in I4M"), (axes[1], dn, BLUE, f"Top {n} lower in I4M")]:
        y = np.arange(len(d))
        ax.barh(y, d.logFC, color=col, height=0.7)
        labs = [f"* {g}" if g in key else g for g in d.index]
        ax.set_yticks(y, labs, fontsize=10, style="italic")
        for t, g in zip(ax.get_yticklabels(), d.index):
            if g in key:
                t.set_weight("bold"); t.set_color("#7A4A00")
        ax.set_ylim(len(d) - 0.5, -0.5)
        lim = max(d.logFC.abs().max(), 0.5) * 1.15
        ax.set_xlim((0, lim) if col == RED else (-lim, 0))
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_title(lab, fontsize=11.5, weight="bold", color=col, loc="left")
        ax.tick_params(axis="x", labelsize=9); ax.tick_params(axis="y", length=0)
    axes[1].set_xlabel("log2FC (I4M vs WT)", fontsize=10)
    fig.text(0.02, 0.005, "ranked by P; * bold = key early-change gene (next panel)", fontsize=9, color="#555555")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out / name, dpi=220)
    plt.close(fig)
    print("top up:", list(up.index)); print("top down:", list(dn.index))


def short(db, s, n=48):
    s = nice(db, s).rsplit("  [", 1)[0]
    s = re.sub(r"(?i)\bmrna\b", "mRNA", s)
    s = re.sub(r"(?i)\brna\b", "RNA", s)
    for a, b in [("Esc ", "ESC "), ("fgfr1", "FGFR1"), ("Uv ", "UV "), ("Kras", "KRAS"), ("Myc ", "MYC "), ("pre mRNA", "pre-mRNA")]:
        s = s.replace(a, b)
    return s if len(s) <= n else s[:n - 1] + "…"


def top_pw_fig(wd, out, method, k=3):
    """Unbiased top CM pathways in each of the 5 databases for one method (GSEA: top k up + k down by P;
    ORA: top k sets in robust up DEGs + top k in robust down DEGs, by P). Grey = padj > 0.05."""
    dbs = list(DBL)
    cols = [dbs[:3], dbs[3:]]
    fig = plt.figure(figsize=(12.3, 4.7))
    for c, dl in enumerate(cols):
        rows = []
        for db in dl:
            rows.append(("head", db, None, None))
            if method == "GSEA":
                g = pd.read_csv(wd / f"pw_GSEA_{db}_{CT}.csv.gz")
                for d in (g[g.NES > 0].sort_values("p").head(k), g[g.NES < 0].sort_values("p").head(k)):
                    rows += [("row", short(db, r.set), r.NES, r.padj) for r in d.itertuples()]
            else:
                for dr, sign in (("Up", 1), ("Down", -1)):
                    o = pd.read_csv(wd / f"pw_ORA_{db}_{CT}_{dr}.csv.gz").sort_values("p").head(k)
                    rows += [("row", short(db, r.set), sign * -np.log10(r.p), r.padj) for r in o.itertuples()]
        ax = fig.add_axes([0.255 + c * 0.5, 0.11, 0.2, 0.83])
        y = 0
        for kind, lab, v, pa in rows:
            if kind == "head":
                ax.text(-0.02, y, {"GO_BP": "GO", "WikiPathways": "WikiPathways"}.get(lab, lab), transform=ax.get_yaxis_transform(),
                        ha="right", va="center", fontsize=10, weight="bold", color="#333333")
                y += 1
                continue
            sig = pa <= 0.05
            ax.barh(y, v, color=(RED if v > 0 else BLUE) if sig else GREY, height=0.7)
            ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5,
                    color="#222222" if sig else "#8A8A8A")
            ax.text(v + (0.05 if v > 0 else -0.05) * (3 if method == "GSEA" else 6), y,
                    f"{pa:.1e}" if pa < 0.01 else f"{pa:.2f}", va="center", ha="left" if v > 0 else "right",
                    fontsize=7.5, color="#444444")
            y += 1
        ax.set_ylim(y - 0.5, -0.6)
        if method == "GSEA":
            ax.set_xlim(-3.6, 3.6); ax.set_xlabel("GSEA NES (label = padj)", fontsize=9)
        else:
            m = max(abs(r[2]) for r in rows if r[0] == "row")
            ax.set_xlim(-m * 1.5, m * 1.5); ax.set_xlabel("−log10 P, + = in up DEGs, − = in down DEGs (label = padj)", fontsize=8.5)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.set_yticks([]); ax.spines["left"].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
    title = ("GSEA (all genes ranked): top 3 up and top 3 down per database" if method == "GSEA" else
             "ORA (robust DEGs): top 3 sets among up DEGs and top 3 among down DEGs per database")
    fig.text(0.01, 0.975, f"Cardiomyocytes, {title}. Coloured = padj <= 0.05, grey = n.s.", fontsize=10.5, va="top")
    fig.savefig(out / f"hf_cm_top_{method.lower()}.png", dpi=220)
    plt.close(fig)


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    de = pd.read_csv(wd / f"de_{CT}.csv.gz").set_index("gene")
    top_fig(de, out)
    top_pw_fig(wd, out, "GSEA")
    top_pw_fig(wd, out, "ORA")
    gene_fig(de, out)
    for ct, ab, nm in [("Endocardial", "endo", "Endocardial cells"), ("Fibroblast", "fib", "Fibroblasts")]:
        d2 = pd.read_csv(wd / f"de_{ct}.csv.gz").set_index("gene")
        top_fig(d2, out, groups=KEY[ct], name=f"hf_top_{ab}.png")
        gene_fig(d2, out, groups=KEY[ct], name=f"hf_key_{ab}.png", title=f"{nm}: HF-related robust DEGs")
    T = path_fig(wd, de, out)
    pd.set_option("display.width", 220)
    print(T.drop(columns="genes").round(4).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:3])
