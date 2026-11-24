"""WT vs TNNC1-I4M: result figures redrawn in the styles of the example figures in the collaborators' slides
(Seurat / ggplot look: split UMAP, feature plots, box + points composition, stacked composition bar, up / down DEG
bars, labelled volcano plots, GO dot plots, violin plots with test labels, sample x gene heatmaps with a pathway bar).

Inputs: final.h5ad (flag == "ok" nuclei), sample_sheet.csv, 12/13 outputs (composition_test.csv, de_summary.csv,
go_top.csv, work_dir/genotype/de_*.csv.gz, pb_*.csv.gz).
Outputs: results_dir/genotype/styled/*.png

Usage: python 15_genotype_figures_styled.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text
except ImportError:  # labels are still drawn, just not repelled
    adjust_text = None

DE_FDR = 0.1  # DE cut-off, as in 13_pseudobulk_de.R
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
ORDER = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
         "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"]
NAME = {"Cardiomyocyte": "Cardiomyocyte", "Fibroblast": "Fibroblast", "Endothelial": "Endothelial",
        "Endocardial": "Endocardial", "Lymphatic_EC": "Lymphatic endothelial", "Pericyte": "Pericyte",
        "SMC": "VSMC", "Macrophage": "Macrophage", "Lymphocyte": "Lymphocyte", "Dendritic": "Dendritic (cDC1)",
        "Epicardial": "Epicardial", "Neural_Schwann": "Neural / Schwann"}
ABBR = {"Cardiomyocyte": "CM", "Fibroblast": "FB", "Endothelial": "EC", "Endocardial": "Endo",
        "Lymphatic_EC": "LEC", "Pericyte": "PC", "SMC": "SMC", "Macrophage": "Mac", "Lymphocyte": "Lym",
        "Dendritic": "DC", "Epicardial": "Epi", "Neural_Schwann": "NS"}
# ggplot2 default hue palette for 12 groups (as in the example Seurat figures)
HUE = ["#F8766D", "#DE8C00", "#B79F00", "#7CAE00", "#00BA38", "#00C08B", "#00BFC4", "#00B4F0", "#619CFF",
       "#C77CFF", "#F564E3", "#FF64B0"]
TCOL = dict(zip(ORDER, HUE))
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}       # example: brown control, green disease
UPC, DNC = "#8B1A2B", "#3434C8"                    # example DEG bars: dark red up, blue down
GLAB = {"WT": "WT", "I4M": "I4M"}
SEED = 0


def save(fig, out, name, dpi=220):
    fig.savefig(out / name, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel()


# ---------------- UMAP ----------------
def umap_split(a, out, basis="X_umap", xlab="umap_1", ylab="umap_2", fname="umap_split_genotype.png"):
    um = a.obsm[basis]
    ct = a.obs["cell_type"].values
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
    for ax, g in zip(axes, ["WT", "I4M"]):
        m = rng.permutation(np.where((a.obs["genotype"] == g).values)[0])
        ax.scatter(um[m, 0], um[m, 1], s=0.6, c=[TCOL[t] for t in ct[m]], rasterized=True, linewidths=0)
        ax.set_title(GLAB[g], fontsize=14, weight="bold")
        ax.set_xlabel(xlab)
        for t in ORDER:
            mm = (ct == t)
            if mm.sum() >= 300:
                ax.text(*np.median(um[mm], 0), ABBR[t], fontsize=9, ha="center", va="center", weight="bold")
    axes[0].set_ylabel(ylab)
    h = [Line2D([], [], marker="o", ls="", color=TCOL[t], markersize=7, label=NAME[t]) for t in ORDER]
    fig.legend(handles=h, loc="center left", bbox_to_anchor=(0.91, 0.5), frameon=False, fontsize=10)
    fig.subplots_adjust(wspace=0.04, right=0.9)
    save(fig, out, fname)


def feature_plots(a, out, genes):
    um = a.obsm["X_umap"]
    cmap = LinearSegmentedColormap.from_list("parula", ["#352A87", "#0F5CDD", "#1481D6", "#06A4CA", "#2EB7A4",
                                                        "#87BF77", "#D1BB59", "#FEC832", "#F9FB0E"])
    fig, axes = plt.subplots(1, len(genes), figsize=(4.2 * len(genes), 4))
    for ax, g in zip(axes, genes):
        v = dense(a[:, g].X)
        o = np.argsort(v)
        sc_ = ax.scatter(um[o, 0], um[o, 1], c=v[o], s=0.4, cmap=cmap, rasterized=True, linewidths=0)
        ax.set_title(g, fontsize=14)
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(True)
    cb = fig.colorbar(sc_, ax=axes, fraction=0.015)
    cb.set_label("Expression")
    save(fig, out, "feature_umap_markers.png")


# ---------------- composition ----------------
def composition(a, meta, ct_test, out):
    cnt = pd.crosstab(a.obs["sample"], a.obs["cell_type"]).reindex(columns=ORDER, fill_value=0)
    prop = cnt.div(cnt.sum(1), axis=0)
    fig, ax = plt.subplots(figsize=(11, 4.2))
    w = 0.34
    for i, t in enumerate(ORDER):
        for j, g in enumerate(["WT", "I4M"]):
            v = prop.loc[meta.index[meta.genotype == g], t].values
            x = i + (j - 0.5) * 0.38
            ax.boxplot(v, positions=[x], widths=w, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=GCOL[g], color="#333333", lw=0.8),
                       medianprops=dict(color="#333333", lw=1), whiskerprops=dict(color="#333333", lw=0.8),
                       capprops=dict(color="#333333", lw=0.8))
            ax.scatter(np.full(len(v), x) + np.linspace(-0.06, 0.06, len(v)), v, s=10, color="black", zorder=3)
        f = ct_test.loc[t, "fdr"]
        if f < 0.05:
            top = prop[t].max() + 0.02
            ax.plot([i - 0.19, i + 0.19], [top, top], color="red", lw=1)
            ax.text(i, top, "*", color="red", ha="center", fontsize=14)
        if i:
            ax.axvline(i - 0.5, color="#d9d9d9", lw=0.8)
    ax.set_xticks(range(len(ORDER)), [NAME[t] for t in ORDER], rotation=40, ha="right")
    ax.set_ylabel("Proportion")
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.legend([Patch(facecolor=GCOL[g], edgecolor="#333333") for g in ["WT", "I4M"]], ["WT", "I4M"],
              loc="upper right", frameon=True, fontsize=10)
    save(fig, out, "composition_box.png")

    # stacked bar by genotype (pooled nuclei, as in the example)
    pc = pd.crosstab(a.obs["genotype"], a.obs["cell_type"], normalize="index").reindex(
        index=["WT", "I4M"], columns=ORDER, fill_value=0)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    bottom = np.zeros(2)
    for t in ORDER[::-1]:
        ax.bar([0, 1], pc[t].values, bottom=bottom, color=TCOL[t], width=0.8, label=NAME[t])
        bottom += pc[t].values
    ax.set_xticks([0, 1], ["WT", "I4M"])
    ax.set_xlabel("Genotype")
    ax.set_ylabel("Cell proportion")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"])
    hd, lb = ax.get_legend_handles_labels()
    ax.legend(hd[::-1], lb[::-1], title="Cell type", loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    save(fig, out, "composition_stacked.png")


# ---------------- DEG bars ----------------
def deg_bars(summ, out):
    s = summ.reindex([t for t in ORDER if t in summ.index])
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    x = np.arange(len(s))
    ax.bar(x, s["up_robust"], color=UPC, width=0.8, label="Up-regulated in I4M", edgecolor="black", lw=0.3)
    ax.bar(x, s["down_robust"], bottom=s["up_robust"], color=DNC, width=0.8, label="Down-regulated in I4M",
           edgecolor="black", lw=0.3)
    ax.set_xticks(x, [NAME[t] for t in s.index], rotation=40, ha="right")
    ax.set_ylabel("# differentially\nexpressed genes")
    ax.legend(frameon=True, loc="upper right")
    save(fig, out, "deg_counts_bar.png")


# ---------------- volcano ----------------
def volcano(de, out, types):
    n = len(types)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.6))
    for ax, t in zip(np.atleast_1d(axes), types):
        d = de[t].copy()
        d["y"] = -np.log10(d["FDR"].clip(1e-300))
        up = d["robust"] & (d["logFC"] > 0)
        dn = d["robust"] & (d["logFC"] < 0)
        flag = (d["FDR"] <= DE_FDR) & (d["logFC"].abs() >= 0.25) & ~d["robust"]
        rest = ~(up | dn | flag)
        ax.scatter(d.loc[rest & (d.logFC > 0), "logFC"], d.loc[rest & (d.logFC > 0), "y"], s=4, c="#F2A3A8", lw=0,
                   rasterized=True)
        ax.scatter(d.loc[rest & (d.logFC <= 0), "logFC"], d.loc[rest & (d.logFC <= 0), "y"], s=4, c="#A3A8E0",
                   lw=0, rasterized=True)
        ax.scatter(d.loc[flag, "logFC"], d.loc[flag, "y"], s=7, facecolors="none", edgecolors="#9A9A9A", lw=0.6)
        ax.scatter(d.loc[up, "logFC"], d.loc[up, "y"], s=7, c="#D7263D", lw=0)
        ax.scatter(d.loc[dn, "logFC"], d.loc[dn, "y"], s=7, c="#1B2A7A", lw=0)
        ax.axhline(-np.log10(DE_FDR), color="#BBBBBB", ls="--", lw=0.8)
        for xv in (-0.25, 0.25):  # |log2FC| cut-off
            ax.axvline(xv, color="#BBBBBB", ls="--", lw=0.8)
        # label the strongest robust genes (|log2FC| x -log10 P), 6 per direction; many FDR values are tied
        d["score"] = d["logFC"].abs() * -np.log10(d["PValue"].clip(1e-300))
        ok = ~d["gene"].str.match(r"^(Gm\d|.*Rik$)")
        lab = pd.concat([d[up & ok].nlargest(6, "score"), d[dn & ok].nlargest(6, "score")])
        texts = [ax.text(r.logFC, r.y, r.gene, fontsize=7,
                         bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#D7263D" if r.logFC > 0 else "#1B2A7A",
                                   lw=0.6)) for r in lab.itertuples()]
        if adjust_text and texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4))
        ax.set_title(NAME[t], fontsize=15)
        ax.set_xlabel("log2FoldChange")
        ax.set_ylabel("-log10(p.adj)")
        ax.grid(color="#EEEEEE", lw=0.6)
        ax.set_axisbelow(True)
    save(fig, out, "volcano_major_types.png")


# ---------------- GO dot plots ----------------
def go_dots(go, out):
    go = go.copy()
    go["FE"] = (go["DE"] / go["n_deg_entrez"]) / (go["N"] / go["universe"])
    go["nl"] = -np.log10(go["FDR"].clip(1e-10))
    blocks = [("Cardiomyocyte", "Up"), ("Cardiomyocyte", "Down"), ("Endocardial", "Up")]
    cmap = LinearSegmentedColormap.from_list("pink", ["#D9D9D9", "#F08A9B"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.2), gridspec_kw={"wspace": 2.6})
    fe_max = min(20.0, float(go["FE"].max()))  # cap so one extreme term does not wash out the scale
    for ax, (t, d) in zip(axes, blocks):
        g = go[(go.cell_type == t) & (go.direction == d)].head(12).iloc[::-1]
        y = np.arange(len(g))
        sc_ = ax.scatter(np.zeros(len(g)), y, s=40 + g["nl"] * 60, c=g["FE"], cmap=cmap, vmin=0, vmax=fe_max,
                         edgecolors="none")
        ax.set_yticks(y, ["GOBP " + s_ if len(s_) < 55 else "GOBP " + s_[:52] + "..." for s_ in g["term"]],
                      fontsize=8.5)
        ax.set_xticks([0], ["robust DEGs"])
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(-0.7, len(g) - 0.3)
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(True)
        ax.set_title(f"{NAME[t]}\n{'up' if d == 'Up' else 'down'} in I4M ({int(g['n_deg'].iloc[0])} genes)", fontsize=12)
    cb = fig.colorbar(sc_, ax=axes, fraction=0.012, pad=0.02, extend="max")
    cb.set_label("Fold Enrichment")
    hs = [plt.scatter([], [], s=40 + v * 60, c="black") for v in [1, 2, 3]]
    fig.legend(hs, ["1", "2", "3"], title="-log10(padj)", loc="upper right", bbox_to_anchor=(0.99, 0.98),
               frameon=False, labelspacing=1.2)
    save(fig, out, "go_dotplot.png")


# ---------------- violins ----------------
def violins(a, de, t, genes, out, name, ncol=4):
    b = a[(a.obs["cell_type"] == t).values]
    d = de[t].set_index("gene")
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 3.2 * nrow))
    rng = np.random.default_rng(SEED)
    for ax, g in zip(axes.ravel(), genes):
        v = dense(b[:, g].X)
        for j, grp in enumerate(["WT", "I4M"]):
            vv = v[(b.obs["genotype"] == grp).values]
            if vv.max() > 0:
                parts = ax.violinplot(vv, positions=[j], widths=0.85, showextrema=False, bw_method=0.3)
                for pc in parts["bodies"]:
                    pc.set_facecolor(GCOL[grp]); pc.set_edgecolor("#333333"); pc.set_alpha(0.85); pc.set_lw(0.6)
            k = rng.choice(len(vv), min(1500, len(vv)), replace=False)
            ax.scatter(j + rng.uniform(-0.25, 0.25, len(k)), vv[k], s=1.2, c="black", lw=0, rasterized=True)
        f = d.loc[g, "FDR"] if g in d.index else np.nan
        lab = "not tested" if f != f else ("ns" if f > DE_FDR else f"FDR = {f:.2g}")
        ax.text(0.04, 0.96, lab, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_title(g, fontsize=13, weight="bold")
        ax.set_xticks([0, 1], ["WT", "I4M"], rotation=45)
        ax.set_ylabel("Expression Level")
    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")
    fig.suptitle(f"{NAME[t]} nuclei (pseudobulk FDR, ~ sex + genotype)", fontsize=12, y=1.0)
    fig.tight_layout()
    save(fig, out, name)


# ---------------- heatmaps (sample x gene, pathway bar) ----------------
def heatmap(pb, de, meta, groups, out, name, title):
    d = de.set_index("gene")
    cpm = pb / pb.sum() * 1e6
    genes, gcol, glab = [], [], []
    pal = ["#A6C82B", "#F28BD7", "#2ED1B9", "#20C452", "#F4A27E", "#F5A623", "#8DA0CB", "#E78AC3"]
    for k, (gname, gl) in enumerate(groups.items()):
        gl = [g for g in gl if g in d.index and g in cpm.index]
        genes += gl; gcol += [pal[k % len(pal)]] * len(gl); glab.append((gname, pal[k % len(pal)]))
    lc = np.log2(cpm.loc[genes] + 1)
    z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, 1), axis=0).clip(-2, 2)
    order = list(meta.index[meta.genotype == "WT"]) + list(meta.index[meta.genotype == "I4M"])
    z = z[order].T
    fig = plt.figure(figsize=(max(8, 0.32 * len(genes) + 3), 4.6))
    ax = fig.add_axes([0.1, 0.22, 0.75, 0.52])
    axb = fig.add_axes([0.1, 0.17, 0.75, 0.035])
    cmap = LinearSegmentedColormap.from_list("bwr2", ["#1A1AA6", "#7B7BD0", "#E8E8E8", "#F07070", "#E8141E"])
    im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-2, vmax=2)
    ax.set_yticks(range(len(order)), [f"{s} ({meta.loc[s, 'genotype']}, {meta.loc[s, 'sex']})" for s in order],
                  fontsize=9)
    ax.axhline(3.5, color="white", lw=3)
    marks = ["*" if (d.loc[g, "FDR"] <= DE_FDR and abs(d.loc[g, "logFC"]) >= 0.25) else ("·" if d.loc[g, "FDR"] <= 0.2 else "") for g in genes]
    ax.set_xticks(range(len(genes)), [f"{g} {m}".strip() for g, m in zip(genes, marks)], rotation=90, fontsize=9,
                  style="italic")
    ax.xaxis.tick_top()
    for s_ in ax.spines.values():
        s_.set_visible(False)
    axb.imshow(np.array([[matplotlib.colors.to_rgb(c) for c in gcol]]), aspect="auto")
    axb.set_xticks([]); axb.set_yticks([])
    for s_ in axb.spines.values():
        s_.set_visible(False)
    cax = fig.add_axes([0.87, 0.3, 0.012, 0.35])
    cax.set_title("z-score", fontsize=8)
    fig.colorbar(im, cax=cax)
    fig.legend([Patch(color=c) for _, c in glab], [n_ for n_, _ in glab], title="Pathway", loc="center left",
               bbox_to_anchor=(0.92, 0.45), frameon=False, fontsize=9)
    fig.text(0.1, 0.93, title, fontsize=12, weight="bold")
    fig.text(0.1, 0.06, "* FDR <= 0.1 and |log2FC| >= 0.25   · FDR <= 0.2 (pseudobulk edgeR, ~ sex + genotype)", fontsize=8, color="#555555")
    save(fig, out, name)


# ---------------- Seurat-style dot plot ----------------
def dotplot(a, groupkey, groups, genes, labels, out, name, figsize):
    X = a[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    gv = a.obs[groupkey].astype(str).values
    pct = np.array([(X[gv == g] > 0).mean(0) for g in groups]) * 100
    avg = np.array([np.expm1(X[gv == g]).mean(0) for g in groups])
    avg = np.log1p(avg)
    sc_ = (avg - avg.mean(0)) / avg.std(0).clip(1e-9)
    sc_ = sc_.clip(-2.5, 2.5)
    fig, ax = plt.subplots(figsize=figsize)
    yy, xx = np.meshgrid(np.arange(len(groups)), np.arange(len(genes)), indexing="ij")
    cm = LinearSegmentedColormap.from_list("seurat", ["#D3D3D3", "#0000FF"])
    p = ax.scatter(xx.ravel(), yy.ravel(), s=pct.ravel() * 1.6, c=sc_.ravel(), cmap=cm, vmin=-2.5, vmax=2.5,
                   edgecolors="none")
    ax.set_xticks(range(len(genes)), genes, rotation=90, style="italic")
    ax.set_yticks(range(len(groups)), labels)
    ax.set_xlim(-0.7, len(genes) - 0.3)
    ax.set_ylim(len(groups) - 0.5, -0.5)
    ax.set_xlabel("Features"); ax.set_ylabel("Identity")
    cb = fig.colorbar(p, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("Average Expression")
    hs = [plt.scatter([], [], s=v * 1.6, c="#555555") for v in [25, 50, 75, 100]]
    ax.legend(hs, ["25", "50", "75", "100"], title="Percent Expressed", loc="upper left", bbox_to_anchor=(1.08, 1),
              frameon=False, labelspacing=1.1)
    save(fig, out, name)


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    out = res_dir / "genotype" / "styled"
    out.mkdir(parents=True, exist_ok=True)
    wd = work_dir / "genotype"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad")
    a = a[(a.obs["flag"] == "ok").values].copy()
    a.obs["cell_type"] = a.obs["cell_type_final"].astype(str)
    a.obs["genotype"] = a.obs["sample"].map(meta["genotype"]).astype(str)
    print(a.n_obs, flush=True)
    de = {t: pd.read_csv(wd / f"de_{t}.csv.gz") for t in ORDER if (wd / f"de_{t}.csv.gz").exists()}
    summ = pd.read_csv(res_dir / "genotype" / "de_summary.csv").set_index("cell_type")
    ct_test = pd.read_csv(res_dir / "genotype" / "composition_test.csv").set_index("cell_type")
    go = pd.read_csv(res_dir / "genotype" / "go_top.csv")

    umap_split(a, out)
    feature_plots(a, out, ["Tnnt2", "Col1a2", "Pecam1", "F13a1"])
    composition(a, meta, ct_test, out)
    deg_bars(summ, out)
    volcano(de, out, ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"])
    go_dots(go, out)
    violins(a, de, "Fibroblast", ["Col1a1", "Col3a1", "Fn1", "Ltbp2", "Fap", "Acta2", "Comp", "Postn"], out,
            "violin_fibroblast_activation.png")
    violins(a, de, "Cardiomyocyte", ["Tnnc1", "Tnni3", "Tnnt2", "Ttn", "Myh6", "Atp2a2", "Nppa", "Cpt1b"], out,
            "violin_cardiomyocyte.png")
    pbf = pd.read_csv(wd / "pb_Fibroblast.csv.gz", index_col=0)[meta.index]
    heatmap(pbf, de["Fibroblast"], meta, {
        "fibril forming collagens": ["Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col5a3", "Col11a1", "Col24a1",
                                     "Col27a1"],
        "FACITs": ["Col9a1", "Col9a3", "Col12a1", "Col14a1", "Col16a1", "Col19a1", "Col21a1"],
        "proteoglycans": ["Cspg4", "Dcn", "Bgn", "Vcan", "Hspg2", "Spock1", "Fmod", "Lum"],
        "multiplexin collagens": ["Col15a1", "Col18a1"],
        "beaded collagens": ["Col6a1", "Col6a2", "Col6a3", "Col6a6", "Col7a1", "Col28a1"],
        "MMPs and TIMPs": ["Mmp2", "Mmp11", "Mmp14", "Mmp15", "Mmp16", "Mmp17", "Mmp19", "Timp1", "Timp2", "Timp3",
                           "Timp4"],
        "activation": ["Postn", "Cthrc1", "Fn1", "Acta2", "Lox", "Fap", "Comp", "Ltbp2"]},
        out, "heatmap_fibroblast_ecm.png", "Fibroblasts: extracellular matrix genes (per-sample pseudobulk)")
    pbc = pd.read_csv(wd / "pb_Cardiomyocyte.csv.gz", index_col=0)[meta.index]
    heatmap(pbc, de["Cardiomyocyte"], meta, {
        "thin filament / sarcomere": ["Tnnc1", "Tnni3", "Tnnt2", "Tpm1", "Actc1", "Myh6", "Myh7", "Myl2", "Myl3",
                                      "Mybpc3", "Ttn", "Actn2"],
        "calcium handling": ["Atp2a2", "Pln", "Ryr2", "Cacna1c", "Slc8a1", "Casq2", "Camk2d", "Trdn"],
        "stress / hypertrophy": ["Nppa", "Nppb", "Acta1", "Ankrd1", "Rcan1", "Fhl1", "Myh7b"],
        "metabolism": ["Ppargc1a", "Ppara", "Cpt1b", "Acadm", "Hadha", "Pdk4", "Ucp3", "Slc2a4", "Hk2", "Ldha"]},
        out, "heatmap_cardiomyocyte_panel.png", "Cardiomyocytes: sarcomere, calcium, stress and metabolism genes")

    # marker dot plot (cell types) and EC subclusters
    mk = pd.read_csv(res_dir / "provider_compare" / "cmp_markers_ours.csv").set_index("cell_type")
    canon = {"Cardiomyocyte": ["Tnnt2", "Myh6", "Ryr2"], "Fibroblast": ["Dcn", "Pdgfra", "Col1a2"],
             "Endothelial": ["Pecam1", "Flt1", "Ptprb"], "Endocardial": ["Npr3", "Tmem108", "Cdh11"],
             "Lymphatic_EC": ["Prox1", "Mmrn1", "Reln"], "Pericyte": ["Pdgfrb", "Abcc9", "Kcnj8"],
             "SMC": ["Myh11", "Acta2", "Kcnab1"], "Macrophage": ["F13a1", "Mrc1", "Cd163"],
             "Lymphocyte": ["Ptprc", "Skap1", "Pax5"], "Dendritic": ["Clec9a", "Flt3", "Wdfy4"],
             "Epicardial": ["Msln", "Upk3b", "Muc16"], "Neural_Schwann": ["Plp1", "Cdh19", "Scn7a"]}
    genes = []
    for t in ORDER:
        for g in canon[t]:
            if g in a.var_names and g not in genes:
                genes.append(g)
    dotplot(a, "cell_type", ORDER, genes, [NAME[t] for t in ORDER], out, "dotplot_markers.png", (15, 5.2))
    ec = a[a.obs["state"].astype(str).isin(["EC_arterial", "EC_capillary_high_UMI", "EC_capillary_low_UMI",
                                            "Endocardial_1", "Endocardial_2", "Lymphatic_EC"]).values].copy()
    ec.obs["st"] = ec.obs["state"].astype(str).replace({"EC_capillary_high_UMI": "EC_capillary",
                                                        "EC_capillary_low_UMI": "EC_capillary", "Endocardial_1": "Endocardial",
                                                        "Endocardial_2": "Endocardial"})  # the two endocardial states differ by RNA quality
    sts = ["EC_capillary", "EC_arterial", "Endocardial", "Lymphatic_EC"]
    eg = ["Ptprg", "Plcb4", "Kdr", "Car4", "Arl15", "Eln", "Prdm16", "Gja5", "Npr3", "Tmem108", "Cdh11", "Hmcn1",
          "Chrm3", "Cfh", "Pde3a", "Prox1", "Reln", "Mmrn1", "Flt4"]
    eg = [g for g in eg if g in ec.var_names]
    dotplot(ec, "st", sts, eg, ["Capillary EC", "Arterial EC", "Endocardial", "Lymphatic EC"], out,
            "dotplot_ec_subclusters.png", (10.5, 3.1))
    print("done", sorted(p.name for p in out.glob("*.png")), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
