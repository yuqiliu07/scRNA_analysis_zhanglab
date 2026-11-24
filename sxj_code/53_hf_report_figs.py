"""HF-focused figures for the collaborator report (cardiomyocytes, endocardium, fibroblasts).

  hf_main_violin.png        classic HF markers per mouse (CM stress / fetal and Ca2+ handling; fibroblast fibrosis):
                            none is a DEG -> the main HF programme is unchanged
  hf_deg_heatmap_<ab>.png   HF-related robust DEGs (genes in the HF categories of 47_hf_deg_classes), per-mouse z-score
                            heatmap + log2FC; CM top 30 up + 30 down by P, endocardium / fibroblasts all
  hf_deg_violin_<ab>.png    single-nucleus violins of the most significant HF-related DEGs (6 up + 6 down, detected)
  hf_same_*.png             only HF-related DEGs / pathways changing in the same direction as in HF (--same)
  hf_pw_<ab>.png            significant HF pathways (77 curated sets, GSEA or ORA padj <= 0.05) with the direction they
                            usually take in heart failure -> same as / opposite to HF
  ora_<ab>.png / .csv        ORA of HF-related robust DEGs, up and down separately, vs all genes tested,
                            top sets by P over 5 databases (--ora)
  sub_markers_<ct>.csv      top marker genes per subcluster (--markers; no expression data needed)
  sub_markerhm_<ct>.png     single-nucleus heatmap of subcluster marker genes (cells grouped by subcluster)
Usage: python 53_hf_report_figs.py <work_dir> <res_dir>
"""
import importlib.util
import re
import textwrap
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))


def _load(name, file):
    sp = importlib.util.spec_from_file_location(name, HERE / file)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


cls = _load("cls", "47_hf_deg_classes.py")
subf = _load("subf", "32_subcluster_figures.py")
plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
RED, BLUE, GREY = "#B2182B", "#2166AC", "#BDBDBD"
CMAP = LinearSegmentedColormap.from_list("bwr2", ["#1A1AA6", "#7B7BD0", "#E8E8E8", "#F07070", "#E8141E"])
CMAP_SC = LinearSegmentedColormap.from_list("purple_yellow", ["#FF00FF", "#000000", "#FFFF00"])
AB = {"Cardiomyocyte": "cm", "Endocardial": "endo", "Fibroblast": "fib"}
MINDET = {"Cardiomyocyte": 0.25, "Endocardial": 0.15, "Fibroblast": 0.08}


def meta_order():
    meta = pd.read_csv(HERE / "sample_sheet.csv").set_index("sample")
    return meta, list(meta.index[meta.genotype == "WT"]) + list(meta.index[meta.genotype == "I4M"])


def dense(X):
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)


# ---------------------------------------------------------------- violins
def violin_grid(a, wd, out, rows, name, title, mindet=0.0, show_rank=None):
    """rows = [(row label, cell type, [genes])]; one violin per mouse."""
    meta, order = meta_order()
    ncol = max(len(g) for _, _, g in rows)
    fig, axes = plt.subplots(len(rows), ncol, figsize=(2.35 * ncol, 2.75 * len(rows) + 0.6), squeeze=False)
    rng = np.random.default_rng(0)
    for r, (lab, ct, genes) in enumerate(rows):
        b = a[(a.obs["cell_type_final"].astype(str) == ct).values]
        de = pd.read_csv(wd / f"de_{ct}.csv.gz").set_index("gene")
        smp = b.obs["sample"].astype(str).values
        fig.text(0.006, 1 - (r + 0.5) / len(rows) * (1 - 0.1) - 0.03, lab, rotation=90, fontsize=11.5, weight="bold",
                 va="center")
        for c in range(ncol):
            ax = axes[r][c]
            if c >= len(genes):
                ax.axis("off")
                continue
            g = genes[c]
            v = dense(b[:, g].X).ravel()
            for j, s in enumerate(order):
                vv = v[smp == s]
                if vv.max() > 0:
                    parts = ax.violinplot(vv, positions=[j], widths=0.85, showextrema=False, bw_method=0.3)
                    for pc in parts["bodies"]:
                        pc.set_facecolor(GCOL[meta.loc[s, "genotype"]]); pc.set_edgecolor("#333333")
                        pc.set_alpha(0.85); pc.set_lw(0.5)
                k = rng.choice(len(vv), min(250, len(vv)), replace=False)
                ax.scatter(j + rng.uniform(-0.22, 0.22, len(k)), vv[k], s=0.7, c="black", lw=0, rasterized=True)
                ax.scatter(j, vv.mean(), s=24, c="white", edgecolor="black", lw=0.9, zorder=5)
            d = de.loc[g]
            sig = bool(d.robust)
            ax.set_title(g + (f" (#{show_rank[ct][g]})" if show_rank and g in show_rank.get(ct, {}) else ""),
                         fontsize=11, weight="bold", style="italic")
            ax.text(0.03, 0.97, f"{d.logFC:+.2f}, FDR {d.FDR:.2g}" + ("" if sig else "  n.s."), transform=ax.transAxes,
                    va="top", fontsize=8, color=(RED if d.logFC > 0 else BLUE) if sig else "#555555")
            ax.set_xticks(range(8), order, fontsize=8)
            ax.axvline(3.5, color="#999999", lw=0.8, ls="--")
            ax.tick_params(axis="y", labelsize=7.5)
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=GCOL[k]) for k in GCOL] +
               [plt.scatter([], [], s=24, c="white", edgecolor="black")], ["WT (A-D)", "I4M (E-H)", "mean"],
               loc="lower center", ncol=3, frameon=False, fontsize=9.5)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0.018, 0.05, 1, 0.95))
    fig.savefig(out / name, dpi=200)
    plt.close(fig)


def hf_genes(wd):
    members, via, classic = cls.load_sets(wd)
    R, _ = cls.classify(wd, members, via, classic)
    return R[R.n_cat > 0]


def detected(a, ct, genes, mindet):
    b = a[(a.obs["cell_type_final"].astype(str) == ct).values]
    genes = [g for g in genes if g in b.var_names]
    det = np.asarray((b[:, genes].X > 0).mean(0)).ravel()
    return [g for g, d in zip(genes, det) if d >= mindet]


# ---------------------------------------------------------------- heatmap of HF-related DEGs
def deg_heatmap(wd, R, ct, out, n):
    meta, order = meta_order()
    de = pd.read_csv(wd / f"de_{ct}.csv.gz").set_index("gene")
    pb = pd.read_csv(wd / f"pb_{ct}.csv.gz", index_col=0)
    cpm = pb[order] / pb[order].sum() * 1e6
    hf = R[R.cell_type == ct].sort_values("PValue")
    up, dn = hf[hf.logFC > 0], hf[hf.logFC < 0]
    nup, ndn = len(up), len(dn)
    up, dn = up.head(n), dn.head(n)
    nr = max(len(up), len(dn))
    H = 0.205 * nr + 1.45
    fig = plt.figure(figsize=(9.6, H))
    bot, top = 0.8 / H, 1 - 0.95 / H
    for k, (d, lab, col, ntot) in enumerate([(up, "higher in I4M", RED, nup), (dn, "lower in I4M", BLUE, ndn)]):
        x0 = 0.1 + k * 0.5
        hh = (top - bot) * len(d) / nr
        hm = fig.add_axes([x0, top - hh, 0.17, hh])
        br = fig.add_axes([x0 + 0.18, top - hh, 0.19, hh])
        lc = np.log2(cpm.loc[d.gene] + 1)
        z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, 1), axis=0).clip(-2, 2)
        hm.imshow(z.values, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
        hm.set_yticks(range(len(d)), d.gene, fontsize=9.5, style="italic")
        hm.set_xticks(range(8), order, fontsize=8.5)
        hm.xaxis.tick_top()
        hm.axvline(3.5, color="white", lw=2.5)
        hm.text(1.5, -1.6, "WT", ha="center", fontsize=9, weight="bold", color=GCOL["WT"])
        hm.text(5.5, -1.6, "I4M", ha="center", fontsize=9, weight="bold", color=GCOL["I4M"])
        for s in hm.spines.values():
            s.set_visible(False)
        hm.tick_params(length=0)
        y = np.arange(len(d))
        br.barh(y, d.logFC, color=col, height=0.7)
        lim = max(hf.logFC.abs().clip(upper=3).max(), 0.5)
        for yy, r in zip(y, d.itertuples()):
            br.text(r.logFC + (0.04 if r.logFC > 0 else -0.04) * lim, yy, f"{r.logFC:+.2f}", va="center",
                    ha="left" if r.logFC > 0 else "right", fontsize=8, color="#333333")
        br.set_ylim(len(d) - 0.5, -0.5)
        br.set_xlim((0, lim * 1.35) if k == 0 else (-lim * 1.35, 0))
        br.axvline(0, color="#333333", lw=0.8)
        br.set_yticks([]); br.spines["left"].set_visible(False)
        br.tick_params(axis="x", labelsize=8)
        br.set_xlabel("log2FC", fontsize=8.5)
        shown = f"top {len(d)} of {ntot}" if ntot > len(d) else f"all {ntot}"
        fig.text(x0 - 0.08, 1 - 0.3 / H, f"HF-related, {lab} ({shown})", fontsize=11, weight="bold", color=col)
    fig.text(0.5, 0.01, f"{ct}: HF-related robust DEGs ranked by P. Heatmap = per-mouse z-score of log2(CPM + 1); "
             "red = high, blue = low", ha="center", fontsize=8.5, color="#555555")
    fig.savefig(out / f"hf_deg_heatmap_{AB[ct]}.png", dpi=220)
    plt.close(fig)
    return hf


# ---------------------------------------------------------------- HF pathways: same / opposite to HF
DIR_THEME = {"Stress / hypertrophy / cardiomyopathy": 1, "Contraction / sarcomere": -1, "Ca2+ handling": -1,
             "Electrical remodelling": -1, "Energy metabolism": -1, "Protein stress / hypoxia": 1, "Fibrosis": 1,
             "Inflammation": 1}
DIR_SET = {"Glycolysis / Gluconeogenesis [KEGG]": 1, "Glycolysis [Hallmark]": 1, "CGMP-PKG signaling pathway [KEGG]": -1,
           "Apelin signaling pathway [KEGG]": -1}
# themes whose direction in HF refers to cardiomyocytes only
CM_ONLY = {"Contraction / sarcomere", "Ca2+ handling", "Electrical remodelling", "Energy metabolism",
           "Stress / hypertrophy / cardiomyopathy"}


def hf_pw(res_dir, ct, out):
    h = pd.read_csv(Path(res_dir) / "genotype" / "hf_pathways.csv")
    gs, og_u, og_d = h[f"padj_{ct}"] <= 0.05, h[f"ORApadj_{ct}_Up"] <= 0.05, h[f"ORApadj_{ct}_Down"] <= 0.05
    d = h[gs | og_u | og_d].copy()
    rows = []
    for r in d.itertuples(index=False):
        nes = getattr(r, f"NES_{ct}")
        obs = np.sign(nes) if getattr(r, f"padj_{ct}") <= 0.05 else (1 if getattr(r, f"ORApadj_{ct}_Up") <= 0.05 else -1)
        exp = DIR_SET.get(r.label, DIR_THEME[r.theme])
        if ct != "Cardiomyocyte" and r.theme in CM_ONLY and r.label not in DIR_SET:
            verdict = "not applicable (CM programme)"
        else:
            verdict = "same as HF" if obs == exp else "opposite to HF"
        rows.append((r.theme, r.label, nes, getattr(r, f"padj_{ct}"), getattr(r, f"ORApadj_{ct}_Up"),
                     getattr(r, f"ORApadj_{ct}_Down"), exp, obs, verdict))
    T = pd.DataFrame(rows, columns=["theme", "label", "NES", "padj", "ora_up", "ora_down", "exp", "obs", "verdict"])
    T = T.sort_values(["theme", "NES"])
    n = len(T)
    H = 0.42 * n + 1.5
    fig = plt.figure(figsize=(11.0, H))
    ax = fig.add_axes([0.38, 0.95 / H, 0.25, (H - 1.55) / H])
    y = np.arange(n)
    for yy, r in zip(y, T.itertuples()):
        sig = r.padj <= 0.05
        ax.barh(yy, r.NES, color=(RED if r.NES > 0 else BLUE) if sig else GREY, height=0.6)
        if r.ora_up <= 0.05:
            ax.scatter(3.1, yy, marker="^", s=70, color=RED, zorder=3)
        if r.ora_down <= 0.05:
            ax.scatter(-3.1, yy, marker="v", s=70, color=BLUE, zorder=3)
        ax.text(-3.55, yy, r.label, ha="right", va="center", fontsize=10, transform=ax.transData)
        ax.text(4.0, yy, "↑" if r.exp > 0 else "↓", ha="center", va="center", fontsize=13, weight="bold")
        ax.text(5.0, yy, "↑" if r.obs > 0 else "↓", ha="center", va="center", fontsize=13, weight="bold",
                color=RED if r.obs > 0 else BLUE)
        vc = {"same as HF": "#B2182B", "opposite to HF": "#1B7837"}.get(r.verdict, "#777777")
        ax.text(5.7, yy, r.verdict, ha="left", va="center", fontsize=10, weight="bold", color=vc)
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(n - 0.5, -0.9)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("GSEA NES (grey = n.s.)   ▲/▼ = ORA significant in up / down DEGs", fontsize=9)
    ax.tick_params(axis="x", labelsize=8.5)
    for x, t in [(4.0, "in HF"), (5.0, "I4M"), (5.7, "verdict")]:
        ax.text(x, -0.85, t, ha="center" if x < 5.5 else "left", va="bottom", fontsize=9.5, weight="bold")
    fig.text(0.02, 1 - 0.35 / H, f"{ct}: {n} of 77 curated HF pathways significant (GSEA or ORA padj ≤ 0.05); the rest unchanged", fontsize=11, weight="bold")
    fig.savefig(out / f"hf_pw_{AB[ct]}.png", dpi=220)
    plt.close(fig)
    T.to_csv(out / f"hf_pw_{AB[ct]}.csv", index=False)


# ---------------------------------------------------------------- subcluster marker heatmap
BAD = re.compile(r"^(mt-|Rp[ls]\d|Gm\d|.*Rik$|Malat1$|Hb[ab]-)")


def marker_table(work_dir, ct, out, k=5):
    """Top k marker genes per subcluster: naming genes first, then the most specific (pct in group - pct in rest)."""
    wd = Path(work_dir) / "subclusters"
    o = pd.read_csv(wd / f"sub_{ct}.csv.gz", index_col=0, dtype={"sub": str})
    m = pd.read_csv(wd / f"markers_{ct}.csv.gz", dtype={"group": str})
    nm = subf.names_for(ct, m)
    rows = []
    for c in sorted(nm, key=int):
        n = int((o["sub"] == c).sum())
        mm = m[(m.group == c) & (m.pvals_adj < 0.05) & ~m.names.str.match(BAD)].head(40).copy()
        mm["spec"] = mm.pct_nz_group - mm.pct_nz_reference
        req = [g for g in subf.NAMES[ct].get(c, ("", []))[1] if g in set(mm.names)]
        sel = req[:2] + [g for g in mm.sort_values("spec", ascending=False).names if g not in req][: k - len(req[:2])]
        rows.append((c, nm[c], n, ", ".join(sel)))
    T = pd.DataFrame(rows, columns=["cluster", "name", "n_nuclei", "markers"])
    T.to_csv(out / f"sub_markers_{ct}.csv", index=False)
    print(T.to_string())


def marker_heatmap(a, work_dir, ct, out, k=4, per=250):
    wd = Path(work_dir) / "subclusters"
    o = pd.read_csv(wd / f"sub_{ct}.csv.gz", index_col=0, dtype={"sub": str})
    m = pd.read_csv(wd / f"markers_{ct}.csv.gz", dtype={"group": str})
    nm = subf.names_for(ct, m)
    keep = [c for c in sorted(nm, key=int) if (o["sub"] == c).sum() >= 50]
    rng = np.random.default_rng(0)
    cells, groups, genes, gstart = [], [], [], []
    for c in keep:
        idx = o.index[o["sub"] == c].to_numpy()
        cells += list(rng.choice(idx, min(per, len(idx)), replace=False))
        groups.append((nm[c], min(per, len(idx))))
        mm = m[(m.group == c) & (m.pvals_adj < 0.05) & ~m.names.str.match(BAD)].head(40).copy()
        mm["spec"] = mm.pct_nz_group - mm.pct_nz_reference
        req = [g for g in subf.NAMES[ct].get(c, ("", []))[1] if g in set(mm.names)]
        sel = req[:2] + [g for g in mm.sort_values("spec", ascending=False).names if g not in req][: k - len(req[:2])]
        sel = [g for g in sel if g not in genes]
        gstart.append((nm[c], len(genes), len(sel)))
        genes += sel
    b = a[cells]
    X = dense(b[:, genes].X).T
    Z = ((X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True).clip(1e-9)).clip(-1.5, 2)
    fig = plt.figure(figsize=(11, 0.2 * len(genes) + 1.6))
    ax = fig.add_axes([0.13, 0.04, 0.72, 0.78])
    axg = fig.add_axes([0.13, 0.83, 0.72, 0.03])
    ax.imshow(Z, aspect="auto", cmap=CMAP_SC, vmin=-1.5, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(genes)), genes, fontsize=8.5, style="italic")
    ax.set_xticks([])
    hue = subf.HUE
    x = 0
    for i, (n_, w) in enumerate(groups):
        axg.add_patch(plt.Rectangle((x, 0), w, 1, color=hue[i % len(hue)]))
        if x > 0:
            ax.axvline(x - 0.5, color="white", lw=1.2)
        x += w
    axg.set_xlim(0, x); axg.set_ylim(0, 1); axg.axis("off")
    for i, (n_, w) in enumerate(groups):
        x0 = sum(g[1] for g in groups[:i])
        fig.text(0.13 + 0.72 * (x0 + w / 2) / x, 0.87, n_.replace(" ", chr(10), 1), ha="center", va="bottom", fontsize=8.5)
    for _, s0, ln in gstart[1:]:
        ax.axhline(s0 - 0.5, color="white", lw=1.2)
    cax = fig.add_axes([0.87, 0.3, 0.012, 0.3])
    plt.colorbar(ax.images[0], cax=cax).ax.tick_params(labelsize=7.5)
    cax.set_title("z", fontsize=8)
    fig.savefig(out / f"sub_markerhm_{ct}.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- direction of each HF-related DEG relative to HF
# literature corrections where pathway membership gives the wrong HF direction
# (anti-fibrotic Ccn5 / Timp3 fall in HF; Hcn1 (If) and Mef2c (hypertrophy) rise; endothelial NO / flow genes fall)
GENE_DIR = {"Ccn5": -1, "Timp3": -1, "Hcn1": 1, "Mef2c": 1, "Nos3": -1, "Klf4": -1, "Klf2": -1, "Kdr": -1,
            "Piezo1": -1, "Bcl6b": -1}


def hf_direction(wd, res_dir):
    """Expected HF direction of every HF-related robust DEG: classic HF genes use their listed direction; other genes
    take the HF direction of the curated HF pathways that contain them (theme default, set-level exceptions in DIR_SET;
    the broad cardiomyopathy theme is ignored when a gene also sits in a more specific theme; CM-only themes are
    ignored outside cardiomyocytes). Genes with conflicting directions are 'ambiguous'."""
    from hf_sets import GENES, THEMES
    h = pd.read_csv(Path(res_dir) / "genotype" / "hf_pathways.csv")
    lab = {(r.db, r.set): r.label for r in h.itertuples()}
    g = {}
    for db in ["GO_BP", "KEGG", "Reactome", "WikiPathways", "Hallmark"]:
        for line in open(wd / f"gmt_{db}.gmt", encoding="utf-8"):
            x = line.rstrip("\n").split("\t")
            g[(db, x[0])] = set(x[2:])
    classic = {(gene, ct): e for L in GENES.values() for gene, e, ct in L}
    R = hf_genes(wd)
    rows = []
    for r in R[R.cell_type.isin(list(AB))].itertuples():
        ct = r.cell_type
        if r.gene in GENE_DIR:
            exp = {GENE_DIR[r.gene]}
        elif (r.gene, ct) in classic:
            exp = {classic[(r.gene, ct)]}
        else:
            d = {}
            for th, L in THEMES.items():
                if ct != "Cardiomyocyte" and th in CM_ONLY:
                    continue
                for k in L:
                    if r.gene in g[k]:
                        l_ = lab.get(k, "")
                        key = "override" if l_ in DIR_SET else th
                        d.setdefault(key, set()).add(DIR_SET.get(l_, DIR_THEME[th]))
            if len(d) > 1:
                d.pop("Stress / hypertrophy / cardiomyopathy", None)
            exp = set().union(*d.values()) if d else set()
        obs = 1 if r.logFC > 0 else -1
        v = "ambiguous" if len(exp) != 1 else ("same" if obs in exp else "opposite")
        rows.append((ct, r.gene, r.logFC, r.PValue, r.FDR, v))
    D = pd.DataFrame(rows, columns=["cell_type", "gene", "logFC", "PValue", "FDR", "verdict"])
    D.to_csv(Path(res_dir) / "genotype" / "hf_deg_direction.csv", index=False)
    return D


def heat_panels(wd, ct, panels, out, name, foot):
    """panels = [(title, colour, DataFrame with gene / logFC)]; per-mouse z-score heatmap + log2FC bar per panel."""
    meta, order = meta_order()
    pb = pd.read_csv(wd / f"pb_{ct}.csv.gz", index_col=0)
    cpm = pb[order] / pb[order].sum() * 1e6
    nr = max(len(d) for _, _, d in panels)
    lim = max(max(d.logFC.abs().clip(upper=3).max() for _, _, d in panels), 0.5)
    H = 0.205 * nr + 1.45
    W = max(4.9 * len(panels), 6.5)
    fig = plt.figure(figsize=(W, H))
    bot, top = 0.8 / H, 1 - 0.95 / H
    for k, (ttl, col, d) in enumerate(panels):
        x0 = (k + 0.2) / len(panels)
        hh = (top - bot) * len(d) / nr
        hm = fig.add_axes([x0, top - hh, 0.34 / len(panels), hh])
        br = fig.add_axes([x0 + 0.36 / len(panels), top - hh, 0.38 / len(panels), hh])
        lc = np.log2(cpm.loc[d.gene] + 1)
        z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, 1), axis=0).clip(-2, 2)
        hm.imshow(z.values, aspect="auto", cmap=CMAP, vmin=-2, vmax=2)
        hm.set_yticks(range(len(d)), d.gene, fontsize=9.5, style="italic")
        hm.set_xticks(range(8), order, fontsize=8.5)
        hm.xaxis.tick_top()
        hm.axvline(3.5, color="white", lw=2.5)
        hm.text(1.5, -1.6, "WT", ha="center", fontsize=9, weight="bold", color=GCOL["WT"])
        hm.text(5.5, -1.6, "I4M", ha="center", fontsize=9, weight="bold", color=GCOL["I4M"])
        for sp_ in hm.spines.values():
            sp_.set_visible(False)
        hm.tick_params(length=0)
        y = np.arange(len(d))
        br.barh(y, d.logFC, color=[RED if v > 0 else BLUE for v in d.logFC], height=0.7)
        for yy, v in zip(y, d.logFC):
            br.text(v + (0.04 if v > 0 else -0.04) * lim, yy, f"{v:+.2f}", va="center", ha="left" if v > 0 else "right",
                    fontsize=8, color="#333333")
        br.set_ylim(len(d) - 0.5, -0.5)
        br.set_xlim((0, lim * 1.35) if d.logFC.iloc[0] > 0 else (-lim * 1.35, 0))
        br.axvline(0, color="#333333", lw=0.8)
        br.set_yticks([]); br.spines["left"].set_visible(False)
        br.tick_params(axis="x", labelsize=8); br.set_xlabel("log2FC", fontsize=8.5)
        fig.text(x0 - 0.15 / len(panels), 1 - 0.3 / H, ttl, fontsize=11, weight="bold", color=col)
    import textwrap
    fig.text(0.5, 0.01, textwrap.fill(foot, int(W * 15)), ha="center", fontsize=8.5, color="#555555")
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def same_figs(a, wd, res_dir, out):
    D = hf_direction(wd, res_dir)
    S = D[D.verdict == "same"].sort_values("PValue")
    foot = ("HF-related robust DEGs, ranked by P. Heatmap = per-mouse z-score of log2(CPM + 1); red = high, blue = low")
    for ct in AB:
        s = S[S.cell_type == ct]
        up, dn = s[s.logFC > 0], s[s.logFC < 0]
        if ct == "Cardiomyocyte":
            for lab, d, col, nm in [("higher in I4M", up, RED, "up"), ("lower in I4M", dn, BLUE, "down")]:
                h2 = (len(d) + 1) // 2
                heat_panels(wd, ct, [(lab[0].upper() + lab[1:], col, d.iloc[:h2]),
                                     ("(continued)", col, d.iloc[h2:])], out, f"hf_same_heatmap_cm_{nm}.png",
                            f"Cardiomyocyte: {foot}")
        else:
            panels = [("Higher in I4M", RED, up)] if len(up) else []
            if len(dn):
                panels.append(("Lower in I4M", BLUE, dn))
            heat_panels(wd, ct, panels, out, f"hf_same_heatmap_{AB[ct]}.png", f"{ct}: {foot}")
        rank = {g: i + 1 for i, g in enumerate(s.gene)}
        nv = 8 if ct == "Cardiomyocyte" else 6
        upv = detected(a, ct, list(up.gene), MINDET[ct])[:nv]
        dnv = detected(a, ct, list(dn.gene), MINDET[ct])[:nv]
        rows = [("Up in I4M", ct, upv)] if upv else []
        rows += [("Down in I4M", ct, dnv)] if dnv else []
        if not dnv:
            more = detected(a, ct, list(up.gene), MINDET[ct])[nv:2 * nv]
            rows += [("Up in I4M (cont.)", ct, more)] if more else []
        violin_grid(a, wd, out, rows, f"hf_same_violin_{AB[ct]}.png",
                    f"{ct}: HF-related DEGs detected in >= {int(MINDET[ct] * 100)}% "
                    "of nuclei (rank in brackets); log2FC / FDR from pseudobulk edgeR", show_rank={ct: rank})
        print(ct, len(up), len(dn), "violin", upv, dnv)
    # significant HF pathways in the same direction as HF, all three cell types
    T = pd.concat([pd.read_csv(out / f"hf_pw_{AB[ct]}.csv").assign(cell_type=ct) for ct in AB])
    T = T[T.verdict == "same as HF"]
    fig, ax = plt.subplots(figsize=(8.5, 0.55 * len(T) + 1.3))
    y = np.arange(len(T))
    ax.barh(y, T.NES, color=[RED if v > 0 else BLUE for v in T.NES], height=0.55)
    for yy, r in zip(y, T.itertuples()):
        ax.text(r.NES + (0.08 if r.NES > 0 else -0.08), yy, f"NES {r.NES:+.2f}, padj {r.padj:.2g}", va="center",
                ha="left" if r.NES > 0 else "right", fontsize=9)
    ax.set_yticks(y, [f"{r.cell_type}: {r.label}" for r in T.itertuples()], fontsize=10)
    ax.set_xlim(-3.5, 3.5); ax.axvline(0, color="#333333", lw=0.8); ax.set_ylim(len(T) - 0.5, -0.5)
    ax.set_xlabel("GSEA NES (+ higher, - lower in I4M)", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(out / "hf_pw_same.png", dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------- ORA of the robust DEGs (all, not only HF-related)
ORA_DB = {"GO_BP": "GO", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WP", "Hallmark": "Hallmark"}


def ora_label(s):
    s = re.sub(r"\s*\(GO:\d+\)$", "", s)
    if re.match(r"^(REACTOME|WP|HALLMARK)_", s):
        s = re.sub(r"^(REACTOME|WP|HALLMARK)_", "", s).replace("_", " ").lower()
        s = s[0].upper() + s[1:]
        for a_, b_ in [("Phase i ", "Phase I "), ("Uv response dn", "UV response down"), ("Uv response up", "UV response up"), ("Rna ", "RNA "), (" ii ", " II "), ("Ncam ", "NCAM "), ("Hsf1 ", "HSF1 "), ("pi3kaktmtor", "PI3K-AKT-mTOR"), ("Mitf m ", "MITF-M "), ("Ncam ", "NCAM ")]:
            s = s.replace(a_, b_)
    return s if len(s) <= 52 else s[:50] + "..."


def ora_all(wd, ct, keep=None, universe_keep=False, direction=None):
    """ORA of the robust DEGs of one cell type (up and down together) against the genes tested in that cell type,
    each of the 5 databases separately (sets with 10-500 genes in the universe; BH within database).
    keep = gene set restricting the DEGs (HF-related genes); universe_keep = also restrict the universe to it;
    direction = 'Up' / 'Down' to use only up- / down-regulated DEGs."""
    from scipy.stats import hypergeom
    de = pd.read_csv(wd / f"de_{ct}.csv.gz")
    if keep is not None and universe_keep:
        de = de[de.gene.isin(keep)]
    U = set(de.gene)
    deg = de[de.robust & (de.gene.isin(keep) if keep is not None else True)]
    up, dn = set(deg.gene[deg.logFC > 0]), set(deg.gene[deg.logFC < 0])
    g = up if direction == "Up" else dn if direction == "Down" else up | dn
    rows = []
    for db, ab in ORA_DB.items():
        R = []
        for line in open(wd / f"gmt_{db}.gmt", encoding="utf-8"):
            x = line.rstrip("\n").split("\t")
            S = set(x[2:]) & U
            if not 10 <= len(S) <= 500:
                continue
            hit = S & g
            p = hypergeom.sf(len(hit) - 1, len(U), len(S), len(g)) if hit else 1.0
            R.append([x[0], ab, len(S), len(hit), len(hit & up), len(hit & dn), p,
                      ",".join(sorted(hit, key=lambda v: deg.set_index("gene").PValue[v]))])
        R = pd.DataFrame(R, columns=["set", "db", "size", "overlap", "n_up", "n_down", "p", "genes"])
        q = R.p.values * len(R) / (R.p.rank(method="first").values)
        o = np.argsort(R.p.values)
        qq = np.minimum.accumulate(q[o][::-1])[::-1]
        R["padj"] = 0.0
        R.loc[R.index[o], "padj"] = np.minimum(qq, 1)
        rows.append(R)
    return pd.concat(rows).sort_values("p"), len(g), len(up), len(dn)


def ora_figs(wd, out, top=10, pmax=0.01, hf=True):
    """ora_<ab>.png / .csv: ORA of HF-related (hf=True) or all robust DEGs, up- and down-regulated genes separately,
    against all genes tested; top sets by p over the 5 databases, skipping sets whose genes mostly repeat an already
    chosen set; bar = -log10 p, solid = FDR <= 0.05."""
    from matplotlib.transforms import blended_transform_factory
    keep = None
    if hf:
        members, _, _ = cls.load_sets(wd)
        keep = set().union(*members.values())
    for ct, ab in AB.items():
        parts, full, ns = [], [], {}
        for d in ["Up", "Down"]:
            T, n, nu, nd = ora_all(wd, ct, keep, direction=d)
            ns[d] = n
            full.append(T.assign(direction=d))
            T = T[(T.p < pmax) & (T.overlap >= 3)]
            nsig = int((T.padj <= 0.05).sum())
            pick, seen = [], []
            for r in T.itertuples():
                g = set(r.genes.split(","))
                if any(len(g & s) / len(g | s) > 0.5 for s in seen):
                    continue
                seen.append(g)
                pick.append(r)
                if len(pick) == top:
                    break
            t = pd.DataFrame(pick).drop(columns="Index") if pick else T.head(0)
            parts.append((d, t.assign(direction=d), nsig))
            print(ct, d, n, "genes; sets FDR<=0.05:", nsig)
            print(t[["set", "db", "overlap", "p", "padj", "genes"]].to_string(max_colwidth=50))
        pd.concat(full).to_csv(out / f"ora_all_{ab}.csv.gz", index=False)
        pd.concat([t for _, t, _ in parts]).to_csv(out / f"ora_{ab}.csv", index=False)
        nrow = sum(max(len(t), 1) + 1 for _, t, _ in parts)
        fig, ax = plt.subplots(figsize=(12.5, 0.33 * nrow + 1.4))
        tr = blended_transform_factory(fig.transFigure, ax.transData)
        y, ticks, labs = 0, [], []
        xmax = max([(-np.log10(t.p)).max() for _, t, _ in parts if len(t)] + [2.5])
        for d, t, nsig in parts:
            col = RED if d == "Up" else BLUE
            ax.text(0.33 - 0.005, y, f"{'Higher' if d == 'Up' else 'Lower'} in I4M ({ns[d]} HF-related DEGs)",
                    transform=tr, ha="right", va="center", fontsize=10.5, weight="bold", color=col)
            if y == 0:
                for xx, h in [(0.595, "DEGs / set"), (0.665, "FDR"), (0.725, "genes (most significant first)")]:
                    ax.text(xx, y, h, transform=tr, fontsize=8, color="#666666", weight="bold", va="center")
            y += 1
            if not len(t):
                ax.text(0.01 * xmax, y, "no set with P < 0.01", va="center", fontsize=9, color="#777777", style="italic")
                y += 1
                continue
            for r in t.itertuples():
                v = -np.log10(r.p)
                ax.barh(y, v, color=col if r.padj <= 0.05 else "white", edgecolor=col, lw=1.2, height=0.65)
                ax.text(0.595, y, f"{r.overlap}/{r.size}", transform=tr, va="center", fontsize=8.3, color="#333333")
                ax.text(0.665, y, f"{r.padj:.2g}", transform=tr, va="center", fontsize=8.3, color="#333333",
                        weight="bold" if r.padj <= 0.05 else "normal")
                ax.text(0.725, y, ", ".join(r.genes.split(",")[:5]), transform=tr, va="center", fontsize=8.3,
                        color="#333333", style="italic")
                ticks.append(y); labs.append(f"{ora_label(r.set)} [{r.db}]")
                y += 1
        ax.set_yticks(ticks, labs, fontsize=9)
        ax.set_ylim(y - 0.4, -0.6)
        ax.set_xlim(0, xmax * 1.05)
        ax.axvline(2, color="#999999", lw=0.6, ls=":")
        ax.set_xlabel("-log10 P (ORA, hypergeometric)", fontsize=9.5)
        for k in ["top", "right"]:
            ax.spines[k].set_visible(False)
        fig.text(0.01, 0.005, textwrap.fill(
            f"{ct}: {'HF-related' if hf else 'all'} robust DEGs (FDR <= 0.1, |log2FC| >= 0.25), up and down separately, "
            "vs all genes tested"
            + (" (HF-related = in the 77 curated HF sets or classic HF genes; pre-selected from HF pathways, so HF-type "
               "sets are expected to rank high)" if hf else "")
            + "; GO BP, KEGG, Reactome, WikiPathways, Hallmark; top 10 sets by P (P < 0.01); solid bar = FDR <= 0.05, "
            "open = FDR > 0.05; numbers = DEGs in set / set size", 190), fontsize=7.8, color="#555555")
        fig.subplots_adjust(left=0.33, right=0.575, top=1 - 0.25 / fig.get_figheight(),
                            bottom=0.9 / fig.get_figheight())
        fig.savefig(out / f"ora_{ab}.png", dpi=200)
        plt.close(fig)


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype"
    if "--ora" in sys.argv:
        ora_figs(wd, out)
        return
    if "--markers" in sys.argv:
        for ct in AB:
            marker_table(work_dir, ct, out)
        return
    a = ad.read_h5ad(Path(work_dir) / "merged" / "final.h5ad", backed="r")
    a = a[a.obs["cell_type_final"].astype(str).isin(list(AB)).values].to_memory()
    if "--same" in sys.argv:
        same_figs(a, wd, res_dir, out)
        return
    R = hf_genes(wd)
    for ct in AB:
        hf_pw(res_dir, ct, out)
        hf = deg_heatmap(wd, R, ct, out, 30)
        rank = {g: i + 1 for i, g in enumerate(hf.gene)}
        up = detected(a, ct, list(hf[hf.logFC > 0].gene), MINDET[ct])[:6]
        dn = detected(a, ct, list(hf[hf.logFC < 0].gene), MINDET[ct])[:6]
        violin_grid(a, wd, out, [("Up in I4M", ct, up), ("Down in I4M", ct, dn)], f"hf_deg_violin_{AB[ct]}.png",
                    f"{ct}: most significant HF-related robust DEGs detected in ≥ {int(MINDET[ct] * 100)}% of nuclei "
                    "(rank among HF-related DEGs in brackets); log2FC / FDR from pseudobulk edgeR",
                    show_rank={ct: rank})
        print(ct, "violin up", up, "down", dn)
    violin_grid(a, wd, out, [("CM: stress / fetal", "Cardiomyocyte", ["Nppa", "Nppb", "Myh7", "Acta1", "Ankrd1", "Myh6"]),
                             ("CM: Ca2+ handling", "Cardiomyocyte", ["Atp2a2", "Pln", "Ryr2", "Slc8a1", "Camk2d", "Cacna1c"]),
                             ("Fibroblast: fibrosis", "Fibroblast", ["Postn", "Col1a1", "Col3a1", "Fn1", "Tgfb1", "Lox"])],
                "hf_main_violin.png", "Classic heart-failure markers per mouse: none is a DEG (n.s. = not significant)")
    for ct in AB:
        marker_heatmap(a, work_dir, ct, out)
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
