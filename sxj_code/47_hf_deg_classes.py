"""Classify every robust DEG (main analysis) into the heart-failure (HF) categories of hf_sets.py.

A gene is put in a category when it belongs to any of that category's HF gene sets (hf_sets.THEMES, five databases)
or is one of that category's classic HF genes (hf_sets.GENES). A gene can sit in several categories; genes in none
are "not HF-related". The eight categories are the same as on the HF pathway slides.

Outputs
  results/genotype/styled/hf_deg_class_overview.png  categories x cell types: robust DEGs up / down per category,
                                                      tile colour = enrichment (share of robust DEGs in the category
                                                      vs share of tested genes), * = BH-adjusted P < 0.05
  results/genotype/styled/hf_deg_class_cm.png        cardiomyocytes: per category, number up / down and the most
                                                      significant genes
  results/genotype/styled/hf_deg_class_endo.png      same for endocardial cells
  results/genotype/styled/hf_deg_class_fib.png       same for fibroblasts
  results/genotype/hf_deg_class_summary.csv          counts and enrichment per cell type and category (no gene lists)
  <xlsx>                                             per-gene table (unpublished: keep outside the repository)
Enrichment: one-sided hypergeometric test of robust DEGs (both directions) among the category's genes, universe =
genes tested in that cell type; BH over all tiles. Only cell types with >= 30 robust DEGs are tested (CM, Endo, Fib, as for ORA).
Usage: python 47_hf_deg_classes.py <work_dir> <results_dir> <out.xlsx>
"""
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from hf_sets import GENES, THEMES  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "SMC"]
LAB = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC",
       "Lymphatic_EC": "LEC", "SMC": "SMC"}
ZH = {"Cardiomyocyte": "心肌", "Endocardial": "心内膜", "Fibroblast": "成纤维", "Endothelial": "血管内皮",
      "Lymphatic_EC": "淋巴管内皮", "SMC": "平滑肌"}
CATS = list(THEMES)
ANY = "any HF category"
CZH = {"Stress / hypertrophy / cardiomyopathy": "压力 / 肥大 / 心肌病", "Contraction / sarcomere": "收缩 / 肌节",
       "Ca2+ handling": "钙处理", "Electrical remodelling": "电生理重构", "Energy metabolism": "能量代谢",
       "Protein stress / hypoxia": "蛋白应激 / 缺氧", "Fibrosis": "纤维化", "Inflammation": "炎症",
       ANY: "任一心衰类别"}
# classic-gene groups of hf_sets.GENES -> category
GMAP = {"Stress / fetal genes": CATS[0], "Hypertrophy signalling": CATS[0], "Ca2+ handling": CATS[2],
        "Electrical (K+ / Na+ / gap junction)": CATS[3], "Energy metabolism": CATS[4], "Sarcomere": CATS[1],
        "Fibrosis (fibroblasts)": CATS[6], "Inflammation (macrophages)": CATS[7]}
DBS = ["GO_BP", "KEGG", "Reactome", "WikiPathways", "Hallmark"]
RED, BLUE = "#B2182B", "#2166AC"
MIN = 30  # robust DEGs needed for the enrichment test (CM, Endo, Fib: the cell types that also had ORA)


def load_sets(wd):
    g = {}
    for db in DBS:
        for line in open(wd / f"gmt_{db}.gmt", encoding="utf-8"):
            x = line.rstrip("\n").split("\t")
            g[(db, x[0])] = set(x[2:])
    members, via = {c: set() for c in CATS}, {}
    for c, L in THEMES.items():
        for k in L:
            for gene in g[k]:
                members[c].add(gene)
                via.setdefault((gene, c), []).append(f"{k[0]}: {k[1]}")
    classic = {}
    for grp, L in GENES.items():
        for gene, exp, ct in L:
            members[GMAP[grp]].add(gene)
            classic[gene] = ("↑" if exp > 0 else "↓")
            via.setdefault((gene, GMAP[grp]), []).insert(0, "经典心衰基因")
    return members, via, classic


def classify(wd, members, via, classic):
    rows, summ = [], []
    for ct in CTS:
        de = pd.read_csv(wd / f"de_{ct}.csv.gz")
        tested = set(de.gene)
        rb = de[de.robust].sort_values("FDR")
        nr = len(rb)
        anyset = set().union(*members.values())
        for c in CATS + [ANY]:
            m = (anyset if c == ANY else members[c]) & tested
            hit = rb[rb.gene.isin(m)]
            k, up = len(hit), int((hit.logFC > 0).sum())
            p = hypergeom.sf(k - 1, len(tested), len(m), nr) if nr >= MIN else np.nan
            exp_ = nr * len(m) / len(tested)
            summ.append(dict(cell_type=ct, category=c, genes_tested_in_category=len(m), robust_total=nr,
                             robust_in_category=k, up=up, down=k - up, expected=round(exp_, 2),
                             fold=round(k / exp_, 2) if exp_ > 0 else np.nan, p=p))
        for r in rb.itertuples():
            cats = [c for c in CATS if r.gene in members[c]]
            row = {"cell_type": ct, "gene": r.gene, "dir": "上调" if r.logFC > 0 else "下调", "logFC": r.logFC,
                   "FDR": r.FDR, "PValue": r.PValue, "n_cat": len(cats), "cats": cats,
                   "classic": classic.get(r.gene, ""),
                   "via": "；".join(f"[{CZH[c]}] " + "、".join(via.get((r.gene, c), [])[:4]) for c in cats)}
            rows.append(row)
    S = pd.DataFrame(summ)
    ok = S.p.notna()
    S.loc[ok, "padj"] = multipletests(S.loc[ok, "p"], method="fdr_bh")[1]
    return pd.DataFrame(rows), S


def overview(S, R, out):
    cols = CTS
    rows = CATS + [ANY, "not HF-related"]
    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    cmap = plt.get_cmap("RdBu_r")
    for j, ct in enumerate(cols):
        r = R[R.cell_type == ct]
        for i, c in enumerate(rows):
            if i <= len(CATS):
                s = S[(S.cell_type == ct) & (S.category == c)].iloc[0]
                up, dn, fold, padj = s.up, s.down, s.fold, s.get("padj", np.nan)
            else:
                sub = r[r.n_cat == 0]
                up, dn = int((sub.logFC > 0).sum()), int((sub.logFC < 0).sum())
                fold, padj = np.nan, np.nan
            if pd.notna(padj) and pd.notna(fold):
                col = cmap(0.5 + np.clip(np.log2(max(fold, 1e-3)), -1, 1) / 2 * 0.8)
            else:
                col = "#F2F2F2" if i < len(CATS) else "#E4E4E4"
            ax.add_patch(plt.Rectangle((j, i), 1, 1, color=col, ec="white", lw=2))
            star = " *" if pd.notna(padj) and padj < 0.05 else ""
            ax.text(j + 0.5, i + 0.5, f"↑{up}  ↓{dn}{star}", ha="center", va="center", fontsize=11,
                    color="#222222", weight="bold" if (up + dn) > 0 else "normal")
        ax.text(j + 0.5, -0.15, f"{LAB[ct]}\n({len(r)} robust)", ha="center", va="bottom", fontsize=11)
    ax.axhline(len(CATS), color="#333333", lw=1.2)
    ax.set_xlim(0, len(cols)); ax.set_ylim(len(rows), -0.75)
    ax.set_yticks(np.arange(len(rows)) + 0.5, rows, fontsize=11)
    ax.set_xticks([]); ax.tick_params(axis="y", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-1.25, 1.25))
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("enrichment, log2(observed / expected)\n(grey = not tested, < 30 robust DEGs)", fontsize=9.5)
    fig.suptitle("Robust DEGs (I4M vs WT) sorted into heart-failure categories\n"
                 "↑ / ↓ = number up / down in I4M; a gene can be in several categories; "
                 "* enriched, hypergeometric BH P < 0.05", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "hf_deg_class_overview.png", dpi=220)
    plt.close(fig)


def wrap(genes, n):
    s = ", ".join(genes[:n]) + (f" … (+{len(genes) - n})" if len(genes) > n else "")
    return "\n".join(textwrap.wrap(s, 38))


def detail(R, cts, out, name, ntop):
    fig, axes = plt.subplots(len(cts), 1, figsize=(12, 1.1 * len(CATS) * len(cts) + 1.4), squeeze=False)
    for ax, ct in zip(axes[:, 0], cts):
        r = R[R.cell_type == ct]
        mx = 0
        for i, c in enumerate(CATS):
            sub = r[r.cats.map(lambda x: c in x)]
            u, d = sub[sub.logFC > 0], sub[sub.logFC < 0]
            mx = max(mx, len(u), len(d))
        lim = mx * 1.0 + 1
        for i, c in enumerate(CATS):
            sub = r[r.cats.map(lambda x: c in x)].sort_values("FDR")
            u, d = sub[sub.logFC > 0], sub[sub.logFC < 0]
            ax.barh(i, -len(d), color=BLUE, height=0.6)
            ax.barh(i, len(u), color=RED, height=0.6)
            ax.text(-len(d) - lim * 0.02, i, str(len(d)), ha="right", va="center", fontsize=11.5)
            ax.text(len(u) + lim * 0.02, i, str(len(u)), ha="left", va="center", fontsize=11.5)
            ax.text(-lim * 1.22, i, wrap(list(d.gene), ntop), ha="right", va="center", fontsize=11.5,
                    style="italic", color=BLUE)
            ax.text(lim * 1.22, i, wrap(list(u.gene), ntop), ha="left", va="center", fontsize=11.5,
                    style="italic", color=RED)
        ax.set_yticks(range(len(CATS)), [""] * len(CATS))
        for i, c in enumerate(CATS):
            ax.text(0, i - 0.42, c, ha="center", va="bottom", fontsize=11.5, weight="bold", color="#333333",
                    bbox=dict(fc="white", ec="none", pad=1.5))
        ax.set_ylim(len(CATS) - 0.5, -0.8)
        ax.set_xlim(-lim * 1.22, lim * 1.22)
        ax.axvline(0, color="#333333", lw=0.8)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.set_title(f"{LAB[ct]} ({ct}): {len(r)} robust DEGs, {int((r.n_cat > 0).sum())} in at least one category"
                     f"      ← down in I4M  |  up in I4M →", fontsize=11.5)
    fig.suptitle(f"Robust DEGs per heart-failure category — bars = number of genes, names = most significant first "
                 f"(up to {ntop})", fontsize=12)
    fig.tight_layout(rect=(0.2, 0, 0.8, 0.97))
    fig.subplots_adjust(left=0.33, right=0.67)
    fig.savefig(out / name, dpi=220)
    plt.close(fig)


def excel(R, S, xlsx):
    readme = [
        ["TNNC1-I4M vs WT：全部稳健差异基因按心衰类别分类", ""],
        ["", ""],
        ["稳健 DEG", "主分析（edgeR ~ 性别 + 基因型）中 FDR ≤ 0.1、|log2FC| ≥ 0.25 且无技术 / 背景 / 多细胞共有标记的基因"],
        ["心衰类别", "8 类，与 PPT 心衰通路页相同：" + "、".join(CZH[c] for c in CATS)],
        ["归类规则", "基因属于该类别的任一心衰相关基因集（GO、KEGG、Reactome、WikiPathways、Hallmark 五个库，"
                   "共 " + str(sum(len(v) for v in THEMES.values())) + " 个，见 hf_sets.py），或是该类的经典心衰基因，"
                   "即归入该类。一个基因可属于多类；都不属于的记为“非心衰相关”"],
        ["注意", "基因集按通路功能定义，不分细胞类型；例如心肌里的“纤维化”类基因多为细胞外基质 / 黏附相关基因，"
               "不代表心肌发生纤维化。类别只说明基因的功能归属，不说明是否朝心衰方向变化"],
        ["富集", "超几何检验：该类基因中稳健 DEG 的比例是否高于全部检测基因的平均比例（上下调合计）；"
               "背景 = 该细胞检测到的全部基因；BH 校正跨所有细胞 × 类别；只检验稳健 DEG ≥ 30 的细胞（心肌、心内膜、成纤维，与 ORA 相同）"],
        ["经典心衰基因", "列出心衰时通常的方向（↑ / ↓）；与本数据方向相同即“与心衰同向”"],
        ["代码", "FibrosisAgent/src/mouse_qc/47_hf_deg_classes.py"],
    ]
    with pd.ExcelWriter(xlsx, engine="xlsxwriter") as xw:
        wb = xw.book
        hdr = wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "text_wrap": True, "valign": "top"})
        ttl = wb.add_format({"bold": True, "font_size": 12})
        wrapf = wb.add_format({"text_wrap": True, "valign": "top"})
        sci = wb.add_format({"num_format": "0.00E+00"})
        pd.DataFrame(readme).to_excel(xw, sheet_name="说明", index=False, header=False)
        ws = xw.sheets["说明"]
        ws.set_column(0, 0, 16); ws.set_column(1, 1, 110, wrapf); ws.write(0, 0, readme[0][0], ttl)

        T = S.copy()
        T.insert(0, "细胞类型", T.cell_type.map(lambda c: f"{ZH[c]} ({c})"))
        T.insert(1, "心衰类别", T.category.map(CZH))
        T = T.drop(columns=["cell_type", "category"]).rename(columns={
            "genes_tested_in_category": "该类检测基因数", "robust_total": "该细胞稳健DEG总数",
            "robust_in_category": "该类稳健DEG数", "up": "上调", "down": "下调", "expected": "期望数",
            "fold": "富集倍数", "p": "P", "padj": "BH校正P"})
        T.to_excel(xw, sheet_name="汇总", index=False, startrow=1)
        ws = xw.sheets["汇总"]
        ws.write(0, 0, "各细胞 × 心衰类别的稳健 DEG 数和富集检验", ttl)
        for j, col in enumerate(T.columns):
            ws.write(1, j, col, hdr)
        ws.set_column(0, 0, 26); ws.set_column(1, 1, 20); ws.set_column(2, 8, 11); ws.set_column(9, 10, 11, sci)
        ws.set_row(1, 32); ws.freeze_panes(2, 2)

        for ct in CTS:
            r = R[R.cell_type == ct]
            if r.empty:
                continue
            o = pd.DataFrame({"基因": r.gene.values, "方向": r.dir.values, "log2FC": r.logFC.round(3).values,
                              "P": r.PValue.values, "FDR": r.FDR.values,
                              "心衰类别数": r.n_cat.values,
                              "心衰类别": r.cats.map(lambda x: "、".join(CZH[c] for c in x) or "非心衰相关").values})
            for c in CATS:
                o[CZH[c]] = r.cats.map(lambda x: "●" if c in x else "").values
            o["经典心衰基因（心衰时方向）"] = r.classic.values
            o["与心衰同向"] = [("是" if (k == "↑") == (lf > 0) else "否") if k else ""
                          for k, lf in zip(r.classic, r.logFC)]
            o["归类依据（基因集，每类最多列 4 个）"] = r.via.values
            o = o.sort_values(["心衰类别数", "FDR"], ascending=[False, True], key=lambda s: s if s.name == "FDR"
                              else (s > 0).astype(int))
            name = f"{ZH[ct]}"
            o.to_excel(xw, sheet_name=name, index=False, startrow=1)
            ws = xw.sheets[name]
            nh = int((r.n_cat > 0).sum())
            ws.write(0, 0, f"{ZH[ct]} ({ct})：稳健 DEG {len(r)} 个，其中 {nh} 个属于至少一个心衰类别", ttl)
            for j, col in enumerate(o.columns):
                ws.write(1, j, col, hdr)
            ws.set_row(1, 45)
            ws.set_column(0, 1, 10); ws.set_column(2, 2, 9); ws.set_column(3, 4, 11, sci); ws.set_column(5, 5, 8)
            ws.set_column(6, 6, 30); ws.set_column(7, 14, 7); ws.set_column(15, 16, 10)
            ws.set_column(17, 17, 90)
            ws.freeze_panes(2, 1)
            ws.autofilter(1, 0, len(o) + 1, len(o.columns) - 1)


def main(work_dir, res_dir, xlsx):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    members, via, classic = load_sets(wd)
    R, S = classify(wd, members, via, classic)
    S.to_csv(Path(res_dir) / "genotype" / "hf_deg_class_summary.csv", index=False)
    overview(S, R, out)
    detail(R, ["Cardiomyocyte"], out, "hf_deg_class_cm.png", 12)
    detail(R, ["Endocardial"], out, "hf_deg_class_endo.png", 10)
    detail(R, ["Fibroblast"], out, "hf_deg_class_fib.png", 10)
    excel(R, S, xlsx)
    pd.set_option("display.width", 200)
    print(S.round(3).to_string(index=False))
    print(R.groupby("cell_type").n_cat.apply(lambda x: f"{(x > 0).sum()} / {len(x)}"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
