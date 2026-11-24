"""Excel workbook of the main-analysis DEGs (edgeR QL ~ sex + genotype) split by cell type and direction.

Sheets
  说明              cut-offs, flags, column dictionary
  汇总              per cell type: genes tested, significant up / down, robust up / down, strict (FDR <= 0.05) up / down
  <细胞>_上调 / _下调  every significant gene (FDR <= 0.1, |log2FC| >= 0.25) of that cell type and direction, robust genes
                    first, then by FDR; statistics, flags and per-mouse CPM (A-D WT, E-H I4M)
Cell types without any significant gene in a direction get no sheet for it (they appear in 汇总 with 0).
Inputs: work_dir/genotype/de_<type>.csv.gz, pb_<type>.csv.gz.
The workbook holds unpublished per-gene results: write it outside the repository (reports/), do not commit it.
English version (sheet names Summary, README, <Type>_up / _down): add 'en' as third argument.
Usage: python 45_deg_updown_excel.py <work_dir> <out.xlsx> [zh|en]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ORDER = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "Pericyte", "SMC", "Macrophage",
         "Lymphocyte", "Neural_Schwann"]
ZH = {"Cardiomyocyte": "心肌", "Endocardial": "心内膜", "Fibroblast": "成纤维", "Endothelial": "血管内皮",
      "Lymphatic_EC": "淋巴管内皮", "Pericyte": "周细胞", "SMC": "平滑肌", "Macrophage": "巨噬细胞",
      "Lymphocyte": "淋巴细胞", "Neural_Schwann": "神经施万细胞"}
EN = {"Cardiomyocyte": "Cardiomyocyte", "Endocardial": "Endocardial", "Fibroblast": "Fibroblast",
      "Endothelial": "Endothelial", "Lymphatic_EC": "Lymphatic EC", "Pericyte": "Pericyte", "SMC": "SMC",
      "Macrophage": "Macrophage", "Lymphocyte": "Lymphocyte", "Neural_Schwann": "Neural Schwann"}
SHEET_EN = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC",
            "Lymphatic_EC": "LEC", "Pericyte": "Pericyte", "SMC": "SMC", "Macrophage": "Mac",
            "Lymphocyte": "Lympho", "Neural_Schwann": "Schwann"}
FDR, LFC = 0.1, 0.25
WT, MUT = list("ABCD"), list("EFGH")


LANG = "zh"


def yn(s):
    return s.map({True: "yes", False: "no"} if LANG == "en" else {True: "是", False: "否"})


def main(work_dir, out, lang="zh"):
    global LANG
    LANG = lang
    en = lang == "en"
    wd = Path(work_dir) / "genotype"
    summary, sheets = [], []
    for t in ORDER:
        f = wd / f"de_{t}.csv.gz"
        if not f.exists():
            continue
        de = pd.read_csv(f)
        pb = pd.read_csv(wd / f"pb_{t}.csv.gz", index_col=0)
        cpm = pb / pb.sum() * 1e6
        sig = de[(de.FDR <= FDR) & (de.logFC.abs() >= LFC)].copy()
        sig["strict"] = sig.robust & (sig.FDR <= 0.05)
        up, dn = sig[sig.logFC > 0], sig[sig.logFC < 0]
        vals = [len(de), len(up), len(dn), int(up.robust.sum()), int(dn.robust.sum()), int(sig.robust.sum()),
                int(up.strict.sum()), int(dn.strict.sum())]
        keys = (["Cell type", "Genes tested", "Significant: up", "Significant: down", "Robust: up", "Robust: down",
                 "Robust: total", "Strict (FDR ≤ 0.05) robust: up", "Strict (FDR ≤ 0.05) robust: down"] if en else
                ["细胞类型", "检测基因数", "显著：上调", "显著：下调", "稳健：上调", "稳健：下调", "稳健：合计",
                 "严格 (FDR ≤ 0.05) 稳健：上调", "严格 (FDR ≤ 0.05) 稳健：下调"])
        summary.append(dict(zip(keys, [EN[t] if en else f"{ZH[t]} ({t})"] + vals)))
        for lab, d in ([("up", up), ("down", dn)] if en else [("上调", up), ("下调", dn)]):
            if d.empty:
                continue
            d = d.sort_values(["robust", "FDR"], ascending=[False, True])
            c = cpm.reindex(d.gene)
            cols = (["Gene", "Robust DEG", "log2FC (I4M vs WT)", "Fold change (FC)", "P", "FDR",
                     "Strict robust (FDR ≤ 0.05)", "Mean expression (log2 CPM)", "WT mean CPM", "I4M mean CPM",
                     "Flag: technical", "Flag: ambient RNA", "Flag: shared (≥ 3 cell types)"] if en else
                    ["基因", "稳健 DEG", "log2FC (I4M vs WT)", "倍数变化 (FC)", "P", "FDR", "严格标准 FDR ≤ 0.05 稳健",
                     "平均表达 (log2 CPM)", "WT 平均 CPM", "I4M 平均 CPM", "标记：技术基因", "标记：背景 RNA",
                     "标记：多细胞共有"])
            vals = [d.gene.values, yn(d.robust).values, d.logFC.round(3).values, (2 ** d.logFC).round(2).values,
                    d.PValue.values, d.FDR.values, yn(d.strict).values, d.logCPM.round(2).values,
                    c[WT].mean(1).round(1).values, c[MUT].mean(1).round(1).values, yn(d.technical).values,
                    yn(d.ambient).values, yn(d.shared).values]
            o = pd.DataFrame(dict(zip(cols, vals)))
            for m in WT + MUT:
                o[f"CPM {m} ({'WT' if m in WT else 'I4M'})"] = c[m].round(1).values
            if en:
                sheets.append((f"{SHEET_EN[t]}_{lab}", f"{EN[t]}, {lab} in I4M: {len(d)} significant genes, "
                               f"{int(d.robust.sum())} robust", o))
            else:
                sheets.append((f"{ZH[t]}_{lab}", f"{ZH[t]} ({t}) {lab}基因：{len(d)} 个显著，其中稳健 {int(d.robust.sum())} 个", o))

    readme_en = [
        ["TNNC1-I4M vs WT mouse heart snRNA-seq: differentially expressed genes per cell type (up / down)", ""],
        ["", ""],
        ["Samples", "4-week-old mice: WT n = 4 (A–D), TNNC1-I4M (+/−) n = 4 (E–H); 2 males + 2 females per group"],
        ["Direction", "I4M relative to WT: log2FC > 0 = higher in I4M (up), < 0 = lower (down); FC = 2^log2FC"],
        ["Method", "Pseudobulk counts per mouse and cell type; edgeR quasi-likelihood test, ~ sex + genotype; TMM normalisation"],
        ["Significance", f"FDR ≤ {FDR} and |log2FC| ≥ {LFC} (FDR = Benjamini–Hochberg adjusted P)"],
        ["Robust DEG", "Significant and none of the flags below; DEG counts in the slides are robust DEGs. Robust genes are listed first in each sheet"],
        ["Flag: technical", "Mitochondrial (mt-), ribosomal protein (Rpl / Rps) and ambient / low-quality signature genes"],
        ["Flag: ambient RNA", "Mean expression in another cell type ≥ 5× that in this cell type: may come from free RNA of other cells"],
        ["Flag: shared", "Significant in the same direction in ≥ 3 cell types: more likely a whole-animal / batch-level change"],
        ["Strict", "Robust and FDR ≤ 0.05"],
        ["CPM", "Counts per million of this cell type's pseudobulk (not TMM), for comparing mice"],
        ["Note", "Genotype coincides with sample-processing batch; some differences may come from batch. Key genes need experimental validation."],
        ["Sheets", "Summary: gene counts per cell type; <Type>_up / <Type>_down: all significant genes of that direction (no sheet if none). "
                   "CM = cardiomyocyte, Endo = endocardial, Fib = fibroblast, EC = endothelial, LEC = lymphatic endothelial"],
        ["Code", "FibrosisAgent/src/mouse_qc/45_deg_updown_excel.py (statistics from 13_pseudobulk_de.R)"],
    ]
    readme = [
        ["TNNC1-I4M vs WT 心脏单核测序：各类细胞的差异基因（上调 / 下调）", ""],
        ["", ""],
        ["样本", "4 周龄小鼠，WT 4 只（A–D）、TNNC1-I4M (+/−) 4 只（E–H），每组 2 雄 2 雌"],
        ["方向", "I4M 相对 WT：log2FC > 0 = 突变组升高（上调），< 0 = 降低（下调）；FC = 2^log2FC"],
        ["方法", "每只小鼠、每类细胞的计数加和（pseudobulk），edgeR 拟似然检验 ~ 性别 + 基因型；TMM 标准化"],
        ["显著标准", f"FDR ≤ {FDR} 且 |log2FC| ≥ {LFC}（FDR = BH 校正后的 P 值）"],
        ["稳健 DEG", "显著且没有以下任何标记；PPT 中的差异基因数均指稳健 DEG。每张表稳健基因排在前面"],
        ["标记：技术基因", "线粒体（mt-）、核糖体蛋白（Rpl / Rps）和背景 / 低质量特征基因"],
        ["标记：背景 RNA", "该基因在另一类细胞中的平均表达 ≥ 本细胞的 5 倍，可能来自其他细胞的游离 RNA"],
        ["标记：多细胞共有", "在 ≥ 3 类细胞中同向显著，更像整只小鼠 / 批次层面的变化"],
        ["严格标准", "稳健且 FDR ≤ 0.05"],
        ["CPM", "每百万计数（本细胞类型的 pseudobulk，未经 TMM），用于直观比较每只小鼠的表达量"],
        ["注意", "基因型与样本处理批次重合，部分差异可能来自批次；关键基因需要实验验证。"],
        ["工作表", "汇总：各细胞的基因数；<细胞>_上调 / <细胞>_下调：该方向全部显著基因。某方向没有显著基因则不建表"],
        ["代码", "FibrosisAgent/src/mouse_qc/45_deg_updown_excel.py（统计来自 13_pseudobulk_de.R）"],
    ]
    with pd.ExcelWriter(out, engine="xlsxwriter") as xw:
        wb = xw.book
        hdr = wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "text_wrap": True, "valign": "top"})
        ttl = wb.add_format({"bold": True, "font_size": 12})
        wrap = wb.add_format({"text_wrap": True, "valign": "top"})
        sci = wb.add_format({"num_format": "0.00E+00"})
        bold = wb.add_format({"bold": True})
        upf = wb.add_format({"font_color": "#B2182B"})
        dnf = wb.add_format({"font_color": "#2166AC"})

        rd, rname, sname = (readme_en, "README", "Summary") if en else (readme, "说明", "汇总")
        pd.DataFrame(rd).to_excel(xw, sheet_name=rname, index=False, header=False)
        ws = xw.sheets[rname]
        ws.set_column(0, 0, 18, bold); ws.set_column(1, 1, 110, wrap); ws.write(0, 0, rd[0][0], ttl)

        S = pd.DataFrame(summary)
        S.to_excel(xw, sheet_name=sname, index=False, startrow=1)
        ws = xw.sheets[sname]
        ws.write(0, 0, (f"DEGs per cell type (FDR ≤ {FDR}, |log2FC| ≥ {LFC})" if en else
                        f"各类细胞差异基因数（FDR ≤ {FDR}，|log2FC| ≥ {LFC}）"), ttl)
        for j, col in enumerate(S.columns):
            ws.write(1, j, col, hdr)
        ws.set_column(0, 0, 28); ws.set_column(1, len(S.columns) - 1, 13); ws.set_row(1, 32); ws.freeze_panes(2, 1)

        for name, title, o in sheets:
            o.to_excel(xw, sheet_name=name, index=False, startrow=1)
            ws = xw.sheets[name]
            ws.write(0, 0, title, ttl)
            for j, col in enumerate(o.columns):
                ws.write(1, j, col, hdr)
            ws.set_row(1, 45)
            ws.set_column(0, 0, 14, upf if (name.endswith("上调") or name.endswith("_up")) else dnf)
            ws.set_column(1, 3, 10); ws.set_column(4, 5, 11, sci); ws.set_column(6, len(o.columns) - 1, 10)
            ws.freeze_panes(2, 1)
            ws.autofilter(1, 0, len(o) + 1, len(o.columns) - 1)
    print(S.to_string(index=False))
    print("sheets:", len(sheets), "->", out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
