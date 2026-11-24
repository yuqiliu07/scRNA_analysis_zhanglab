"""Supplementary Excel workbook: DEGs of every DE method and cell type with statistics (WT vs TNNC1-I4M).

Sheets
  说明            methods, models, cut-offs, flags, column dictionary
  汇总_基因数      DEG counts per cell type and method (as Table 1 / 2 of the report)
  DEG_<type>      per cell type: every gene significant in at least one method (FDR <= 0.1, |log2FC| >= 0.25), with
                  logFC / P / FDR of edgeR QL (main), edgeR LRT, limma-voom, voom + quality weights, flags, robust,
                  strict (FDR <= 0.05) and consensus (robust in >= 3 of 4 methods); cell types with no such gene are
                  listed in the summary only
  性别分析         robust DEGs and sex-interaction genes: male / female / interaction contrasts
  内部对照         robust DEGs of CM and endocardium: target-vs-reference interaction tests and category
  亚群内差异 / 亚群汇总  robust DEGs within each subcluster, specificity vs the rest of the cell type
  通路ORA          significant ORA terms (padj <= 0.05), cell types and subclusters
  细胞比例         composition tests, cell types and subclusters
Inputs: work_dir/genotype/{de,dem,sex,ic}_<type>.csv.gz, pw_ORA_*; work_dir/subclusters/sde_*, sora_*;
        results_dir/genotype/*.csv, results_dir/subclusters/*.csv.
The workbook contains unpublished per-gene results: write it outside the repository (reports/), do not commit it.
Usage: python 39_supplementary_excel.py <work_dir> <results_dir> <out.xlsx>
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ORDER = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "Pericyte", "SMC", "Macrophage",
         "Lymphocyte", "Neural_Schwann"]
ZH = {"Cardiomyocyte": "心肌", "Endocardial": "心内膜", "Fibroblast": "成纤维", "Endothelial": "血管内皮",
      "Lymphatic_EC": "淋巴管内皮", "Pericyte": "周细胞", "SMC": "平滑肌", "Macrophage": "巨噬细胞",
      "Lymphocyte": "淋巴细胞", "Neural_Schwann": "神经 / 施万细胞", "Dendritic": "树突状细胞", "Epicardial": "心外膜"}
METH = [("edgeR_QL", "edgeR QL（主分析）"), ("edgeR_LRT", "edgeR LRT"), ("voom", "limma-voom"),
        ("voom_qw", "voom + 样本质量权重")]
FDR, LFC = 0.1, 0.25

spec = importlib.util.spec_from_file_location("icfig", Path(__file__).with_name("20_internal_control_figures.py"))
icfig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(icfig)


def yn(x):
    return x.map({True: "是", False: "否"}) if x.dtype == bool else x


def deg_sheet(wd, t):
    m = pd.read_csv(wd / f"dem_{t}.csv.gz")
    de = pd.read_csv(wd / f"de_{t}.csv.gz")[["gene", "shared", "robust"]]
    m = m.merge(de, on="gene", how="left")
    anysig = m[[f"{k}_sig" for k, _ in METH]].any(axis=1)
    m = m[anysig].copy()
    if m.empty:
        return None
    m["strict"] = m.robust & (m.edgeR_QL_FDR <= 0.05)
    m["direction"] = np.where(m.edgeR_QL_logFC > 0, "上调", "下调")
    cols = {"gene": "基因", "logCPM": "平均表达 (log2 CPM)", "direction": "方向（主分析）"}
    for k, lab in METH:
        cols.update({f"{k}_logFC": f"log2FC | {lab}", f"{k}_P": f"P | {lab}", f"{k}_FDR": f"FDR | {lab}",
                     f"{k}_sig": f"显著 | {lab}", f"{k}_robust": f"稳健 | {lab}"})
        if k == "edgeR_QL":
            cols.update({"technical": "标记：技术基因", "ambient": "标记：背景 RNA", "shared": "标记：多细胞共有",
                         "robust": "主分析稳健 DEG", "strict": "严格标准 (FDR ≤ 0.05) 稳健"})
    cols.update({"n_methods_robust": "稳健的方法数 (0–4)", "consensus": "共识 DEG (≥ 3 种方法)"})
    m = m.sort_values(["robust", "edgeR_QL_FDR"], ascending=[False, True])
    m = m[list(cols)].rename(columns=cols)
    return m.apply(yn)


def main(work_dir, res_dir, out):
    wd, sd = Path(work_dir) / "genotype", Path(work_dir) / "subclusters"
    rg, rs = Path(res_dir) / "genotype", Path(res_dir) / "subclusters"
    sheets = {}

    # ---- summary ----
    ms = pd.read_csv(rg / "de_methods_summary.csv")
    ds = pd.read_csv(rg / "de_summary.csv").set_index("cell_type")
    rows = []
    for t in ORDER:
        if t not in ds.index:
            continue
        g = ms[ms.cell_type == t].set_index("method")
        rows.append({"细胞类型": f"{ZH[t]} ({t})", "检测基因数": int(ds.loc[t, "n_genes_tested"]),
                     "主分析显著 (FDR ≤ 0.1)": int(ds.loc[t, "up"] + ds.loc[t, "down"]),
                     "主分析稳健：上调": int(ds.loc[t, "up_robust"]), "主分析稳健：下调": int(ds.loc[t, "down_robust"]),
                     "主分析稳健：合计": int(g.loc["edgeR_QL", "robust"]),
                     "严格标准 FDR ≤ 0.05 稳健": int(g.loc["edgeR_QL_FDR0.05", "robust"]),
                     "edgeR LRT 稳健": int(g.loc["edgeR_LRT", "robust"]), "limma-voom 稳健": int(g.loc["voom", "robust"]),
                     "voom + 质量权重 稳健": int(g.loc["voom_qw", "robust"]),
                     "共识 (≥ 3 种方法)": int(g.loc["consensus_3of4", "robust"]),
                     "被标记：技术": int(ds.loc[t, "flag_technical"]), "被标记：背景 RNA": int(ds.loc[t, "flag_ambient"]),
                     "被标记：多细胞共有": int(ds.loc[t, "flag_shared"])})
    sheets["汇总_基因数"] = pd.DataFrame(rows)

    # ---- per cell type ----
    no_deg = []
    for t in ORDER:
        if not (wd / f"dem_{t}.csv.gz").exists():
            continue
        d = deg_sheet(wd, t)
        if d is None:
            no_deg.append(ZH[t])
        else:
            sheets[f"DEG_{ZH[t].replace(' / ', '_')}"] = d

    # ---- sex ----
    sx = []
    for t in ORDER:
        f = wd / f"sex_{t}.csv.gz"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        d = d[d.overall_robust | ((d.FDR_interaction <= FDR) & (d.logFC_interaction.abs() >= LFC))].copy()
        d.insert(0, "细胞类型", ZH[t])
        d["sex_interaction"] = (d.FDR_interaction <= FDR) & (d.logFC_interaction.abs() >= LFC)
        sx.append(d)
    sx = pd.concat(sx)
    sx = sx[["细胞类型", "gene", "overall_logFC", "overall_robust", "logFC_male", "P_male", "FDR_male", "logFC_female",
             "P_female", "FDR_female", "logFC_interaction", "P_interaction", "FDR_interaction", "sex_interaction", "flagged"]]
    sx.columns = ["细胞类型", "基因", "log2FC | 主分析（两性合并）", "主分析稳健 DEG", "log2FC | 雄", "P | 雄", "FDR | 雄",
                  "log2FC | 雌", "P | 雌", "FDR | 雌", "log2FC | 交互作用（雌 − 雄）", "P | 交互作用", "FDR | 交互作用",
                  "性别交互基因", "被标记（技术 / 背景 / 共有）"]
    sheets["性别分析"] = sx.apply(yn)

    # ---- internal control ----
    ic = []
    zcat = {"target-specific evidence": "本细胞特有的证据", "same direction in references": "对照细胞同向变化（像批次）",
            "other / inconclusive": "不确定", "not testable": "无法检验（对照细胞不表达）"}
    for t in ["Cardiomyocyte", "Endocardial"]:
        r = icfig.categorise(pd.read_csv(wd / f"ic_{t}.csv.gz"))
        r.insert(0, "细胞类型", ZH[t])
        r["ic_category"] = r.ic_category.map(zcat)
        ic.append(r)
    ic = pd.concat(ic).sort_values(["细胞类型", "FDR_target"], ascending=[False, True])
    icc = {"细胞类型": "细胞类型", "gene": "基因", "logFC_target": "log2FC | 本细胞", "FDR_target": "FDR | 本细胞",
           "logFC_in_Fibroblast": "log2FC | 成纤维", "logFC_int_vs_Fibroblast": "差值 log2FC | vs 成纤维",
           "P_int_vs_Fibroblast": "P | vs 成纤维", "FDR_int_vs_Fibroblast": "FDR | vs 成纤维",
           "logFC_in_Endothelial": "log2FC | 血管内皮", "logFC_int_vs_Endothelial": "差值 log2FC | vs 血管内皮",
           "P_int_vs_Endothelial": "P | vs 血管内皮", "FDR_int_vs_Endothelial": "FDR | vs 血管内皮", "ic_category": "分类"}
    sheets["内部对照"] = ic[list(icc)].rename(columns=icc)

    # ---- subclusters ----
    S = pd.read_csv(rs / "subcluster_de_summary.csv", dtype={"cluster": str})
    sub = []
    for r in S.itertuples():
        d = pd.read_csv(sd / f"sde_{r.cell_type}_{r.cluster}.csv.gz")
        d = d[d.robust].sort_values("PValue")
        d.insert(0, "亚群", r.subcluster); d.insert(0, "细胞类型", ZH[r.cell_type])
        sub.append(d)
    sub = pd.concat(sub)
    spz = {"subcluster-specific": "亚群特有", "shared in cell type": "与同类其他亚群一致", "not_testable": "无法检验"}
    sub["specificity"] = sub.specificity.map(spz)
    sc = {"细胞类型": "细胞类型", "亚群": "亚群", "gene": "基因", "logCPM": "平均表达 (log2 CPM)", "logFC": "log2FC | 亚群内",
          "PValue": "P | 亚群内", "FDR": "FDR | 亚群内", "logFC_celltype": "log2FC | 整类细胞", "FDR_celltype": "FDR | 整类细胞",
          "robust_celltype": "整类分析也是稳健 DEG", "logFC_vs_rest": "差值 log2FC | 亚群 vs 同类其他",
          "FDR_vs_rest": "FDR | 亚群 vs 同类其他", "specificity": "特异性"}
    sheets["亚群内差异"] = sub[list(sc)].rename(columns=sc).apply(yn)
    sz = {"cell_type": "细胞类型", "subcluster": "亚群", "min_nuclei": "每只小鼠最少细胞核", "total_nuclei": "细胞核总数",
          "n_tested": "检测基因数", "deg": "显著 (FDR ≤ 0.1)", "robust": "稳健 DEG", "up": "上调", "down": "下调",
          "also_robust_in_celltype": "整类分析也找到", "only_in_subcluster": "只在亚群中找到",
          "celltype_robust_tested": "整类稳健 DEG 中在亚群检测的", "celltype_robust_same_sign": "其中方向相同的",
          "logFC_cor_celltype": "与整类 log2FC 相关 (全部基因)", "logFC_cor_celltype_robust": "与整类 log2FC 相关 (整类 DEG)",
          "specific": "亚群特有 DEG"}
    S2 = S[list(sz)].rename(columns=sz)
    S2["细胞类型"] = S2["细胞类型"].map(ZH)
    sheets["亚群汇总"] = S2

    # ---- ORA ----
    ora = []
    for db in ["KEGG", "Reactome", "WikiPathways", "Hallmark", "GO_BP"]:
        for t in ORDER:
            for dr in ["Up", "Down"]:
                f = wd / f"pw_ORA_{db}_{t}_{dr}.csv.gz"
                if f.exists():
                    o = pd.read_csv(f); o = o[o.padj <= 0.05]
                    ora.append(o.assign(level="整类细胞", unit=ZH[t], database=db, direction=dr))
    for r in S.itertuples():
        for db in ["KEGG", "Reactome", "WikiPathways", "GO_BP"]:
            for dr in ["Up", "Down"]:
                f = sd / f"sora_{r.cell_type}_{r.cluster}_{db}_{dr}.csv.gz"
                if f.exists():
                    o = pd.read_csv(f); o = o[o.padj <= 0.05]
                    ora.append(o.assign(level="亚群", unit=f"{ZH[r.cell_type]}：{r.subcluster}", database=db, direction=dr))
    ora = pd.concat(ora)
    ora["direction"] = ora.direction.map({"Up": "上调基因", "Down": "下调基因"})
    oc = {"level": "层面", "unit": "细胞 / 亚群", "direction": "输入基因", "database": "数据库", "set": "通路",
          "size": "通路基因数（检测到的）", "overlap": "重合基因数", "expected": "随机预期", "fold_enrichment": "富集倍数",
          "p": "P", "padj": "padj (BH)", "genes": "重合基因（最多 15 个）"}
    sheets["通路ORA"] = ora[list(oc)].rename(columns=oc)

    # ---- composition ----
    c1 = pd.read_csv(rg / "composition_test.csv")
    c1 = c1.assign(层面="整类细胞", 名称=c1.cell_type.map(lambda x: ZH.get(x, x)))
    c2 = pd.read_csv(rg / "subcluster_composition_test.csv")
    c2 = c2.assign(层面="亚群", 名称=[f"{ZH[a]}：{b}" for a, b in zip(c2.cell_type, c2.subcluster)])
    pcol = "p_value" if "p_value" in c1 else "p"
    fcol1 = [c for c in c1.columns if c.lower() in ("fdr", "padj", "p_adj")][0]
    comp = pd.concat([c1[["层面", "名称", "mean_WT", "mean_I4M", pcol, fcol1]].set_axis(["层面", "名称", "WT 平均比例", "I4M 平均比例", "P", "FDR"], axis=1),
                      c2[["层面", "名称", "mean_WT", "mean_I4M", "p", "fdr"]].set_axis(["层面", "名称", "WT 平均比例", "I4M 平均比例", "P", "FDR"], axis=1)])
    sheets["细胞比例"] = comp

    readme = readme_rows(no_deg)
    write(out, readme, sheets)
    for k, v in sheets.items():
        print(f"{k}: {len(v)} rows", flush=True)


def readme_rows(no_deg):
    R = [("TNNC1-I4M vs WT 心脏单核测序：补充表（差异基因与统计值）", ""),
         ("样本", "4 周龄小鼠，WT 4 只（A–D）、TNNC1-I4M (+/−) 4 只（E–H），每组 2 雄 2 雌；比较方向：I4M 相对 WT（log2FC > 0 = I4M 升高）"),
         ("", ""),
         ("一、主分析", ""),
         ("数据", "每只小鼠、每类细胞的 CellBender 计数加和（pseudobulk）；只分析每只小鼠 ≥ 20 个细胞核的细胞类型"),
         ("模型", "edgeR 拟似然（QL）检验，~ 性别 + 基因型；filterByExpr 过滤低表达基因；TMM 标准化"),
         ("差异基因标准", "FDR ≤ 0.1（Benjamini–Hochberg）且 |log2FC| ≥ 0.25（主标准，2026-11-15 起；约 10% 可能为假阳性）"),
         ("严格标准", "同时 FDR ≤ 0.05 的稳健基因（原主标准），列“严格标准 (FDR ≤ 0.05) 稳健”"),
         ("标记：技术基因", "线粒体（mt-）、核糖体蛋白（Rpl / Rps）和背景 / 低质量特征基因（Gm42418、Cmss1、Lars2、Camk1d、Gphn、Gm19951、Cdk8、AY036118、Gm26917、Malat1）"),
         ("标记：背景 RNA", "该基因在另一类细胞中的平均表达 ≥ 本细胞的 5 倍，可能来自其他细胞的游离 RNA"),
         ("标记：多细胞共有", "在 ≥ 3 类细胞中同向显著，更像整只小鼠 / 批次层面的变化"),
         ("稳健 DEG", "显著且没有以上任何标记；PPT 中的差异基因数均指稳健 DEG"),
         ("", ""),
         ("二、方法比较（同样的数据、过滤、标准化和模型，只换检验）", ""),
         ("edgeR QL", "主分析；考虑离散度估计的不确定性，较保守"),
         ("edgeR LRT", "似然比检验；不考虑离散度不确定性，较宽松"),
         ("limma-voom", "log-CPM + 精度权重，稳健经验贝叶斯"),
         ("voom + 样本质量权重", "voomWithQualityWeights，降低噪声较大小鼠的权重"),
         ("共识 DEG", "在 4 种方法中 ≥ 3 种为稳健 DEG（每种方法各自计算“多细胞共有”标记）"),
         ("", ""),
         ("三、其他分析", ""),
         ("性别分析", "~ 0 + 组（WT_M、WT_F、I4M_M、I4M_F），雄、雌各自的 I4M vs WT 和交互作用（雌 − 雄）；每组每种性别只有 2 只，统计力低。表中列出主分析稳健 DEG 和交互作用显著（FDR ≤ 0.1 且 |log2FC| ≥ 0.25）的基因"),
         ("内部对照", "每只小鼠计算 log-CPM（本细胞）− log-CPM（成纤维或血管内皮），limma ~ 性别 + 基因型；批次影响所有细胞，差值中会抵消。只检验对照细胞自身表达的基因。分类：本细胞特有的证据 = 对每个可检验的对照差值 P < 0.05 且同向（无基因达 FDR < 0.05）；对照细胞同向变化 = 对照细胞中同向且幅度 ≥ 本细胞一半"),
         ("亚群内差异", "每个亚群（每只小鼠 ≥ 20 个细胞核，不含低质量 / 疑似双细胞组）单独做与主分析相同的检验；特异性 = 每只小鼠“亚群 − 同类其他细胞”的差值检验（limma），FDR ≤ 0.1、|log2FC| ≥ 0.25 且同向为“亚群特有”"),
         ("通路 ORA", "稳健 DEG（某方向 ≥ 15 个）的超几何检验，背景为检测基因，通路 10–500 个基因，每组内 BH 校正；表中只列 padj ≤ 0.05。基因集：GO BP（org.Mm.eg.db）、KEGG（mmu）、MSigDB 小鼠 Reactome / WikiPathways / Hallmark"),
         ("细胞比例", "每只小鼠的比例，反正弦平方根变换，OLS ~ 性别 + 基因型，BH 校正（整类细胞之间；亚群在同一细胞类型内）"),
         ("", ""),
         ("四、各工作表", ""),
         ("汇总_基因数", "每类细胞、每种方法的差异基因数"),
         ("DEG_<细胞>", "每类细胞一张：至少在一种方法中显著的所有基因，含 4 种方法的 log2FC / P / FDR、各项标记和是否稳健；按“主分析稳健”优先、FDR 从小到大排序。筛选“主分析稳健 DEG = 是”即得主结果"),
         ("没有差异基因的细胞类型", "、".join(no_deg) + "：4 种方法均无显著基因，未单独建表" if no_deg else "无"),
         ("性别分析 / 内部对照 / 亚群内差异 / 亚群汇总 / 通路ORA / 细胞比例", "见上文"),
         ("", ""),
         ("注意", "基因型与样本处理批次重合：WT 样本的背景 / 低质量 RNA 比 I4M 多 2–3 倍。内部对照显示心肌 74–79%、心内膜 87–91% 的 DEG 在对照细胞中也同向变化，部分差异可能来自批次；关键基因需要在新的、同批处理的样本中验证。"),
         ("代码", "src/mouse_qc/13_pseudobulk_de.R、17、19、20、23、29、37 和本表 39_supplementary_excel.py")]
    return pd.DataFrame(R, columns=["项目", "说明"])


def write(out, readme, sheets):
    with pd.ExcelWriter(out, engine="xlsxwriter") as xw:
        wb = xw.book
        hdr = wb.add_format({"bold": True, "bg_color": "#4A1520", "font_color": "#FFFFFF", "text_wrap": True,
                             "valign": "vcenter", "border": 1, "font_name": "Microsoft YaHei", "font_size": 10})
        sci = wb.add_format({"num_format": "0.00E+00", "font_size": 10})
        f3 = wb.add_format({"num_format": "0.000", "font_size": 10})
        pct = wb.add_format({"num_format": "0.0%", "font_size": 10})
        txt = wb.add_format({"font_size": 10})
        wrap = wb.add_format({"text_wrap": True, "valign": "top", "font_name": "Microsoft YaHei", "font_size": 10})
        bold = wb.add_format({"bold": True, "valign": "top", "font_name": "Microsoft YaHei", "font_size": 10})
        title = wb.add_format({"bold": True, "font_size": 14, "font_name": "Microsoft YaHei"})
        readme.to_excel(xw, sheet_name="说明", index=False, header=False)
        ws = xw.sheets["说明"]
        ws.set_column(0, 0, 30, bold); ws.set_column(1, 1, 120, wrap)
        ws.write(0, 0, readme.iloc[0, 0], title)
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name, index=False, startrow=1, header=False)
            ws = xw.sheets[name]
            for j, c in enumerate(df.columns):
                ws.write(0, j, c, hdr)
                cl = str(c)
                if cl.startswith(("P |", "FDR", "P", "padj")) and df[c].dtype.kind == "f":
                    fmt, w = sci, 11
                elif "比例" in cl and df[c].dtype.kind == "f":
                    fmt, w = pct, 11
                elif df[c].dtype.kind == "f":
                    fmt, w = f3, 11
                else:
                    w = min(max(10, int(df[c].astype(str).str.len().quantile(0.95)) + 2 if len(df) else 10, len(cl) // 2 + 2), 60)
                    fmt = txt
                ws.set_column(j, j, w, fmt)
            ws.set_row(0, 45)
            ws.freeze_panes(1, 1 if name.startswith("DEG_") or name == "汇总_基因数" else 2)
            ws.autofilter(0, 0, len(df), len(df.columns) - 1)


if __name__ == "__main__":
    main(*sys.argv[1:4])
