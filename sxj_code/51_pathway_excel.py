"""English Excel workbook of the main pathway analysis: every cell type, 5 databases x 2 methods.

Sheets
  README          methods, cut-offs, column dictionary
  Summary         per cell type x database: GSEA sets tested / significant up / down; ORA sets tested and significant
                  among robust up and down DEGs (ORA only where a list has >= 15 robust DEGs)
  <Type>_GSEA     all tested sets of the 5 databases (GO BP, KEGG, Reactome, WikiPathways, MSigDB Hallmark), sorted by
                  database then padj
  <Type>_ORA      all sets overlapping the robust up / down DEGs, 5 databases, sorted by database, list, padj
Inputs: work_dir/genotype/pw_GSEA_<db>_<type>.csv.gz, pw_ORA_<db>_<type>_<Up|Down>.csv.gz (main analysis,
28_pathway_methods; nothing is re-tested here).
The workbook holds unpublished results: write it outside the repository (reports/), do not commit it.
Usage: python 51_pathway_excel.py <work_dir> <out.xlsx>
"""
import sys
from pathlib import Path

import pandas as pd

ORDER = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "Pericyte", "SMC", "Macrophage",
         "Lymphocyte", "Neural_Schwann"]
SHEET = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "Fib", "Endothelial": "EC",
         "Lymphatic_EC": "LEC", "Pericyte": "Pericyte", "SMC": "SMC", "Macrophage": "Mac", "Lymphocyte": "Lympho",
         "Neural_Schwann": "Schwann"}
NAME = {"Lymphatic_EC": "Lymphatic EC", "Neural_Schwann": "Neural Schwann"}
DBS = {"GO_BP": "GO BP", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WikiPathways", "Hallmark": "Hallmark"}
PADJ = 0.05


def nice(db, s):
    if db in ("Reactome", "WikiPathways", "Hallmark") and "_" in s:
        return s.split("_", 1)[1].replace("_", " ").capitalize()
    return s


def main(work_dir, out):
    wd = Path(work_dir) / "genotype"
    summ, sheets = [], []
    for ct in ORDER:
        g_all, o_all = [], []
        for db, dbl in DBS.items():
            f = wd / f"pw_GSEA_{db}_{ct}.csv.gz"
            row = {"Cell type": NAME.get(ct, ct), "Database": dbl}
            if f.exists():
                g = pd.read_csv(f)
                row.update({"GSEA sets tested": len(g), "GSEA sig. up": int(((g.padj <= PADJ) & (g.NES > 0)).sum()),
                            "GSEA sig. down": int(((g.padj <= PADJ) & (g.NES < 0)).sum())})
                g = g.sort_values("padj")
                g_all.append(pd.DataFrame({
                    "Database": dbl, "Pathway": [nice(db, x) for x in g.set], "Set ID": g.set.values,
                    "Size": g["size"].values, "NES": g.NES.round(3).values,
                    "Direction in I4M": ["up" if v > 0 else "down" for v in g.NES],
                    "P": g.p.values, "padj": g.padj.values,
                    "Significant (padj ≤ 0.05)": ["yes" if v <= PADJ else "no" for v in g.padj],
                    "Leading-edge genes": g.leading_edge.fillna("").values}))
            for dr in ("Up", "Down"):
                f = wd / f"pw_ORA_{db}_{ct}_{dr}.csv.gz"
                if f.exists():
                    o = pd.read_csv(f)
                    o = o[o.overlap > 0].sort_values("padj")
                    row[f"ORA ({dr.lower()} DEGs) sets with overlap"] = len(o)
                    row[f"ORA ({dr.lower()} DEGs) sig."] = int((o.padj <= PADJ).sum())
                    o_all.append(pd.DataFrame({
                        "Database": dbl, "DEG list": f"robust {dr.lower()} in I4M",
                        "Pathway": [nice(db, x) for x in o.set], "Set ID": o.set.values, "Set size": o["size"].values,
                        "Overlap": o.overlap.values, "Expected": o.expected.round(2).values,
                        "Fold enrichment": o.fold_enrichment.round(2).values, "P": o.p.values, "padj": o.padj.values,
                        "Significant (padj ≤ 0.05)": ["yes" if v <= PADJ else "no" for v in o.padj],
                        "Overlapping genes": o.genes.fillna("").values}))
                else:
                    row[f"ORA ({dr.lower()} DEGs) sets with overlap"] = "not run"
                    row[f"ORA ({dr.lower()} DEGs) sig."] = "not run"
            summ.append(row)
        nm = NAME.get(ct, ct)
        if g_all:
            sheets.append((f"{SHEET[ct]}_GSEA", f"{nm}: GSEA, all genes ranked by sign(log2FC) × −log10 P; 5 databases",
                           pd.concat(g_all, ignore_index=True)))
        if o_all:
            sheets.append((f"{SHEET[ct]}_ORA", f"{nm}: ORA on robust DEGs (up and down separately); 5 databases",
                           pd.concat(o_all, ignore_index=True)))
    S = pd.DataFrame(summ)
    readme = [
        ["TNNC1-I4M vs WT mouse heart snRNA-seq: pathway analysis per cell type (5 databases × 2 methods)", ""],
        ["", ""],
        ["Samples", "4-week-old mice: WT n = 4 (A–D), TNNC1-I4M (+/−) n = 4 (E–H); 2 males + 2 females per group"],
        ["Gene statistics", "Pseudobulk edgeR quasi-likelihood test per cell type, ~ sex + genotype"],
        ["Databases", "GO Biological Process, KEGG, Reactome, WikiPathways, MSigDB Hallmark (mouse gene sets)"],
        ["GSEA", "fgsea on all tested genes ranked by sign(log2FC) × −log10 P, technical / ambient-RNA genes removed; "
                 "sets of 15–500 genes; NES > 0 = pathway higher in I4M; BH correction within each database and cell type"],
        ["ORA", "Hypergeometric test of robust DEGs (FDR ≤ 0.1, |log2FC| ≥ 0.25, not flagged), up and down separately; "
                "only lists with ≥ 15 genes (cardiomyocyte, endocardial, fibroblast); sets of 10–500 genes; background = "
                "genes tested in that cell type; BH correction within each list and database. Only sets with ≥ 1 "
                "overlapping gene are listed"],
        ["Significance", "padj ≤ 0.05"],
        ["Columns", "NES = normalised enrichment score; Leading-edge genes = genes driving the GSEA score; Expected = "
                    "overlap expected by chance; Fold enrichment = Overlap / Expected"],
        ["Note", "Genotype coincides with sample-processing batch; pathways shared by many cell types (e.g. RNA "
                 "splicing, ribosome biogenesis) may partly reflect batch. Results need experimental validation."],
        ["Sheets", "Summary: counts per cell type and database; <Type>_GSEA and <Type>_ORA: full results. CM = "
                   "cardiomyocyte, Endo = endocardial, Fib = fibroblast, EC = endothelial, LEC = lymphatic endothelial, "
                   "Mac = macrophage, Lympho = lymphocyte"],
        ["Code", "FibrosisAgent/src/mouse_qc/51_pathway_excel.py (results from 28_pathway_methods)"],
    ]
    with pd.ExcelWriter(out, engine="xlsxwriter") as xw:
        wb = xw.book
        hdr = wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "text_wrap": True, "valign": "top"})
        ttl = wb.add_format({"bold": True, "font_size": 12})
        wrap = wb.add_format({"text_wrap": True, "valign": "top"})
        bold = wb.add_format({"bold": True})
        sci = wb.add_format({"num_format": "0.00E+00"})
        pd.DataFrame(readme).to_excel(xw, sheet_name="README", index=False, header=False)
        ws = xw.sheets["README"]
        ws.set_column(0, 0, 18, bold); ws.set_column(1, 1, 110, wrap); ws.write(0, 0, readme[0][0], ttl)
        S.to_excel(xw, sheet_name="Summary", index=False, startrow=1)
        ws = xw.sheets["Summary"]
        ws.write(0, 0, "Number of pathways per cell type and database (significant = padj ≤ 0.05)", ttl)
        for j, col in enumerate(S.columns):
            ws.write(1, j, col, hdr)
        ws.set_column(0, 1, 15); ws.set_column(2, len(S.columns) - 1, 12); ws.set_row(1, 45); ws.freeze_panes(2, 2)
        for name, title, d in sheets:
            d.to_excel(xw, sheet_name=name, index=False, startrow=1)
            ws = xw.sheets[name]
            ws.write(0, 0, title, ttl)
            for j, col in enumerate(d.columns):
                ws.write(1, j, col, hdr)
            ws.set_row(1, 32)
            ws.set_column(0, 1, 13); ws.set_column(2 if "GSEA" in name else 2, 2, 48)
            if name.endswith("GSEA"):
                ws.set_column(3, 3, 30); ws.set_column(4, 7, 10); ws.set_column(6, 7, 11, sci)
                ws.set_column(8, 8, 12); ws.set_column(9, 9, 70)
            else:
                ws.set_column(3, 3, 30); ws.set_column(4, 7, 10); ws.set_column(8, 9, 11, sci)
                ws.set_column(10, 10, 12); ws.set_column(11, 11, 70)
            ws.freeze_panes(2, 3)
            ws.autofilter(1, 0, len(d) + 1, len(d.columns) - 1)
    print(S.to_string(index=False))
    print("sheets:", len(sheets), "->", out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
