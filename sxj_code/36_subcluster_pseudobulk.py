"""Pseudobulk counts per subcluster (31_subclusters.py), for the within-subcluster WT vs I4M test (37_subcluster_de.R).

For each of CM / endocardial / fibroblast / endothelial, every named subcluster with >= MIN_NUC nuclei in every mouse
(as in 12_genotype_prep.py) gets summed CellBender counts per mouse. Subclusters named low quality, doublet-like or
"Minor" are skipped. The rest of the same cell type (all its other nuclei) is summed too, for the specificity test.

Outputs (local): <work_dir>/subclusters/pbsub_<type>_<k>.csv.gz       genes x mice, the subcluster
                 <work_dir>/subclusters/pbrest_<type>_<k>.csv.gz      genes x mice, the rest of that cell type
                 <work_dir>/subclusters/pbsub_index.csv              type, cluster, name, min / total nuclei
Usage: python 36_subcluster_pseudobulk.py <work_dir>
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

MIN_NUC = 20
TYPES = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial"]
SKIP = ("Low quality", "doublet", "Minor")

spec = importlib.util.spec_from_file_location("subfig", Path(__file__).with_name("32_subcluster_figures.py"))
subfig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subfig)


def main(work_dir):
    wd = Path(work_dir)
    sd = wd / "subclusters"
    samples = pd.read_csv(Path(__file__).with_name("sample_sheet.csv"))["sample"].tolist()
    a = ad.read_h5ad(wd / "merged" / "final.h5ad", backed="r")
    idx = []
    for t in TYPES:
        s = pd.read_csv(sd / f"sub_{t}.csv.gz", index_col=0)
        s["sub"] = s["sub"].astype(str)
        b = a[s.index].to_memory()
        X = b.layers["counts"].tocsr()
        smp = s["sample"].values
        cnt = pd.crosstab(s["sub"], s["sample"]).reindex(columns=samples, fill_value=0)
        for k in sorted(cnt.index, key=int):
            name = subfig.NAMES[t].get(k, (f"{t[:3]}_{k}", []))[0]
            ok = cnt.loc[k].min() >= MIN_NUC and not any(w in name for w in SKIP)
            idx.append({"cell_type": t, "cluster": k, "name": name, "min_nuclei": int(cnt.loc[k].min()),
                        "total_nuclei": int(cnt.loc[k].sum()), "tested": ok})
            if not ok:
                continue
            ins = (s["sub"] == k).values
            for tag, m0 in [("pbsub", ins), ("pbrest", ~ins)]:
                cols = {p: np.asarray(X[m0 & (smp == p)].sum(0)).ravel() for p in samples}
                pd.DataFrame(cols, index=b.var_names).astype(int).to_csv(sd / f"{tag}_{t}_{k}.csv.gz")
            print(t, k, name, "min", cnt.loc[k].min(), flush=True)
    pd.DataFrame(idx).to_csv(sd / "pbsub_index.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
