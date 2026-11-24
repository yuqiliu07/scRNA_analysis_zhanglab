"""Compare CellBender runs with 150 vs 60 epochs on the same samples (A-C were run with both).

For the QC-pass nuclei of each sample (02_per_sample_qc.py): fraction of counts removed, per-gene counts removed
(correlation between runs), per-nucleus denoised UMI (correlation), CellBender cell probability (> 0.5 and > 0.9),
and the final training / test ELBO.

Usage: python compare_cellbender_epochs.py <work_dir> <dir_150> <dir_60> <out_csv> [samples]
"""
import re
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc


def load(f, names):
    a = sc.read_10x_h5(f)
    a.var_names_make_unique()
    with h5py.File(f, "r") as h:
        bc = h["matrix/barcodes"][:].astype(str)
        idx = h["droplet_latents/barcode_indices_for_latents"][:]
        prob = pd.Series(h["droplet_latents/cell_probability"][:], index=bc[idx])
    names = names[names.isin(a.obs_names)]
    return a[names].X, prob.reindex(names), names


def elbo(log):
    txt = Path(log).read_text(errors="ignore")
    tr = re.findall(r"\[epoch (\d+)\]\s+average training loss: ([0-9.]+)", txt)
    te = re.findall(r"\[epoch (\d+)\] average test loss: ([0-9.]+)", txt)
    best = min(float(v) for _, v in te)
    return float(tr[-1][1]), float(te[-1][1]), best


def main(work_dir, d150, d60, out, samples="ABC"):
    work_dir, d150, d60 = Path(work_dir), Path(d150), Path(d60)
    rows = []
    for s in samples:
        r = ad.read_h5ad(work_dir / "per_sample" / f"{s}.h5ad", backed="r")
        q = r.obs_names[r.obs["qc_pass"].values]
        raw = r[q.values].to_memory().X
        r.file.close()
        x1, p1, n1 = load(d150 / s / f"{s}_cellbender.h5", q)
        x2, p2, n2 = load(d60 / s / f"{s}_cellbender.h5", q)
        assert (n1 == q).all() and (n2 == q).all()
        raw_g = np.asarray(raw.sum(0)).ravel()
        g1, g2 = np.asarray(x1.sum(0)).ravel(), np.asarray(x2.sum(0)).ravel()
        rm1, rm2 = raw_g - g1, raw_g - g2
        top = raw_g >= np.percentile(raw_g[raw_g > 0], 90)
        c1, c2 = np.asarray(x1.sum(1)).ravel(), np.asarray(x2.sum(1)).ravel()
        tr1, te1, b1 = elbo(d150 / s / f"{s}_cellbender.log")
        tr2, te2, b2 = elbo(d60 / s / f"{s}_cellbender.log")
        rows.append({
            "sample": s, "n_qc_pass": len(q),
            "frac_removed_150": 1 - g1.sum() / raw_g.sum(), "frac_removed_60": 1 - g2.sum() / raw_g.sum(),
            "r_gene_removed": np.corrcoef(rm1, rm2)[0, 1],
            "r_gene_removed_top10pct_genes": np.corrcoef(rm1[top], rm2[top])[0, 1],
            "r_cell_umi": np.corrcoef(c1, c2)[0, 1],
            "median_cell_umi_150": np.median(c1), "median_cell_umi_60": np.median(c2),
            "median_abs_rel_diff_cell_umi": np.median(np.abs(c1 - c2) / np.maximum(c1, 1)),
            "frac_p_gt_0.5_150": float((p1 > 0.5).mean()), "frac_p_gt_0.5_60": float((p2 > 0.5).mean()),
            "frac_p_gt_0.9_150": float((p1 > 0.9).mean()), "frac_p_gt_0.9_60": float((p2 > 0.9).mean()),
            "cell_call_agreement": float(((p1 > 0.5) == (p2 > 0.5)).mean()),
            "final_train_150": tr1, "final_test_150": te1, "best_test_150": b1,
            "final_train_60": tr2, "final_test_60": te2, "best_test_60": b2,
        })
        print(rows[-1], flush=True)
    df = pd.DataFrame(rows).set_index("sample")
    df.to_csv(out)
    print(df.T.round(4).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:5], *(sys.argv[5:6] or ["ABC"]))
