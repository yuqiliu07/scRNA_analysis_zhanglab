"""Barcode-rank curves and ambient-plateau estimate for each sample (raw_feature_bc_matrix.h5).

In this batch every partition carries a few hundred UMIs of ambient RNA, so the barcode-rank curve has a long
plateau (empty droplets, ~10^5 barcodes) followed by a cliff to ~0. EmptyDrops' default ambient set
(barcodes <= 100 UMIs) only sees the tail beyond the cliff, which does not represent the empty droplets.
Here we find the cliff (n_partitions = barcodes with > CLIFF_UMI UMIs) and summarise the plateau
(ranks 0.6-1.0 x n_partitions); 01_emptydrops.R then uses the plateau to set `lower`.

Usage: python 00_barcode_rank.py <data_dir> <work_dir> <results_dir>
Output: <work_dir>/emptydrops/plateau.csv, <results_dir>/figures/barcode_rank.png
"""
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SAMPLES = list("ABCDEFGH")
CLIFF_UMI = 50


def totals(h5):
    with h5py.File(h5, "r") as f:
        g = f["matrix"]
        data = g["data"][:]
        indptr = g["indptr"][:]
    t = np.add.reduceat(data, indptr[:-1].clip(max=len(data) - 1)).astype(np.int64)
    t[np.diff(indptr) == 0] = 0
    return t


def main(data_dir, work_dir, res_dir):
    data_dir, work_dir, res_dir = Path(data_dir), Path(work_dir), Path(res_dir)
    (work_dir / "emptydrops").mkdir(parents=True, exist_ok=True)
    (res_dir / "figures").mkdir(parents=True, exist_ok=True)
    rows = []
    fig, ax = plt.subplots(figsize=(7, 5))
    for s in SAMPLES:
        t = np.sort(totals(data_dir / s / "raw_feature_bc_matrix.h5"))[::-1]
        n_cr = h5py.File(data_dir / s / "filtered_feature_bc_matrix.h5", "r")["matrix/barcodes"].shape[0]
        n_part = int((t > CLIFF_UMI).sum())
        plateau = t[int(0.6 * n_part):n_part]
        rows.append(dict(sample=s, barcodes_nonzero=int((t > 0).sum()), n_partitions=n_part,
                         plateau_q25=int(np.percentile(plateau, 25)), plateau_median=int(np.median(plateau)),
                         plateau_q75=int(np.percentile(plateau, 75)), cellranger_cells=n_cr,
                         umi_at_cellranger_edge=int(t[n_cr - 1]),
                         frac_umi_in_plateau=round(plateau.sum() / t.sum(), 3)))
        r = np.arange(1, len(t) + 1)
        k = t > 0
        ax.loglog(r[k], t[k], lw=1, label=s)
    ax.set_xlabel("barcode rank")
    ax.set_ylabel("UMI count")
    ax.set_title("Barcode rank (raw matrix)")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(res_dir / "figures" / "barcode_rank.png", dpi=150)
    df = pd.DataFrame(rows)
    df.to_csv(work_dir / "emptydrops" / "plateau.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:4])
