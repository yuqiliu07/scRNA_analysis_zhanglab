#!/usr/bin/env python3
"""
Barcode rank / ambient plateau analysis for FSU_seq samples A-H.

For each sample:
  1. Sum UMIs per barcode from Cell Ranger raw_feature_bc_matrix.h5, sort descending.
  2. Cliff: n_partitions = number of barcodes with UMI > 50.
  3. Plateau = barcodes at ranks [0.6, 1.0] x n_partitions; compute UMI q25 / median / q75.
  4. Log-log barcode rank plot marking the Cell Ranger knee (n cells called in
     filtered_feature_bc_matrix.h5) and the plateau range.
  5. Write per-sample summary (incl. plateau_q75) to <work_dir>/emptydrops/plateau.csv.

Usage (on Hoffman2, login node is fine: ~5 min, <2 GB RAM):
  module load python/3.12.13
  python3 02_plateau_barcode_rank.py
"""
import argparse
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH", "/u/scratch/y/yuqi")
DEFAULT_CR_DIR = f"{SCRATCH}/FSU_seq/02_cellranger"
DEFAULT_OUT_DIR = f"{SCRATCH}/FSU_seq/work_dir/emptydrops"
DEFAULT_SAMPLE_SHEET = f"{SCRATCH}/FSU_seq/sample_sheet.csv"

CLIFF_UMI = 50          # barcodes with UMI > this define n_partitions
PLATEAU_LO, PLATEAU_HI = 0.6, 1.0


def umi_per_barcode(h5_path):
    """Total UMI per barcode from a 10x CSC matrix (columns = barcodes)."""
    with h5py.File(h5_path, "r") as f:
        data = f["matrix/data"][:]
        indptr = f["matrix/indptr"][:]
    cs = np.concatenate(([0], np.cumsum(data, dtype=np.int64)))
    return cs[indptr[1:]] - cs[indptr[:-1]]


def n_barcodes(h5_path):
    with h5py.File(h5_path, "r") as f:
        return int(f["matrix/shape"][1])


def analyze(sample, cr_dir):
    outs = f"{cr_dir}/{sample}/outs"
    umi = np.sort(umi_per_barcode(f"{outs}/raw_feature_bc_matrix.h5"))[::-1]
    umi = umi[umi > 0]
    n_cells_cr = n_barcodes(f"{outs}/filtered_feature_bc_matrix.h5")

    n_part = int((umi > CLIFF_UMI).sum())
    lo, hi = int(PLATEAU_LO * n_part), int(PLATEAU_HI * n_part)
    plateau = umi[lo:hi]
    q25, med, q75 = np.percentile(plateau, [25, 50, 75])

    stats = dict(
        sample=sample,
        n_barcodes_raw=len(umi),
        cellranger_cells=n_cells_cr,
        cellranger_knee_umi=int(umi[n_cells_cr - 1]),
        n_partitions=n_part,
        plateau_rank_start=lo + 1,
        plateau_rank_end=hi,
        plateau_q25=float(q25),
        plateau_median=float(med),
        plateau_q75=float(q75),
    )
    return umi, stats


def plot_panel(ax, umi, s, meta):
    ranks = np.arange(1, len(umi) + 1)
    ax.plot(ranks, umi, color="#2b6cb0", lw=1.2)
    # plateau range
    ax.axvspan(s["plateau_rank_start"], s["plateau_rank_end"],
               color="#dd6b20", alpha=0.18, lw=0,
               label=f"plateau {PLATEAU_LO}-{PLATEAU_HI}×n_part")
    ax.axhline(s["plateau_q75"], color="#dd6b20", ls="--", lw=1,
               label=f"plateau q75 = {s['plateau_q75']:.0f}")
    # cliff
    ax.axvline(s["n_partitions"], color="#718096", ls=":", lw=1,
               label=f"cliff n_part = {s['n_partitions']:,}")
    # Cell Ranger knee
    ax.scatter([s["cellranger_cells"]], [s["cellranger_knee_umi"]], s=40,
               color="#c53030", zorder=5,
               label=f"CR cells = {s['cellranger_cells']:,} (UMI {s['cellranger_knee_umi']})")
    ax.axvline(s["cellranger_cells"], color="#c53030", lw=0.8, alpha=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, len(umi) * 1.2)
    ax.set_ylim(0.8, umi[0] * 1.5)
    title = s["sample"] if meta is None else f"{s['sample']}  ({meta})"
    ax.set_title(title, fontsize=11)
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(fontsize=7, loc="lower left", frameon=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cr-dir", default=DEFAULT_CR_DIR)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--sample-sheet", default=DEFAULT_SAMPLE_SHEET)
    ap.add_argument("--samples", default="A,B,C,D,E,F,G,H")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    samples = args.samples.split(",")

    meta = {}
    if os.path.exists(args.sample_sheet):
        ss = pd.read_csv(args.sample_sheet)
        meta = {r["sample"]: f"{r['genotype']}, {r['sex']}" for _, r in ss.iterrows()}

    ncol = 4
    nrow = int(np.ceil(len(samples) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow),
                             sharex=True, sharey=True, squeeze=False)
    rows = []
    for i, smp in enumerate(samples):
        umi, s = analyze(smp, args.cr_dir)
        rows.append(s)
        plot_panel(axes.flat[i], umi, s, meta.get(smp))
        print(f"{smp}: n_part={s['n_partitions']:,}  CR cells={s['cellranger_cells']:,}  "
              f"plateau q25/med/q75 = {s['plateau_q25']:.0f}/{s['plateau_median']:.0f}/{s['plateau_q75']:.0f}")
    for ax in axes.flat[len(samples):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Barcode rank")
    for ax in axes[:, 0]:
        ax.set_ylabel("UMI count")
    fig.suptitle(f"Barcode rank plots (cliff: UMI > {CLIFF_UMI}; "
                 f"plateau: ranks {PLATEAU_LO}-{PLATEAU_HI} × n_partitions)", fontsize=12)
    fig.tight_layout()
    png = f"{args.out_dir}/barcode_rank_plateau.png"
    fig.savefig(png, dpi=150)

    df = pd.DataFrame(rows)
    csv = f"{args.out_dir}/plateau.csv"
    df.to_csv(csv, index=False)
    print(f"\nSaved: {png}\nSaved: {csv}")


if __name__ == "__main__":
    main()
