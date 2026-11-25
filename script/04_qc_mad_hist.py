#!/usr/bin/env python3
"""
Per-sample QC distributions of EmptyDrops-called barcodes (candidate nuclei/cells).

For each sample:
  - take barcodes with is_cell_ed == TRUE from <ED_DIR>/<sample>_emptydrops_full.csv
  - compute nUMI, nGene, pct_mt (genes starting with "mt-") from the raw matrix
  - MAD thresholds: median +/- 3 and 5 x MAD (scaled MAD, = R mad() / scater)
      nUMI, nGene : on log10 scale (skewed counts)
      pct_mt      : on linear scale
  - histogram of each metric with median, +/-3 MAD, +/-5 MAD and pct_mt = 5% marked
Figure: 2 x 4 grid of samples, each sample panel = 3 stacked histograms.

Outputs (in <ED_DIR>/):
  qc_mad_hist.png          the figure
  qc_mad_thresholds.csv    per sample: medians, MADs, thresholds, % removed
  <sample>_qc_metrics.csv  per-cell nUMI / nGene / pct_mt

Usage (Hoffman2 login node is fine; keep BLAS single-threaded):
  module load python/3.12.13
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MALLOC_ARENA_MAX=2
  python3 04_qc_mad_hist.py [--cr-dir ...] [--ed-dir ...]
"""
import argparse
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH", "/u/scratch/y/yuqi")
DEFAULT_CR_DIR = f"{SCRATCH}/FSU_seq/02_cellranger/shared_h5"
DEFAULT_ED_DIR = f"{SCRATCH}/FSU_seq/work_dir/emptydrops_mm10"
DEFAULT_SAMPLE_SHEET = f"{SCRATCH}/FSU_seq/sample_sheet.csv"
MT_CUTOFF = 5.0
MAD_K = (3, 5)
MAD_SCALE = 1.4826


def qc_metrics(h5_path, keep_barcodes, block=200_000):
    """nUMI, nGene, pct_mt for the selected barcodes, reading the CSC matrix in blocks."""
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        names = np.array([x.decode() for x in g["features/name"][:]])
        is_mt = np.char.startswith(np.char.lower(names.astype(str)), "mt-")
        barcodes = np.array([x.decode() for x in g["barcodes"][:]])
        sel = np.isin(barcodes, list(keep_barcodes))
        indptr = g["indptr"][:].astype(np.int64)
        data_ds, idx_ds = g["data"], g["indices"]
        n = len(barcodes)
        umi = np.zeros(n, np.int64)
        ngene = np.zeros(n, np.int64)
        mt = np.zeros(n, np.int64)
        for a in range(0, n, block):
            b = min(a + block, n)
            if not sel[a:b].any():
                continue
            s, e = indptr[a], indptr[b]
            d = data_ds[s:e].astype(np.int64)
            mt_d = np.where(is_mt[idx_ds[s:e]], d, 0)
            loc = indptr[a:b + 1] - s
            cs = np.concatenate(([0], np.cumsum(d)))
            cm = np.concatenate(([0], np.cumsum(mt_d)))
            umi[a:b] = cs[loc[1:]] - cs[loc[:-1]]
            mt[a:b] = cm[loc[1:]] - cm[loc[:-1]]
            ngene[a:b] = np.diff(loc)
    df = pd.DataFrame({"barcode": barcodes[sel], "nUMI": umi[sel],
                       "nGene": ngene[sel], "mt_umi": mt[sel]})
    df["pct_mt"] = 100.0 * df["mt_umi"] / df["nUMI"].clip(lower=1)
    return df, int(is_mt.sum())


def mad_bounds(x, log):
    v = np.log10(x) if log else np.asarray(x, float)
    med = np.median(v)
    mad = MAD_SCALE * np.median(np.abs(v - med))
    out = {"median": 10 ** med if log else med, "mad": mad}
    for k in MAD_K:
        lo, hi = med - k * mad, med + k * mad
        out[f"lo{k}"], out[f"hi{k}"] = (10 ** lo, 10 ** hi) if log else (lo, hi)
    return out


def draw_hist(ax, x, st, label, log, color):
    x = np.asarray(x, float)
    if log:
        # widen the axis so all +/-3 and +/-5 MAD lines are visible
        lo = min(max(x.min(), 1), st["lo5"]) / 1.2
        hi = max(x.max(), st["hi5"]) * 1.2
        bins = np.logspace(np.log10(lo), np.log10(hi), 70)
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        bins = np.linspace(0, max(np.percentile(x, 99.5), MT_CUTOFF * 1.5, st["hi5"]), 60)
    ax.hist(np.clip(x, bins[0], bins[-1]), bins=bins, color=color, alpha=0.75, lw=0)
    ax.axvline(st["median"], color="black", lw=1)
    for k, c in ((3, "#dd6b20"), (5, "#c53030")):
        for side in ("lo", "hi"):
            v = st[f"{side}{k}"]
            if bins[0] <= v <= bins[-1]:
                ax.axvline(v, color=c, ls="--", lw=1)
    if label == "pct_mt":
        ax.axvline(MT_CUTOFF, color="#2f855a", lw=1.3)
    ax.set_ylabel(label, fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_yticks([])
    out5 = ((x < st["lo5"]) | (x > st["hi5"])).mean() * 100
    txt = f"±5MAD out: {out5:.1f}%"
    if label == "pct_mt":
        txt += f"\n>5% mt: {(x > MT_CUTOFF).mean() * 100:.1f}%"
    ax.text(0.98, 0.95, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cr-dir", default=DEFAULT_CR_DIR)
    ap.add_argument("--ed-dir", default=DEFAULT_ED_DIR)
    ap.add_argument("--sample-sheet", default=DEFAULT_SAMPLE_SHEET)
    ap.add_argument("--samples", default="A,B,C,D,E,F,G,H")
    args = ap.parse_args()
    samples = args.samples.split(",")

    meta = {}
    if os.path.exists(args.sample_sheet):
        ss = pd.read_csv(args.sample_sheet)
        meta = {r["sample"]: f"{r['genotype']}, {r['sex']}" for _, r in ss.iterrows()}

    fig = plt.figure(figsize=(20, 11))
    outer = fig.add_gridspec(2, 4, hspace=0.28, wspace=0.12)
    rows = []
    specs = [("nUMI", True, "#2b6cb0"), ("nGene", True, "#6b46c1"), ("pct_mt", False, "#718096")]
    for i, smp in enumerate(samples):
        ed = pd.read_csv(f"{args.ed_dir}/{smp}_emptydrops_full.csv")
        keep = set(ed.loc[ed["is_cell_ed"].astype(str).str.upper() == "TRUE", "barcode"])
        df, n_mt = qc_metrics(f"{args.cr_dir}/{smp}/outs/raw_feature_bc_matrix.h5", keep)
        df.to_csv(f"{args.ed_dir}/{smp}_qc_metrics.csv", index=False)

        inner = GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[i // 4, i % 4], hspace=0.35)
        row = {"sample": smp, "n_ed_cells": len(df), "n_mt_genes": n_mt}
        fail5 = np.zeros(len(df), bool)
        fail3 = np.zeros(len(df), bool)
        for j, (m, log, col) in enumerate(specs):
            st = mad_bounds(df[m].clip(lower=1) if log else df[m], log)
            ax = fig.add_subplot(inner[j])
            draw_hist(ax, df[m], st, m, log, col)
            if j == 0:
                ax.set_title(f"{smp} ({meta.get(smp, '')})  n={len(df):,}", fontsize=10)
            x = df[m].to_numpy(float)
            f5 = (x < st["lo5"]) | (x > st["hi5"])
            f3 = (x < st["lo3"]) | (x > st["hi3"])
            fail5 |= f5
            fail3 |= f3
            row.update({f"{m}_median": round(st["median"], 3),
                        f"{m}_mad{'_log10' if log else ''}": round(st["mad"], 4),
                        f"{m}_lo5": round(st["lo5"], 3), f"{m}_hi5": round(st["hi5"], 3),
                        f"{m}_lo3": round(st["lo3"], 3), f"{m}_hi3": round(st["hi3"], 3),
                        f"{m}_pct_out_5mad": round(f5.mean() * 100, 2),
                        f"{m}_pct_out_3mad": round(f3.mean() * 100, 2)})
        row["pct_mt_gt5"] = round((df["pct_mt"] > MT_CUTOFF).mean() * 100, 2)
        row["any_pct_out_5mad"] = round(fail5.mean() * 100, 2)
        row["any_n_out_5mad"] = int(fail5.sum())
        row["any_pct_out_3mad"] = round(fail3.mean() * 100, 2)
        rows.append(row)
        print(f"{smp}: n={len(df):,}  ±5MAD any={row['any_pct_out_5mad']}%  "
              f"(UMI {row['nUMI_pct_out_5mad']}%, gene {row['nGene_pct_out_5mad']}%, "
              f"mt {row['pct_mt_pct_out_5mad']}%)  mt>5%: {row['pct_mt_gt5']}%")

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="black", lw=1, label="median"),
               Line2D([], [], color="#dd6b20", ls="--", label="median ± 3 MAD"),
               Line2D([], [], color="#c53030", ls="--", label="median ± 5 MAD"),
               Line2D([], [], color="#2f855a", lw=1.3, label="pct_mt = 5%")]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=10, frameon=False,
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("EmptyDrops-called barcodes: nUMI / nGene (log10, MAD on log scale) / pct_mt (linear)",
                 fontsize=13, y=0.995)
    png = f"{args.ed_dir}/qc_mad_hist.png"
    fig.savefig(png, dpi=130, bbox_inches="tight")
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.ed_dir}/qc_mad_thresholds.csv", index=False)
    print(f"\nSaved: {png}\nSaved: {args.ed_dir}/qc_mad_thresholds.csv")


if __name__ == "__main__":
    main()
