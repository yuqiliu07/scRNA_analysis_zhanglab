#!/usr/bin/env python3
"""
Summary figure + combined table for 05_qc_doublet.R outputs.

Left  : per sample, fraction of EmptyDrops nuclei that are
        pass-QC singlets / pass-QC doublets / fail QC.
Right : nUMI (log10) of pass-QC singlets vs predicted doublets per sample
        (sanity check: doublets should carry more UMIs).

Outputs (in <QC_DIR>/):
  qc_doublet_summary.png
  qc_doublet_summary_all_samples.csv

Usage (Hoffman2 login node is fine):
  module load python/3.12.13
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MALLOC_ARENA_MAX=2
  python3 06_qc_doublet_summary_plot.py [--qc-dir ...]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH", "/u/scratch/y/yuqi")
DEFAULT_QC_DIR = f"{SCRATCH}/FSU_seq/work_dir/qc_mm10"
DEFAULT_SAMPLE_SHEET = f"{SCRATCH}/FSU_seq/sample_sheet.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", default=DEFAULT_QC_DIR)
    ap.add_argument("--sample-sheet", default=DEFAULT_SAMPLE_SHEET)
    ap.add_argument("--samples", default="A,B,C,D,E,F,G,H")
    args = ap.parse_args()
    samples = args.samples.split(",")

    meta = {}
    if os.path.exists(args.sample_sheet):
        ss = pd.read_csv(args.sample_sheet)
        meta = {r["sample"]: f"{r['genotype']}" for _, r in ss.iterrows()}

    summ = pd.concat([pd.read_csv(f"{args.qc_dir}/{s}_qc_doublet_summary.csv") for s in samples])
    summ.to_csv(f"{args.qc_dir}/qc_doublet_summary_all_samples.csv", index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios": [1, 1.3]})
    labels = [f"{s}\n{meta.get(s, '')}" for s in samples]
    x = np.arange(len(samples))

    n = summ["n_nuclei"].to_numpy()
    singlet = summ["n_pass_qc_singlet"].to_numpy() / n * 100
    dbl_pass = summ["n_doublet_in_pass"].to_numpy() / n * 100
    fail = (n - summ["n_pass_qc"].to_numpy()) / n * 100
    ax1.bar(x, singlet, color="#2b6cb0", label="pass QC, singlet")
    ax1.bar(x, dbl_pass, bottom=singlet, color="#dd6b20", label="pass QC, predicted doublet")
    ax1.bar(x, fail, bottom=singlet + dbl_pass, color="#a0aec0", label="fail QC")
    for i in range(len(x)):
        ax1.text(x[i], singlet[i] / 2, f"{singlet[i]:.1f}%\n{int(summ['n_pass_qc_singlet'].iloc[i]):,}",
                 ha="center", va="center", color="white", fontsize=8)
    ax1.set_xticks(x, labels, fontsize=9)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("% of EmptyDrops nuclei")
    ax1.set_title("QC (scuttle) + scDblFinder outcome")
    ax1.legend(loc="lower left", fontsize=8, frameon=True)

    data, pos, colors = [], [], []
    for i, s in enumerate(samples):
        fl = pd.read_csv(f"{args.qc_dir}/{s}_qc_doublet_flags.csv",
                         usecols=["qc_nUMI", "pass_qc", "predicted_doublet"])
        fl = fl[fl["pass_qc"]]
        for j, (sel, col) in enumerate(((~fl["predicted_doublet"], "#2b6cb0"),
                                        (fl["predicted_doublet"], "#dd6b20"))):
            data.append(np.log10(fl.loc[sel, "qc_nUMI"].clip(lower=1)))
            pos.append(i * 3 + j)
            colors.append(col)
    bp = ax2.boxplot(data, positions=pos, widths=0.8, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    for m in bp["medians"]:
        m.set_color("black")
    ax2.set_xticks([i * 3 + 0.5 for i in range(len(samples))], labels, fontsize=9)
    ax2.set_ylabel("log10 nUMI (pass-QC nuclei)")
    ax2.set_title("nUMI: singlet (blue) vs predicted doublet (orange)")
    ax2.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    png = f"{args.qc_dir}/qc_doublet_summary.png"
    fig.savefig(png, dpi=140)
    print(f"Saved: {png}\nSaved: {args.qc_dir}/qc_doublet_summary_all_samples.csv")


if __name__ == "__main__":
    main()
