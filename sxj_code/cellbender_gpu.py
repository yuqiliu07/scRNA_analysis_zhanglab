"""CellBender remove-background for the 8 mouse samples (workflow sec. 1.2 / 3.3-1), CUDA GPU.

In this batch every partition carries ~240-400 ambient UMIs: the barcode-rank curve is a long empty-droplet
plateau followed by a cliff to ~0 (00_barcode_rank.py). CellBender's heuristic priors fail on that shape (it
put the empty prior at ~3 UMI, i.e. in the post-cliff tail, and then found no empty droplets to train on), so
the priors and the droplet window are set from this run's own tables, consistently with EmptyDrops:

  --force-empty-umi-prior    = plateau median UMI (plateau.csv; ranks 0.6-1.0 x n_partitions)
  --force-cell-umi-prior     = median UMI of EmptyDrops-called barcodes (<S>_summary.csv, median_umi_ed)
  --expected-cells           = nuclei passing basic QC (02_per_sample_qc.py, n_pass_basic)
  --total-droplets-included  = 0.6 x n_partitions, so the plateau (the same droplets EmptyDrops used as
                               ambient) lies beyond it and serves as CellBender's "surely empty" set
  --low-count-threshold      = 50 (the cliff, CLIFF_UMI in 00_barcode_rank.py): post-cliff tail excluded

Check each sample's learning curve and cell-probability histogram in <S>_cellbender_report.html before use.

Epochs: A-C were run with 150. Their test ELBO was best at epochs 45-75 and did not improve afterwards, so D-H
use --epochs 60: a full OneCycleLR cycle over 60 epochs, not a 150-epoch run stopped early (epochs is part of
CellBender's checkpoint hash, so a 150-epoch checkpoint is never resumed by a 60-epoch run). Whether 60 and 150
epochs give equivalent output is still to be checked (e.g. rerun A with --epochs 60 --out-subdir cellbender_e60).

Setup: Python 3.10 venv with torch (CUDA) and cellbender==0.3.2.
Usage: python cellbender_gpu.py <data_dir> <work_dir> [sample ...] [--epochs N] [--out-subdir NAME]
Output: <work_dir>/<out-subdir>/<S>/<S>_cellbender.h5, *_filtered.h5, *_cell_barcodes.csv, *_report.html,
        <work_dir>/<out-subdir>/settings.csv (updated after each sample)
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

SAMPLES = list("ABCDEFGH")
LOW_COUNT = 50
WINDOW_FRAC = 0.6
EPOCHS = 150
FPR = 0.01


def settings(work_dir, s, epochs):
    p = pd.read_csv(work_dir / "emptydrops" / "plateau.csv").set_index("sample").loc[s]
    ed = pd.read_csv(work_dir / "emptydrops" / f"{s}_summary.csv").iloc[0]
    thr = pd.read_csv(work_dir / "per_sample" / f"{s}_thresholds.csv").iloc[0]
    return dict(sample=s,
                expected_cells=int(thr["n_pass_basic"]),
                total_droplets=int(round(WINDOW_FRAC * p["n_partitions"])),
                empty_umi_prior=float(p["plateau_median"]),
                cell_umi_prior=float(ed["median_umi_ed"]),
                low_count_threshold=LOW_COUNT, epochs=epochs, fpr=FPR)


def run(data_dir, work_dir, s, exe, epochs, subdir):
    st = settings(work_dir, s, epochs)
    out = work_dir / subdir / s
    out.mkdir(parents=True, exist_ok=True)
    cmd = exe + ["remove-background", "--cuda",
           "--input", str((data_dir / s / "raw_feature_bc_matrix.h5").resolve()),
           "--output", f"{s}_cellbender.h5",
           "--expected-cells", str(st["expected_cells"]),
           "--total-droplets-included", str(st["total_droplets"]),
           "--force-empty-umi-prior", str(st["empty_umi_prior"]),
           "--force-cell-umi-prior", str(st["cell_umi_prior"]),
           "--low-count-threshold", str(LOW_COUNT),
           "--fpr", str(FPR), "--epochs", str(epochs), "--checkpoint-mins", "10"]
    print(f"=== {s}: {st}", flush=True)
    with open(out / f"{s}_run.log", "w") as log:
        rc = subprocess.call(cmd, cwd=out, stdout=log, stderr=subprocess.STDOUT)
    st["exit_code"] = rc
    print(f"    exit {rc}", flush=True)
    return st


def save_settings(f, row):
    old = pd.read_csv(f) if f.exists() else pd.DataFrame()
    if len(old):
        old = old[old["sample"] != row["sample"]]
    pd.concat([old, pd.DataFrame([row])]).sort_values("sample").to_csv(f, index=False)


def main(data_dir, work_dir, samples, epochs=EPOCHS, subdir="cellbender"):
    data_dir, work_dir = Path(data_dir), Path(work_dir)
    # cellbender through _cellbender_run.py, which fixes checkpoint saving under torch 2.x (see that file)
    exe = [sys.executable, str(Path(__file__).with_name("_cellbender_run.py"))]
    for s in samples:
        save_settings(work_dir / subdir / "settings.csv", run(data_dir, work_dir, s, exe, epochs, subdir))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir")
    ap.add_argument("work_dir")
    ap.add_argument("samples", nargs="*", default=SAMPLES)
    ap.add_argument("--epochs", type=int, default=EPOCHS)
    ap.add_argument("--out-subdir", default="cellbender")
    a = ap.parse_args()
    main(a.data_dir, a.work_dir, a.samples or SAMPLES, a.epochs, a.out_subdir)
