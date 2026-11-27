#!/usr/bin/env python3
"""
CellBender learning curves + convergence summary for every finished sample.

Reads <OUT_ROOT>/<s>/<s>_cellbender.log ("[epoch N] average training/test loss") and
<s>_cellbender_metrics.csv, then writes to <OUT_ROOT>:
  cellbender_learning_curves.png   one panel per sample: train loss (line) + test loss (points)
  cellbender_convergence.csv       per sample: final losses, plateau drift, metrics
Plateau check: |mean train loss of last 10 epochs - previous 10| / loss range < 1% -> "plateau".

Run on Hoffman2 login node with the CellBender env (has matplotlib):
  source $HOME/envs/cellbender_t113/bin/activate   # after: module load python/3.10.20
  python3 13_cellbender_curves.py [--upload]
--upload copies the PNG to plots/cellbender and the CSV to results/cellbender on Drive (rclone copy + check).
"""
import argparse, csv, os, re, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRATCH = os.environ.get("SCRATCH", "/u/scratch/y/yuqi")
ap = argparse.ArgumentParser()
ap.add_argument("--out-root", default=f"{SCRATCH}/FSU_seq/work_dir/cellbender_mm10")
ap.add_argument("--upload", action="store_true")
ap.add_argument("--gdrive", default=os.environ.get("GDRIVE_ROOT", "gdrive:FSU_seq_results"))
a = ap.parse_args()

pat = re.compile(r"\[epoch (\d+)\]\s+average (training|test) loss: ([-\d.]+)")
done = []
for s in "ABCDEFGH":
    log = os.path.join(a.out_root, s, f"{s}_cellbender.log")
    met = os.path.join(a.out_root, s, f"{s}_cellbender_metrics.csv")
    if not (os.path.exists(log) and os.path.exists(met)):
        continue
    cur = {"training": {}, "test": {}}; lr = "1e-4"
    for line in open(log):
        if "--learning-rate" in line:
            lr = line.split("--learning-rate")[1].split()[0]
        m = pat.search(line)
        if m:
            cur[m[2]][int(m[1])] = float(m[3])   # later runs (resume) overwrite earlier epochs
    metrics = {r[0]: r[1] for r in csv.reader(open(met))}
    done.append((s, cur, metrics, lr))
if not done:
    raise SystemExit("no finished CellBender sample found")

TRAIN, TEST = "#2a6fbb", "#d9822b"
n = len(done); ncol = min(4, n); nrow = (n + ncol - 1) // ncol
fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.9 * nrow), squeeze=False)
rows = []
for ax, (s, cur, met, lr) in zip(axes.flat, done):
    tr = sorted(cur["training"].items()); te = sorted(cur["test"].items())
    ax.plot([e for e, _ in tr], [v for _, v in tr], color=TRAIN, lw=2, label="train")
    ax.plot([e for e, _ in te], [v for _, v in te], "o", color=TEST, ms=5, label="test")
    vals = [v for _, v in tr]
    last, prev = vals[-10:], vals[-20:-10]
    rng = max(vals) - min(vals)
    drift = abs(sum(last) / len(last) - sum(prev) / len(prev)) / rng if prev and rng else float("nan")
    lo = min(vals[len(vals) // 3:] + [v for _, v in te]); hi = max(vals[len(vals) // 3:] + [v for _, v in te])
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.6 * (hi - lo))   # zoom on the plateau; early spikes clipped
    ax.set_title(f"{s}  lr={lr}  drift={drift:.1%}", fontsize=10)
    ax.set_xlabel("epoch"); ax.set_ylabel("loss (-ELBO)")
    ax.grid(color="#e6e6e6", lw=0.6); ax.spines[["top", "right"]].set_visible(False)
    rows.append(dict(sample=s, learning_rate=lr, epochs=tr[-1][0], train_loss_final=round(vals[-1], 1),
                     test_loss_final=round(te[-1][1], 1) if te else "",
                     train_drift_last10_vs_prev10=round(drift, 4),
                     plateau="yes" if drift < 0.01 else "no",
                     convergence_indicator=met.get("convergence_indicator"),
                     found_cells=met.get("found_cells"), expected_cells=met.get("expected_cells"),
                     fraction_counts_removed=met.get("fraction_counts_removed")))
for ax in list(axes.flat)[n:]:
    ax.axis("off")
axes.flat[0].legend(frameon=False, fontsize=9)
fig.suptitle("CellBender learning curves", fontsize=12)
fig.tight_layout()
png = os.path.join(a.out_root, "cellbender_learning_curves.png")
fig.savefig(png, dpi=130)
out = os.path.join(a.out_root, "cellbender_convergence.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
print("->", png, out)

if a.upload:
    for name, dest in [("cellbender_learning_curves.png", "plots/cellbender"), ("cellbender_convergence.csv", "results/cellbender")]:
        for cmd in ("copy", "check"):
            extra = ["--one-way"] if cmd == "check" else []
            subprocess.run(["rclone", cmd, a.out_root, f"{a.gdrive}/{dest}", "--include", name] + extra, check=True)
