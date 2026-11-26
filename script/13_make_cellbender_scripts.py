#!/usr/bin/env python3
"""
Generate one CellBender (remove-background, GPU) Slurm script per FSU_seq sample.

Parameters per sample:
  --expected-cells           = emptydrops_cells   (<ED_DIR>/<s>_emptydrops_summary.csv)
  --total-droplets-included  = plateau_rank_start + PLATEAU_PAD (5000), capped at plateau_rank_end
                               (<ED_DIR>/plateau.csv; the plateau starts at 0.6 x n_partitions,
                               so this reaches a few thousand droplets into the empty plateau,
                               as CellBender recommends)
  --epochs 60, --fpr 0.01, --checkpoint-mins 10
  --learning-rate 1e-4 (CellBender default); per-sample override with --lr, e.g. --lr D=5e-5
  (D hit "NaN in mu" at the start of training with 1e-4; CellBender suggests halving the rate)

Environment: $HOME/envs/cellbender_t113 (CellBender 0.3.2, torch 1.13.1+cu117, pyro-ppl 1.8.6).
With torch 2.5 + pyro 1.9 the checkpoint cannot be saved ("cannot pickle weakref") and the
posterior step then fails; see 13_install_cellbender.sh.

Writes (next to this script):
  13_cellbender_<s>.sh            one GPU job per sample
  13_submit_all_cellbender.sh     sbatch all 8
and <OUT_ROOT>/cellbender_params.csv (the parameters used).

Run on Hoffman2 (login node is fine, standard library only):
  python3 13_make_cellbender_scripts.py
Options: --ed-dir, --out-root, --raw-dir, --epochs, --pad, --time, --lr S=RATE (repeatable)
"""
import argparse
import csv
import os

SCRATCH = os.environ.get("SCRATCH", "/u/scratch/y/yuqi")
HERE = os.path.dirname(os.path.abspath(__file__))

TEMPLATE = """#!/bin/bash
#SBATCH --job-name=cellbender_{s}
#SBATCH --account=pi_hzh
#SBATCH --partition=campus
#SBATCH --gres=gpu:1
#SBATCH --constraint="V100|A100|RTX2080Ti"
#SBATCH --time={time}
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output={scratch}/FSU_seq/work_dir/logs/cellbender_{s}.%j.out
#SBATCH --error={scratch}/FSU_seq/work_dir/logs/cellbender_{s}.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# CellBender remove-background：样本 {s}（由 13_make_cellbender_scripts.py 生成，请勿手改，改参数后重新生成）
#   expected-cells          = {exp_cells}  (EmptyDrops 细胞数)
#   total-droplets-included = {total}  (plateau_rank_start {pstart} + {pad}，plateau_rank_end = {pend})
#   epochs = {epochs}；learning-rate = {lr}；每 10 分钟写 checkpoint，超时后重新 sbatch 本脚本会从 ckpt.tar.gz 接着跑
# GPU 限定 V100 / A100 / RTX2080Ti（torch 1.13.1+cu117 不支持 H200、Blackwell 和 P4）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load python/3.10.20
source $HOME/envs/cellbender_t113/bin/activate

OUT={out_root}/{s}
GDRIVE_ROOT=${{GDRIVE_ROOT:-gdrive:FSU_seq_results}}
mkdir -p $OUT
cd $OUT                      # CellBender 在当前目录找 ckpt.tar.gz 续跑
nvidia-smi --query-gpu=name,memory.total --format=csv

cellbender remove-background --cuda \\
  --input {raw_dir}/{s}/outs/raw_feature_bc_matrix.h5 \\
  --output $OUT/{s}_cellbender.h5 \\
  --expected-cells {exp_cells} \\
  --total-droplets-included {total} \\
  --epochs {epochs} \\
  --learning-rate {lr} \\
  --fpr 0.01 \\
  --checkpoint-mins 10

if [ "$GDRIVE_ROOT" = "none" ]; then echo "[{s}] 跳过上传；结果在 $OUT"; exit 0; fi
# 只上传报告 PDF 和 metrics CSV（h5 计数矩阵留在 Hoffman）
rclone copy $OUT $GDRIVE_ROOT/plots/cellbender --include "{s}_cellbender.pdf"
rclone check $OUT $GDRIVE_ROOT/plots/cellbender --one-way --include "{s}_cellbender.pdf"
rclone copy $OUT $GDRIVE_ROOT/results/cellbender --include "{s}_cellbender_metrics.csv"
rclone check $OUT $GDRIVE_ROOT/results/cellbender --one-way --include "{s}_cellbender_metrics.csv"
echo "[{s}] 完成"
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ed-dir", default=f"{SCRATCH}/FSU_seq/work_dir/emptydrops_mm10")
    ap.add_argument("--out-root", default=f"{SCRATCH}/FSU_seq/work_dir/cellbender_mm10")
    ap.add_argument("--raw-dir", default=f"{SCRATCH}/FSU_seq/02_cellranger/shared_h5")
    ap.add_argument("--script-dir", default=HERE)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--pad", type=int, default=5000)
    ap.add_argument("--time", default="02:00:00")
    ap.add_argument("--lr", action="append", default=["D=5e-5"],
                    help="per-sample learning rate, S=RATE (default D=5e-5; others 1e-4)")
    a = ap.parse_args()
    lr_over = dict(x.split("=", 1) for x in a.lr)

    with open(os.path.join(a.ed_dir, "plateau.csv")) as f:
        plateau = {r["sample"]: r for r in csv.DictReader(f)}

    rows, names = [], []
    for s in sorted(plateau):
        with open(os.path.join(a.ed_dir, f"{s}_emptydrops_summary.csv")) as f:
            ed = next(csv.DictReader(f))
        p = plateau[s]
        exp_cells = int(ed["emptydrops_cells"])
        pstart, pend = int(p["plateau_rank_start"]), int(p["plateau_rank_end"])
        total = min(pstart + a.pad, pend)
        lr = lr_over.get(s, "1e-4")
        if total <= exp_cells:
            raise SystemExit(f"{s}: total_droplets {total} <= expected_cells {exp_cells}")
        name = f"13_cellbender_{s}.sh"
        with open(os.path.join(a.script_dir, name), "w") as f:
            f.write(TEMPLATE.format(s=s, time=a.time, scratch=SCRATCH, exp_cells=exp_cells, total=total,
                                    pstart=pstart, pend=pend, pad=a.pad, epochs=a.epochs, lr=lr,
                                    out_root=a.out_root, raw_dir=a.raw_dir))
        os.chmod(os.path.join(a.script_dir, name), 0o755)
        names.append(name)
        rows.append(dict(sample=s, expected_cells=exp_cells, total_droplets_included=total,
                         plateau_rank_start=pstart, plateau_rank_end=pend, plateau_q75=p["plateau_q75"],
                         cellranger_cells=ed["cellranger_cells"], epochs=a.epochs, learning_rate=lr, fpr=0.01))
        print(f"{s}: expected_cells={exp_cells:,}  total_droplets_included={total:,}  "
              f"(plateau ranks {pstart:,}-{pend:,})  lr={lr}  -> {name}")

    with open(os.path.join(a.script_dir, "13_submit_all_cellbender.sh"), "w") as f:
        f.write("#!/bin/bash\n# 提交 8 个 CellBender GPU 作业（由 13_make_cellbender_scripts.py 生成）\n"
                "for s in " + " ".join(names) + "; do\n"
                "  sbatch $HOME/scrna_project/script/$s\ndone\n")
    os.chmod(os.path.join(a.script_dir, "13_submit_all_cellbender.sh"), 0o755)

    os.makedirs(a.out_root, exist_ok=True)
    with open(os.path.join(a.out_root, "cellbender_params.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"params -> {a.out_root}/cellbender_params.csv")


if __name__ == "__main__":
    main()
