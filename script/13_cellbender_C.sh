#!/bin/bash
#SBATCH --job-name=cellbender_C
#SBATCH --account=pi_hzh
#SBATCH --partition=campus
#SBATCH --gres=gpu:1
#SBATCH --constraint="V100|A100|H200|RTX2080Ti"
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellbender_C.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellbender_C.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# CellBender remove-background：样本 C（由 13_make_cellbender_scripts.py 生成，请勿手改，改参数后重新生成）
#   expected-cells          = 32610  (EmptyDrops 细胞数)
#   total-droplets-included = 86821  (plateau_rank_start 81821 + 5000，plateau_rank_end = 136367)
#   epochs = 60；每 10 分钟写 checkpoint，超时后重新 sbatch 本脚本会从 ckpt.tar.gz 接着跑
# GPU 限定 V100 / A100 / H200 / RTX2080Ti（torch 2.5.1+cu124 不支持 Blackwell 和 P4）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load python/3.10.20
source $HOME/envs/cellbender/bin/activate

OUT=/u/scratch/y/yuqi/FSU_seq/work_dir/cellbender_mm10/C
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}
mkdir -p $OUT
cd $OUT                      # CellBender 在当前目录找 ckpt.tar.gz 续跑
nvidia-smi --query-gpu=name,memory.total --format=csv

cellbender remove-background --cuda \
  --input /u/scratch/y/yuqi/FSU_seq/02_cellranger/shared_h5/C/outs/raw_feature_bc_matrix.h5 \
  --output $OUT/C_cellbender.h5 \
  --expected-cells 32610 \
  --total-droplets-included 86821 \
  --epochs 60 \
  --fpr 0.01 \
  --checkpoint-mins 10

if [ "$GDRIVE_ROOT" = "none" ]; then echo "[C] 跳过上传；结果在 $OUT"; exit 0; fi
# 只上传报告 PDF 和 metrics CSV（h5 计数矩阵留在 Hoffman）
rclone copy $OUT $GDRIVE_ROOT/plots/cellbender --include "C_cellbender.pdf"
rclone check $OUT $GDRIVE_ROOT/plots/cellbender --one-way --include "C_cellbender.pdf"
rclone copy $OUT $GDRIVE_ROOT/results/cellbender --include "C_cellbender_metrics.csv"
rclone check $OUT $GDRIVE_ROOT/results/cellbender --one-way --include "C_cellbender_metrics.csv"
echo "[C] 完成"
