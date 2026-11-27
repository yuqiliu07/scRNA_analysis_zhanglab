#!/bin/bash
#SBATCH --job-name=cellbender_B
#SBATCH --account=pi_hzh
#SBATCH --partition=campus
#SBATCH --gres=gpu:1
#SBATCH --constraint="V100|A100|RTX2080Ti"
#SBATCH --exclude=g12
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellbender_B.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellbender_B.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# CellBender remove-background：样本 B（由 13_make_cellbender_scripts.py 生成，请勿手改，改参数后重新生成）
#   expected-cells          = 49986  (EmptyDrops 细胞数)
#   total-droplets-included = 87300  (plateau_rank_start 82300 + 5000，plateau_rank_end = 137165)
#   epochs = 150；learning-rate = 1e-4；RTX2080Ti 上 60 epochs 约 55 分钟（150 epochs 约 2.3 小时）；每 10 分钟写 checkpoint，超时后重新 sbatch 本脚本会从 ckpt.tar.gz 接着跑
# GPU 限定 V100 / A100 / RTX2080Ti（torch 1.13.1+cu117 不支持 H200、Blackwell 和 P4）；排除 g12（GPU 故障，第 1 个 epoch 即 CUDA 报错）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load python/3.10.20
source $HOME/envs/cellbender_t113/bin/activate

OUT=/u/scratch/y/yuqi/FSU_seq/work_dir/cellbender_mm10/B
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}
mkdir -p $OUT
cd $OUT                      # CellBender 在当前目录找 ckpt.tar.gz 续跑
nvidia-smi --query-gpu=name,memory.total --format=csv

cellbender remove-background --cuda \
  --input /u/scratch/y/yuqi/FSU_seq/02_cellranger/shared_h5/B/outs/raw_feature_bc_matrix.h5 \
  --output $OUT/B_cellbender.h5 \
  --expected-cells 49986 \
  --total-droplets-included 87300 \
  --epochs 150 \
  --learning-rate 1e-4 \
  --fpr 0.01 \
  --checkpoint-mins 10

if [ "$GDRIVE_ROOT" = "none" ]; then echo "[B] 跳过上传；结果在 $OUT"; exit 0; fi
# 只上传报告 PDF 和 metrics CSV（h5 计数矩阵留在 Hoffman）
rclone copy $OUT $GDRIVE_ROOT/plots/cellbender --include "B_cellbender.pdf"
rclone check $OUT $GDRIVE_ROOT/plots/cellbender --one-way --include "B_cellbender.pdf"
rclone copy $OUT $GDRIVE_ROOT/results/cellbender --include "B_cellbender_metrics.csv"
rclone check $OUT $GDRIVE_ROOT/results/cellbender --one-way --include "B_cellbender_metrics.csv"
echo "[B] 完成"
