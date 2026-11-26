#!/bin/bash
#SBATCH --job-name=soupx
#SBATCH --account=pi_hzh
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/soupx.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/soupx.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# SoupX：8 个样本依次运行
#   soup profile = 原始 h5 中 UMI <= plateau_q75 的 barcode；cluster = 第 11 步清理后的 leiden_clean
#   rho：autoEstCont；失败时用心肌基因（CM panel）在非心肌 cluster 核中的计数作 fallback
#   每个样本的 rho 打印在 .out 日志末尾，并写入 soupx_rho.csv；分布图 soupx_rho_*.png
# 可用 sbatch --export=ALL,RAW_DIR=...,PLATEAU_CSV=...,ANNOT_DIR=...,OUT_DIR=...,CM_SCORE_MAX=...,GDRIVE_ROOT=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export RAW_DIR=${RAW_DIR:-$SCRATCH/FSU_seq/02_cellranger/shared_h5}
export PLATEAU_CSV=${PLATEAU_CSV:-$SCRATCH/FSU_seq/work_dir/emptydrops_mm10/plateau.csv}
export ANNOT_DIR=${ANNOT_DIR:-$SCRATCH/FSU_seq/work_dir/annot_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/soupx_mm10}
export CM_SCORE_MAX=${CM_SCORE_MAX:-1.8}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript -e 'if (!requireNamespace("SoupX", quietly = TRUE)) install.packages("SoupX", repos = "https://cloud.r-project.org")'
Rscript $HOME/scrna_project/script/12_soupx.R

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传；结果在 $OUT_DIR"; exit 0; fi
# 只上传图片和汇总 CSV（校正后的计数 .rds 留在 Hoffman）
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/soupx --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/soupx --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/soupx --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/soupx --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/soupx；CSV → $GDRIVE_ROOT/results/soupx"
