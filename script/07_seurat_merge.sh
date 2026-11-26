#!/bin/bash
#SBATCH --job-name=seurat_merge
#SBATCH --account=pi_hzh
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/seurat_merge.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/seurat_merge.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# Seurat v5：读取 8 个样本 → 保留 pass_qc_singlet → merge → 基因过滤（≥10 核）
#            → NormalizeData → 按样本 FindVariableFeatures(3000) 取并集
# 可用 sbatch --export=ALL,CR_DIR=...,QC_DIR=...,OUT_DIR=...,MIN_CELLS=...,N_HVG=...,GDRIVE_DIR=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export CR_DIR=${CR_DIR:-$SCRATCH/FSU_seq/02_cellranger/shared_h5}
export QC_DIR=${QC_DIR:-$SCRATCH/FSU_seq/work_dir/qc_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/seurat_mm10}
export MIN_CELLS=${MIN_CELLS:-10}
export N_HVG=${N_HVG:-3000}
GDRIVE_DIR=${GDRIVE_DIR:-gdrive:FSU_seq_results/results/seurat}

Rscript $HOME/scrna_project/script/07_seurat_merge.R

if [ "$GDRIVE_DIR" = "none" ]; then
    echo "GDRIVE_DIR=none，跳过上传；结果在 $OUT_DIR"
    exit 0
fi
# 只上传汇总 CSV（不上传 .rds 数据文件）
rclone copy $OUT_DIR $GDRIVE_DIR --include "seurat_*.csv"
rclone check $OUT_DIR $GDRIVE_DIR --one-way --include "seurat_*.csv"
echo "CSV 已上传到 $GDRIVE_DIR"
