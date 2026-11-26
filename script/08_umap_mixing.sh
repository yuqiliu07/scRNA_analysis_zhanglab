#!/bin/bash
#SBATCH --job-name=umap_mixing
#SBATCH --account=pi_hzh
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/umap_mixing.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/umap_mixing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 未整合（无 Harmony）的 PCA → UMAP（PC 1-20），按样本着色；
# 每个细胞 30 近邻中同样本比例 vs 该样本占比 = mixing score
# 可用 sbatch --export=ALL,SEURAT_DIR=...,OUT_DIR=...,N_PCS=...,K_NN=...,GDRIVE_ROOT=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export SEURAT_DIR=${SEURAT_DIR:-$SCRATCH/FSU_seq/work_dir/seurat_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/umap_mm10}
export N_PCS=${N_PCS:-20}
export K_NN=${K_NN:-30}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript $HOME/scrna_project/script/08_umap_mixing.R

if [ "$GDRIVE_ROOT" = "none" ]; then
    echo "GDRIVE_ROOT=none，跳过上传；结果在 $OUT_DIR"
    exit 0
fi
# 只上传图片和汇总 CSV（.rds 和逐细胞 csv.gz 是数据文件，留在 Hoffman）
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/umap_mixing --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/umap_mixing --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/umap_mixing --include "mixing_*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/umap_mixing --one-way --include "mixing_*.csv"
echo "图片 → $GDRIVE_ROOT/plots/umap_mixing；CSV → $GDRIVE_ROOT/results/umap_mixing"
