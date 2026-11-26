#!/bin/bash
#SBATCH --job-name=harmony_cluster
#SBATCH --account=pi_hzh
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/harmony_cluster.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/harmony_cluster.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# Harmony（按样本，max_iter 20）→ UMAP（harmony 1-30，n_neighbors 20）
# → FindNeighbors(k 20) → Leiden（igraph，modularity，resolution 1.0，n_iterations 2）
# PCA 计算 50 个 PC，使用前 30 个；同时生成整合前的 UMAP（PC 1-30）用于对比
# 可用 sbatch --export=ALL,N_PCS=...,RESOLUTION=...,OUT_DIR=...,GDRIVE_ROOT=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

# 首次运行时安装 harmony（CRAN）到用户库；已安装则跳过
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-4}"
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
if (!requireNamespace("harmony", quietly = TRUE))
  install.packages("harmony", lib = lib, repos = "https://cloud.r-project.org")
cat("harmony", as.character(packageVersion("harmony")), "\n")'

export SEURAT_DIR=${SEURAT_DIR:-$SCRATCH/FSU_seq/work_dir/seurat_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/harmony_mm10}
export N_PCS_COMPUTE=${N_PCS_COMPUTE:-50}
export N_PCS=${N_PCS:-30}
export MAX_ITER_HARMONY=${MAX_ITER_HARMONY:-20}
export N_NEIGHBORS=${N_NEIGHBORS:-20}
export RESOLUTION=${RESOLUTION:-1.0}
export LEIDEN_ITER=${LEIDEN_ITER:-2}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript $HOME/scrna_project/script/09_harmony_cluster.R

if [ "$GDRIVE_ROOT" = "none" ]; then
    echo "GDRIVE_ROOT=none，跳过上传；结果在 $OUT_DIR"
    exit 0
fi
# 只上传图片和汇总 CSV（.rds 和逐细胞 csv.gz 是数据文件，留在 Hoffman）
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/harmony --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/harmony --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/harmony --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/harmony --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/harmony；CSV → $GDRIVE_ROOT/results/harmony"
