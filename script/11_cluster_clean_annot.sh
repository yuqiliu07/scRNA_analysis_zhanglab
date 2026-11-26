#!/bin/bash
#SBATCH --job-name=cluster_annot
#SBATCH --account=pi_hzh
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cluster_annot.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cluster_annot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 清理 Leiden cluster（< 10 个核的 cluster 按 harmony 近邻多数票并入），
# UMAP（cluster / doublet score / pct_mt）、各 cluster nGene 中位数柱状图，
# 原始代码（sxj_code/04_merge_cluster.py PANELS）marker panel 的 cluster × 细胞类型热图
# 可用 sbatch --export=ALL,HARMONY_DIR=...,OUT_DIR=...,MIN_CLUSTER_SIZE=...,GDRIVE_ROOT=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export HARMONY_DIR=${HARMONY_DIR:-$SCRATCH/FSU_seq/work_dir/harmony_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/annot_mm10}
export MIN_CLUSTER_SIZE=${MIN_CLUSTER_SIZE:-10}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript $HOME/scrna_project/script/11_cluster_clean_annot.R

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传；结果在 $OUT_DIR"; exit 0; fi
# 只上传图片和汇总 CSV（.rds 和逐细胞 csv.gz 留在 Hoffman）
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/annotation --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/annotation --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/annotation --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/annotation --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/annotation；CSV → $GDRIVE_ROOT/results/annotation"
