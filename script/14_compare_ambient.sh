#!/bin/bash
#SBATCH --job-name=ambient_cmp
#SBATCH --account=pi_hzh
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/ambient_cmp.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/ambient_cmp.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 比较 SoupX（第 12 步）与 CellBender（第 13 步）的 ambient RNA 去除，以及 CellBender 与 EmptyDrops 的细胞判定
# 结果在 work_dir/ambient_compare_mm10/；图片 → plots/ambient_compare，CSV → results/ambient_compare
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/ambient_compare_mm10}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}
Rscript $HOME/scrna_project/script/14_compare_ambient.R

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传；结果在 $OUT_DIR"; exit 0; fi
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/ambient_compare --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/ambient_compare --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/ambient_compare --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/ambient_compare --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/ambient_compare；CSV → $GDRIVE_ROOT/results/ambient_compare"
