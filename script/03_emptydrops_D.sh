#!/bin/bash
#SBATCH --job-name=emptydrops_D
#SBATCH --account=pi_hzh
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/emptydrops_D.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/emptydrops_D.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# EmptyDrops：样本 D（lower = plateau_q75，niters = 10000，FDR < 0.001）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

SAMPLE=D
# 可用 sbatch --export=ALL,CR_DIR=...,OUT_DIR=...,GDRIVE_DIR=... 覆盖；默认值 = 原 GRCm39 流程
#   CR_DIR      Cell Ranger 输出根目录（读取 <CR_DIR>/<sample>/outs/）
#   OUT_DIR     plateau.csv 所在目录，也是结果输出目录
#   GDRIVE_DIR  结果上传位置；设为 none 则不上传
export CR_DIR=${CR_DIR:-$SCRATCH/FSU_seq/02_cellranger}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/emptydrops}
GDRIVE_DIR=${GDRIVE_DIR:-gdrive:FSU_seq_results/results/emptydrops}

Rscript $HOME/scrna_project/script/03_emptydrops.R $SAMPLE ${SLURM_CPUS_PER_TASK:-8}

if [ "$GDRIVE_DIR" = "none" ]; then
    echo "[$SAMPLE] GDRIVE_DIR=none，跳过上传；结果在 $OUT_DIR"
    exit 0
fi

# 结果 CSV 推到 Google Drive（只传 CSV，不传数据和代码）
rclone copy $OUT_DIR $GDRIVE_DIR \
  --include "${SAMPLE}_emptydrops_full.csv" \
  --include "${SAMPLE}_emptydrops_summary.csv"
rclone check $OUT_DIR $GDRIVE_DIR --one-way \
  --include "${SAMPLE}_emptydrops_*.csv"
echo "[$SAMPLE] CSV 已上传到 $GDRIVE_DIR"
