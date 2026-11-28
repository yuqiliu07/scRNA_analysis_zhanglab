#!/bin/bash
#SBATCH --job-name=rescue_5c
#SBATCH --account=pi_hzh
#SBATCH --partition=campus
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/rescue_5c.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/rescue_5c.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 第 17 步（step 5C）：双细胞 cluster 中单谱系核的 KNN 救回 + 淋巴细胞识别 + 二次审阅，输出最终注释对象
# 参数见 17_rescue_doublets.R 文件头；可用 --export=ALL,DOUBLET_CLUSTERS=...,LOWQ_CLUSTERS=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

R_SCRIPT=$HOME/scrna_project/script/17_rescue_doublets.R
TABLE_DIR=${TABLE_DIR:-$SCRATCH/FSU_seq/results/tables}
FIG_DIR=${FIG_DIR:-$SCRATCH/FSU_seq/results/figures}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}
export TABLE_DIR FIG_DIR
Rscript $R_SCRIPT

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传"; exit 0; fi
# 只传小表和 PNG；rds、逐核表和 PDF 留在 Hoffman
T="--include rescue_summary.csv --include new_cluster_check.csv --include final_celltype_counts.csv"
P="--include rescue_umap_celltype.png --include rescue_umap_status.png"
rclone copy  $TABLE_DIR $GDRIVE_ROOT/results/rescue $T
rclone check $TABLE_DIR $GDRIVE_ROOT/results/rescue --one-way $T
rclone copy  $FIG_DIR   $GDRIVE_ROOT/plots/rescue $P
rclone check $FIG_DIR   $GDRIVE_ROOT/plots/rescue --one-way $P
echo "表 → $GDRIVE_ROOT/results/rescue；图 → $GDRIVE_ROOT/plots/rescue"
