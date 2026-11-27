#!/bin/bash
#SBATCH --job-name=cb_harmony
#SBATCH --account=pi_hzh
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cb_harmony.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cb_harmony.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 用 CellBender 校正后的 counts 重做第 07/09/11 步：
#   核 = pass_qc_singlet ∩ CellBender 判为细胞 ∩ 校正后 nGene >= 200
#   → 基因过滤（>= 10 核）→ LogNormalize → 每个样本 HVG(vst, 3000) 取并集
#   → PCA(50，用 1-30) → Harmony（按样本）→ UMAP + SNN(k 20) → Leiden(resolution 1.0)
#   → 小 cluster(<10) 并入近邻 → 17 个 marker panel 打分 → 与第 11 步（原始 counts）的 cluster 对比
# 可用 sbatch --export=ALL,MIN_GENES=...,N_PCS=...,RESOLUTION=...,OUT_DIR=...,GDRIVE_ROOT=... 覆盖
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/harmony_cb_mm10}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript $HOME/scrna_project/script/15_cellbender_harmony.R

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传；结果在 $OUT_DIR"; exit 0; fi
# 只上传图片和汇总 CSV（.rds 和逐细胞 csv.gz 是数据文件，留在 Hoffman）
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/harmony_cellbender --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/harmony_cellbender --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/harmony_cellbender --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/harmony_cellbender --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/harmony_cellbender；CSV → $GDRIVE_ROOT/results/harmony_cellbender"
