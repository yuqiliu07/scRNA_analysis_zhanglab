#!/bin/bash
#SBATCH --job-name=compare_novogene
#SBATCH --account=pi_hzh
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/compare_novogene.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/compare_novogene.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 我们的 UMAP/cluster（09）vs Novogene 交付的 UMAP/cluster（压缩包 3.Identification_subpopulation）
# 需要先把 Novogene 的 UMAP/Clustering 文件解压到 $NOVO_DIR（已在 $SCRATCH/compare/novo_umap）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export HARMONY_DIR=${HARMONY_DIR:-$SCRATCH/FSU_seq/work_dir/harmony_mm10}
export NOVO_DIR=${NOVO_DIR:-$SCRATCH/compare/novo_umap}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/compare_novogene}
export K_NN=${K_NN:-30}
GDRIVE_ROOT=${GDRIVE_ROOT:-gdrive:FSU_seq_results}

Rscript $HOME/scrna_project/script/10_compare_novogene_umap.R

if [ "$GDRIVE_ROOT" = "none" ]; then echo "跳过上传；结果在 $OUT_DIR"; exit 0; fi
rclone copy $OUT_DIR $GDRIVE_ROOT/plots/compare_novogene --include "*.png"
rclone check $OUT_DIR $GDRIVE_ROOT/plots/compare_novogene --one-way --include "*.png"
rclone copy $OUT_DIR $GDRIVE_ROOT/results/compare_novogene --include "*.csv"
rclone check $OUT_DIR $GDRIVE_ROOT/results/compare_novogene --one-way --include "*.csv"
echo "图片 → $GDRIVE_ROOT/plots/compare_novogene；CSV → $GDRIVE_ROOT/results/compare_novogene"
