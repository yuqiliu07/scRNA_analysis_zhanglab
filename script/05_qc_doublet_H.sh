#!/bin/bash
#SBATCH --job-name=qc_doublet_H
#SBATCH --account=pi_hzh
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/qc_doublet_H.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/qc_doublet_H.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# QC flags（scuttle perCellQCFilters）+ scDblFinder：样本 H
#   nUMI / nGene: median ± 5 MAD (log)；nGene ≥ 200；pct_mt ≤ median + 3 MAD 且 ≤ 5%
# 可用 sbatch --export=ALL,CR_DIR=...,ED_DIR=...,OUT_DIR=...,DBR_PER1K=... 覆盖；默认 = Novogene mm10 数据
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export CR_DIR=${CR_DIR:-$SCRATCH/FSU_seq/02_cellranger/shared_h5}
export ED_DIR=${ED_DIR:-$SCRATCH/FSU_seq/work_dir/emptydrops_mm10}
export OUT_DIR=${OUT_DIR:-$SCRATCH/FSU_seq/work_dir/qc_mm10}
export DBR_PER1K=${DBR_PER1K:-0.004}

Rscript $HOME/scrna_project/script/05_qc_doublet.R H ${SLURM_CPUS_PER_TASK:-8}
