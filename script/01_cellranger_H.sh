#!/bin/bash
#SBATCH --job-name=cellranger_H
#SBATCH --account=pi_hzh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellranger_H.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellranger_H.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# Cell Ranger count：样本 H
# FASTQ 前缀由文件名 H_S*_L008_R{1,2}_001.fastq.gz 识别为 "H"
set -euo pipefail

# 确保非登录 shell 中 module 命令可用
. /u/local/Modules/default/init/modules.sh
module purge
module load cellranger

REF=$SCRATCH/FSU_seq/reference/refdata-gex-GRCm39-2024-A
if [ ! -d "$REF" ]; then
    echo "错误：参考基因组不存在 $REF，请先运行 00_download_ref.sh" >&2
    exit 1
fi

mkdir -p $SCRATCH/FSU_seq/work_dir/logs
mkdir -p $SCRATCH/FSU_seq/02_cellranger
cd $SCRATCH/FSU_seq/02_cellranger

cellranger count \
  --id=H \
  --transcriptome=$REF \
  --fastqs=$SCRATCH/FSU_seq/01_rawdata/H/ \
  --sample=H \
  --include-introns true \
  --create-bam false \
  --localcores 16 \
  --localmem 64
