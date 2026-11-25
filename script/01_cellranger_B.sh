#!/bin/bash
#SBATCH --job-name=cellranger_B
#SBATCH --account=pi_hzh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellranger_B.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/cellranger_B.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# ---------------------------------------------------------------------------
# 参考基因组 / Cell Ranger 版本说明（与 Novogene 交付结果的区别）
#   本流程：Cell Ranger 9.0.0 + refdata-gex-GRCm39-2024-A
#           (GRCm39, GENCODE vM33 / Ensembl 110)
#   Novogene：Cell Ranger 8.0.1 + refdata-gex-mm10-2020-A
#           (GRCm38/mm10, GENCODE vM23 / Ensembl 98)
#   两者均为 3' v4 化学、include-introns=True。参考基因组与软件版本不同，
#   比对率、细胞数、barcode 列表不能与 Novogene 结果直接逐一对应。
# ---------------------------------------------------------------------------

# Cell Ranger count：样本 B
# FASTQ 前缀由文件名 B_S*_L008_R{1,2}_001.fastq.gz 识别为 "B"
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
  --id=B \
  --transcriptome=$REF \
  --fastqs=$SCRATCH/FSU_seq/01_rawdata/B/ \
  --sample=B \
  --include-introns true \
  --create-bam false \
  --localcores 16 \
  --localmem 64
