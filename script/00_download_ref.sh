#!/bin/bash
#SBATCH --job-name=download_ref
#SBATCH --account=pi_hzh
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/download_ref.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/download_ref.%j.err
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

# 下载 10x 小鼠参考基因组 refdata-gex-GRCm39-2024-A（约 10.3 GB）
set -euo pipefail

REF_DIR=$SCRATCH/FSU_seq/reference
REF_NAME=refdata-gex-GRCm39-2024-A

mkdir -p "$REF_DIR"
cd "$REF_DIR"

if [ -d "$REF_NAME" ]; then
    echo "参考基因组已存在，跳过下载：$REF_DIR/$REF_NAME"
    exit 0
fi

# -c 支持断点续传
wget -c https://cf.10xgenomics.com/supp/cell-exp/${REF_NAME}.tar.gz
tar -xzf ${REF_NAME}.tar.gz
rm ${REF_NAME}.tar.gz
echo "参考基因组下载完成：$REF_DIR/$REF_NAME"
