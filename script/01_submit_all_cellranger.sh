#!/bin/bash

# ---------------------------------------------------------------------------
# 参考基因组 / Cell Ranger 版本说明（与 Novogene 交付结果的区别）
#   本流程：Cell Ranger 9.0.0 + refdata-gex-GRCm39-2024-A
#           (GRCm39, GENCODE vM33 / Ensembl 110)
#   Novogene：Cell Ranger 8.0.1 + refdata-gex-mm10-2020-A
#           (GRCm38/mm10, GENCODE vM23 / Ensembl 98)
#   两者均为 3' v4 化学、include-introns=True。参考基因组与软件版本不同，
#   比对率、细胞数、barcode 列表不能与 Novogene 结果直接逐一对应。
# ---------------------------------------------------------------------------
# 一键提交 8 个样本的 Cell Ranger 作业
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REF=$SCRATCH/FSU_seq/reference/refdata-gex-GRCm39-2024-A

if [ ! -d "$REF" ]; then
    echo "错误：参考基因组不存在 $REF，请先 sbatch $SCRIPT_DIR/00_download_ref.sh 并等待完成" >&2
    exit 1
fi

# Slurm 日志目录必须在提交前存在
mkdir -p $SCRATCH/FSU_seq/work_dir/logs

for sample in A B C D E F G H; do
    JOB=$(sbatch "$SCRIPT_DIR/01_cellranger_${sample}.sh" | awk '{print $4}')
    echo "样本 $sample 提交成功，Job ID: $JOB"
done
