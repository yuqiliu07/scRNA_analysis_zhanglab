#!/bin/bash
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
