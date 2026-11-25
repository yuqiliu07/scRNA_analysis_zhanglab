#!/bin/bash
# 一键提交 EmptyDrops：如 R 包未安装，先提交安装作业，8 个样本作业在其成功后自动开始
# 用法（在登录节点用 bash 运行，不要 sbatch）：
#   bash $SCRATCH/FSU_seq/script/03_submit_all_emptydrops.sh
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p $SCRATCH/FSU_seq/work_dir/logs

if [ ! -f $SCRATCH/FSU_seq/work_dir/emptydrops/plateau.csv ]; then
    echo "错误：缺少 plateau.csv，请先运行 02_plateau_barcode_rank.py" >&2
    exit 1
fi

. /u/local/Modules/default/init/modules.sh
module load R/4.6.1 >/dev/null 2>&1
DEP=""
if Rscript -e 'quit(status = !all(sapply(c("DropletUtils","BiocParallel","data.table"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    echo "R 包已安装，直接提交 EmptyDrops 作业"
else
    INSTALL_JOB=$(sbatch --parsable "$SCRIPT_DIR/03_install_R_packages.sh")
    echo "R 包安装作业提交成功，Job ID: $INSTALL_JOB"
    DEP="--dependency=afterok:$INSTALL_JOB"
fi

for sample in A B C D E F G H; do
    JOB=$(sbatch --parsable $DEP "$SCRIPT_DIR/03_emptydrops_${sample}.sh")
    echo "样本 $sample 提交成功，Job ID: $JOB"
done
