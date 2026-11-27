#!/bin/bash
# 提交 8 个 CellBender GPU 作业（由 13_make_cellbender_scripts.py 生成）
# 用法: bash 13_submit_all_cellbender.sh [样本...]，如 bash 13_submit_all_cellbender.sh C D E F G H
S=${@:-A B C D E F G H}
for s in $S; do
  sbatch $HOME/scrna_project/script/13_cellbender_$s.sh
done
