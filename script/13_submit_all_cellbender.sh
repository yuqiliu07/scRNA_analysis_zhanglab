#!/bin/bash
# 提交 8 个 CellBender GPU 作业（由 13_make_cellbender_scripts.py 生成）
for s in 13_cellbender_A.sh 13_cellbender_B.sh 13_cellbender_C.sh 13_cellbender_D.sh 13_cellbender_E.sh 13_cellbender_F.sh 13_cellbender_G.sh 13_cellbender_H.sh; do
  sbatch $HOME/scrna_project/script/$s
done
