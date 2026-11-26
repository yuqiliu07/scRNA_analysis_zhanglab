#!/bin/bash
# CellBender 0.3.2 环境：$HOME/envs/cellbender_t113（2026-11-26 已在登录节点装好，约 4 GB）
# 重新安装时在登录节点运行：bash 13_install_cellbender.sh
#   torch 1.13.1+cu117 + pyro-ppl 1.8.6：用 torch 2.5 + pyro 1.9 时 checkpoint 无法保存
#     （"cannot pickle 'weakref.ReferenceType' object"，CellBender 自带的 test_checkpoint 同样失败），
#     之后的 posterior 步骤因此报错；此组合下 test_save_and_load_cellbender_checkpoint 通过
#   支持 V100 / RTX2080Ti / A100（不支持 H200、Blackwell RTXPRO6000、P4）
#   numpy<2、scipy<1.14、pandas<2.2：CellBender 0.3.2 与更新版本不兼容
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module load python/3.10.20
python3 -m venv $HOME/envs/cellbender_t113
source $HOME/envs/cellbender_t113/bin/activate
pip install --upgrade pip
pip install --no-cache-dir torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
pip install --no-cache-dir "numpy<2" "scipy<1.14" "pandas<2.2" "pyro-ppl==1.8.6" cellbender==0.3.2
cellbender --version
