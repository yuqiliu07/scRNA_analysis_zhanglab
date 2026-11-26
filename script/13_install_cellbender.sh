#!/bin/bash
# CellBender 0.3.2 环境（已在 2026-11-26 于登录节点装好：$HOME/envs/cellbender，约 6 GB）
# 重新安装时在登录节点运行：bash 13_install_cellbender.sh
#   torch 2.5.1+cu124：支持 V100 / RTX2080Ti / A100 / H200（不支持 Blackwell RTXPRO6000 和 P4）
#   numpy<2、scipy<1.14、pandas<2.2：CellBender 0.3.2 与更新版本不兼容
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module load python/3.10.20
python3 -m venv $HOME/envs/cellbender
source $HOME/envs/cellbender/bin/activate
pip install --upgrade pip
pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cu124
pip install --no-cache-dir "numpy<2" "scipy<1.14" "pandas<2.2" cellbender==0.3.2
cellbender --version
