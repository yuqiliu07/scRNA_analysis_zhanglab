#!/bin/bash
#SBATCH --job-name=install_Seurat
#SBATCH --account=pi_hzh
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/install_Seurat.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/install_Seurat.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 一次性安装 Seurat v5（CRAN）到用户库（R/4.6.1）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-8}"
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
pkgs <- c("SeuratObject", "Seurat")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, lib = lib)
for (p in pkgs) library(p, character.only = TRUE)
stopifnot(packageVersion("Seurat") >= "5.0.0")
cat("Seurat", as.character(packageVersion("Seurat")), "SeuratObject", as.character(packageVersion("SeuratObject")), "OK\n")
'
