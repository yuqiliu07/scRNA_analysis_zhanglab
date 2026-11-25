#!/bin/bash
#SBATCH --job-name=install_scDblFinder
#SBATCH --account=pi_hzh
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/install_scDblFinder.%j.out
#SBATCH --error=/u/scratch/y/yuqi/FSU_seq/work_dir/logs/install_scDblFinder.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuqiliu07@g.ucla.edu

# 一次性安装 05_qc_doublet.R 所需 R 包（scuttle 已随 DropletUtils 安装）到用户库（R/4.6.1）
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module purge
module load R/4.6.1

export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-4}"
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
pkgs <- c("scuttle", "scDblFinder", "DropletUtils", "BiocParallel", "data.table")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) BiocManager::install(need, lib = lib, update = FALSE, ask = FALSE)
for (p in pkgs) library(p, character.only = TRUE)
cat("All packages OK:", pkgs, "| scDblFinder", as.character(packageVersion("scDblFinder")), "\n")
'
