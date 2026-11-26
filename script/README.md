# FSU_seq scripts（10x snRNA-seq，小鼠，样本 A–H）

在 hoffman2 上运行。路径都以 `$SCRATCH/FSU_seq` 为根目录，sbatch 一律由用户自己提交。
从步骤 02 起，下游使用 Novogene 提供的 mm10 h5（Cell Ranger 8.0.1，`refdata-gex-mm10-2020-A`），
放在 `02_cellranger/shared_h5/<s>/outs/`。

| 步骤 | 脚本 | 运行方式 | 输出（`work_dir/`） |
|---|---|---|---|
| 00 | `00_download_ref.sh` | bash | 下载 GRCm39 参考基因组（自己跑 Cell Ranger 时用） |
| 01 | `01_cellranger_<s>.sh`、`01_submit_all_cellranger.sh` | sbatch | Cell Ranger 9.0.0 count（GRCm39；已被 Novogene 的 mm10 结果取代） |
| 02 | `02_plateau_barcode_rank.py` | 登录节点运行 python | `emptydrops_mm10/plateau.csv`，barcode rank 图 |
| 03 | `03_emptydrops.R`、`03_emptydrops_<s>.sh` | sbatch（设置 `CR_DIR`/`OUT_DIR`） | `emptydrops_mm10/<s>_emptydrops_{full,summary}.csv`；参数 lower = plateau_q75，FDR < 0.001 |
| 04 | `04_qc_mad_hist.py` | 登录节点运行 python | `qc_mad_hist.png`、`qc_mad_thresholds.csv` |
| 05 | `05_qc_doublet.R`、`05_qc_doublet_<s>.sh` | sbatch | `qc_mm10/<s>_qc_doublet_flags.csv`（scuttle QC + scDblFinder） |
| 06 | `06_qc_doublet_summary_plot.py` | 登录节点运行 python | `qc_mm10/qc_doublet_summary.png` |
| 07 | `07_install_Seurat.sh`、`07_seurat_merge.R/.sh` | sbatch | `seurat_mm10/fsu_merged_mm10.rds`（只保留 pass_qc_singlet 的细胞；layer 为 counts.A..H；HVG 取各样本并集） |
| 08 | `08_umap_mixing.R/.sh` | sbatch | `umap_mm10/`：PCA → UMAP（PC 1–20，未做整合），并计算 kNN mixing score（可选，尚未运行） |
| 09 | `09_harmony_cluster.R/.sh` | sbatch | `harmony_mm10/`：PCA（计算 50 个 PC，使用前 30 个）→ Harmony（按样本整合）→ UMAP + SNN（k=20）→ Leiden（igraph，resolution=1.0）；输出整合前后的 UMAP 对比图和 cluster 图 |
| 10 | `10_compare_novogene_umap.R/.sh` | sbatch | `compare_novogene/`：与 Novogene 交付的 UMAP/cluster 对比（细胞重叠、ARI、2D kNN 混合度、交叉着色 UMAP） |

需要安装的 R 包见 `03_install_R_packages.sh`、`05_install_R_packages.sh`、`07_install_Seurat.sh`，
R 版本为 R/4.6.1。

上传 Google Drive（`gdrive:FSU_seq_results`）的规则：只传图片、汇总 CSV、参数和 Slurm 日志，
不传数据文件（h5/rds/逐细胞矩阵）和代码。
