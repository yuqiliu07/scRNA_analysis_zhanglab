# Step A (workflow docs/mouse_sc_qc_workflow_zh.md, sec. 1.2 / 3.3-1):
# re-call cells from the raw droplet matrix with EmptyDrops and compare with Cell Ranger's call.
#
# Usage: Rscript 01_emptydrops.R <data_dir> <work_dir> <sample> [<sample> ...]
#   data_dir: folder holding <sample>/raw_feature_bc_matrix.h5 and filtered_feature_bc_matrix.h5
#   work_dir: output folder for per-sample tables (kept outside the git repo)
# Output per sample (work_dir/emptydrops/):
#   <S>_emptydrops.csv   every barcode with >= LOWER UMIs: total, EmptyDrops LogProb/PValue/FDR/Limited,
#                        called (FDR < 0.001), in_cellranger
#   <S>_rank.csv         barcode-rank curve (log-spaced subsample) + knee/inflection
#   <S>_summary.csv      one-row summary
suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
})

# Ambient set: barcodes with <= LOWER UMIs. In this batch every partition carries a few hundred ambient UMIs,
# so the default (100) only sees the tail beyond the partition cliff. LOWER is set per sample to the
# 75th percentile of the empty-droplet plateau (work_dir/emptydrops/plateau.csv from 00_barcode_rank.py).
FDR_CUT <- 0.001
NITERS <- 10000
set.seed(20261108)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]; work_dir <- args[2]; samples <- args[-(1:2)]
out_dir <- file.path(work_dir, "emptydrops")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_h5 <- function(path) {
  sce <- read10xCounts(path, col.names = TRUE)
  m <- as(counts(sce), "CsparseMatrix")
  rownames(m) <- rowData(sce)$Symbol
  m
}

plateau <- read.csv(file.path(out_dir, "plateau.csv"))

for (s in samples) {
  t0 <- Sys.time()
  LOWER <- plateau$plateau_q75[plateau$sample == s]
  cat(sprintf("  lower = %d (plateau q75)\n", LOWER))
  cat(sprintf("[%s] %s: reading raw matrix\n", format(t0, "%H:%M:%S"), s))
  raw <- read_h5(file.path(data_dir, s, "raw_feature_bc_matrix.h5"))
  filt_bc <- colnames(read10xCounts(file.path(data_dir, s, "filtered_feature_bc_matrix.h5"), col.names = TRUE))
  tot <- Matrix::colSums(raw)
  keep <- tot > 0
  raw <- raw[, keep]; tot <- tot[keep]
  cat(sprintf("  %d barcodes with >0 UMI, %d in Cell Ranger filtered\n", ncol(raw), length(filt_bc)))

  br <- barcodeRanks(raw, lower = LOWER)
  knee <- metadata(br)$knee; infl <- metadata(br)$inflection

  ed <- emptyDrops(raw, lower = LOWER, niters = NITERS, test.ambient = FALSE)
  ed <- ed[!is.na(ed$FDR) | ed$Total > LOWER, ]
  res <- data.frame(barcode = rownames(ed), total = ed$Total, LogProb = ed$LogProb,
                    PValue = ed$PValue, FDR = ed$FDR, Limited = ed$Limited)
  res <- res[res$total > LOWER, ]
  res$called <- !is.na(res$FDR) & res$FDR < FDR_CUT
  res$in_cellranger <- res$barcode %in% filt_bc
  write.csv(res, file.path(out_dir, paste0(s, "_emptydrops.csv")), row.names = FALSE)

  # rank curve: log-spaced subsample so the file stays small
  o <- order(tot, decreasing = TRUE); r <- seq_along(o)
  idx <- unique(round(10^seq(0, log10(length(o)), length.out = 3000)))
  rk <- data.frame(rank = r[idx], total = tot[o][idx], knee = knee, inflection = infl)
  write.csv(rk, file.path(out_dir, paste0(s, "_rank.csv")), row.names = FALSE)

  limited_nonsig <- sum(res$Limited & !res$called)
  sm <- data.frame(
    sample = s, lower = LOWER, barcodes_nonzero = ncol(raw), cellranger_cells = length(filt_bc),
    emptydrops_cells = sum(res$called), both = sum(res$called & res$in_cellranger),
    cellranger_only = sum(!res$called & res$in_cellranger),
    emptydrops_only = sum(res$called & !res$in_cellranger),
    knee = knee, inflection = infl,
    median_umi_ed = median(res$total[res$called]),
    min_umi_ed = min(res$total[res$called]),
    limited_nonsig = limited_nonsig,
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  write.csv(sm, file.path(out_dir, paste0(s, "_summary.csv")), row.names = FALSE)
  print(sm)
  rm(raw, ed, res); gc()
}
