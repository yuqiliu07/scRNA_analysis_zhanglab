#!/usr/bin/env Rscript
# EmptyDrops cell calling for one FSU_seq sample.
#
# lower  = this sample's plateau_q75 from <OUT_DIR>/plateau.csv
#          (barcodes with total UMI <= lower form the ambient profile)
# niters = 10000, FDR threshold = 0.001
#
# Outputs (in <OUT_DIR>/):
#   <sample>_emptydrops_full.csv    every barcode with total UMI > lower
#   <sample>_emptydrops_summary.csv EmptyDrops vs Cell Ranger cell calls
#
# Usage: Rscript 03_emptydrops.R <sample> [n_cores]
#
# Optional environment variables (defaults = original GRCm39 run):
#   CR_DIR   Cell Ranger output root, reads <CR_DIR>/<sample>/outs/
#            default: $SCRATCH/FSU_seq/02_cellranger
#            Novogene mm10: $SCRATCH/FSU_seq/02_cellranger/shared_h5
#   OUT_DIR  directory holding plateau.csv and receiving outputs
#            default: $SCRATCH/FSU_seq/work_dir/emptydrops
#            Novogene mm10: $SCRATCH/FSU_seq/work_dir/emptydrops_mm10

suppressPackageStartupMessages({
  library(DropletUtils)
  library(BiocParallel)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 03_emptydrops.R <sample> [n_cores]")
sample  <- args[1]
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

NITERS  <- 10000
FDR_THR <- 0.001
SEED    <- 100

scratch <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
cr_dir  <- Sys.getenv("CR_DIR",  file.path(scratch, "FSU_seq/02_cellranger"))
out_dir <- Sys.getenv("OUT_DIR", file.path(scratch, "FSU_seq/work_dir/emptydrops"))
cr_outs <- file.path(cr_dir, sample, "outs")
message(sprintf("[%s] CR_DIR = %s ; OUT_DIR = %s", sample, cr_dir, out_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- lower = plateau_q75 ----
plateau <- fread(file.path(out_dir, "plateau.csv"))
lower <- plateau$plateau_q75[plateau$sample == sample]
if (length(lower) != 1 || is.na(lower)) stop("plateau_q75 not found for sample ", sample)
message(sprintf("[%s] lower = plateau_q75 = %g, niters = %d, FDR < %g, cores = %d",
                sample, lower, NITERS, FDR_THR, n_cores))

# ---- raw counts ----
sce <- read10xCounts(file.path(cr_outs, "raw_feature_bc_matrix.h5"), col.names = TRUE)
m <- as(counts(sce), "dgCMatrix")
rm(sce); invisible(gc())
totals <- colSums(m)
message(sprintf("[%s] raw barcodes: %d; barcodes with UMI > lower: %d",
                sample, ncol(m), sum(totals > lower)))

# ---- Cell Ranger calls ----
cr_bc <- fread(file.path(cr_outs, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
               header = FALSE)$V1

# ---- EmptyDrops ----
set.seed(SEED)
bp <- if (n_cores > 1) MulticoreParam(workers = n_cores, RNGseed = SEED) else SerialParam(RNGseed = SEED)
ed <- emptyDrops(m, lower = lower, niters = NITERS, BPPARAM = bp)

# ---- (1) full results: all barcodes with UMI > lower ----
tested <- !is.na(ed$PValue)               # == totals > lower
full <- data.table(
  barcode       = rownames(ed)[tested],
  total_umi     = ed$Total[tested],
  log_prob      = signif(ed$LogProb[tested], 6),
  p_value       = signif(ed$PValue[tested], 4),
  fdr           = signif(ed$FDR[tested], 4),
  limited       = ed$Limited[tested],
  is_cell_ed    = ed$FDR[tested] < FDR_THR,
  is_cell_cr    = rownames(ed)[tested] %in% cr_bc
)
setorder(full, -total_umi)
fwrite(full, file.path(out_dir, paste0(sample, "_emptydrops_full.csv")))

# ---- (2) summary ----
ed_bc <- full[is_cell_ed == TRUE, barcode]
summary_dt <- data.table(
  sample                 = sample,
  lower_plateau_q75      = lower,
  niters                 = NITERS,
  fdr_threshold          = FDR_THR,
  n_tested               = nrow(full),
  emptydrops_cells       = length(ed_bc),
  cellranger_cells       = length(cr_bc),
  both                   = length(intersect(ed_bc, cr_bc)),
  emptydrops_only        = length(setdiff(ed_bc, cr_bc)),
  cellranger_only        = length(setdiff(cr_bc, ed_bc)),
  # CR cells with UMI <= lower are never tested by EmptyDrops
  cellranger_only_below_lower = sum(cr_bc %in% names(totals)[totals <= lower]),
  # non-significant but Limited: more iterations could lower their p-value
  limited_not_significant = sum(full$limited & !full$is_cell_ed),
  ed_cells_median_umi    = median(full[is_cell_ed == TRUE, total_umi]),
  ed_only_median_umi     = median(full[barcode %in% setdiff(ed_bc, cr_bc), total_umi])
)
fwrite(summary_dt, file.path(out_dir, paste0(sample, "_emptydrops_summary.csv")))

print(summary_dt)
message(sprintf("[%s] done", sample))
