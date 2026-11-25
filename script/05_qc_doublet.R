#!/usr/bin/env Rscript
# Per-nucleus QC flags + scDblFinder doublet calls for one FSU_seq sample.
#
# Input : EmptyDrops-called barcodes (is_cell_ed == TRUE) from
#         <ED_DIR>/<sample>_emptydrops_full.csv, counts from
#         <CR_DIR>/<sample>/outs/raw_feature_bc_matrix.h5
#
# QC (scuttle; MAD = scaled MAD as in mad(); nUMI/nGene MADs on log scale):
#   nUMI   : outside median +/- 5 MAD (log)   -> fail
#            lower side via perCellQCFilters(nmads = 5), upper side via isOutlier
#   nGene  : outside median +/- 5 MAD (log)   -> fail   (same as nUMI)
#   nGene  : < 200                            -> fail
#   pct_mt : > median + 3 MAD                 -> fail   (perCellQCFilters(nmads = 3))
#   pct_mt : > 5 %                            -> fail
#   pass_qc = no fail flag
#
# Doublets: scDblFinder on all EmptyDrops-called nuclei of the sample
#   (before QC filtering, so every nucleus gets a call).
#   Expected doublet rate: DBR_PER1K per 1000 cells (default 0.004 = 10x GEM-X /
#   3' v4 multiplet rate; scDblFinder's own default 0.008 is for older chemistry).
#
# Output (<OUT_DIR>/):
#   <sample>_qc_doublet_flags.csv  one row per nucleus, first column = barcode.
#       Seurat:  md <- read.csv(f, row.names = 1); obj <- AddMetaData(obj, md)
#       (barcode format "XXXX-1" as in Cell Ranger; prefix it if the Seurat
#        object uses sample-prefixed cell names)
#   <sample>_qc_doublet_summary.csv thresholds and counts
#
# Usage: Rscript 05_qc_doublet.R <sample> [n_cores]
# Env (optional): CR_DIR, ED_DIR, OUT_DIR, DBR_PER1K

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 05_qc_doublet.R <sample> [n_cores]")
sample  <- args[1]
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

UMI_GENE_NMADS <- 5
MT_NMADS       <- 3
MT_MAX         <- 5
MIN_GENES      <- 200
SEED           <- 100

scratch   <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
cr_dir    <- Sys.getenv("CR_DIR",  file.path(scratch, "FSU_seq/02_cellranger/shared_h5"))
ed_dir    <- Sys.getenv("ED_DIR",  file.path(scratch, "FSU_seq/work_dir/emptydrops_mm10"))
out_dir   <- Sys.getenv("OUT_DIR", file.path(scratch, "FSU_seq/work_dir/qc_mm10"))
dbr_per1k <- as.numeric(Sys.getenv("DBR_PER1K", "0.004"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message(sprintf("[%s] CR_DIR=%s ED_DIR=%s OUT_DIR=%s DBR_PER1K=%g cores=%d",
                sample, cr_dir, ed_dir, out_dir, dbr_per1k, n_cores))

# ---- EmptyDrops-called nuclei ----
ed <- fread(file.path(ed_dir, paste0(sample, "_emptydrops_full.csv")))
ed_cells <- ed[toupper(as.character(is_cell_ed)) == "TRUE", barcode]

sce <- read10xCounts(file.path(cr_dir, sample, "outs/raw_feature_bc_matrix.h5"),
                     col.names = TRUE)
sce <- sce[, colnames(sce) %in% ed_cells]
counts(sce) <- as(counts(sce), "dgCMatrix")
stopifnot(ncol(sce) == length(ed_cells))
message(sprintf("[%s] EmptyDrops nuclei: %d", sample, ncol(sce)))

# ---- QC metrics ----
is_mt <- grepl("^mt-", rowData(sce)$Symbol, ignore.case = TRUE)
message(sprintf("[%s] mitochondrial genes: %d", sample, sum(is_mt)))
qc <- perCellQCMetrics(sce, subsets = list(Mito = is_mt))

f5 <- perCellQCFilters(qc, sub.fields = "subsets_Mito_percent", nmads = UMI_GENE_NMADS)
f3 <- perCellQCFilters(qc, sub.fields = "subsets_Mito_percent", nmads = MT_NMADS)
high_umi  <- isOutlier(qc$sum,      nmads = UMI_GENE_NMADS, type = "higher", log = TRUE)
high_gene <- isOutlier(qc$detected, nmads = UMI_GENE_NMADS, type = "higher", log = TRUE)

flags <- data.table(
  barcode          = colnames(sce),
  sample           = sample,
  qc_nUMI          = qc$sum,
  qc_nGene         = qc$detected,
  qc_pct_mt        = round(qc$subsets_Mito_percent, 4),
  fail_low_umi     = as.logical(f5$low_lib_size),
  fail_high_umi    = as.logical(high_umi),
  fail_low_gene    = as.logical(f5$low_n_features),
  fail_high_gene   = as.logical(high_gene),
  fail_gene_lt200  = qc$detected < MIN_GENES,
  fail_mt_3mad     = as.logical(f3$high_subsets_Mito_percent),
  fail_mt_gt5      = qc$subsets_Mito_percent > MT_MAX
)
fail_cols <- grep("^fail_", names(flags), value = TRUE)
flags[, pass_qc := !Reduce(`|`, .SD), .SDcols = fail_cols]

# ---- scDblFinder ----
set.seed(SEED)
bp <- if (n_cores > 1) MulticoreParam(workers = n_cores, RNGseed = SEED) else SerialParam(RNGseed = SEED)
dbl_args <- list(sce, BPPARAM = bp)
if ("dbr.per1k" %in% names(formals(scDblFinder))) {
  dbl_args$dbr.per1k <- dbr_per1k
} else {
  dbl_args$dbr <- dbr_per1k * ncol(sce) / 1000
}
sce <- do.call(scDblFinder, dbl_args)
stopifnot(identical(colnames(sce), flags$barcode))
flags[, scDblFinder_score := round(sce$scDblFinder.score, 5)]
flags[, scDblFinder_class := as.character(sce$scDblFinder.class)]
flags[, predicted_doublet := scDblFinder_class == "doublet"]
flags[, pass_qc_singlet := pass_qc & !predicted_doublet]

fwrite(flags, file.path(out_dir, paste0(sample, "_qc_doublet_flags.csv")))

# ---- summary ----
thr <- function(f, name) as.numeric(attr(f[[name]], "thresholds"))
hu <- attr(high_umi, "thresholds"); hg <- attr(high_gene, "thresholds")
summary_dt <- data.table(
  sample              = sample,
  n_nuclei            = nrow(flags),
  umi_lower_5mad      = round(thr(f5, "low_lib_size")[1], 1),
  umi_upper_5mad      = round(hu[["higher"]], 1),
  gene_lower_5mad     = round(thr(f5, "low_n_features")[1], 1),
  gene_upper_5mad     = round(hg[["higher"]], 1),
  gene_min            = MIN_GENES,
  mt_upper_3mad       = round(thr(f3, "high_subsets_Mito_percent")[2], 3),
  mt_max              = MT_MAX,
  n_fail_low_umi      = sum(flags$fail_low_umi),
  n_fail_high_umi     = sum(flags$fail_high_umi),
  n_fail_low_gene     = sum(flags$fail_low_gene),
  n_fail_high_gene    = sum(flags$fail_high_gene),
  n_fail_gene_lt200   = sum(flags$fail_gene_lt200),
  n_fail_mt_3mad      = sum(flags$fail_mt_3mad),
  n_fail_mt_gt5       = sum(flags$fail_mt_gt5),
  n_pass_qc           = sum(flags$pass_qc),
  pct_fail_qc         = round(100 * mean(!flags$pass_qc), 2),
  dbr_per1k           = dbr_per1k,
  n_doublet           = sum(flags$predicted_doublet),
  pct_doublet         = round(100 * mean(flags$predicted_doublet), 2),
  n_doublet_in_pass   = sum(flags$pass_qc & flags$predicted_doublet),
  n_pass_qc_singlet   = sum(flags$pass_qc_singlet),
  pct_pass_qc_singlet = round(100 * mean(flags$pass_qc_singlet), 2)
)
fwrite(summary_dt, file.path(out_dir, paste0(sample, "_qc_doublet_summary.csv")))
print(t(summary_dt))
message(sprintf("[%s] done", sample))
