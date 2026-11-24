# Step C (workflow sec. 1.2): ambient-RNA correction with SoupX, per sample.
# CellBender is the preferred method (needs a GPU; see cellbender_gpu.py) -- SoupX is the CPU stand-in and
# the cross-check. Primary analyses keep the raw counts; SoupX counts are only for robustness checks
# (discovery_qc_strategy_zh.md sec. 7c-2).
#
# Soup profile: all partitions in the empty-droplet plateau (UMI <= plateau q75 from 00_barcode_rank.py),
#   i.e. the same droplets EmptyDrops used as ambient -- not SoupX's default (<= 100 UMI), which here would
#   only contain the post-cliff tail.
# Cells: barcodes passing basic QC in 02_per_sample_qc.py, with their quick clusters.
# rho: autoEstCont; if it fails or gives an implausible value, a fallback estimate from cardiomyocyte genes
#   in clusters that do not express them (estimateNonExpressingCells) is recorded alongside.
#
# Usage: Rscript 03_soupx.R <data_dir> <work_dir> <sample> [<sample> ...]
# Output: <work_dir>/soupx/<S>_soupx_counts.h5 (10x h5, integer counts, QC-pass cells),
#         <work_dir>/soupx/<S>_soupx.csv (rho estimates), <work_dir>/soupx/<S>_soup_top.csv
suppressPackageStartupMessages({
  library(DropletUtils)
  library(SoupX)
  library(Matrix)
})
set.seed(20261108)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]; work_dir <- args[2]; samples <- args[-(1:2)]
out_dir <- file.path(work_dir, "soupx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plateau <- read.csv(file.path(work_dir, "emptydrops", "plateau.csv"))
CM_GENES <- c("Myh6", "Tnnt2", "Ttn", "Ryr2", "Myl2", "Myl3", "Actc1", "Tnni3", "Nppa", "Mb")

for (s in samples) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), s))
  sce <- read10xCounts(file.path(data_dir, s, "raw_feature_bc_matrix.h5"), col.names = TRUE)
  tod <- as(counts(sce), "CsparseMatrix")
  rownames(tod) <- make.unique(rowData(sce)$Symbol)
  gene_ids <- rowData(sce)$ID
  rm(sce); gc()
  tot <- colSums(tod)
  tod <- tod[, tot > 0]; tot <- tot[tot > 0]

  cl <- read.csv(file.path(work_dir, "per_sample", paste0(s, "_clusters.csv")), colClasses = "character")
  toc <- tod[, cl$barcode]
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  soup_upper <- plateau$plateau_q75[plateau$sample == s]
  amb <- tot > 50 & tot <= soup_upper
  soup <- rowSums(tod[, amb])
  sc <- setSoupProfile(sc, data.frame(est = soup / sum(soup), counts = soup, row.names = rownames(tod)))
  sc <- setClusters(sc, setNames(cl$quick_cluster, cl$barcode))

  rho_auto <- tryCatch({
    pdf(NULL); x <- autoEstCont(sc, doPlot = FALSE, verbose = FALSE); dev.off()
    x$metaData$rho[1]
  }, error = function(e) { cat("  autoEstCont failed:", conditionMessage(e), "\n"); NA })

  cm <- intersect(CM_GENES, rownames(toc))
  rho_cm <- tryCatch({
    ute <- estimateNonExpressingCells(sc, nonExpressedGeneList = list(CM = cm), clusters = setNames(cl$quick_cluster, cl$barcode))
    x <- calculateContaminationFraction(sc, list(CM = cm), useToEst = ute, verbose = FALSE)
    x$metaData$rho[1]
  }, error = function(e) { cat("  CM-gene estimate failed:", conditionMessage(e), "\n"); NA })

  rho_used <- if (!is.na(rho_auto) && rho_auto > 0.01 && rho_auto < 0.8) rho_auto else rho_cm
  rho_src <- if (identical(rho_used, rho_auto)) "autoEstCont" else "cardiomyocyte genes"
  sc <- setContaminationFraction(sc, rho_used, forceAccept = TRUE)
  adj <- adjustCounts(sc, roundToInt = TRUE)

  write10xCounts(file.path(out_dir, paste0(s, "_soupx_counts.h5")), adj, gene.id = gene_ids,
                 gene.symbol = rownames(adj), type = "HDF5", overwrite = TRUE, version = "3")
  top <- head(sort(soup / sum(soup), decreasing = TRUE), 50)
  write.csv(data.frame(gene = names(top), frac = as.numeric(top)),
            file.path(out_dir, paste0(s, "_soup_top.csv")), row.names = FALSE)
  res <- data.frame(sample = s, n_cells = ncol(toc), n_soup_droplets = sum(amb),
                    rho_auto = rho_auto, rho_cm_genes = rho_cm, rho_used = rho_used, rho_source = rho_src,
                    frac_counts_removed = 1 - sum(adj) / sum(toc))
  write.csv(res, file.path(out_dir, paste0(s, "_soupx.csv")), row.names = FALSE)
  print(res)
  rm(tod, toc, sc, adj); gc()
}
