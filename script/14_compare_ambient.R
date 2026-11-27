#!/usr/bin/env Rscript
# Compare ambient-RNA removal: SoupX (step 12) vs CellBender (step 13), plus cell calls
# CellBender vs EmptyDrops, for the 8 FSU_seq samples.
#
# Common nuclei = the QC-passed, clustered nuclei of step 11 (cluster_clean_per_cell.csv.gz).
# Per sample, on those nuclei:
#   raw        raw_feature_bc_matrix.h5
#   soupx      <SOUPX_DIR>/<s>_soupx_counts.rds   (adjustCounts with rho_final = autoEstCont)
#   cellbender <CB_DIR>/<s>/<s>_cellbender.h5     (all analysed droplets, FPR 0.01)
# Metrics:
#   1. fraction of UMI removed (SoupX, CellBender) vs the SoupX rho estimates (auto, CM fallback)
#   2. fraction removed per cell type (cluster -> cell type map below, from the step-11 annotation),
#      on nuclei that CellBender also calls as cells (it zeroes the others)
#   3. CM-gene UMI left in nuclei of non-CM clusters (CM panel score < CM_SCORE_MAX), and the
#      fraction of those nuclei still carrying >= 1 CM-gene UMI
#   4. per-nucleus removed fraction, SoupX vs CellBender (Spearman rho)
#   5. top genes removed by each method
#   6. cell calls: CellBender cell_barcodes.csv vs EmptyDrops (is_cell_ed) vs the step-11 nuclei
#
# Outputs (<OUT_DIR>/, PNG + CSV are uploaded by 14_compare_ambient.sh):
#   ambient_summary.csv            one row per sample (metrics 1, 3, 4, 6)
#   ambient_by_celltype.csv        sample x cell type x method removed fraction
#   ambient_top_removed_genes.csv  top 20 removed genes per sample and method
#   cellcall_by_celltype.csv       step-11 nuclei not called by CellBender, per cell type
#   ambient_removed_by_sample.png  rho_auto, rho_fallback, SoupX and CellBender removed fraction
#   ambient_removed_by_celltype.png
#   ambient_cm_genes_nonCM.png     CM-gene UMI left in non-CM nuclei
#   ambient_per_cell_scatter.png   per-nucleus removed fraction, SoupX vs CellBender
#   cellcall_overlap.png           CellBender vs EmptyDrops cell calls
#
# Usage: Rscript 14_compare_ambient.R [samples, default A..H]
# Env (optional): TEST_MAX_NUCLEI, RAW_DIR, ANNOT_DIR, ED_DIR, SOUPX_DIR, CB_DIR, OUT_DIR, CM_SCORE_MAX (1.8)

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(data.table)
  library(ggplot2)
})

scratch   <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
wd        <- file.path(scratch, "FSU_seq/work_dir")
raw_dir   <- Sys.getenv("RAW_DIR",   file.path(scratch, "FSU_seq/02_cellranger/shared_h5"))
annot_dir <- Sys.getenv("ANNOT_DIR", file.path(wd, "annot_mm10"))
ed_dir    <- Sys.getenv("ED_DIR",    file.path(wd, "emptydrops_mm10"))
soupx_dir <- Sys.getenv("SOUPX_DIR", file.path(wd, "soupx_mm10"))
cb_dir    <- Sys.getenv("CB_DIR",    file.path(wd, "cellbender_mm10"))
out_dir   <- Sys.getenv("OUT_DIR",   file.path(wd, "ambient_compare_mm10"))
cm_max    <- as.numeric(Sys.getenv("CM_SCORE_MAX", "1.8"))
test_n    <- as.integer(Sys.getenv("TEST_MAX_NUCLEI", "0"))   # >0: first n nuclei per sample (login-node test)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args    <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else LETTERS[1:8]

CM_GENES <- c("Myh6", "Tnnt2", "Ryr2", "Ttn", "Ankrd1", "Mhrt", "Rbm20", "Myl2", "Myl4", "Nppa")
# cluster -> cell type (step-11 annotation; mixed clusters go with their main type)
CELLTYPE <- c(
  setNames(rep("CM", 3), c(0, 4, 6)),
  setNames(rep("CM+X (doublet?)", 9), c(10, 11, 14, 18, 19, 30, 31, 33, 34)),
  setNames(rep("Fibroblast", 5), c(1, 5, 23, 13, 24)),
  setNames(rep("Endothelial", 4), c(3, 7, 20, 2)),
  setNames(rep("Endocardial", 2), c(12, 22)),
  "16" = "Lymphatic EC",
  setNames(rep("Pericyte", 2), c(8, 21)),
  "25" = "SMC",
  setNames(rep("Macrophage", 4), c(9, 15, 17, 27)),
  "28" = "T/NK + B", "29" = "Schwann", "35" = "Epicardial",
  "26" = "Other (26, 32)", "32" = "Other (26, 32)")
CT_ORDER <- c("CM", "CM+X (doublet?)", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic EC",
              "Pericyte", "SMC", "Macrophage", "T/NK + B", "Schwann", "Epicardial", "Other (26, 32)")

cells  <- fread(file.path(annot_dir, "cluster_clean_per_cell.csv.gz"))
pscore <- fread(file.path(annot_dir, "cluster_panel_scores.csv"))
nonCM_clusters <- as.character(pscore[Cardiomyocyte < cm_max, cluster])
rho_tab <- fread(file.path(soupx_dir, "soupx_rho.csv"))

# read a 10x-format h5 (raw or CellBender output), keep barcodes `bc` in that order; rows = Ensembl IDs.
# Streams the CSC arrays with rhdf5 in blocks and keeps only the wanted columns, so memory stays
# ~ size of the result (DelayedArray/TENxMatrix fails on the CellBender h5).
read_h5_cols <- function(path, bc, block = 4e6) {
  g <- function(x, ...) rhdf5::h5read(path, paste0("matrix/", x), bit64conversion = "double", ...)
  on.exit(rhdf5::h5closeAll())
  bcs <- as.character(g("barcodes"))
  pos <- match(bcs, bc)                         # column of the result for each h5 column (NA = drop)
  if (sum(!is.na(pos)) != length(bc)) stop(sprintf("%s: %d barcodes missing", path, length(setdiff(bc, bcs))))
  ids <- as.character(g("features/id")); sym <- as.character(g("features/name"))
  indptr <- as.numeric(g("indptr")); nnz <- indptr[length(indptr)]
  starts <- indptr[-length(indptr)]
  I <- list(); J <- list(); X <- list(); k <- 0
  for (s0 in seq(0, nnz - 1, by = block)) {
    n <- min(block, nnz - s0)
    col <- findInterval(s0 + seq_len(n) - 1, starts)   # 1-based h5 column of each element
    jj  <- pos[col]; keep <- !is.na(jj)
    if (!any(keep)) next
    k <- k + 1
    I[[k]] <- as.integer(g("indices", start = s0 + 1, count = n))[keep] + 1L
    X[[k]] <- as.numeric(g("data", start = s0 + 1, count = n))[keep]
    J[[k]] <- jj[keep]
  }
  m <- sparseMatrix(i = unlist(I), j = unlist(J), x = unlist(X),
                    dims = c(length(ids), length(bc)), dimnames = list(ids, bc))
  list(m = m, symbol = setNames(sym, ids))
}
colsum_nz <- function(m) as.numeric(Matrix::colSums(m))

summ <- list(); by_ct <- list(); topg <- list(); cc_ct <- list(); percell <- list()
for (s in samples) {
  t0 <- Sys.time()
  cs  <- cells[sample == s]
  if (test_n > 0) cs <- head(cs, test_n)
  bc  <- sub(paste0("^", s, "_"), "", cs$cell)
  ct  <- CELLTYPE[as.character(cs$cluster_clean)]; ct[is.na(ct)] <- "Other (26, 32)"
  nonCM <- as.character(cs$cluster_clean) %in% nonCM_clusters

  sx  <- readRDS(file.path(soupx_dir, paste0(s, "_soupx_counts.rds")))
  # SoupX rows are make.unique(Symbol) in raw-h5 order, columns are "<s>_<barcode>"
  sx  <- sx[, paste0(s, "_", bc), drop = FALSE]
  colnames(sx) <- bc
  invisible(gc())
  r   <- read_h5_cols(file.path(raw_dir, s, "outs/raw_feature_bc_matrix.h5"), bc)
  raw <- r$m; sym <- r$symbol
  cb  <- read_h5_cols(file.path(cb_dir, s, paste0(s, "_cellbender.h5")), bc)$m
  cb  <- cb[rownames(raw), , drop = FALSE]
  if (nrow(sx) != nrow(raw)) stop(s, ": SoupX gene count differs from raw h5")
  rownames(sx) <- rownames(raw)

  n_raw <- colsum_nz(raw); n_sx <- colsum_nz(sx); n_cb <- colsum_nz(cb)
  cm_id  <- names(sym)[sym %in% CM_GENES]
  cm_raw <- colsum_nz(raw[cm_id, , drop = FALSE])
  cm_sx  <- colsum_nz(sx[cm_id, , drop = FALSE])
  cm_cb  <- colsum_nz(cb[cm_id, , drop = FALSE])

  # cell calls
  ed  <- fread(file.path(ed_dir, paste0(s, "_emptydrops_full.csv")))
  ed_cells <- ed[is_cell_ed == TRUE, barcode]
  cb_cells <- fread(file.path(cb_dir, s, paste0(s, "_cellbender_cell_barcodes.csv")), header = FALSE)$V1
  umi <- setNames(ed$total_umi, ed$barcode)
  ed_only <- setdiff(ed_cells, cb_cells); cb_only <- setdiff(cb_cells, ed_cells)
  annot_not_cb <- !(bc %in% cb_cells)

  rr <- rho_tab[sample == s]
  pc_sx <- 1 - n_sx / n_raw; pc_cb <- 1 - n_cb / n_raw
  summ[[s]] <- data.table(
    sample = s, n_nuclei = length(bc),
    rho_auto = rr$rho_auto, rho_cm_fallback = rr$rho_cm_fallback,
    frac_removed_soupx = 1 - sum(n_sx) / sum(n_raw),
    frac_removed_cellbender = 1 - sum(n_cb) / sum(n_raw),
    # CellBender zeroes droplets it does not call as cells; restrict to nuclei it did call
    frac_removed_soupx_cbcalled = 1 - sum(n_sx[!annot_not_cb]) / sum(n_raw[!annot_not_cb]),
    frac_removed_cellbender_cbcalled = 1 - sum(n_cb[!annot_not_cb]) / sum(n_raw[!annot_not_cb]),
    cb_metrics_frac_removed = fread(file.path(cb_dir, s, paste0(s, "_cellbender_metrics.csv")),
                                    header = FALSE)[V1 == "fraction_counts_removed", V2],
    percell_spearman = suppressWarnings(cor(pc_sx, pc_cb, method = "spearman")),
    median_percell_removed_soupx = median(pc_sx), median_percell_removed_cb = median(pc_cb),
    cm_umi_left_nonCM_soupx = sum(cm_sx[nonCM]) / sum(cm_raw[nonCM]),
    cm_umi_left_nonCM_cb    = sum(cm_cb[nonCM]) / sum(cm_raw[nonCM]),
    cm_detect_nonCM_raw   = mean(cm_raw[nonCM] > 0),
    cm_detect_nonCM_soupx = mean(cm_sx[nonCM] > 0),
    cm_detect_nonCM_cb    = mean(cm_cb[nonCM] > 0),
    cm_umi_left_CM_soupx = sum(cm_sx[ct == "CM"]) / sum(cm_raw[ct == "CM"]),
    cm_umi_left_CM_cb    = sum(cm_cb[ct == "CM"]) / sum(cm_raw[ct == "CM"]),
    ed_cells = length(ed_cells), cb_cells = length(cb_cells),
    both = length(intersect(ed_cells, cb_cells)), ed_only = length(ed_only), cb_only = length(cb_only),
    median_umi_both = median(umi[intersect(ed_cells, cb_cells)], na.rm = TRUE),
    median_umi_ed_only = median(umi[ed_only], na.rm = TRUE),
    median_umi_cb_only = median(umi[cb_only], na.rm = TRUE),
    nuclei_not_cb_called = sum(annot_not_cb), frac_nuclei_not_cb_called = mean(annot_not_cb))

  for (meth in c("SoupX", "CellBender")) {
    n_new <- if (meth == "SoupX") n_sx else n_cb
    by_ct[[paste(s, meth)]] <- data.table(sample = s, method = meth, celltype = ct,
                                          raw = n_raw, new = n_new)[!annot_not_cb][
      , .(n_nuclei = .N, frac_removed = 1 - sum(new) / sum(raw)), by = .(sample, method, celltype)]
    m_new <- if (meth == "SoupX") sx else cb
    rem <- Matrix::rowSums(raw) - Matrix::rowSums(m_new)
    o <- head(order(rem, decreasing = TRUE), 20)
    topg[[paste(s, meth)]] <- data.table(sample = s, method = meth, rank = seq_along(o),
                                         gene = sym[rownames(raw)[o]], umi_removed = rem[o],
                                         frac_of_gene_removed = rem[o] / Matrix::rowSums(raw)[o],
                                         frac_of_all_removed = rem[o] / sum(rem))
  }
  cc_ct[[s]] <- data.table(sample = s, celltype = ct, not_cb = annot_not_cb, umi = n_raw)[
    , .(n_nuclei = .N, n_not_cb_called = sum(not_cb), frac_not_cb_called = mean(not_cb),
        median_umi_not_called = as.numeric(median(umi[not_cb]))), by = .(sample, celltype)]
  set.seed(1)
  keep <- sample(length(bc), min(5000, length(bc)))
  percell[[s]] <- data.table(sample = s, soupx = pc_sx[keep], cellbender = pc_cb[keep])
  message(sprintf("[%s] %d nuclei | removed: SoupX %.3f, CellBender %.3f (rho auto %.3f, fallback %.3f) | %.1f min",
                  s, length(bc), summ[[s]]$frac_removed_soupx, summ[[s]]$frac_removed_cellbender,
                  rr$rho_auto, rr$rho_cm_fallback, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(raw, sx, cb); invisible(gc())
}

summ <- rbindlist(summ); by_ct <- rbindlist(by_ct); topg <- rbindlist(topg)
cc_ct <- rbindlist(cc_ct); percell <- rbindlist(percell)
fwrite(summ, file.path(out_dir, "ambient_summary.csv"))
fwrite(by_ct, file.path(out_dir, "ambient_by_celltype.csv"))
fwrite(topg, file.path(out_dir, "ambient_top_removed_genes.csv"))
fwrite(cc_ct, file.path(out_dir, "cellcall_by_celltype.csv"))
print(summ[, .(sample, rho_auto, rho_cm_fallback, frac_removed_soupx, frac_removed_cellbender,
               frac_removed_cellbender_cbcalled,
               percell_spearman, cm_umi_left_nonCM_soupx, cm_umi_left_nonCM_cb,
               cm_detect_nonCM_soupx, cm_detect_nonCM_cb)], digits = 3)
print(summ[, .(sample, ed_cells, cb_cells, both, ed_only, cb_only, median_umi_both,
               median_umi_ed_only, median_umi_cb_only, nuclei_not_cb_called)])

# ---- plots ----
COL <- c("SoupX rho (autoEstCont)" = "#257de0", "SoupX rho (CM fallback)" = "#6a4ce2",
         "SoupX removed" = "#00b894", "CellBender removed" = "#f47230")
theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                plot.background = element_rect(fill = "white", colour = NA)))

d1 <- melt(summ[, .(sample, `SoupX rho (autoEstCont)` = rho_auto, `SoupX rho (CM fallback)` = rho_cm_fallback,
                    `SoupX removed` = frac_removed_soupx, `CellBender removed` = frac_removed_cellbender)],
           id.vars = "sample", variable.name = "measure", value.name = "fraction")
d1[, measure := factor(measure, levels = names(COL))]
p <- ggplot(d1, aes(sample, fraction, colour = measure)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(x = "sample", y = "fraction of UMI in QC nuclei",
       title = "Ambient RNA: SoupX rho vs UMI actually removed") +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "ambient_removed_by_sample.png"), p, width = 8, height = 4.5, dpi = 150)

by_ct[, celltype := factor(celltype, levels = rev(CT_ORDER))]
p <- ggplot(by_ct, aes(frac_removed, celltype, colour = method)) +
  geom_point(size = 2, alpha = 0.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(SoupX = "#00b894", CellBender = "#f47230"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "fraction of UMI removed (one point per sample; nuclei called by both)", y = NULL,
       title = "UMI removed per cell type") +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "#eeeeee"))
ggsave(file.path(out_dir, "ambient_removed_by_celltype.png"), p, width = 7, height = 5.5, dpi = 150)

d3 <- melt(summ[, .(sample, SoupX = cm_umi_left_nonCM_soupx, CellBender = cm_umi_left_nonCM_cb)],
           id.vars = "sample", variable.name = "method", value.name = "left")
p <- ggplot(d3, aes(sample, left, colour = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(SoupX = "#00b894", CellBender = "#f47230"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(x = "sample", y = "CM-gene UMI left (raw = 100%)",
       title = "CM genes (Myh6, Tnnt2, Ryr2, Ttn ...) left in non-CM nuclei") +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "ambient_cm_genes_nonCM.png"), p, width = 7, height = 4.2, dpi = 150)

p <- ggplot(percell, aes(soupx, cellbender)) +
  geom_bin2d(bins = 60) + scale_fill_gradient(low = "#d6e6f9", high = "#0b3d78", trans = "log10", name = "nuclei") +
  geom_abline(linetype = 2, colour = "#666666") +
  facet_wrap(~sample, ncol = 4) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  geom_text(data = summ, aes(x = 0.02, y = 0.97, label = sprintf("Spearman %.2f", percell_spearman)),
            hjust = 0, size = 3, colour = "#191919", inherit.aes = FALSE) +
  labs(x = "SoupX: fraction removed per nucleus", y = "CellBender: fraction removed per nucleus",
       title = "Per-nucleus removed fraction (5,000 random nuclei per sample)")
ggsave(file.path(out_dir, "ambient_per_cell_scatter.png"), p, width = 10, height = 5.8, dpi = 150)

d5 <- melt(summ[, .(sample, `both` = both, `EmptyDrops only` = ed_only, `CellBender only` = cb_only)],
           id.vars = "sample", variable.name = "call", value.name = "n")
d5[, call := factor(call, levels = c("CellBender only", "EmptyDrops only", "both"))]
p <- ggplot(d5, aes(sample, n, fill = call)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(both = "#257de0", `EmptyDrops only` = "#f2a600", `CellBender only` = "#f47230"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "sample", y = "barcodes called as cells", title = "Cell calls: CellBender vs EmptyDrops") +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "cellcall_overlap.png"), p, width = 7, height = 4.5, dpi = 150)
message("done -> ", out_dir)
