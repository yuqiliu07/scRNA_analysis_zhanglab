#!/usr/bin/env Rscript
# SoupX ambient-RNA estimation for the 8 FSU_seq samples (Novogene mm10 raw h5).
#
# Per sample:
#  1. tod = raw_feature_bc_matrix.h5 (all barcodes); toc = raw counts of the nuclei kept after QC
#     and clustered (cluster_clean_per_cell.csv.gz of 11_cluster_clean_annot.R).
#  2. Soup profile = summed counts of barcodes with 0 < total UMI <= plateau_q75
#     (plateau.csv of 02_plateau_barcode_rank.py; same cut-off as the EmptyDrops "lower").
#  3. Clusters = cleaned Leiden clusters (leiden_clean, from step 11).
#  4. rho:
#     a. autoEstCont (SoupX defaults). It fails with an error when no usable marker genes are
#        found or when rho > 0.5; then the fallback is used.
#     b. Fallback, always computed for comparison: cardiomyocyte panel genes (sxj_code PANELS)
#        in nuclei of non-cardiomyocyte clusters (cluster mean Cardiomyocyte panel score
#        < CM_SCORE_MAX in cluster_panel_scores.csv). These nuclei should not express CM genes,
#        so every CM-gene UMI there is soup:
#            rho = sum(CM-gene UMI) / sum(nUMI x soup fraction of CM genes)
#        (the SoupX marker estimate, = calculateContaminationFraction with useToEst = these
#        nuclei). The CM-gene detection rate (fraction of these nuclei with >= 1 CM-gene UMI)
#        is reported next to it.
#     rho_final = autoEstCont if it succeeded, else the fallback.
#  5. adjustCounts(roundToInt = TRUE) with rho_final -> <sample>_soupx_counts.rds (stays on Hoffman)
#
# Outputs (<OUT_DIR>/):
#   soupx_rho.csv                rho per sample (auto, its FWHM interval, fallback, final, method)
#   soupx_top_soup_genes.csv     top 20 soup genes per sample
#   soupx_rho_by_sample.png      rho per sample, autoEstCont vs CM-gene fallback
#   soupx_rho_posterior.png      autoEstCont posterior of rho per sample
#   <sample>_soupx_counts.rds    corrected counts (dgCMatrix, genes x nuclei)
#
# Usage: Rscript 12_soupx.R [samples, default A..H]
# Env (optional): RAW_DIR (<RAW_DIR>/<sample>/outs/raw_feature_bc_matrix.h5), PLATEAU_CSV,
#                 ANNOT_DIR, OUT_DIR, CM_SCORE_MAX (1.8)

suppressPackageStartupMessages({
  library(SoupX)
  library(DropletUtils)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

scratch     <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
raw_dir     <- Sys.getenv("RAW_DIR",     file.path(scratch, "FSU_seq/02_cellranger/shared_h5"))
plateau_csv <- Sys.getenv("PLATEAU_CSV", file.path(scratch, "FSU_seq/work_dir/emptydrops_mm10/plateau.csv"))
annot_dir   <- Sys.getenv("ANNOT_DIR",   file.path(scratch, "FSU_seq/work_dir/annot_mm10"))
out_dir     <- Sys.getenv("OUT_DIR",     file.path(scratch, "FSU_seq/work_dir/soupx_mm10"))
cm_max      <- as.numeric(Sys.getenv("CM_SCORE_MAX", "1.8"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args    <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else LETTERS[1:8]

CM_GENES <- c("Myh6", "Tnnt2", "Ryr2", "Ttn", "Ankrd1", "Mhrt", "Rbm20", "Myl2", "Myl4", "Nppa")

plateau <- fread(plateau_csv)
cells   <- fread(file.path(annot_dir, "cluster_clean_per_cell.csv.gz"))
pscore  <- fread(file.path(annot_dir, "cluster_panel_scores.csv"))
nonCM_clusters <- as.character(pscore[Cardiomyocyte < cm_max, cluster])
message(sprintf("non-CM clusters (CM panel score < %g): %s", cm_max, paste(nonCM_clusters, collapse = ",")))
message(sprintf("CM clusters: %s", paste(setdiff(as.character(pscore$cluster), nonCM_clusters), collapse = ",")))

res <- list(); post <- list(); topsoup <- list()
for (s in samples) {
  t0 <- Sys.time()
  q75 <- plateau[sample == s, plateau_q75]
  if (length(q75) != 1) stop("plateau_q75 missing for ", s)

  sce <- read10xCounts(file.path(raw_dir, s, "outs", "raw_feature_bc_matrix.h5"), col.names = TRUE)
  tod <- as(counts(sce), "dgCMatrix")
  rownames(tod) <- make.unique(as.character(rowData(sce)$Symbol))   # same gene names as 07_seurat_merge.R
  rm(sce); invisible(gc())
  tot <- colSums(tod)
  tod <- tod[, tot > 0]; tot <- tot[tot > 0]

  cs <- cells[sample == s]
  cs[, bc := sub(paste0("^", s, "_"), "", cell)]
  miss <- sum(!cs$bc %in% colnames(tod))
  if (miss > 0) stop(sprintf("[%s] %d clustered barcodes not in raw h5", s, miss))
  toc <- tod[, cs$bc]

  # soup profile from barcodes with total UMI <= plateau_q75
  soup_bc <- tot <= q75
  soup_cnt <- rowSums(tod[, soup_bc, drop = FALSE])
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- setSoupProfile(sc, data.frame(est = soup_cnt / sum(soup_cnt), counts = soup_cnt,
                                      row.names = rownames(tod)))
  rm(tod); invisible(gc())
  sc <- setClusters(sc, setNames(as.character(cs$cluster_clean), cs$bc))
  message(sprintf("[%s] nuclei %d; soup barcodes (UMI <= %g) %d, soup UMI %d",
                  s, ncol(toc), q75, sum(soup_bc), sum(soup_cnt)))

  # ---- a. autoEstCont ----
  auto <- tryCatch(autoEstCont(sc, doPlot = FALSE, verbose = TRUE), error = function(e) e)
  auto_ok <- !inherits(auto, "error")
  if (auto_ok) {
    fit <- auto$fit
    rho_auto <- fit$rhoEst; fwhm <- fit$rhoFWHM; n_mark <- nrow(fit$markersUsed)
    if (!is.null(fit$posterior))
      post[[s]] <- data.table(sample = s, rho = seq(0, 1, 0.001), posterior = fit$posterior / max(fit$posterior))
    auto_msg <- "ok"
  } else {
    rho_auto <- NA_real_; fwhm <- c(NA_real_, NA_real_); n_mark <- NA_integer_
    auto_msg <- gsub("[\r\n,]+", " ", conditionMessage(auto))
  }

  # ---- b. fallback: CM genes in non-CM clusters ----
  g <- intersect(CM_GENES, rownames(toc))
  use <- cs$bc[as.character(cs$cluster_clean) %in% nonCM_clusters]
  cm_umi <- colSums(toc[g, use, drop = FALSE])
  soup_frac_cm <- sum(sc$soupProfile[g, "est"])
  rho_fb <- sum(cm_umi) / sum(sc$metaData[use, "nUMIs"] * soup_frac_cm)
  det_rate <- mean(cm_umi > 0)

  rho_final <- if (auto_ok) rho_auto else rho_fb
  method <- if (auto_ok) "autoEstCont" else "CM-gene fallback"
  message(sprintf("[%s] rho autoEstCont = %s (FWHM %s-%s, %s markers; %s) | CM-gene fallback = %.4f (non-CM nuclei %d, CM-gene detection %.3f) | final = %.4f [%s]",
                  s, format(rho_auto, digits = 3), format(fwhm[1], digits = 3), format(fwhm[2], digits = 3),
                  n_mark, auto_msg, rho_fb, length(use), det_rate, rho_final, method))

  # ---- adjust counts with rho_final ----
  sc <- setContaminationFraction(sc, rho_final, forceAccept = TRUE)
  adj <- adjustCounts(sc, roundToInt = TRUE, verbose = 0)
  colnames(adj) <- paste0(s, "_", colnames(adj))
  saveRDS(adj, file.path(out_dir, paste0(s, "_soupx_counts.rds")))
  frac_removed <- 1 - sum(adj) / sum(toc)

  sp <- sc$soupProfile[order(-sc$soupProfile$est), ][1:20, ]
  topsoup[[s]] <- data.table(sample = s, rank = 1:20, gene = rownames(sp), soup_frac = signif(sp$est, 4))

  res[[s]] <- data.table(sample = s, n_nuclei = ncol(toc), plateau_q75 = q75,
                         n_soup_barcodes = sum(soup_bc), soup_umi = sum(soup_cnt),
                         rho_auto = round(rho_auto, 4), rho_auto_fwhm_low = fwhm[1], rho_auto_fwhm_high = fwhm[2],
                         n_markers_auto = n_mark, auto_status = auto_msg,
                         rho_cm_fallback = round(rho_fb, 4), n_nonCM_nuclei = length(use),
                         cm_gene_detect_rate_nonCM = round(det_rate, 4), cm_soup_frac = signif(soup_frac_cm, 4),
                         rho_final = round(rho_final, 4), method = method,
                         frac_umi_removed = round(frac_removed, 4))
  rm(sc, toc, adj); invisible(gc())
  message(sprintf("[%s] done in %.1f min", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

rho <- rbindlist(res)
fwrite(rho, file.path(out_dir, "soupx_rho.csv"))
fwrite(rbindlist(topsoup), file.path(out_dir, "soupx_top_soup_genes.csv"))
cat("\n==== SoupX rho per sample ====\n")
print(rho[, .(sample, n_nuclei, rho_auto, rho_auto_fwhm_low, rho_auto_fwhm_high, rho_cm_fallback,
              cm_gene_detect_rate_nonCM, rho_final, method, frac_umi_removed)])

# ---- plots ----
theme_set(theme_minimal(base_size = 13) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                  plot.title = element_text(size = 14), plot.subtitle = element_text(size = 11, colour = "#555555")))
COL <- c("autoEstCont" = "#257de0", "CM-gene fallback" = "#f47230")

pd <- rbind(rho[, .(sample, estimate = "autoEstCont", rho = rho_auto, lo = rho_auto_fwhm_low, hi = rho_auto_fwhm_high)],
            rho[, .(sample, estimate = "CM-gene fallback", rho = rho_cm_fallback, lo = NA_real_, hi = NA_real_)])
meth <- setNames(rho$method, rho$sample)          # (pd has its own column "rho")
pd[, used := estimate == meth[sample]]
pd <- pd[!is.na(rho)]
xpos <- setNames(seq_len(nrow(rho)), rho$sample)
pd[, x := xpos[sample] + ifelse(estimate == "autoEstCont", -0.12, 0.12)]
p1 <- ggplot(pd, aes(x, rho, colour = estimate)) +
  geom_hline(yintercept = median(rho$rho_final), colour = "#999999", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.6, na.rm = TRUE) +
  geom_point(aes(shape = used), size = 3.2, stroke = 1.1) +
  scale_x_continuous(breaks = seq_len(nrow(rho)), labels = rho$sample) +
  scale_colour_manual(values = COL, name = "Estimate") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "used (rho_final)", `FALSE` = "not used"),
                     name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "SoupX contamination fraction (rho) per sample",
       subtitle = sprintf("Error bar = autoEstCont FWHM; dashed line = median rho_final (%.3f); soup = barcodes with UMI <= plateau_q75",
                          median(rho$rho_final)),
       x = "Sample", y = "rho")
ggsave(file.path(out_dir, "soupx_rho_by_sample.png"), p1, width = 9, height = 5.5, dpi = 200, bg = "white")

if (length(post)) {
  pp <- rbindlist(post)
  xmax <- min(1, max(0.2, max(rho$rho_final, na.rm = TRUE) * 2.5))
  p2 <- ggplot(pp[rho <= xmax], aes(rho, posterior)) +
    geom_area(fill = "#257de0", alpha = 0.18) +
    geom_line(colour = "#257de0", linewidth = 0.7) +
    geom_vline(data = rho, aes(xintercept = rho_final), colour = "#0d2f63", linewidth = 0.5) +
    geom_vline(data = rho, aes(xintercept = rho_cm_fallback), colour = "#f47230", linewidth = 0.5, linetype = "dashed") +
    facet_wrap(~sample, ncol = 4) +
    labs(title = "autoEstCont posterior of rho per sample",
         subtitle = "Curve = posterior (scaled to max 1); solid line = rho_final; dashed orange = CM-gene fallback",
         x = "rho", y = "Posterior (scaled)") +
    theme(panel.grid.major.x = element_line(colour = "#eeeeee"))
  ggsave(file.path(out_dir, "soupx_rho_posterior.png"), p2, width = 11, height = 5.5, dpi = 200, bg = "white")
}
message("All done.")
