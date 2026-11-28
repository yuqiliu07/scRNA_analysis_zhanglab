#!/usr/bin/env Rscript
# Step 16: drop low-quality clusters from the CellBender clustering (step 15) and report, for every
# kept cluster, the top two marker panels (panel_1, panel_2).
#
# Panel score = mean over the panel genes of the log-normalised expression (Seurat data layer,
# LogNormalize, scale 1e4), averaged over the nuclei of the cluster. This is exactly what step 15
# wrote to cluster_panel_scores.csv with the same 17 PANELS. LogNormalize works on each nucleus
# on its own, so removing other clusters does not change the scores of the kept clusters, and
# this step only needs the step-15 tables (it runs on the login node, no rds is loaded).
#
# Low-quality clusters (default DROP_CLUSTERS = 14,26,30), from step-15 cluster_summary.csv:
#   14  median nGene 193 (< 200), 60% of UMI removed by CellBender: CM fragments left with ambient only
#   26  median pct_mt 1.48 (other clusters <= 0.84), median nGene 271, mixed origin
#   30  median pct_mt 1.90, 59% of nuclei from sample D (+ C): sample-specific technical cluster
# CM+X clusters (likely doublets) are kept here and handled later.
#
# Inputs : <CB_DIR>/cluster_panel_scores.csv, <CB_DIR>/cluster_summary.csv
# Outputs: <CB_DIR>/clean/  cluster_panels_clean.csv, cluster_panel_heatmap_clean.png
#          (uploaded to <GDRIVE_ROOT>/results/harmony_cellbender and plots/harmony_cellbender)
# Usage  : Rscript 16_cb_clean_panels.R      Env: CB_DIR, DROP_CLUSTERS, GDRIVE_ROOT (none = no upload)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

cb_dir <- Sys.getenv("CB_DIR", file.path(Sys.getenv("SCRATCH"), "FSU_seq/work_dir/harmony_cb_mm10"))
drop <- strsplit(Sys.getenv("DROP_CLUSTERS", "14,26,30"), ",")[[1]]
gdrive <- Sys.getenv("GDRIVE_ROOT", "gdrive:FSU_seq_results")
out_dir <- file.path(cb_dir, "clean")
dir.create(out_dir, showWarnings = FALSE)

ps <- fread(file.path(cb_dir, "cluster_panel_scores.csv"), colClasses = list(character = "cluster"))
sm <- fread(file.path(cb_dir, "cluster_summary.csv"), colClasses = list(character = "cluster"))
stopifnot(all(drop %in% ps$cluster))
pan <- setdiff(names(ps), "cluster")
cat("dropped clusters:", paste(drop, collapse = ", "), "|", sum(sm[cluster %in% drop, n_nuclei]), "nuclei\n")

ps <- ps[!cluster %in% drop]
m <- as.matrix(ps[, ..pan]); rownames(m) <- ps$cluster
o <- t(apply(m, 1, order, decreasing = TRUE))
res <- data.table(cluster = ps$cluster,
                  panel_1 = pan[o[, 1]], panel_1_score = round(m[cbind(seq_len(nrow(m)), o[, 1])], 3),
                  panel_2 = pan[o[, 2]], panel_2_score = round(m[cbind(seq_len(nrow(m)), o[, 2])], 3))
res[, panel_ratio := round(panel_1_score / pmax(panel_2_score, 1e-6), 2)]
res <- merge(sm[, .(cluster, n_nuclei, median_nGene_cb, median_pct_mt, mean_doublet_score, main_old_celltype)],
             res, by = "cluster")
res <- res[order(match(panel_1, pan), -n_nuclei)]
fwrite(res, file.path(out_dir, "cluster_panels_clean.csv"))
print(res[, .(cluster, n_nuclei, panel_1, panel_1_score, panel_2, panel_2_score, panel_ratio, main_old_celltype)], nrows = 100)
cat("\nnuclei kept:", sum(res$n_nuclei), "in", nrow(res), "clusters\n")
print(res[, .(clusters = .N, nuclei = sum(n_nuclei)), by = panel_1])

# heatmap: numbers = mean log-expression; fill = score / panel max (max floored at 0.5 so that
# near-zero panels such as Neutrophil stay pale); panel_1 boxed
hm <- melt(ps, id.vars = "cluster", variable.name = "panel", value.name = "score")
hm[, scaled := pmin(score / max(max(score), 0.5), 1), by = panel]
hm[, `:=`(cluster = factor(cluster, levels = rev(res$cluster)), panel = factor(panel, levels = pan))]
top <- merge(hm, res[, .(cluster, panel = panel_1)], by = c("cluster", "panel"))
p <- ggplot(hm, aes(panel, cluster, fill = scaled)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.1, colour = ifelse(hm$scaled > 0.6, "white", "grey20")) +
  geom_tile(data = top, fill = NA, colour = "black", linewidth = 0.6) +
  scale_fill_gradient(low = "#f4f8fd", high = "#0b3d78", name = "score /\npanel max") +
  labs(x = NULL, y = "cluster (CellBender)",
       title = sprintf("Marker panel score per cluster, low-quality clusters %s removed", paste(drop, collapse = "/")),
       subtitle = "numbers = mean log-normalised expression of the panel genes; box = panel_1") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(out_dir, "cluster_panel_heatmap_clean.png"), p, width = 10, height = 0.22 * nrow(res) + 2, dpi = 150)

if (gdrive != "none") {
  run <- function(...) stopifnot(system2("rclone", c(...)) == 0)
  run("copy", out_dir, file.path(gdrive, "results/harmony_cellbender"), "--include", "cluster_panels_clean.csv")
  run("check", out_dir, file.path(gdrive, "results/harmony_cellbender"), "--one-way", "--include", "cluster_panels_clean.csv")
  run("copy", out_dir, file.path(gdrive, "plots/harmony_cellbender"), "--include", "cluster_panel_heatmap_clean.png")
  run("check", out_dir, file.path(gdrive, "plots/harmony_cellbender"), "--one-way", "--include", "cluster_panel_heatmap_clean.png")
}
cat("done ->", out_dir, "\n")
