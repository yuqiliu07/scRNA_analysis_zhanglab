#!/usr/bin/env Rscript
# Clean up the Leiden clusters of 09_harmony_cluster.R and score cell-type marker panels per cluster.
#
#  1. Read <HARMONY_DIR>/fsu_merged_mm10_harmony.rds (obj$leiden_res1, harmony kNN graph "harmony_nn")
#  2. Cluster cleanup: clusters with < MIN_CLUSTER_SIZE nuclei (default 10; the singleton clusters
#     that igraph Leiden leaves for nuclei without SNN edges) are dissolved; each of their nuclei
#     gets the most common cluster among its k = 20 harmony neighbours (kNN graph, ties -> larger
#     cluster). Kept cluster ids are unchanged, so they match earlier tables. -> obj$leiden_clean
#  3. Marker panels = PANELS of sxj_code/04_merge_cluster.py (original pipeline):
#     per nucleus, panel score = mean log-normalised expression ("data" layers) of the panel genes
#     present after the gene filter; per cluster, mean of the nucleus scores (same as the original
#     cluster_review). panel_1 / panel_2 = highest / second-highest panel per cluster.
#  4. Figures
#     umap_cluster_clean.png        UMAP (harmony) coloured + labelled by cleaned Leiden cluster
#     umap_doublet_score.png        UMAP coloured by scDblFinder score (05_qc_doublet.R)
#     umap_pct_mt.png               UMAP coloured by % mitochondrial UMI (05_qc_doublet.R)
#     cluster_median_ngene.png      median genes per nucleus per cluster (bar plot)
#     cluster_celltype_heatmap.png  cluster x panel, mean panel score (top panel per cluster outlined)
#     cluster_celltype_heatmap_scaled.png  same, each panel min-max scaled across clusters
#  Tables: cluster_cleanup.csv, cluster_summary.csv (size, QC medians, sample spread, panel_1/2),
#          cluster_panel_scores.csv, marker_panels_used.csv
#  Data (stays on Hoffman): fsu_merged_mm10_harmony_clean.rds, cluster_clean_per_cell.csv.gz
#
# Usage: Rscript 11_cluster_clean_annot.R
# Env (optional): HARMONY_DIR, OUT_DIR, MIN_CLUSTER_SIZE (10)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})
options(Seurat.object.assay.version = "v5")

scratch     <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
harmony_dir <- Sys.getenv("HARMONY_DIR", file.path(scratch, "FSU_seq/work_dir/harmony_mm10"))
out_dir     <- Sys.getenv("OUT_DIR",     file.path(scratch, "FSU_seq/work_dir/annot_mm10"))
min_size    <- as.integer(Sys.getenv("MIN_CLUSTER_SIZE", "10"))
SEED        <- 42
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# PANELS from sxj_code/04_merge_cluster.py (17 panels)
PANELS <- list(
  Cardiomyocyte  = c("Myh6", "Tnnt2", "Ryr2", "Ttn", "Ankrd1", "Mhrt", "Rbm20", "Myl2", "Myl4", "Nppa"),
  Fibroblast     = c("Dcn", "Col1a1", "Col1a2", "Pdgfra", "Gsn", "Lama2", "Abca8a", "Bicc1"),
  Endothelial    = c("Pecam1", "Cdh5", "Flt1", "Egfl7", "Kdr", "Ptprb", "Emcn"),
  Endocardial    = c("Npr3", "Cdh11", "Vwf"),
  Lymphatic_EC   = c("Prox1", "Lyve1", "Mmrn1", "Reln"),
  Pericyte       = c("Rgs5", "Kcnj8", "Abcc9", "Pdgfrb", "Vtn"),
  SMC            = c("Myh11", "Acta2", "Cnn1", "Tagln", "Lmod1"),
  Macrophage     = c("C1qa", "C1qb", "F13a1", "Mrc1", "Cd163", "Adgre1", "Csf1r"),
  Neutrophil     = c("S100a8", "S100a9", "Retnlg"),
  T_NK           = c("Cd3e", "Cd3g", "Skap1", "Nkg7", "Il7r"),
  B              = c("Cd79a", "Ms4a1", "Pax5", "Bank1"),
  Epicardial     = c("Msln", "Upk3b", "Wt1", "Tcf21"),
  Adipocyte      = c("Adipoq", "Plin1", "Lep", "Cidec"),
  Neural_Schwann = c("Plp1", "Nrxn1", "Scn7a", "Kcna1"),
  Hepatocyte     = c("Alb", "Apoa1", "Cyp2e1"),
  Kidney_epi     = c("Slc34a1", "Lrp2", "Umod", "Nphs1"),
  Lung_epi       = c("Sftpc", "Ager", "Scgb1a1")
)

obj <- readRDS(file.path(harmony_dir, "fsu_merged_mm10_harmony.rds"))
message(sprintf("%d genes x %d nuclei; layers: %s", nrow(obj), ncol(obj), paste(Layers(obj), collapse = ",")))
num_lv <- function(x) as.character(sort(as.integer(unique(as.character(x)))))

# ---- 2. cluster cleanup ----
cl <- as.character(obj$leiden_res1)
sizes <- table(cl)
small <- names(sizes)[sizes < min_size]
nn <- obj[["harmony_nn"]]                               # kNN graph (cells x cells, 1 = neighbour)
clean <- cl
reassign <- list()
for (i in which(cl %in% small)) {
  nb <- setdiff(which(nn[i, ] > 0), i)
  nb_cl <- cl[nb][!cl[nb] %in% small]
  if (!length(nb_cl)) stop("nucleus ", colnames(obj)[i], " has no neighbour in a kept cluster")
  tt <- table(nb_cl)
  best <- names(tt)[tt == max(tt)]
  best <- best[which.max(sizes[best])]                  # tie -> larger cluster
  clean[i] <- best
  reassign[[length(reassign) + 1]] <- data.table(cell = colnames(obj)[i], sample = obj$orig.ident[i],
                                                  old_cluster = cl[i], old_cluster_size = as.integer(sizes[cl[i]]),
                                                  new_cluster = best, n_neighbors = length(nb),
                                                  frac_neighbors_new = round(max(tt) / length(nb), 3))
}
reassign <- rbindlist(reassign)
fwrite(reassign, file.path(out_dir, "cluster_cleanup.csv"))
obj$leiden_clean <- factor(clean, levels = num_lv(clean))
Idents(obj) <- "leiden_clean"
message(sprintf("cleanup: %d clusters < %d nuclei dissolved (%d nuclei reassigned); %d -> %d clusters",
                length(small), min_size, nrow(reassign), length(sizes), nlevels(obj$leiden_clean)))
if (nrow(reassign)) print(reassign)

# ---- 3. panel scores ----
genes_used <- lapply(PANELS, function(g) intersect(g, rownames(obj)))
pan_tab <- rbindlist(lapply(names(PANELS), function(p)
  data.table(panel = p, n_genes_panel = length(PANELS[[p]]), n_genes_used = length(genes_used[[p]]),
             genes_used = paste(genes_used[[p]], collapse = ";"),
             genes_missing = paste(setdiff(PANELS[[p]], genes_used[[p]]), collapse = ";"))))
fwrite(pan_tab, file.path(out_dir, "marker_panels_used.csv"))
print(pan_tab[, .(panel, n_genes_panel, n_genes_used, genes_missing)])
genes_used <- genes_used[lengths(genes_used) > 0]
all_g <- unique(unlist(genes_used))

score <- matrix(NA_real_, ncol(obj), length(genes_used), dimnames = list(colnames(obj), names(genes_used)))
for (l in Layers(obj, search = "data")) {
  x <- LayerData(obj, layer = l, features = all_g)       # genes x cells of this sample
  x <- as(x, "CsparseMatrix")
  for (p in names(genes_used)) score[colnames(x), p] <- Matrix::colMeans(x[genes_used[[p]], , drop = FALSE])
}
stopifnot(!anyNA(score))
for (p in colnames(score)) obj[[paste0("panel_", p)]] <- score[, p]

clv <- levels(obj$leiden_clean)
ps <- t(sapply(clv, function(k) colMeans(score[obj$leiden_clean == k, , drop = FALSE])))
dimnames(ps) <- list(cluster = clv, panel = colnames(score))
ps_dt <- data.table(cluster = clv, round(ps, 4))
fwrite(ps_dt, file.path(out_dir, "cluster_panel_scores.csv"))

md <- data.table(cell = colnames(obj), sample = as.character(obj$orig.ident),
                 cluster = as.character(obj$leiden_clean),
                 nGene = obj$qc_nGene, nUMI = obj$qc_nUMI, pct_mt = obj$qc_pct_mt,
                 dbl = obj$scDblFinder_score)
top2 <- t(apply(ps, 1, function(r) { o <- order(-r); c(names(r)[o[1]], r[o[1]], names(r)[o[2]], r[o[2]]) }))
summ <- md[, .(n_nuclei = .N, median_nGene = median(nGene), median_nUMI = median(nUMI),
               median_pct_mt = round(median(pct_mt), 3), mean_doublet_score = round(mean(dbl), 4),
               n_samples_ge1pct = sum(table(factor(sample, levels = LETTERS[1:8])) / .N >= 0.01),
               top_sample = names(which.max(table(sample))),
               top_sample_frac = round(max(table(sample)) / .N, 3)), by = cluster]
summ <- merge(summ, data.table(cluster = rownames(top2), panel_1 = top2[, 1],
                               panel_1_score = round(as.numeric(top2[, 2]), 3), panel_2 = top2[, 3],
                               panel_2_score = round(as.numeric(top2[, 4]), 3)), by = "cluster")
summ[, panel_ratio := round(panel_1_score / pmax(panel_2_score, 1e-6), 2)]
summ <- summ[order(as.integer(cluster))]
fwrite(summ, file.path(out_dir, "cluster_summary.csv"))

saveRDS(obj, file.path(out_dir, "fsu_merged_mm10_harmony_clean.rds"))
um <- Embeddings(obj, "umap")
pc <- cbind(md, UMAP_1 = um[, 1], UMAP_2 = um[, 2])
fwrite(pc[, .(cell, sample, cluster_raw = cl, cluster_clean = cluster,
              UMAP_1 = round(UMAP_1, 4), UMAP_2 = round(UMAP_2, 4))],
       file.path(out_dir, "cluster_clean_per_cell.csv.gz"))

# ---- 4. figures ----
theme_set(theme_classic(base_size = 11) +
          theme(axis.line = element_line(colour = "grey60", linewidth = 0.3),
                axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
                axis.text = element_text(colour = "grey30")))
n_lab <- format(nrow(pc), big.mark = ",")
k <- length(clv)
cl_pal <- setNames(grDevices::hcl(h = seq(15, 375, length.out = k + 1)[1:k], c = 70, l = 62), clv)
set.seed(SEED)
pd <- pc[sample(.N)]
cent <- pc[, .(x = median(UMAP_1), y = median(UMAP_2)), by = cluster]
p1 <- ggplot(pd, aes(UMAP_1, UMAP_2, colour = factor(cluster, levels = clv))) +
  geom_point(size = 0.06, alpha = 0.5, stroke = 0, show.legend = FALSE) +
  geom_label(data = cent, aes(x, y, label = cluster), inherit.aes = FALSE, size = 3,
             linewidth = 0, fill = alpha("white", 0.8), label.padding = unit(0.12, "lines")) +
  scale_colour_manual(values = cl_pal) + coord_equal() +
  labs(title = sprintf("Leiden clusters after cleanup (%d clusters, %s nuclei)", k, n_lab),
       subtitle = sprintf("Harmony UMAP; clusters < %d nuclei merged into their neighbours' cluster", min_size))
ggsave(file.path(out_dir, "umap_cluster_clean.png"), p1, width = 9, height = 8, dpi = 200)

seq_umap <- function(var, title, legend, cap = NULL) {
  d <- copy(pc)[, v := get(var)]
  if (!is.null(cap)) d[, v := pmin(v, cap)]
  d <- d[order(v)]                                      # high values drawn on top
  ggplot(d, aes(UMAP_1, UMAP_2, colour = v)) +
    geom_point(size = 0.06, stroke = 0) +
    scale_colour_gradientn(colours = c("#e6e6e6", "#b6d5f8", "#6aabf1", "#257de0", "#1757ac", "#0d2f63"),
                           name = legend) +
    geom_label(data = cent, aes(x, y, label = cluster), inherit.aes = FALSE, size = 2.6,
               linewidth = 0, fill = alpha("white", 0.7), label.padding = unit(0.1, "lines")) +
    coord_equal() + labs(title = title)
}
ggsave(file.path(out_dir, "umap_doublet_score.png"),
       seq_umap("dbl", "scDblFinder doublet score (all nuclei here were called singlets)", "score"),
       width = 9.5, height = 8, dpi = 200)
mt_cap <- quantile(pc$pct_mt, 0.99)
ggsave(file.path(out_dir, "umap_pct_mt.png"),
       seq_umap("pct_mt", sprintf("%% mitochondrial UMI (colour capped at 99th percentile = %.2f%%)", mt_cap),
                "% mt", cap = mt_cap),
       width = 9.5, height = 8, dpi = 200)

bd <- summ[, .(cluster = factor(cluster, levels = clv), median_nGene, panel_1, n_nuclei)]
xlab_map <- setNames(sprintf("%s  %s", bd$cluster, bd$panel_1), as.character(bd$cluster))
p4 <- ggplot(bd, aes(cluster, median_nGene)) +
  geom_col(fill = "#257de0", width = 0.75) +
  geom_hline(yintercept = median(pc$nGene), linetype = 2, colour = "grey40") +
  geom_text(aes(label = median_nGene), vjust = -0.3, size = 2.4, colour = "grey25") +
  annotate("text", x = k + 0.4, y = median(pc$nGene), label = sprintf("all nuclei: %d", median(pc$nGene)),
           hjust = 1, vjust = -0.5, size = 3, colour = "grey30") +
  scale_x_discrete(labels = xlab_map) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Leiden cluster (cleaned) and its top marker panel", y = "median genes per nucleus",
       title = "Median nGene per cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
ggsave(file.path(out_dir, "cluster_median_ngene.png"), p4, width = 14, height = 6, dpi = 170)

# heatmaps: clusters ordered by top panel (panel order as in PANELS), then by score
pan_lv <- colnames(ps)
ord <- summ[, .(cluster, p1 = factor(panel_1, levels = pan_lv), panel_1_score)][order(p1, -panel_1_score), cluster]
hm <- melt(ps_dt, id.vars = "cluster", variable.name = "panel", value.name = "score")
hm[, scaled := (score - min(score)) / max(max(score) - min(score), 1e-9), by = panel]
hm[, cluster := factor(cluster, levels = rev(ord))]
hm[, panel := factor(panel, levels = pan_lv)]
top <- merge(hm, summ[, .(cluster, panel_1)], by = "cluster")[panel == panel_1]
ylab_map <- setNames(sprintf("%s  (%s nuclei, %s genes)", summ$cluster,
                             format(summ$n_nuclei, big.mark = ","), summ$median_nGene), summ$cluster)
heat <- function(fill_var, legend, title, fmt_var) {
  ggplot(hm, aes(panel, cluster, fill = .data[[fill_var]])) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", .data[[fmt_var]]),
                  colour = .data[[fill_var]] > 0.55 * max(.data[[fill_var]])), size = 2.2) +
    geom_tile(data = top, fill = NA, colour = "black", linewidth = 0.7) +
    scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "white"), guide = "none") +
    scale_fill_gradient(low = "#f4f8fd", high = "#1757ac", name = legend) +
    scale_y_discrete(labels = ylab_map) +
    labs(x = NULL, y = "Leiden cluster (median genes)", title = title,
         subtitle = "Black box = top panel of each cluster; clusters ordered by top panel") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7.5))
}
ggsave(file.path(out_dir, "cluster_celltype_heatmap.png"),
       heat("score", "mean\nlog-expr", "Marker panel score per cluster (mean log-normalised expression)", "score"),
       width = 13, height = 11, dpi = 170)
ggsave(file.path(out_dir, "cluster_celltype_heatmap_scaled.png"),
       heat("scaled", "scaled\n(per panel)",
            "Marker panel score, each panel scaled 0-1 across clusters (numbers = raw score)", "score"),
       width = 13, height = 11, dpi = 170)

cat("\n===== cluster summary =====\n")
print(summ[, .(cluster, n_nuclei, median_nGene, median_pct_mt, mean_doublet_score, top_sample, top_sample_frac,
               panel_1, panel_1_score, panel_2, panel_2_score, panel_ratio)])
message("done")
