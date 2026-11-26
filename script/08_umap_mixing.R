#!/usr/bin/env Rscript
# Uncorrected embedding (no Harmony / integration) of the merged FSU_seq object
# and a per-cell sample-mixing score, to check for sample/batch effects.
#
#  1. Read <SEURAT_DIR>/fsu_merged_mm10.rds (07_seurat_merge.R output;
#     split layers counts.A..H / data.A..H, VariableFeatures = per-sample HVG union)
#  2. ScaleData on the HVG union (all layers stitched -> scaled across all cells)
#  3. RunPCA (npcs = N_PCS_COMPUTE), RunUMAP on PC 1..N_PCS (default 20)
#  4. Mixing score: kNN (k = K_NN, default 30, self excluded) in the same PC 1..N_PCS
#     space; for each cell
#         same_frac = fraction of its k neighbours from the same sample
#         expected  = n_cells(sample) / n_cells(total)
#         ratio     = same_frac / expected   (1 = perfectly mixed, max = 1/expected)
#
# Outputs
#   <OUT_DIR>/ (data, stay on Hoffman)
#     fsu_merged_mm10_pca_umap.rds   object + pca/umap reductions (scale.data dropped)
#     umap_mixing_per_cell.csv.gz    cell, sample, UMAP_1/2, same_frac, expected, ratio
#   <OUT_DIR>/ (summaries + plots, uploaded by 08_umap_mixing.sh)
#     umap_by_sample.png             UMAP coloured by sample
#     umap_by_sample_split.png       one panel per sample (sample highlighted on grey)
#     mixing_score.png               same-sample neighbour fraction vs expected,
#                                    log2 ratio, neighbour-composition heatmap
#     pca_elbow.png                  PC standard deviations
#     mixing_summary_per_sample.csv  per sample: expected, observed median/mean, ratio
#     mixing_neighbor_composition.csv rows = cell's sample, cols = neighbour sample
#
# Usage: Rscript 08_umap_mixing.R
# Env (optional): SEURAT_DIR, OUT_DIR, N_PCS (20), N_PCS_COMPUTE (30), K_NN (30)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
options(Seurat.object.assay.version = "v5")

scratch    <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
seurat_dir <- Sys.getenv("SEURAT_DIR", file.path(scratch, "FSU_seq/work_dir/seurat_mm10"))
out_dir    <- Sys.getenv("OUT_DIR",    file.path(scratch, "FSU_seq/work_dir/umap_mm10"))
n_pcs      <- as.integer(Sys.getenv("N_PCS", "20"))
n_pcs_comp <- max(n_pcs, as.integer(Sys.getenv("N_PCS_COMPUTE", "30")))
k_nn       <- as.integer(Sys.getenv("K_NN", "30"))
SEED       <- 42
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message(sprintf("SEURAT_DIR=%s | OUT_DIR=%s | N_PCS=%d (compute %d) | K_NN=%d",
                seurat_dir, out_dir, n_pcs, n_pcs_comp, k_nn))

obj <- readRDS(file.path(seurat_dir, "fsu_merged_mm10.rds"))
samples <- sort(unique(as.character(obj$orig.ident)))
hvg <- VariableFeatures(obj)
message(sprintf("%d genes x %d nuclei; %d HVGs; layers: %s", nrow(obj), ncol(obj), length(hvg),
                paste(Layers(obj), collapse = ",")))

# ---- PCA / UMAP (no integration) ----
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
set.seed(SEED)
obj <- RunPCA(obj, features = hvg, npcs = n_pcs_comp, seed.use = SEED, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, seed.use = SEED, verbose = FALSE)
obj[["RNA"]]$scale.data <- NULL                     # ~14 GB dense; recompute if needed
obj@misc$umap_params <- list(n_pcs = n_pcs, n_pcs_computed = n_pcs_comp, features = "hvg union",
                             integration = "none", umap_seed = SEED)

# ---- mixing score ----
pcs <- Embeddings(obj, "pca")[, 1:n_pcs]
smp <- as.character(obj$orig.ident)
nn  <- RANN::nn2(pcs, k = k_nn + 1)$nn.idx
# drop self (normally column 1; handle exact duplicates where self is not first)
self_col <- nn == seq_len(nrow(nn))
nn <- t(vapply(seq_len(nrow(nn)), function(i) {
  r <- nn[i, !self_col[i, ]]; r[seq_len(k_nn)]
}, integer(k_nn)))
nb_smp <- matrix(smp[nn], nrow = nrow(nn))
same_frac <- rowMeans(nb_smp == smp)
expected  <- as.numeric(table(smp)[smp]) / length(smp)

um <- Embeddings(obj, "umap")
per_cell <- data.table(cell = colnames(obj), sample = smp,
                       genotype = obj$genotype, sex = obj$sex,
                       UMAP_1 = round(um[, 1], 4), UMAP_2 = round(um[, 2], 4),
                       same_frac = round(same_frac, 4), expected = round(expected, 5),
                       ratio = round(same_frac / expected, 4))
obj$mix_same_frac <- same_frac
obj$mix_ratio <- same_frac / expected
saveRDS(obj, file.path(out_dir, "fsu_merged_mm10_pca_umap.rds"))
fwrite(per_cell, file.path(out_dir, "umap_mixing_per_cell.csv.gz"))

summ <- per_cell[, .(n_cells = .N, expected = round(expected[1], 4),
                     same_frac_median = round(median(same_frac), 4),
                     same_frac_mean = round(mean(same_frac), 4),
                     same_frac_q25 = round(quantile(same_frac, 0.25), 4),
                     same_frac_q75 = round(quantile(same_frac, 0.75), 4),
                     ratio_median = round(median(ratio), 3),
                     ratio_mean = round(mean(ratio), 3),
                     pct_all_neighbors_same = round(100 * mean(same_frac == 1), 2)),
                 by = .(sample, genotype, sex)][order(sample)]
summ <- rbind(summ, per_cell[, .(sample = "all", genotype = "", sex = "", n_cells = .N,
                                 expected = round(sum((table(smp) / .N)^2), 4),  # cell-weighted
                                 same_frac_median = round(median(same_frac), 4),
                                 same_frac_mean = round(mean(same_frac), 4),
                                 same_frac_q25 = round(quantile(same_frac, 0.25), 4),
                                 same_frac_q75 = round(quantile(same_frac, 0.75), 4),
                                 ratio_median = round(median(ratio), 3),
                                 ratio_mean = round(mean(ratio), 3),
                                 pct_all_neighbors_same = round(100 * mean(same_frac == 1), 2))])
fwrite(summ, file.path(out_dir, "mixing_summary_per_sample.csv"))

# neighbour composition: mean fraction of neighbours from each sample, per cell's sample
comp <- t(sapply(samples, function(s) {
  idx <- smp == s
  tabulate(match(nb_smp[idx, ], samples), nbins = length(samples)) / (sum(idx) * k_nn)
}))
dimnames(comp) <- list(sample = samples, neighbor_sample = samples)
comp_dt <- data.table(sample = samples, round(comp, 4))
fwrite(comp_dt, file.path(out_dir, "mixing_neighbor_composition.csv"))

# ---- plots ----
ss_lab <- unique(per_cell[, .(sample, genotype, sex)])[order(sample)]
lab <- setNames(sprintf("%s (%s, %s)", ss_lab$sample, ss_lab$genotype, ss_lab$sex), ss_lab$sample)
pal <- setNames(c("#257de0", "#f47230", "#00b894", "#f2a600",
                  "#ee7fac", "#008a17", "#6a4ce2", "#ec4647")[seq_along(samples)], samples)
theme_set(theme_classic(base_size = 11) +
          theme(axis.line = element_line(colour = "grey60", linewidth = 0.3),
                axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
                axis.text = element_text(colour = "grey30"),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold", hjust = 0)))

set.seed(SEED)
pd <- per_cell[sample(.N)]                           # random draw order, no sample on top
p_umap <- ggplot(pd, aes(UMAP_1, UMAP_2, colour = sample)) +
  geom_point(size = 0.05, alpha = 0.4, stroke = 0) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_equal() +
  labs(title = sprintf("UMAP, PC 1-%d, no integration (%s nuclei)", n_pcs,
                       format(nrow(pd), big.mark = ",")),
       subtitle = "Coloured by sample; see the split figure for each sample on its own")
ggsave(file.path(out_dir, "umap_by_sample.png"), p_umap, width = 9, height = 7.5, dpi = 200)

bg <- per_cell[, .(UMAP_1, UMAP_2)]
fg <- copy(per_cell)[, panel := factor(lab[sample], levels = lab)]
p_split <- ggplot(fg, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey85", size = 0.02, stroke = 0) +
  geom_point(aes(colour = sample), size = 0.08, alpha = 0.6, stroke = 0, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  facet_wrap(~panel, nrow = 2) + coord_equal() +
  labs(title = sprintf("UMAP by sample (PC 1-%d, no integration); grey = all other nuclei", n_pcs))
ggsave(file.path(out_dir, "umap_by_sample_split.png"), p_split, width = 16, height = 8.5, dpi = 170)

exp_dt <- summ[sample %in% samples, .(sample, expected)]
p1 <- ggplot(per_cell, aes(factor(sample, levels = samples), same_frac)) +
  geom_violin(aes(fill = sample), colour = NA, alpha = 0.55, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_errorbar(data = exp_dt, aes(x = sample, ymin = expected, ymax = expected),
                inherit.aes = FALSE, width = 0.7, linewidth = 0.8, colour = "black") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = sub(" \\(", "\n(", lab)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = sprintf("same-sample fraction of %d nearest neighbours", k_nn),
       title = "Observed same-sample neighbour fraction",
       subtitle = "Black bar = expected under perfect mixing (sample share of all nuclei)")
p2 <- ggplot(per_cell, aes(factor(sample, levels = samples), log2(pmax(ratio, 1 / k_nn)))) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_violin(aes(fill = sample), colour = NA, alpha = 0.55, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "log2(observed / expected)",
       title = "Mixing ratio", subtitle = "0 = perfectly mixed; higher = more sample-specific")
cd <- melt(comp_dt, id.vars = "sample", variable.name = "neighbor_sample", value.name = "frac")
p3 <- ggplot(cd, aes(neighbor_sample, factor(sample, levels = rev(samples)), fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", frac), colour = frac > 0.45), size = 3) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient(low = "#eef5fd", high = "#1757ac", limits = c(0, NA), name = "fraction") +
  labs(x = "neighbour's sample", y = "cell's sample",
       title = "Neighbour composition", subtitle = sprintf("Mean fraction of the %d neighbours", k_nn)) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
p_mix <- (p1 | p2) / (p3 | plot_spacer()) +
  plot_annotation(title = sprintf("Sample mixing, kNN k = %d in PC 1-%d space (no integration)", k_nn, n_pcs))
ggsave(file.path(out_dir, "mixing_score.png"), p_mix, width = 15, height = 11, dpi = 160)

sdev <- data.table(PC = seq_along(obj[["pca"]]@stdev), sd = obj[["pca"]]@stdev)
p_elb <- ggplot(sdev, aes(PC, sd)) + geom_point(size = 2, colour = "#257de0") +
  geom_vline(xintercept = n_pcs + 0.5, linetype = 2, colour = "grey40") +
  labs(y = "standard deviation", title = sprintf("PCA elbow (dashed line: PCs used = 1-%d)", n_pcs))
ggsave(file.path(out_dir, "pca_elbow.png"), p_elb, width = 6, height = 4, dpi = 150)

cat("\n===== mixing summary =====\n"); print(summ)
cat("\n===== neighbour composition =====\n"); print(round(comp, 3))
message("done")
