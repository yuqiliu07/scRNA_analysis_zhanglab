#!/usr/bin/env Rscript
# Harmony integration (by sample) + UMAP + Leiden clustering of the merged FSU_seq object.
#
#  1. Read <SEURAT_DIR>/fsu_merged_mm10.rds (07_seurat_merge.R; layers counts/data.A..H,
#     VariableFeatures = per-sample HVG union)
#  2. ScaleData on the HVG union -> RunPCA (N_PCS_COMPUTE = 50); use PC 1..N_PCS (30)
#  3. Before integration: RunUMAP on PCA 1..N_PCS (n.neighbors = N_NEIGHBORS)
#  4. RunHarmony (group.by.vars = orig.ident, dims 1..N_PCS, max_iter = MAX_ITER_HARMONY = 20)
#  5. On the harmony embedding (1..N_PCS):
#       RunUMAP (n.neighbors = N_NEIGHBORS = 20)
#       FindNeighbors (k.param = N_NEIGHBORS) -> SNN graph
#       Leiden on the SNN graph with igraph::cluster_leiden
#         (objective = modularity, resolution = RESOLUTION = 1.0, n_iterations = 2;
#          same settings as scanpy tl.leiden(flavor = "igraph", n_iterations = 2))
#
# Outputs (<OUT_DIR>/)
#   data, stay on Hoffman:
#     fsu_merged_mm10_harmony.rds     object with pca / umap_pca / harmony / umap reductions,
#                                     SNN graph, obj$leiden_res1 (scale.data dropped)
#     harmony_umap_per_cell.csv.gz    cell, sample, UMAP before/after, cluster
#   plots + summaries (uploaded by 09_harmony_cluster.sh):
#     umap_pre_harmony_by_sample.png  UMAP on PCA, coloured by sample
#     umap_harmony_by_sample.png      UMAP on harmony, coloured by sample
#     umap_harmony_by_cluster.png     UMAP on harmony, coloured + labelled by Leiden cluster
#     umap_pre_vs_post_by_sample_split.png  per-sample panels, before (top) vs after (bottom)
#     pca_elbow.png
#     cluster_sample_counts.csv       nuclei per cluster x sample (+ % of each sample)
#     harmony_params.csv
#
# Usage: Rscript 09_harmony_cluster.R
# Env (optional): SEURAT_DIR, OUT_DIR, N_PCS_COMPUTE (50), N_PCS (30), MAX_ITER_HARMONY (20),
#                 N_NEIGHBORS (20), RESOLUTION (1.0), LEIDEN_ITER (2)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(igraph)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})
options(Seurat.object.assay.version = "v5")

scratch     <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
seurat_dir  <- Sys.getenv("SEURAT_DIR", file.path(scratch, "FSU_seq/work_dir/seurat_mm10"))
out_dir     <- Sys.getenv("OUT_DIR",    file.path(scratch, "FSU_seq/work_dir/harmony_mm10"))
n_pcs_comp  <- as.integer(Sys.getenv("N_PCS_COMPUTE", "50"))
n_pcs       <- as.integer(Sys.getenv("N_PCS", "30"))
max_iter_h  <- as.integer(Sys.getenv("MAX_ITER_HARMONY", "20"))
n_nb        <- as.integer(Sys.getenv("N_NEIGHBORS", "20"))
resolution  <- as.numeric(Sys.getenv("RESOLUTION", "1.0"))
leiden_iter <- as.integer(Sys.getenv("LEIDEN_ITER", "2"))
SEED        <- 42
stopifnot(n_pcs <= n_pcs_comp)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
params <- data.table(harmony_version = as.character(packageVersion("harmony")),
                     igraph_version = as.character(packageVersion("igraph")),
                     seurat_version = as.character(packageVersion("Seurat")),
                     features = "per-sample HVG union", n_pcs_computed = n_pcs_comp, n_pcs_used = n_pcs,
                     harmony_group_by = "orig.ident (sample)", max_iter_harmony = max_iter_h,
                     n_neighbors = n_nb, snn_graph = "Seurat FindNeighbors SNN (prune 1/15)",
                     leiden = "igraph::cluster_leiden, objective = modularity",
                     resolution = resolution, leiden_n_iterations = leiden_iter, seed = SEED)
print(t(params))

obj <- readRDS(file.path(seurat_dir, "fsu_merged_mm10.rds"))
samples <- sort(unique(as.character(obj$orig.ident)))
hvg <- VariableFeatures(obj)
message(sprintf("%d genes x %d nuclei; %d HVGs; layers: %s", nrow(obj), ncol(obj), length(hvg),
                paste(Layers(obj), collapse = ",")))

# ---- PCA ----
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = n_pcs_comp, seed.use = SEED, verbose = FALSE)
obj[["RNA"]]$scale.data <- NULL                      # ~14 GB dense; not needed downstream
invisible(gc())

# ---- before integration ----
obj <- RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, n.neighbors = n_nb, seed.use = SEED,
               reduction.name = "umap_pca", reduction.key = "umapPCA_", verbose = FALSE)
message("UMAP (PCA) done")

# ---- Harmony ----
set.seed(SEED)
h_args <- list(object = obj, group.by.vars = "orig.ident", reduction.use = "pca",
               dims.use = 1:n_pcs, reduction.save = "harmony",
               project.dim = FALSE,                  # needs scale.data (dropped above); not used
               verbose = TRUE)
# harmony >= 1.0 uses max_iter; older versions used max.iter.harmony
h_formals <- c(names(formals(getS3method("RunHarmony", "Seurat", envir = asNamespace("harmony")))),
               names(formals(getS3method("RunHarmony", "default", envir = asNamespace("harmony")))))
if ("max_iter" %in% h_formals) h_args$max_iter <- max_iter_h else h_args$max.iter.harmony <- max_iter_h
obj <- do.call(RunHarmony, h_args)
message("Harmony done: ", paste(dim(Embeddings(obj, "harmony")), collapse = " x "))

# ---- after integration: UMAP, neighbours, Leiden ----
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs, n.neighbors = n_nb, seed.use = SEED,
               reduction.name = "umap", reduction.key = "umap_", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs, k.param = n_nb,
                     graph.name = c("harmony_nn", "harmony_snn"), verbose = FALSE)
g <- graph_from_adjacency_matrix(as(obj[["harmony_snn"]], "dgCMatrix"),
                                 mode = "undirected", weighted = TRUE, diag = FALSE)
set.seed(SEED)
cl <- cluster_leiden(g, objective_function = "modularity", weights = E(g)$weight,
                     resolution = resolution, n_iterations = leiden_iter)
memb <- membership(cl)
mod <- igraph::modularity(g, memb, weights = E(g)$weight, resolution = resolution)
# relabel clusters 0..K-1 by decreasing size (scanpy / Seurat convention)
ord <- order(-tabulate(memb))
memb <- match(memb, ord) - 1L
lvl <- as.character(sort(unique(memb)))
obj$leiden_res1 <- factor(as.character(memb), levels = lvl)
Idents(obj) <- "leiden_res1"
rm(g); invisible(gc())
message(sprintf("Leiden: %d clusters, modularity = %.4f", length(lvl), mod))
params[, `:=`(n_clusters = length(lvl), modularity = round(mod, 4))]
fwrite(params, file.path(out_dir, "harmony_params.csv"))

saveRDS(obj, file.path(out_dir, "fsu_merged_mm10_harmony.rds"))

u0 <- Embeddings(obj, "umap_pca"); u1 <- Embeddings(obj, "umap")
pc <- data.table(cell = colnames(obj), sample = as.character(obj$orig.ident),
                 genotype = obj$genotype, sex = obj$sex,
                 UMAP_pca_1 = round(u0[, 1], 4), UMAP_pca_2 = round(u0[, 2], 4),
                 UMAP_harmony_1 = round(u1[, 1], 4), UMAP_harmony_2 = round(u1[, 2], 4),
                 cluster = obj$leiden_res1)
fwrite(pc, file.path(out_dir, "harmony_umap_per_cell.csv.gz"))

ct <- dcast(pc[, .N, by = .(cluster, sample)], cluster ~ sample, value.var = "N", fill = 0)
ct[, total := rowSums(.SD), .SDcols = samples]
pct <- copy(ct)[, (samples) := lapply(.SD, function(v) round(100 * v / sum(v), 3)), .SDcols = samples]
setnames(pct, samples, paste0("pct_", samples))
ct <- cbind(ct, pct[, paste0("pct_", samples), with = FALSE])
ct <- ct[order(as.integer(as.character(cluster)))]
fwrite(ct, file.path(out_dir, "cluster_sample_counts.csv"))

# ---- plots ----
ss_lab <- unique(pc[, .(sample, genotype, sex)])[order(sample)]
lab <- setNames(sprintf("%s (%s, %s)", ss_lab$sample, ss_lab$genotype, ss_lab$sex), ss_lab$sample)
pal <- setNames(c("#257de0", "#f47230", "#00b894", "#f2a600",
                  "#ee7fac", "#008a17", "#6a4ce2", "#ec4647")[seq_along(samples)], samples)
theme_set(theme_classic(base_size = 11) +
          theme(axis.line = element_line(colour = "grey60", linewidth = 0.3),
                axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
                axis.text = element_text(colour = "grey30"),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold", hjust = 0)))
n_lab <- format(nrow(pc), big.mark = ",")
set.seed(SEED)
pd <- pc[sample(.N)]                                  # random draw order

umap_sample <- function(x, y, title) {
  ggplot(pd, aes(.data[[x]], .data[[y]], colour = sample)) +
    geom_point(size = 0.05, alpha = 0.4, stroke = 0) +
    scale_colour_manual(values = pal, labels = lab, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal() + labs(x = "UMAP_1", y = "UMAP_2", title = title)
}
ggsave(file.path(out_dir, "umap_pre_harmony_by_sample.png"),
       umap_sample("UMAP_pca_1", "UMAP_pca_2",
                   sprintf("Before integration: UMAP on PC 1-%d (%s nuclei)", n_pcs, n_lab)),
       width = 9, height = 7.5, dpi = 200)
ggsave(file.path(out_dir, "umap_harmony_by_sample.png"),
       umap_sample("UMAP_harmony_1", "UMAP_harmony_2",
                   sprintf("After Harmony (by sample, max_iter %d): UMAP on harmony 1-%d",
                           max_iter_h, n_pcs)),
       width = 9, height = 7.5, dpi = 200)

cent <- pc[, .(x = median(UMAP_harmony_1), y = median(UMAP_harmony_2)), by = cluster]
k <- length(lvl)
cl_pal <- setNames(grDevices::hcl(h = seq(15, 375, length.out = k + 1)[1:k], c = 70, l = 62), lvl)
p_cl <- ggplot(pd, aes(UMAP_harmony_1, UMAP_harmony_2, colour = cluster)) +
  geom_point(size = 0.05, alpha = 0.4, stroke = 0, show.legend = FALSE) +
  geom_label(data = cent, aes(x, y, label = cluster), inherit.aes = FALSE, size = 3,
             linewidth = 0, fill = alpha("white", 0.75), label.padding = unit(0.12, "lines")) +
  scale_colour_manual(values = cl_pal) + coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2",
       title = sprintf("After Harmony: Leiden clusters (resolution %.1f, %d clusters)", resolution, k),
       subtitle = sprintf("igraph Leiden (modularity, n_iterations = %d) on SNN graph, k = %d, harmony 1-%d",
                          leiden_iter, n_nb, n_pcs))
ggsave(file.path(out_dir, "umap_harmony_by_cluster.png"), p_cl, width = 9, height = 8, dpi = 200)

# per-sample panels, before (row 1) vs after (row 2)
long <- rbind(pc[, .(sample, x = UMAP_pca_1, y = UMAP_pca_2, stage = "Before (PCA)")],
              pc[, .(sample, x = UMAP_harmony_1, y = UMAP_harmony_2, stage = "After (Harmony)")])
long[, stage := factor(stage, levels = c("Before (PCA)", "After (Harmony)"))]
long[, panel := factor(lab[sample], levels = lab)]
bg <- long[, .(x, y, stage)]
p_split <- ggplot(long, aes(x, y)) +
  geom_point(data = bg, colour = "grey85", size = 0.02, stroke = 0) +
  geom_point(aes(colour = sample), size = 0.06, alpha = 0.5, stroke = 0, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  facet_grid(stage ~ panel, scales = "free") +
  labs(x = "UMAP_1", y = "UMAP_2",
       title = "Each sample (coloured) on all other nuclei (grey): before vs after Harmony")
ggsave(file.path(out_dir, "umap_pre_vs_post_by_sample_split.png"), p_split,
       width = 24, height = 7, dpi = 150)

sdev <- data.table(PC = seq_along(obj[["pca"]]@stdev), sd = obj[["pca"]]@stdev)
p_elb <- ggplot(sdev, aes(PC, sd)) + geom_point(size = 2, colour = "#257de0") +
  geom_vline(xintercept = n_pcs + 0.5, linetype = 2, colour = "grey40") +
  labs(y = "standard deviation", title = sprintf("PCA elbow (dashed line: PCs used = 1-%d)", n_pcs))
ggsave(file.path(out_dir, "pca_elbow.png"), p_elb, width = 6, height = 4, dpi = 150)

cat("\n===== clusters x samples =====\n"); print(ct[, c("cluster", samples, "total"), with = FALSE])
message("done")
