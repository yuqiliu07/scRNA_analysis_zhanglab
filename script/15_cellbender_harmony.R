#!/usr/bin/env Rscript
# Re-run the clustering pipeline (steps 07, 09, 11) on CellBender-corrected counts.
#
#  1. Counts = <CB_DIR>/<s>/<s>_cellbender.h5 (13_cellbender_<s>.sh; FPR 0.01), all 8 samples.
#     Nuclei = pass_qc_singlet of step 05 (same set as fsu_merged_mm10.rds)
#              AND called as a cell by CellBender (<s>_cellbender_cell_barcodes.csv)
#              AND >= MIN_GENES (200) genes detected after correction.
#     Gene names = make.unique(Symbol), same as 07_seurat_merge.R.
#  2. As step 07: gene filter (count > 0 in >= MIN_CELLS = 10 nuclei), LogNormalize (1e4) per sample
#     layer, FindVariableFeatures (vst, N_HVG = 3000) per sample -> union of the per-sample HVGs.
#  3. As step 09: ScaleData (HVG union) -> PCA (50, use 1..N_PCS = 30) -> UMAP before integration
#     -> Harmony (by sample, max_iter 20) -> UMAP + SNN (k = 20) -> igraph Leiden
#     (modularity, resolution 1.0, n_iterations 2), clusters relabelled 0..K-1 by size.
#  4. As step 11: clusters < MIN_CLUSTER_SIZE (10) nuclei dissolved into the majority cluster of their
#     k = 20 harmony neighbours (-> cb_cluster); 17 marker panels scored per cluster.
#  5. Comparison with the raw-count clustering (step 11, annot_mm10/cluster_clean_per_cell.csv.gz):
#     ARI, new cluster x old cluster table, UMAP coloured by the old cell type, and the
#     Cardiomyocyte panel score per old cell type, raw (step 11) vs CellBender.
#
# Outputs (<OUT_DIR>/)
#   data, stay on Hoffman: fsu_cellbender_harmony.rds, cb_per_cell.csv.gz
#   plots: umap_pre_harmony_by_sample.png, umap_harmony_by_sample.png, umap_harmony_by_cluster.png,
#          umap_old_celltype.png, pca_elbow.png, cluster_celltype_heatmap.png,
#          old_vs_new_cluster_heatmap.png, cm_panel_by_celltype.png
#   tables: cb_params.csv, cb_cells_per_sample.csv, cb_hvg_union.csv, cluster_sample_counts.csv,
#           cluster_summary.csv, cluster_panel_scores.csv, old_vs_new_crosstab.csv,
#           cm_panel_by_celltype.csv
#
# Usage: Rscript 15_cellbender_harmony.R
# Env (optional): TEST_MAX_NUCLEI, CB_DIR, QC_DIR, ANNOT_DIR, OUT_DIR, SAMPLE_SHEET, MIN_GENES (200), MIN_CELLS (10),
#                 N_HVG (3000), N_PCS_COMPUTE (50), N_PCS (30), MAX_ITER_HARMONY (20),
#                 N_NEIGHBORS (20), RESOLUTION (1.0), LEIDEN_ITER (2), MIN_CLUSTER_SIZE (10)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(igraph)
  library(Matrix)
  library(rhdf5)
  library(data.table)
  library(ggplot2)
})
options(Seurat.object.assay.version = "v5")

scratch      <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
wd           <- file.path(scratch, "FSU_seq/work_dir")
cb_dir       <- Sys.getenv("CB_DIR",    file.path(wd, "cellbender_mm10"))
qc_dir       <- Sys.getenv("QC_DIR",    file.path(wd, "qc_mm10"))
annot_dir    <- Sys.getenv("ANNOT_DIR", file.path(wd, "annot_mm10"))
out_dir      <- Sys.getenv("OUT_DIR",   file.path(wd, "harmony_cb_mm10"))
sample_sheet <- Sys.getenv("SAMPLE_SHEET", file.path(scratch, "FSU_seq/sample_sheet.csv"))
min_genes    <- as.integer(Sys.getenv("MIN_GENES", "200"))
min_cells    <- as.integer(Sys.getenv("MIN_CELLS", "10"))
n_hvg        <- as.integer(Sys.getenv("N_HVG", "3000"))
n_pcs_comp   <- as.integer(Sys.getenv("N_PCS_COMPUTE", "50"))
n_pcs        <- as.integer(Sys.getenv("N_PCS", "30"))
max_iter_h   <- as.integer(Sys.getenv("MAX_ITER_HARMONY", "20"))
n_nb         <- as.integer(Sys.getenv("N_NEIGHBORS", "20"))
resolution   <- as.numeric(Sys.getenv("RESOLUTION", "1.0"))
leiden_iter  <- as.integer(Sys.getenv("LEIDEN_ITER", "2"))
min_size     <- as.integer(Sys.getenv("MIN_CLUSTER_SIZE", "10"))
test_n       <- as.integer(Sys.getenv("TEST_MAX_NUCLEI", "0"))   # >0: first n nuclei per sample (login-node test)
SEED <- 42
samples <- LETTERS[1:8]
stopifnot(n_pcs <= n_pcs_comp)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PANELS <- list(   # sxj_code/04_merge_cluster.py, same as 11_cluster_clean_annot.R
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
  Lung_epi       = c("Sftpc", "Ager", "Scgb1a1"))
# step-11 cluster -> cell type (same map as 14_compare_ambient.R)
CELLTYPE <- c(
  setNames(rep("CM", 3), c(0, 4, 6)),
  setNames(rep("CM+X (doublet?)", 9), c(10, 11, 14, 18, 19, 30, 31, 33, 34)),
  setNames(rep("Fibroblast", 5), c(1, 5, 23, 13, 24)),
  setNames(rep("Endothelial", 4), c(3, 7, 20, 2)),
  setNames(rep("Endocardial", 2), c(12, 22)),
  "16" = "Lymphatic EC", setNames(rep("Pericyte", 2), c(8, 21)), "25" = "SMC",
  setNames(rep("Macrophage", 4), c(9, 15, 17, 27)),
  "28" = "T/NK + B", "29" = "Schwann", "35" = "Epicardial",
  "26" = "Other (26, 32)", "32" = "Other (26, 32)")
CT_ORDER <- c("CM", "CM+X (doublet?)", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic EC",
              "Pericyte", "SMC", "Macrophage", "T/NK + B", "Schwann", "Epicardial", "Other (26, 32)")

params <- data.table(counts = "CellBender remove-background 0.3.2 (FPR 0.01)",
                     nuclei = "pass_qc_singlet & CellBender cell & nGene >= MIN_GENES",
                     min_genes = min_genes, min_cells = min_cells, n_hvg_per_sample = n_hvg,
                     n_pcs_computed = n_pcs_comp, n_pcs_used = n_pcs, max_iter_harmony = max_iter_h,
                     n_neighbors = n_nb, resolution = resolution, leiden_n_iterations = leiden_iter,
                     min_cluster_size = min_size, seed = SEED,
                     seurat_version = as.character(packageVersion("Seurat")),
                     harmony_version = as.character(packageVersion("harmony")))
print(t(params))

# Stream a 10x-format h5 (CellBender output) and keep barcodes `bc` (same reader as 14_compare_ambient.R)
read_h5_cols <- function(path, bc, block = 4e6) {
  g <- function(x, ...) rhdf5::h5read(path, paste0("matrix/", x), bit64conversion = "double", ...)
  on.exit(rhdf5::h5closeAll())
  bcs <- as.character(g("barcodes"))
  pos <- match(bcs, bc)
  if (sum(!is.na(pos)) != length(bc)) stop(sprintf("%s: %d barcodes missing", path, length(setdiff(bc, bcs))))
  ids <- as.character(g("features/id")); sym <- as.character(g("features/name"))
  indptr <- as.numeric(g("indptr")); nnz <- indptr[length(indptr)]
  starts <- indptr[-length(indptr)]
  I <- list(); J <- list(); X <- list(); k <- 0
  for (s0 in seq(0, nnz - 1, by = block)) {
    n <- min(block, nnz - s0)
    jj <- pos[findInterval(s0 + seq_len(n) - 1, starts)]; keep <- !is.na(jj)
    if (!any(keep)) next
    k <- k + 1
    I[[k]] <- as.integer(g("indices", start = s0 + 1, count = n))[keep] + 1L
    X[[k]] <- as.numeric(g("data", start = s0 + 1, count = n))[keep]
    J[[k]] <- jj[keep]
  }
  m <- sparseMatrix(i = unlist(I), j = unlist(J), x = unlist(X),
                    dims = c(length(ids), length(bc)), dimnames = list(make.unique(sym), bc))
  list(m = m, gene_ids = data.table(gene = make.unique(sym), gene_id = ids))
}

# ---- 1. per-sample objects ----
ss <- fread(sample_sheet)
objs <- list(); cell_tab <- list(); gene_ids <- NULL
for (s in samples) {
  flags <- fread(file.path(qc_dir, paste0(s, "_qc_doublet_flags.csv")))
  qc_keep <- flags[toupper(as.character(pass_qc_singlet)) == "TRUE", barcode]
  cb_cells <- fread(file.path(cb_dir, s, paste0(s, "_cellbender_cell_barcodes.csv")), header = FALSE)$V1
  bc <- intersect(qc_keep, cb_cells)
  if (test_n > 0) bc <- head(bc, test_n)
  r <- read_h5_cols(file.path(cb_dir, s, paste0(s, "_cellbender.h5")), bc)
  m <- r$m
  if (is.null(gene_ids)) gene_ids <- r$gene_ids
  stopifnot(identical(gene_ids$gene, rownames(m)))
  ng <- diff(m@p)                                   # genes detected per nucleus after correction
  m <- m[, ng >= min_genes, drop = FALSE]
  colnames(m) <- paste0(s, "_", colnames(m))

  o <- CreateSeuratObject(counts = m, project = "FSU_seq", min.cells = 0, min.features = 0)
  md <- as.data.frame(flags[match(sub("^[A-H]_", "", colnames(o)), barcode)])
  rownames(md) <- colnames(o)
  md <- md[, setdiff(names(md), "barcode"), drop = FALSE]
  sm <- ss[sample == s]
  md$genotype <- sm$genotype; md$sex <- sm$sex; md$mouse_id <- sm$mouse_id; md$age <- sm$age
  o <- AddMetaData(o, md)
  o$orig.ident <- s
  o$cb_frac_removed <- 1 - o$nCount_RNA / o$qc_nUMI
  objs[[s]] <- o
  cell_tab[[s]] <- data.table(sample = s, n_pass_qc_singlet = length(qc_keep), n_cellbender_cells = length(cb_cells),
                              n_both = length(bc), n_qc_not_cb = length(setdiff(qc_keep, cb_cells)),
                              n_lt_min_genes = sum(ng < min_genes), n_kept = ncol(o),
                              median_umi_raw = median(o$qc_nUMI), median_umi_cb = median(o$nCount_RNA),
                              median_gene_cb = median(o$nFeature_RNA))
  message(sprintf("[%s] pass_qc_singlet %d | & CellBender cell %d | & nGene >= %d: %d",
                  s, length(qc_keep), length(bc), min_genes, ncol(o)))
  rm(m, o, r); invisible(gc())
}
cell_tab <- rbindlist(cell_tab)
fwrite(cell_tab, file.path(out_dir, "cb_cells_per_sample.csv"))
print(cell_tab)

obj <- merge(objs[[1]], y = objs[-1], project = "FSU_seq")
rm(objs); invisible(gc())
obj <- JoinLayers(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
stopifnot(setequal(Layers(obj, search = "counts"), paste0("counts.", samples)))

# ---- 2. gene filter, normalise, HVG union ----
count_layers <- Layers(obj, search = "counts")
n_expr <- Reduce(`+`, lapply(count_layers, function(l) {
  x <- drop0(as(LayerData(obj, layer = l), "CsparseMatrix"))
  tabulate(x@i + 1L, nbins = nrow(x))
}))
names(n_expr) <- rownames(obj)
n_genes_before <- nrow(obj)
obj <- subset(obj, features = names(n_expr)[n_expr >= min_cells])
message(sprintf("Gene filter (>= %d nuclei): %d -> %d genes", min_cells, n_genes_before, nrow(obj)))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
hvg_list <- list()
for (l in count_layers) {
  tmp <- CreateSeuratObject(counts = LayerData(obj, layer = l))
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  hvg_list[[sub("^counts\\.", "", l)]] <- VariableFeatures(tmp)
  rm(tmp)
}
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
hvg <- unique(unlist(hvg_list))
VariableFeatures(obj) <- hvg
n_sel <- table(unlist(hvg_list))
old_hvg_f <- file.path(wd, "seurat_mm10/seurat_hvg_union.csv")
old_hvg <- if (file.exists(old_hvg_f)) fread(old_hvg_f)$gene else character(0)
fwrite(data.table(gene = hvg, n_samples_selected = as.integer(n_sel[hvg]), in_raw_hvg_union = hvg %in% old_hvg)[
  order(-n_samples_selected)], file.path(out_dir, "cb_hvg_union.csv"))
message(sprintf("HVG union: %d genes (%d shared with the raw-count HVG union of %d)",
                length(hvg), sum(hvg %in% old_hvg), length(old_hvg)))
params[, `:=`(n_nuclei = ncol(obj), n_genes = nrow(obj), n_hvg_union = length(hvg),
              n_hvg_shared_with_raw = sum(hvg %in% old_hvg))]

# ---- 3. PCA, Harmony, UMAP, Leiden ----
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = n_pcs_comp, seed.use = SEED, verbose = FALSE)
obj[["RNA"]]$scale.data <- NULL
invisible(gc())
obj <- RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, n.neighbors = n_nb, seed.use = SEED,
               reduction.name = "umap_pca", reduction.key = "umapPCA_", verbose = FALSE)
set.seed(SEED)
h_args <- list(object = obj, group.by.vars = "orig.ident", reduction.use = "pca", dims.use = 1:n_pcs,
               reduction.save = "harmony", project.dim = FALSE, verbose = TRUE)
h_formals <- c(names(formals(getS3method("RunHarmony", "Seurat", envir = asNamespace("harmony")))),
               names(formals(getS3method("RunHarmony", "default", envir = asNamespace("harmony")))))
if ("max_iter" %in% h_formals) h_args$max_iter <- max_iter_h else h_args$max.iter.harmony <- max_iter_h
obj <- do.call(RunHarmony, h_args)
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
memb <- match(memb, order(-tabulate(memb))) - 1L
obj$cb_leiden <- factor(as.character(memb), levels = as.character(sort(unique(memb))))
rm(g); invisible(gc())
message(sprintf("Leiden: %d clusters, modularity = %.4f", nlevels(obj$cb_leiden), mod))

# ---- 4. cleanup of tiny clusters (as step 11) ----
clv <- as.character(obj$cb_leiden)
sizes <- table(clv)
small <- names(sizes)[sizes < min_size]
nn <- obj[["harmony_nn"]]
clean <- clv
for (i in which(clv %in% small)) {
  nb <- setdiff(which(nn[i, ] > 0), i)
  nb_cl <- clv[nb][!clv[nb] %in% small]
  if (!length(nb_cl)) stop("nucleus ", colnames(obj)[i], " has no neighbour in a kept cluster")
  tt <- table(nb_cl); best <- names(tt)[tt == max(tt)]
  clean[i] <- best[which.max(sizes[best])]
}
lv <- as.character(sort(as.integer(unique(clean))))
obj$cb_cluster <- factor(clean, levels = lv)
Idents(obj) <- "cb_cluster"
message(sprintf("cleanup: %d clusters < %d nuclei dissolved (%d nuclei); %d clusters",
                length(small), min_size, sum(clv %in% small), length(lv)))
params[, `:=`(n_clusters_leiden = nlevels(obj$cb_leiden), modularity = round(mod, 4), n_clusters_clean = length(lv))]

# panel scores (mean log-normalised expression of the panel genes)
genes_used <- lapply(PANELS, function(x) intersect(x, rownames(obj)))
genes_used <- genes_used[lengths(genes_used) > 0]
all_g <- unique(unlist(genes_used))
score <- matrix(NA_real_, ncol(obj), length(genes_used), dimnames = list(colnames(obj), names(genes_used)))
for (l in Layers(obj, search = "data")) {
  x <- as(LayerData(obj, layer = l, features = all_g), "CsparseMatrix")
  for (p in names(genes_used)) score[colnames(x), p] <- Matrix::colMeans(x[genes_used[[p]], , drop = FALSE])
}
stopifnot(!anyNA(score))
for (p in colnames(score)) obj[[paste0("panel_", p)]] <- score[, p]
ps <- t(sapply(lv, function(k) colMeans(score[obj$cb_cluster == k, , drop = FALSE])))
dimnames(ps) <- list(cluster = lv, panel = colnames(score))
ps_dt <- data.table(cluster = lv, round(ps, 4))
fwrite(ps_dt, file.path(out_dir, "cluster_panel_scores.csv"))

# ---- 5. comparison with the raw-count clustering ----
old <- fread(file.path(annot_dir, "cluster_clean_per_cell.csv.gz"))
old_cl <- setNames(as.character(old$cluster_clean), old$cell)
pc <- data.table(cell = colnames(obj), sample = as.character(obj$orig.ident),
                 cluster = as.character(obj$cb_cluster), old_cluster = old_cl[colnames(obj)],
                 nUMI_raw = obj$qc_nUMI, nUMI_cb = obj$nCount_RNA, nGene_cb = obj$nFeature_RNA,
                 frac_removed = round(obj$cb_frac_removed, 4), pct_mt = obj$qc_pct_mt,
                 dbl = obj$scDblFinder_score)
pc[, old_celltype := ifelse(is.na(old_cluster), "not in step 11", CELLTYPE[old_cluster])]
pc[is.na(old_celltype), old_celltype := "Other (26, 32)"]
u0 <- Embeddings(obj, "umap_pca"); u1 <- Embeddings(obj, "umap")
pc[, `:=`(UMAP_pca_1 = round(u0[, 1], 4), UMAP_pca_2 = round(u0[, 2], 4),
          UMAP_1 = round(u1[, 1], 4), UMAP_2 = round(u1[, 2], 4))]
pc[, CM_panel := round(score[, "Cardiomyocyte"], 4)]
fwrite(pc, file.path(out_dir, "cb_per_cell.csv.gz"))

ari <- function(a, b) {
  t <- table(a, b); n <- sum(t)
  s <- sum(choose(t, 2)); sa <- sum(choose(rowSums(t), 2)); sb <- sum(choose(colSums(t), 2))
  e <- sa * sb / choose(n, 2); (s - e) / ((sa + sb) / 2 - e)
}
both <- pc[!is.na(old_cluster)]
params[, `:=`(n_nuclei_in_step11 = nrow(both), ari_vs_step11_clusters = round(ari(both$cluster, both$old_cluster), 4),
              ari_vs_step11_celltypes = round(ari(both$cluster, both$old_celltype), 4))]
fwrite(params, file.path(out_dir, "cb_params.csv"))
xt <- dcast(both[, .N, by = .(cluster, old_cluster)], cluster ~ old_cluster, value.var = "N", fill = 0)
xt <- xt[order(as.integer(cluster))]
fwrite(xt, file.path(out_dir, "old_vs_new_crosstab.csv"))

# cluster summary: size, QC, sample spread, main old cluster / cell type, top panels
top2 <- t(apply(ps, 1, function(r) { o <- order(-r); c(names(r)[o[1]], r[o[1]], names(r)[o[2]], r[o[2]]) }))
summ <- pc[, .(n_nuclei = .N, median_nGene_cb = as.numeric(median(nGene_cb)),
               median_nUMI_cb = as.numeric(median(nUMI_cb)), median_frac_removed = round(median(frac_removed), 3),
               median_pct_mt = round(as.numeric(median(pct_mt)), 3), mean_doublet_score = round(mean(dbl), 4),
               n_samples_ge1pct = sum(table(factor(sample, levels = samples)) / .N >= 0.01),
               top_sample = names(which.max(table(sample))), top_sample_frac = round(max(table(sample)) / .N, 3),
               main_old_cluster = names(which.max(table(old_cluster))),
               main_old_cluster_frac = round(max(table(old_cluster)) / .N, 3),
               main_old_celltype = names(which.max(table(old_celltype))),
               main_old_celltype_frac = round(max(table(old_celltype)) / .N, 3)), by = cluster]
summ <- merge(summ, data.table(cluster = rownames(top2), panel_1 = top2[, 1],
                               panel_1_score = round(as.numeric(top2[, 2]), 3), panel_2 = top2[, 3],
                               panel_2_score = round(as.numeric(top2[, 4]), 3)), by = "cluster")
summ[, panel_ratio := round(panel_1_score / pmax(panel_2_score, 1e-6), 2)]
summ <- summ[order(as.integer(cluster))]
fwrite(summ, file.path(out_dir, "cluster_summary.csv"))

ct <- dcast(pc[, .N, by = .(cluster, sample)], cluster ~ sample, value.var = "N", fill = 0)
ct[, total := rowSums(.SD), .SDcols = samples]
for (s in samples) ct[[paste0("pct_", s)]] <- round(100 * ct[[s]] / sum(ct[[s]]), 3)
fwrite(ct[order(as.integer(cluster))], file.path(out_dir, "cluster_sample_counts.csv"))

# CM panel score per old cell type: raw (step 11, cluster means weighted by size) vs CellBender
old_ps <- fread(file.path(annot_dir, "cluster_panel_scores.csv"))
old_sz <- fread(file.path(annot_dir, "cluster_summary.csv"))[, .(cluster, n_nuclei)]
old_ps <- merge(old_ps[, .(cluster, Cardiomyocyte)], old_sz, by = "cluster")
old_ps[, celltype := CELLTYPE[as.character(cluster)]]
old_ps[is.na(celltype), celltype := "Other (26, 32)"]
cmp <- merge(old_ps[, .(raw = sum(Cardiomyocyte * n_nuclei) / sum(n_nuclei)), by = celltype],
             both[, .(cellbender = mean(CM_panel), n_nuclei = .N), by = .(celltype = old_celltype)], by = "celltype")
cmp[, celltype := factor(celltype, levels = CT_ORDER)]
cmp <- cmp[order(celltype)]
fwrite(cmp, file.path(out_dir, "cm_panel_by_celltype.csv"))

saveRDS(obj, file.path(out_dir, "fsu_cellbender_harmony.rds"))

# ---- plots ----
ss_lab <- unique(data.table(sample = as.character(obj$orig.ident), genotype = obj$genotype, sex = obj$sex))[order(sample)]
lab <- setNames(sprintf("%s (%s, %s)", ss_lab$sample, ss_lab$genotype, ss_lab$sex), ss_lab$sample)
pal <- setNames(c("#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"), samples)
theme_set(theme_classic(base_size = 11) +
          theme(axis.line = element_line(colour = "grey60", linewidth = 0.3),
                axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
                axis.text = element_text(colour = "grey30"),
                plot.background = element_rect(fill = "white", colour = NA)))
n_lab <- format(nrow(pc), big.mark = ",")
set.seed(SEED)
pd <- pc[sample(.N)]
umap_sample <- function(x, y, title) {
  ggplot(pd, aes(.data[[x]], .data[[y]], colour = sample)) +
    geom_point(size = 0.05, alpha = 0.4, stroke = 0) +
    scale_colour_manual(values = pal, labels = lab, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal() + labs(x = "UMAP_1", y = "UMAP_2", title = title)
}
ggsave(file.path(out_dir, "umap_pre_harmony_by_sample.png"),
       umap_sample("UMAP_pca_1", "UMAP_pca_2", sprintf("CellBender counts, before integration (%s nuclei)", n_lab)),
       width = 9, height = 7.5, dpi = 200)
ggsave(file.path(out_dir, "umap_harmony_by_sample.png"),
       umap_sample("UMAP_1", "UMAP_2", "CellBender counts, after Harmony (by sample)"),
       width = 9, height = 7.5, dpi = 200)

k <- length(lv)
cl_pal <- setNames(grDevices::hcl(h = seq(15, 375, length.out = k + 1)[1:k], c = 70, l = 62), lv)
cent <- pc[, .(x = median(UMAP_1), y = median(UMAP_2)), by = cluster]
p <- ggplot(pd, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = 0.05, alpha = 0.4, stroke = 0, show.legend = FALSE) +
  geom_label(data = cent, aes(x, y, label = cluster), inherit.aes = FALSE, size = 3,
             linewidth = 0, fill = alpha("white", 0.75), label.padding = unit(0.12, "lines")) +
  scale_colour_manual(values = cl_pal) + coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2",
       title = sprintf("CellBender counts: Leiden clusters (resolution %.1f, %d clusters)", resolution, k),
       subtitle = sprintf("ARI vs raw-count clusters (step 11) = %.3f; vs step-11 cell types = %.3f",
                          params$ari_vs_step11_clusters, params$ari_vs_step11_celltypes))
ggsave(file.path(out_dir, "umap_harmony_by_cluster.png"), p, width = 9, height = 8, dpi = 200)

ct_lv <- c(CT_ORDER, "not in step 11")
ct_pal <- setNames(c("#257de0", "#9cc3ef", "#f47230", "#00b894", "#7fdcc9", "#008a17", "#f2a600",
                     "#a87400", "#ee7fac", "#ec4647", "#6a4ce2", "#8b5a2b", "#666666", "#c8c8c8"), ct_lv)
pd[, old_celltype := factor(old_celltype, levels = ct_lv)]
p <- ggplot(pd, aes(UMAP_1, UMAP_2, colour = old_celltype)) +
  geom_point(size = 0.05, alpha = 0.4, stroke = 0) +
  geom_label(data = cent, aes(x, y, label = cluster), inherit.aes = FALSE, size = 2.6,
             linewidth = 0, fill = alpha("white", 0.6), label.padding = unit(0.1, "lines")) +
  scale_colour_manual(values = ct_pal, name = "step-11 cell type", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2", title = "CellBender UMAP coloured by the raw-count (step 11) cell type")
ggsave(file.path(out_dir, "umap_old_celltype.png"), p, width = 10, height = 7.5, dpi = 200)

sdev <- data.table(PC = seq_along(obj[["pca"]]@stdev), sd = obj[["pca"]]@stdev)
p <- ggplot(sdev, aes(PC, sd)) + geom_point(size = 2, colour = "#257de0") +
  geom_vline(xintercept = n_pcs + 0.5, linetype = 2, colour = "grey40") +
  labs(y = "standard deviation", title = sprintf("PCA elbow, CellBender counts (PCs used = 1-%d)", n_pcs))
ggsave(file.path(out_dir, "pca_elbow.png"), p, width = 6, height = 4, dpi = 150)

pan_lv <- colnames(ps)
ord <- summ[, .(cluster, p1 = factor(panel_1, levels = pan_lv), panel_1_score)][order(p1, -panel_1_score), cluster]
hm <- melt(ps_dt, id.vars = "cluster", variable.name = "panel", value.name = "score")
hm[, scaled := (score - min(score)) / max(max(score) - min(score), 1e-9), by = panel]
hm[, `:=`(cluster = factor(cluster, levels = rev(ord)), panel = factor(panel, levels = pan_lv))]
top <- merge(hm, summ[, .(cluster, panel = panel_1)], by = c("cluster", "panel"))
p <- ggplot(hm, aes(panel, cluster, fill = scaled)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.1, colour = ifelse(hm$scaled > 0.6, "white", "grey20")) +
  geom_tile(data = top, fill = NA, colour = "black", linewidth = 0.6) +
  scale_fill_gradient(low = "#f4f8fd", high = "#0b3d78", name = "scaled\n(per panel)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "cluster (CellBender)",
       title = "Marker panel score per cluster, CellBender counts (numbers = mean log-expr)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.line = element_blank())
ggsave(file.path(out_dir, "cluster_celltype_heatmap.png"), p, width = 11, height = 0.25 * k + 2.5, dpi = 170)

xl <- melt(xt, id.vars = "cluster", variable.name = "old_cluster", value.name = "n")
xl[, frac := n / sum(n), by = cluster]
old_lv <- as.character(sort(as.integer(unique(as.character(xl$old_cluster)))))
xl[, `:=`(cluster = factor(cluster, levels = rev(ord)), old_cluster = factor(old_cluster, levels = old_lv))]
p <- ggplot(xl, aes(old_cluster, cluster, fill = frac)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(data = xl[frac >= 0.1], aes(label = round(100 * frac)), size = 2, colour = "white") +
  scale_fill_gradient(low = "#f4f8fd", high = "#0b3d78", labels = scales::percent, name = "% of new\ncluster") +
  labs(x = "raw-count cluster (step 11)", y = "CellBender cluster",
       title = "Where each CellBender cluster comes from (row = 100%)",
       subtitle = sprintf("ARI = %.3f", params$ari_vs_step11_clusters)) +
  theme(axis.line = element_blank())
ggsave(file.path(out_dir, "old_vs_new_cluster_heatmap.png"), p, width = 12, height = 0.25 * k + 2.5, dpi = 170)

cl_long <- melt(cmp[, .(celltype, raw, cellbender)], id.vars = "celltype", variable.name = "counts", value.name = "score")
cl_long[, counts := factor(ifelse(counts == "raw", "raw (step 11)", "CellBender"), levels = c("raw (step 11)", "CellBender"))]
cl_long[, celltype := factor(celltype, levels = rev(CT_ORDER))]
p <- ggplot(cl_long, aes(score, celltype, colour = counts)) +
  geom_line(aes(group = celltype), colour = "grey75") + geom_point(size = 3) +
  scale_colour_manual(values = c(`raw (step 11)` = "#257de0", CellBender = "#f47230"), name = NULL) +
  labs(x = "Cardiomyocyte panel score (mean log-expr of Myh6, Tnnt2, Ryr2, Ttn ...)", y = NULL,
       title = "CM marker signal per cell type: raw vs CellBender counts") +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "cm_panel_by_celltype.png"), p, width = 7.5, height = 5, dpi = 150)

cat("\n===== params =====\n"); print(t(params))
cat("\n===== CM panel score by step-11 cell type =====\n"); print(cmp)
cat("\n===== clusters =====\n")
print(summ[, .(cluster, n_nuclei, median_nGene_cb, median_frac_removed, top_sample, top_sample_frac,
               main_old_cluster, main_old_cluster_frac, main_old_celltype, panel_1, panel_1_score,
               panel_2, panel_2_score)])
message("done")
