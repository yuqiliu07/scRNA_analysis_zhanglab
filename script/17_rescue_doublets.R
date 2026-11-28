#!/usr/bin/env Rscript
# Step 17 (= "step 5C"): rescue single-lineage nuclei from the doublet clusters of the CellBender
# clustering, add lymphocytes, re-embed, re-check, and write the final annotated object.
#
# Input: step-15 Seurat object (<CB_DIR>/fsu_cellbender_harmony.rds; counts = CellBender, RNA assay
# split by sample, metadata incl. scDblFinder_class and cb_cluster) and the step-16 panel table
# (<CB_DIR>/clean/cluster_panels_clean.csv; panel_1 per cluster -> cell type of the kept clusters).
#
# "5B" (cluster triage of step 15/16, defaults below, override with env):
#   LOWQ_CLUSTERS    14,26,30   low quality, dropped, no rescue (14: nGene 193; 26/30: pct_mt 1.5-1.9%)
#   DOUBLET_CLUSTERS 6,7,13,19,23,31,32,33 (CM+X)  12,15,18,20,21,22 (Fib+EC, Fib+Mac, EC+Mac,
#                    Endocardial+Fib, Peri+EC, Fib+Peri): 45-83% of their nuclei carry >= 2 lineages
#   all other clusters are kept ("confirmed"); cell type = panel_1 of step 16 with
#   B / T_NK -> Lymphocyte, Neural_Schwann -> Schwann, and CELLTYPE_OVERRIDE (default 34=Epicardial:
#   panel_1 Endocardial 0.79 vs Epicardial 0.70, step-11 label Epicardial, WT-only like step 11).
#
# Detection = CellBender count > 0. A nucleus "carries" a lineage if >= 2 of its LINEAGE genes are detected.
# STEP 1  lymphocytes, all nuclei of the object: >= 2 lymphocyte-specific genes (candidate list LYM,
#         kept if detected in < 0.5% of all nuclei), Ptprc detected, no non-CM lineage carried,
#         scDblFinder_class != doublet  -> cell type "Lymphocyte" (added back wherever they are).
# STEP 2  rescue candidates, nuclei of DOUBLET_CLUSTERS: exactly one lineage carried (others < 2
#         genes), not lymphocyte, not scDblFinder doublet.
# STEP 3  confirmed + candidates re-embedded (per-sample vst HVG union, PCA 50, Harmony by sample on
#         PCs 1-30); 15 nearest confirmed neighbours of every candidate (harmony 1-30). Rescued if
#         >= 70% of them share one cell type and that type belongs to the candidate's lineage
#         (CM -> Cardiomyocyte; Fib -> Fibroblast; EC -> Endothelial/Endocardial/Lymphatic_EC;
#          Peri -> Pericyte/SMC; SMC -> SMC/Pericyte; Mac -> Macrophage); otherwise discarded.
# STEP 4  kept = confirmed + rescued + lymphocytes; HVG/PCA/Harmony/UMAP/Leiden again (same
#         parameters as step 15). A new cluster is removed (no re-embedding) if > 50% of its nuclei
#         are doublets (carry >= 2 lineages, or scDblFinder doublet) AND > 50% carry its second
#         lineage (most frequent carried lineage after the cluster's main one).
# STEP 5  outputs:
#   <FINAL_RDS>   ($SCRATCH/FSU_seq/work_dir/merged/final_annotated.rds) Seurat object with
#                 cell_type_final, rescue_status (kept/rescued/lymphocyte), cluster_5c, umap/harmony
#   <TABLE_DIR>   ($SCRATCH/FSU_seq/results/tables) rescue_summary.csv, new_cluster_check.csv,
#                 final_celltype_counts.csv
#   <FIG_DIR>     ($SCRATCH/FSU_seq/results/figures) rescue_umap.pdf (+ PNGs for Drive)
#   <STATE_JSON>  ($SCRATCH/FSU_seq/pipeline_state.json) "step_5C" added to "completed"
#   rescue_per_cell.csv.gz next to FINAL_RDS (stays on Hoffman)
#
# Usage: Rscript 17_rescue_doublets.R
# Env: CB_DIR, FINAL_RDS, TABLE_DIR, FIG_DIR, STATE_JSON, LOWQ_CLUSTERS, DOUBLET_CLUSTERS,
#      CELLTYPE_OVERRIDE, N_HVG (3000), N_PCS_COMPUTE (50), N_PCS (30), MAX_ITER_HARMONY (20),
#      N_NEIGHBORS (20), RESOLUTION (1.0), LEIDEN_ITER (2), MIN_CLUSTER_SIZE (10), K_RESCUE (15),
#      RESCUE_FRAC (0.7), LYM_MAX_DETECT (0.005)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(harmony); library(igraph)
  library(Matrix); library(data.table); library(ggplot2); library(RANN); library(jsonlite)
})
options(Seurat.object.assay.version = "v5")

scratch    <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
cb_dir     <- Sys.getenv("CB_DIR", file.path(scratch, "FSU_seq/work_dir/harmony_cb_mm10"))
final_rds  <- Sys.getenv("FINAL_RDS", file.path(scratch, "FSU_seq/work_dir/merged/final_annotated.rds"))
table_dir  <- Sys.getenv("TABLE_DIR", file.path(scratch, "FSU_seq/results/tables"))
fig_dir    <- Sys.getenv("FIG_DIR", file.path(scratch, "FSU_seq/results/figures"))
state_json <- Sys.getenv("STATE_JSON", file.path(scratch, "FSU_seq/pipeline_state.json"))
csv_int    <- function(x) as.character(as.integer(strsplit(x, ",")[[1]]))
lowq       <- csv_int(Sys.getenv("LOWQ_CLUSTERS", "14,26,30"))
dbl_cl     <- csv_int(Sys.getenv("DOUBLET_CLUSTERS", "6,7,13,19,23,31,32,33,12,15,18,20,21,22"))
ct_over    <- strsplit(strsplit(Sys.getenv("CELLTYPE_OVERRIDE", "34=Epicardial"), ",")[[1]], "=")
n_hvg      <- as.integer(Sys.getenv("N_HVG", "3000"))
n_pcs_comp <- as.integer(Sys.getenv("N_PCS_COMPUTE", "50"))
n_pcs      <- as.integer(Sys.getenv("N_PCS", "30"))
max_iter_h <- as.integer(Sys.getenv("MAX_ITER_HARMONY", "20"))
n_nb       <- as.integer(Sys.getenv("N_NEIGHBORS", "20"))
resolution <- as.numeric(Sys.getenv("RESOLUTION", "1.0"))
leiden_it  <- as.integer(Sys.getenv("LEIDEN_ITER", "2"))
min_size   <- as.integer(Sys.getenv("MIN_CLUSTER_SIZE", "10"))
k_resc     <- as.integer(Sys.getenv("K_RESCUE", "15"))
resc_frac  <- as.numeric(Sys.getenv("RESCUE_FRAC", "0.7"))
lym_max    <- as.numeric(Sys.getenv("LYM_MAX_DETECT", "0.005"))
SEED <- 42
for (d in c(dirname(final_rds), table_dir, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

LINEAGE <- list(CM = c("Myh6", "Tnnt2", "Ryr2", "Ttn", "Rbm20"),
                Fib = c("Dcn", "Col1a1", "Pdgfra", "Abca8a", "Gsn"),
                EC = c("Pecam1", "Cdh5", "Flt1", "Kdr", "Ptprb"),
                Peri = c("Rgs5", "Kcnj8", "Abcc9", "Pdgfrb"),
                Mac = c("C1qa", "F13a1", "Mrc1", "Cd163", "Adgre1", "Csf1r"),
                SMC = c("Myh11", "Acta2", "Cnn1"))
LIN2CT <- list(CM = "Cardiomyocyte", Fib = "Fibroblast", EC = c("Endothelial", "Endocardial", "Lymphatic_EC"),
               Peri = c("Pericyte", "SMC"), SMC = c("SMC", "Pericyte"), Mac = "Macrophage")
LYM <- c("Cd3e", "Cd3g", "Cd3d", "Cd2", "Cd247", "Lck", "Itk", "Skap1", "Il7r", "Cd8a", "Cd8b1", "Cd4",
         "Trac", "Trbc2", "Nkg7", "Klrd1", "Klrk1", "Gzma", "Ccl5", "Cd79a", "Cd79b", "Ms4a1", "Pax5",
         "Bank1", "Cd19", "Ighm", "Cd22", "Fcmr")
HEME <- "Ptprc"
CT_LEVELS <- c("Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte",
               "SMC", "Macrophage", "Lymphocyte", "Schwann", "Epicardial")

# ---- load ----
obj <- readRDS(file.path(cb_dir, "fsu_cellbender_harmony.rds"))
obj$cb_cluster <- as.character(obj$cb_cluster)
message(sprintf("step-15 object: %d nuclei, %d genes, %d clusters", ncol(obj), nrow(obj), length(unique(obj$cb_cluster))))
pan <- fread(file.path(cb_dir, "clean/cluster_panels_clean.csv"), colClasses = list(character = "cluster"))
stopifnot(all(dbl_cl %in% obj$cb_cluster), all(lowq %in% obj$cb_cluster), !any(dbl_cl %in% lowq))
rename <- c(B = "Lymphocyte", T_NK = "Lymphocyte", Neural_Schwann = "Schwann")
pan[, cell_type := ifelse(panel_1 %in% names(rename), rename[panel_1], panel_1)]
for (o in ct_over) pan[cluster == o[1], cell_type := o[2]]
kept_cl <- setdiff(unique(obj$cb_cluster), c(lowq, dbl_cl))
stopifnot(all(kept_cl %in% pan$cluster))
ct_map <- setNames(pan$cell_type, pan$cluster)[kept_cl]
stopifnot(all(ct_map %in% CT_LEVELS))
triage <- data.table(cluster = sort(as.integer(unique(obj$cb_cluster))))[, cluster := as.character(cluster)]
triage[, role := fifelse(cluster %in% lowq, "low_quality", fifelse(cluster %in% dbl_cl, "doublet", "kept"))]
triage[, cell_type := ct_map[cluster]]
triage[, n_nuclei := as.integer(table(obj$cb_cluster)[cluster])]
print(triage, nrows = 100)

# ---- detection of lineage / lymphocyte genes (CellBender counts > 0) ----
det_matrix <- function(o, genes) {
  g <- intersect(genes, rownames(o))
  m <- do.call(cbind, lapply(Layers(o, search = "counts"), function(l) {
    x <- as(LayerData(o, layer = l, features = g), "CsparseMatrix")[g, , drop = FALSE]
    x@x <- as.numeric(x@x > 0); drop0(x)
  }))
  m[, colnames(o), drop = FALSE]
}
all_genes <- unique(c(unlist(LINEAGE), LYM, HEME))
miss <- setdiff(all_genes, rownames(obj))
if (length(miss)) message("genes not in the object (detected in < 10 nuclei): ", paste(miss, collapse = ", "))
D <- det_matrix(obj, all_genes)
lin_n <- sapply(LINEAGE, function(g) Matrix::colSums(D[intersect(g, rownames(D)), , drop = FALSE]))
carry <- lin_n >= 2
n_lin <- rowSums(carry)
md <- data.table(cell = colnames(obj), sample = obj$orig.ident, cb_cluster = obj$cb_cluster,
                 scdbl = tolower(as.character(obj$scDblFinder_class)) == "doublet", n_lin = n_lin,
                 lineage = ifelse(n_lin == 1, colnames(carry)[max.col(carry, "first")], ifelse(n_lin == 0, "none", "multi")),
                 combo = apply(carry, 1, function(r) paste(colnames(carry)[r], collapse = "+")))
md[, role := triage$role[match(cb_cluster, triage$cluster)]]
message(sprintf("scDblFinder doublets in the object: %d", sum(md$scdbl)))

# ---- STEP 1: lymphocytes ----
lym_det <- Matrix::rowMeans(D[intersect(LYM, rownames(D)), , drop = FALSE])
lym_genes <- names(lym_det)[lym_det < lym_max]
message(sprintf("lymphocyte-specific genes (detected in < %.1f%% of %d nuclei): %s", 100 * lym_max, ncol(obj),
                paste(lym_genes, collapse = ", ")))
message("  excluded (too widely detected): ", paste(setdiff(names(lym_det), lym_genes), collapse = ", "))
md[, n_lym := Matrix::colSums(D[lym_genes, , drop = FALSE])]
md[, heme := if (HEME %in% rownames(D)) D[HEME, ] > 0 else FALSE]
md[, nonCM_lin := rowSums(carry[, setdiff(colnames(carry), "CM"), drop = FALSE])]
md[, lymphocyte := n_lym >= 2 & heme & nonCM_lin == 0 & !scdbl]
message(sprintf("STEP 1: %d lymphocytes", sum(md$lymphocyte)))
print(md[lymphocyte == TRUE, .N, by = .(cb_cluster, role)][order(-N)])

# ---- STEP 2: candidates ----
md[, candidate := role == "doublet" & !lymphocyte & n_lin == 1 & !scdbl]
cand_tab <- md[role == "doublet", .(n_nuclei = .N, n_lymphocyte = sum(lymphocyte), n_candidate = sum(candidate),
                                    candidate_lineages = paste(names(sort(-table(lineage[candidate]))),
                                                               sort(table(lineage[candidate]), decreasing = TRUE),
                                                               sep = ":", collapse = "; ")), by = cb_cluster]
cand_tab <- cand_tab[order(as.integer(cb_cluster))]
message("STEP 2: rescue candidates per doublet cluster")
print(cand_tab)

# confirmed nuclei and their cell types
md[, cell_type := ifelse(role == "kept", ct_map[cb_cluster], NA_character_)]
md[lymphocyte == TRUE, cell_type := "Lymphocyte"]
md[, confirmed := (role == "kept") | lymphocyte]

# ---- embedding helpers (same parameters as step 15) ----
embed <- function(o) {
  hv <- list()
  for (l in Layers(o, search = "counts")) {
    tmp <- CreateSeuratObject(counts = LayerData(o, layer = l))
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
    hv[[l]] <- VariableFeatures(tmp); rm(tmp)
  }
  hvg <- unique(unlist(hv)); VariableFeatures(o) <- hvg
  message(sprintf("  HVG union %d genes, %d nuclei", length(hvg), ncol(o)))
  o <- ScaleData(o, features = hvg, verbose = FALSE)
  o <- RunPCA(o, features = hvg, npcs = n_pcs_comp, seed.use = SEED, verbose = FALSE)
  o[["RNA"]]$scale.data <- NULL; invisible(gc())
  set.seed(SEED)
  a <- list(object = o, group.by.vars = "orig.ident", reduction.use = "pca", dims.use = 1:n_pcs,
            reduction.save = "harmony", project.dim = FALSE, verbose = FALSE)
  f <- c(names(formals(getS3method("RunHarmony", "Seurat", envir = asNamespace("harmony")))),
         names(formals(getS3method("RunHarmony", "default", envir = asNamespace("harmony")))))
  if ("max_iter" %in% f) a$max_iter <- max_iter_h else a$max.iter.harmony <- max_iter_h
  do.call(RunHarmony, a)
}
leiden <- function(o) {
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:n_pcs, k.param = n_nb,
                     graph.name = c("harmony_nn", "harmony_snn"), verbose = FALSE)
  g <- graph_from_adjacency_matrix(as(o[["harmony_snn"]], "dgCMatrix"), mode = "undirected", weighted = TRUE, diag = FALSE)
  set.seed(SEED)
  cl <- cluster_leiden(g, objective_function = "modularity", weights = E(g)$weight,
                       resolution = resolution, n_iterations = leiden_it)
  memb <- membership(cl); memb <- match(memb, order(-tabulate(memb))) - 1L
  clv <- as.character(memb); sizes <- table(clv); small <- names(sizes)[sizes < min_size]
  nn <- o[["harmony_nn"]]; clean <- clv
  for (i in which(clv %in% small)) {
    nb <- setdiff(which(nn[i, ] > 0), i); nb_cl <- clv[nb][!clv[nb] %in% small]
    if (!length(nb_cl)) stop("nucleus ", colnames(o)[i], " has no neighbour in a kept cluster")
    tt <- table(nb_cl); best <- names(tt)[tt == max(tt)]; clean[i] <- best[which.max(sizes[best])]
  }
  message(sprintf("  Leiden: %d clusters (%d < %d nuclei dissolved)", length(sizes), length(small), min_size))
  o$cluster_5c <- factor(clean, levels = as.character(sort(as.integer(unique(clean)))))
  o
}

# ---- STEP 3: KNN rescue ----
message("STEP 3: re-embedding confirmed + candidate nuclei")
o3 <- subset(obj, cells = md[confirmed | candidate, cell])
rm(obj); invisible(gc())
for (r in setdiff(Reductions(o3), character(0))) o3[[r]] <- NULL   # step-15 embeddings are recomputed
o3@graphs <- list()
o3 <- embed(o3)
H <- Embeddings(o3, "harmony")[, 1:n_pcs]
md[, `:=`(knn_top = NA_character_, knn_frac = NA_real_, rescued = FALSE)]
ref <- md[confirmed == TRUE & cell %in% rownames(H)]
qry <- md[candidate == TRUE & cell %in% rownames(H)]
if (nrow(qry)) {
  nn <- nn2(H[ref$cell, , drop = FALSE], H[qry$cell, , drop = FALSE], k = k_resc)$nn.idx
  vt <- matrix(ref$cell_type[nn], nrow = nrow(nn))
  top <- apply(vt, 1, function(r) { t <- table(r); names(t)[which.max(t)] })
  frac <- rowMeans(vt == top)
  qry[, `:=`(knn_top = top, knn_frac = round(frac, 3))]
  qry[, rescued := knn_frac >= resc_frac & mapply(function(t, l) t %in% LIN2CT[[l]], knn_top, lineage)]
  md[qry, on = "cell", `:=`(knn_top = i.knn_top, knn_frac = i.knn_frac, rescued = i.rescued)]
}
md[rescued == TRUE, cell_type := knn_top]
message(sprintf("STEP 3: %d of %d candidates rescued", sum(md$rescued), sum(md$candidate)))
fail <- md[candidate == TRUE & !rescued, .(n = .N), by = .(reason = fifelse(knn_frac < resc_frac, "neighbours < 70% agree", "type not in own lineage"))]
print(fail)

rs <- md[role == "doublet", .(n_nuclei = .N, n_lymphocyte = sum(lymphocyte), n_candidate = sum(candidate),
                              n_rescued = sum(rescued), n_discarded = sum(candidate & !rescued),
                              rescued_cell_type = {t <- sort(table(cell_type[rescued]), decreasing = TRUE)
                                                   if (length(t)) paste(names(t), t, sep = ":", collapse = "; ") else ""}),
         by = .(cluster = cb_cluster)][order(as.integer(cluster))]
rs[, n_not_candidate := n_nuclei - n_lymphocyte - n_candidate]
rs <- rbind(rs, rs[, .(cluster = "total", n_nuclei = sum(n_nuclei), n_lymphocyte = sum(n_lymphocyte),
                       n_candidate = sum(n_candidate), n_rescued = sum(n_rescued), n_discarded = sum(n_discarded),
                       rescued_cell_type = {t <- sort(table(md[rescued == TRUE, cell_type]), decreasing = TRUE)
                                            paste(names(t), t, sep = ":", collapse = "; ")},
                       n_not_candidate = sum(n_not_candidate))])
setcolorder(rs, c("cluster", "n_candidate", "n_rescued", "n_discarded", "rescued_cell_type"))
fwrite(rs, file.path(table_dir, "rescue_summary.csv"))
print(rs)

# ---- STEP 4: re-embed kept nuclei, re-cluster, check for mixed clusters ----
md[, keep4 := confirmed | rescued]
md[, rescue_status := fifelse(lymphocyte, "lymphocyte", fifelse(rescued, "rescued", "kept"))]
message("STEP 4: re-embedding kept + rescued + lymphocytes")
o4 <- subset(o3, cells = md[keep4 == TRUE, cell])
rm(o3); invisible(gc())
o4 <- embed(o4)
o4 <- RunUMAP(o4, reduction = "harmony", dims = 1:n_pcs, n.neighbors = n_nb, seed.use = SEED,
              reduction.name = "umap", reduction.key = "umap_", verbose = FALSE)
o4 <- leiden(o4)
m4 <- md[match(colnames(o4), cell)]
m4[, cluster_5c := as.character(o4$cluster_5c)]
cw <- carry[m4$cell, , drop = FALSE]
chk <- m4[, {
  fr <- colMeans(cw[.I, , drop = FALSE])
  o <- order(-fr)
  tct <- sort(table(cell_type), decreasing = TRUE)
  list(n_nuclei = .N, main_cell_type = names(tct)[1], main_cell_type_frac = round(tct[[1]] / .N, 3),
       main_lineage = names(fr)[o[1]], main_lineage_frac = round(fr[o[1]], 3),
       second_lineage = names(fr)[o[2]], frac_second_lineage = round(fr[o[2]], 3),
       frac_doublet = round(mean(n_lin >= 2 | scdbl), 3), n_rescued = sum(rescue_status == "rescued"),
       n_lymphocyte = sum(rescue_status == "lymphocyte"))
}, by = cluster_5c]
chk[, remove := frac_doublet > 0.5 & frac_second_lineage > 0.5]
chk <- chk[order(as.integer(cluster_5c))]
fwrite(chk, file.path(table_dir, "new_cluster_check.csv"))
print(chk, nrows = 100)
rm_cl <- chk[remove == TRUE, cluster_5c]
message(sprintf("STEP 4: %d mixed cluster(s) removed: %s (%d nuclei)", length(rm_cl),
                paste(rm_cl, collapse = ", "), sum(chk[remove == TRUE, n_nuclei])))
m4[, removed_step4 := cluster_5c %in% rm_cl]

# ---- STEP 5: final object and outputs ----
o4$cell_type_final <- factor(m4$cell_type, levels = CT_LEVELS)
o4$rescue_status <- factor(m4$rescue_status, levels = c("kept", "rescued", "lymphocyte"))
o4$cb_cluster_step15 <- m4$cb_cluster
o4$lineage_carried <- m4$combo
o4$removed_step4 <- m4$removed_step4
if (length(rm_cl)) o4 <- subset(o4, cells = m4[removed_step4 == FALSE, cell])
o4$cluster_5c <- droplevels(o4$cluster_5c)
stopifnot(!anyNA(o4$cell_type_final))
saveRDS(o4, final_rds)
message("saved ", final_rds, " (", ncol(o4), " nuclei)")

md[, `:=`(in_final = cell %in% colnames(o4))]
md[cell %in% m4$cell, cluster_5c := m4$cluster_5c[match(cell, m4$cell)]]
fwrite(md[, .(cell, sample, cb_cluster, role, lineage, combo, n_lym, heme, lymphocyte, candidate, knn_top,
              knn_frac, rescued, cell_type, rescue_status, cluster_5c, in_final)],
       file.path(dirname(final_rds), "rescue_per_cell.csv.gz"))
fc <- dcast(as.data.table(o4@meta.data)[, .N, by = .(cell_type_final, orig.ident)], cell_type_final ~ orig.ident, value.var = "N", fill = 0)
fc[, total := rowSums(.SD), .SDcols = -1]
fwrite(fc, file.path(table_dir, "final_celltype_counts.csv"))
print(fc)

# UMAP figures
pal <- c(Cardiomyocyte = "#c0392b", Fibroblast = "#2e86c1", Endothelial = "#27ae60", Endocardial = "#8e44ad",
         Lymphatic_EC = "#16a085", Pericyte = "#d68910", SMC = "#7f8c8d", Macrophage = "#e84393",
         Lymphocyte = "#000000", Schwann = "#a0522d", Epicardial = "#f1c40f")
ud <- data.table(Embeddings(o4, "umap"), cell_type = o4$cell_type_final, status = o4$rescue_status,
                 cluster = o4$cluster_5c)
setnames(ud, 1:2, c("u1", "u2"))
lab <- ud[, .(u1 = median(u1), u2 = median(u2)), by = cell_type]
th <- theme_classic(base_size = 10) + theme(legend.position = "right")
p1 <- ggplot(ud, aes(u1, u2)) +
  geom_point(data = ud[status == "kept"], aes(colour = cell_type), size = 0.1, shape = 16, alpha = 0.6) +
  geom_point(data = ud[status != "kept"], aes(colour = cell_type, shape = status), size = 1, stroke = 0.3) +
  scale_colour_manual(values = pal, name = "cell_type_final") +
  scale_shape_manual(values = c(rescued = 17, lymphocyte = 15), name = "rescued / added") +
  geom_label(data = lab, aes(label = cell_type), size = 2.6, linewidth = 0, alpha = 0.7) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)), shape = guide_legend(override.aes = list(size = 3))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = sprintf("Final annotation (%s nuclei)", format(ncol(o4), big.mark = ",")),
       subtitle = sprintf("dots = kept; triangles = rescued (%d); squares = lymphocytes (%d)",
                          sum(ud$status == "rescued"), sum(ud$status == "lymphocyte"))) + th
p2 <- ggplot(ud, aes(u1, u2)) +
  geom_point(data = ud[status == "kept"], colour = "grey80", size = 0.05, shape = 16) +
  geom_point(data = ud[status != "kept"], aes(colour = status), size = 0.4, shape = 16) +
  scale_colour_manual(values = c(rescued = "#d35400", lymphocyte = "#000000"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Rescue status", subtitle = "grey = kept (5B confirmed)") + th
cl_lab <- ud[, .(u1 = median(u1), u2 = median(u2)), by = cluster]
p3 <- ggplot(ud, aes(u1, u2, colour = cluster)) + geom_point(size = 0.05, shape = 16) +
  geom_text(data = cl_lab, aes(label = cluster), colour = "black", size = 2.6) +
  guides(colour = "none") + labs(x = "UMAP 1", y = "UMAP 2", title = "Leiden clusters after step 5C",
                                  subtitle = if (length(rm_cl)) paste("removed mixed clusters:", paste(rm_cl, collapse = ", ")) else "no mixed cluster removed") + th
pdf(file.path(fig_dir, "rescue_umap.pdf"), width = 9, height = 7)
print(p1); print(p2); print(p3)
invisible(dev.off())
ggsave(file.path(fig_dir, "rescue_umap_celltype.png"), p1, width = 9, height = 7, dpi = 150)
ggsave(file.path(fig_dir, "rescue_umap_status.png"), p2, width = 9, height = 7, dpi = 150)

# pipeline state
st <- if (file.exists(state_json)) read_json(state_json, simplifyVector = TRUE) else list(completed = character(0))
st$completed <- unique(c(as.character(st$completed), "step_5C"))
st$updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
st$step_5C <- list(script = "script/17_rescue_doublets.R", final_rds = final_rds, n_nuclei = ncol(o4),
                   n_rescued = sum(o4$rescue_status == "rescued"), n_lymphocyte = sum(o4$rescue_status == "lymphocyte"),
                   doublet_clusters = I(dbl_cl), lowq_clusters = I(lowq), removed_step4 = I(rm_cl))
write_json(st, state_json, auto_unbox = TRUE, pretty = TRUE)
message("done")
