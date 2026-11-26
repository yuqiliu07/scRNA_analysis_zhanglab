#!/usr/bin/env Rscript
# Build a merged Seurat v5 object for FSU_seq samples A-H.
#
#  1. Read each sample's 10x raw matrix (<CR_DIR>/<s>/outs/raw_feature_bc_matrix.h5)
#  2. Keep only nuclei with pass_qc_singlet == TRUE
#     (<QC_DIR>/<s>_qc_doublet_flags.csv: pass_qc & !predicted_doublet)
#  3. CreateSeuratObject per sample (cell names "<sample>_<barcode>"), add QC
#     flags + sample_sheet metadata, merge all samples
#     -> RNA assay keeps one raw "counts.<sample>" layer per sample
#  4. Gene filter: expressed (count > 0) in >= MIN_CELLS nuclei across all samples
#  5. NormalizeData (LogNormalize, 1e4) -> "data.<sample>" layers
#  6. FindVariableFeatures (vst, nfeatures = 3000) separately per sample,
#     VariableFeatures(obj) = union of the per-sample HVGs
#
# Outputs (<OUT_DIR>/):
#   fsu_merged_mm10.rds               Seurat object (data file; stays on Hoffman)
#   seurat_object_dims.csv            genes / nuclei / HVG union size
#   seurat_cells_per_sample.csv       nuclei per sample through each step
#   seurat_hvg_per_sample.csv         per-sample HVG lists (long format)
#   seurat_hvg_union.csv              union HVG list with # samples selecting it
#
# Usage: Rscript 07_seurat_merge.R
# Env (optional): CR_DIR, QC_DIR, OUT_DIR, SAMPLE_SHEET, MIN_CELLS (10), N_HVG (3000)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DropletUtils)
  library(Matrix)
  library(data.table)
})
options(Seurat.object.assay.version = "v5")

scratch      <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
cr_dir       <- Sys.getenv("CR_DIR",  file.path(scratch, "FSU_seq/02_cellranger/shared_h5"))
qc_dir       <- Sys.getenv("QC_DIR",  file.path(scratch, "FSU_seq/work_dir/qc_mm10"))
out_dir      <- Sys.getenv("OUT_DIR", file.path(scratch, "FSU_seq/work_dir/seurat_mm10"))
sample_sheet <- Sys.getenv("SAMPLE_SHEET", file.path(scratch, "FSU_seq/sample_sheet.csv"))
min_cells    <- as.integer(Sys.getenv("MIN_CELLS", "10"))
n_hvg        <- as.integer(Sys.getenv("N_HVG", "3000"))
samples      <- LETTERS[1:8]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message(sprintf("Seurat %s | CR_DIR=%s | QC_DIR=%s | OUT_DIR=%s | MIN_CELLS=%d | N_HVG=%d",
                packageVersion("Seurat"), cr_dir, qc_dir, out_dir, min_cells, n_hvg))

ss <- fread(sample_sheet)

# ---- 1-3. per-sample objects ----
objs <- list(); cell_tab <- list(); gene_ids <- NULL
for (s in samples) {
  flags <- fread(file.path(qc_dir, paste0(s, "_qc_doublet_flags.csv")))
  keep  <- flags[toupper(as.character(pass_qc_singlet)) == "TRUE", barcode]

  sce <- read10xCounts(file.path(cr_dir, s, "outs/raw_feature_bc_matrix.h5"), col.names = TRUE)
  sce <- sce[, colnames(sce) %in% keep]
  stopifnot(ncol(sce) == length(keep))
  m <- as(counts(sce), "dgCMatrix")
  sym <- make.unique(as.character(rowData(sce)$Symbol))
  if (is.null(gene_ids)) gene_ids <- data.table(gene = sym, gene_id = rowData(sce)$ID)
  stopifnot(identical(gene_ids$gene, sym))       # same feature order in all samples
  rownames(m) <- sym
  colnames(m) <- paste0(s, "_", colnames(m))
  rm(sce); invisible(gc())

  o <- CreateSeuratObject(counts = m, project = "FSU_seq", min.cells = 0, min.features = 0)
  md <- as.data.frame(flags[match(sub("^[A-H]_", "", colnames(o)), barcode)])
  rownames(md) <- colnames(o)
  md <- md[, setdiff(names(md), "barcode"), drop = FALSE]
  sm <- ss[sample == s]
  md$genotype <- sm$genotype; md$sex <- sm$sex; md$mouse_id <- sm$mouse_id; md$age <- sm$age
  o <- AddMetaData(o, md)
  o$orig.ident <- s
  objs[[s]] <- o
  cell_tab[[s]] <- data.table(sample = s, genotype = sm$genotype, sex = sm$sex,
                              n_emptydrops = nrow(flags),
                              n_pass_qc = sum(toupper(as.character(flags$pass_qc)) == "TRUE"),
                              n_pass_qc_singlet = length(keep))
  message(sprintf("[%s] %d nuclei kept (of %d EmptyDrops nuclei)", s, ncol(o), nrow(flags)))
  rm(m, o); invisible(gc())
}

obj <- merge(objs[[1]], y = objs[-1], project = "FSU_seq")
rm(objs); invisible(gc())
message("Layers after merge: ", paste(Layers(obj), collapse = ", "))

# ---- 4. gene filter: >= min_cells nuclei across all samples ----
count_layers <- Layers(obj, search = "counts")
n_expr <- Reduce(`+`, lapply(count_layers, function(l) {
  x <- as(LayerData(obj, layer = l), "CsparseMatrix")
  x <- drop0(x)
  tabulate(x@i + 1L, nbins = nrow(x))               # nuclei with count > 0, per gene
}))
names(n_expr) <- rownames(obj)
keep_genes <- names(n_expr)[n_expr >= min_cells]
n_genes_before <- nrow(obj)
obj <- subset(obj, features = keep_genes)
message(sprintf("Gene filter (>= %d nuclei): %d -> %d genes", min_cells, n_genes_before, nrow(obj)))

# ---- 5. normalize (per layer) ----
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)

# ---- 6. per-sample HVGs, union ----
hvg_list <- list()
for (l in count_layers) {
  s <- sub("^counts\\.", "", l)
  tmp <- CreateSeuratObject(counts = LayerData(obj, layer = l))
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  hvg_list[[s]] <- VariableFeatures(tmp)
  rm(tmp)
}
# also run Seurat's per-layer HVG step so HVFInfo is present on the merged object
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
hvg_union <- unique(unlist(hvg_list))
VariableFeatures(obj) <- hvg_union
obj@misc$hvg_per_sample <- hvg_list
obj@misc$gene_ids <- gene_ids[gene %in% rownames(obj)]
obj@misc$params <- list(cr_dir = cr_dir, qc_dir = qc_dir, min_cells = min_cells, n_hvg = n_hvg,
                        hvg_method = "vst per sample, union",
                        reference = "refdata-gex-mm10-2020-A (Novogene, Cell Ranger 8.0.1)")
message(sprintf("HVG union: %d genes", length(hvg_union)))

# ---- outputs ----
saveRDS(obj, file.path(out_dir, "fsu_merged_mm10.rds"))

dims <- data.table(n_genes_raw = n_genes_before, n_genes_filtered = nrow(obj),
                   n_nuclei = ncol(obj), min_cells = min_cells,
                   n_hvg_per_sample = n_hvg, n_hvg_union = length(hvg_union),
                   layers = paste(Layers(obj), collapse = ";"))
fwrite(dims, file.path(out_dir, "seurat_object_dims.csv"))

ct <- rbindlist(cell_tab)
ct[, n_in_object := as.integer(table(obj$orig.ident)[sample])]
ct[, pct_of_emptydrops := round(100 * n_in_object / n_emptydrops, 2)]
ct <- rbind(ct, data.table(sample = "total", genotype = "", sex = "",
                           n_emptydrops = sum(ct$n_emptydrops), n_pass_qc = sum(ct$n_pass_qc),
                           n_pass_qc_singlet = sum(ct$n_pass_qc_singlet),
                           n_in_object = sum(ct$n_in_object),
                           pct_of_emptydrops = round(100 * sum(ct$n_in_object) / sum(ct$n_emptydrops), 2)))
fwrite(ct, file.path(out_dir, "seurat_cells_per_sample.csv"))

fwrite(rbindlist(lapply(names(hvg_list), function(s)
  data.table(sample = s, rank = seq_along(hvg_list[[s]]), gene = hvg_list[[s]]))),
  file.path(out_dir, "seurat_hvg_per_sample.csv"))
n_sel <- table(unlist(hvg_list))
fwrite(data.table(gene = hvg_union, n_samples_selected = as.integer(n_sel[hvg_union]))[order(-n_samples_selected)],
       file.path(out_dir, "seurat_hvg_union.csv"))

cat("\n===== merged object =====\n"); print(obj)
cat("\ndim(obj): ", paste(dim(obj), collapse = " genes x "), " nuclei\n", sep = "")
cat("\n===== nuclei per sample =====\n"); print(ct)
cat("\n===== HVG =====\n")
print(sapply(hvg_list, length)); cat("union:", length(hvg_union), "\n")
message("done")
