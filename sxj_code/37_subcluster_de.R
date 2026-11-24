# WT vs TNNC1-I4M within each subcluster: pseudobulk DE, specificity against the rest of the cell type, ORA.
#
# Input: <work_dir>/subclusters/pbsub_<type>_<k>.csv.gz, pbrest_<type>_<k>.csv.gz, pbsub_index.csv
#        (36_subcluster_pseudobulk.py); <work_dir>/genotype/meta.csv, pb_<type>.csv.gz, de_<type>.csv.gz (13).
# DE: as 13_pseudobulk_de.R: edgeR QL, ~ sex + genotype, filterByExpr, FDR <= 0.1 and |log2FC| >= 0.25.
#     Flags: technical (mt-, ribosomal, ambient / low-quality signature), ambient (>= 5x mean CPM in another major cell
#     type than in this subcluster), shared (flagged shared at the cell-type level in 13). Robust = DEG without flags.
# Comparison with the cell-type level (13): logFC correlation, how many cell-type robust DEGs change in the same
#     direction in the subcluster, and DEGs found only in the subcluster.
# Specificity: per mouse, log-CPM(subcluster) - log-CPM(rest of the cell type), limma ~ sex + genotype (as in
#     19_internal_control.R). A robust DEG is "subcluster-specific" if this difference also changes (FDR <= 0.1,
#     |log2FC| >= 0.25, same sign), i.e. the genotype effect is larger in this subcluster than in its neighbours;
#     otherwise "shared in cell type".
# ORA: hypergeometric test of robust up / down DEGs (>= 15) against tested genes, KEGG / Reactome / WikiPathways /
#     GO BP (GMT files from 27_export_gmt.R), sets of 10-500 tested genes, BH within each comparison.
#
# Outputs: <work_dir>/subclusters/sde_<type>_<k>.csv.gz, sora_<type>_<k>_<db>_<Up|Down>.csv.gz (full tables, local)
#          <results_dir>/subclusters/subcluster_de_summary.csv, subcluster_de_top_genes.csv, subcluster_ora_top.csv
# Usage: Rscript 37_subcluster_de.R <work_dir> <results_dir>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
args <- commandArgs(TRUE)
gd <- file.path(args[1], "genotype"); sd <- file.path(args[1], "subclusters"); od <- file.path(args[2], "subclusters")
dir.create(od, showWarnings = FALSE, recursive = TRUE)
meta <- read.csv(file.path(gd, "meta.csv"), row.names = 1)
meta$genotype <- factor(meta$genotype, levels = c("WT", "I4M")); meta$sex <- factor(meta$sex)
design <- model.matrix(~ sex + genotype, data = meta)
FDR <- 0.1; LFC <- 0.25
TECH <- c("Gm42418", "Cmss1", "Lars2", "Camk1d", "Gphn", "Gm19951", "Cdk8", "AY036118", "Gm26917", "Malat1")
is_tech <- function(g) g %in% TECH | grepl("^mt-", g) | grepl("^Rp[ls][0-9]", g)
rd <- function(f) as.matrix(read.csv(gzfile(f), row.names = 1, check.names = FALSE))[, rownames(meta)]

major <- sub("^pb_(.*)[.]csv[.]gz$", "\\1", list.files(gd, "^pb_.*[.]csv[.]gz$"))
cpm_mean <- sapply(major, function(t) { x <- rd(file.path(gd, paste0("pb_", t, ".csv.gz"))); rowMeans(t(t(x) / colSums(x)) * 1e6) })

read_gmt <- function(f) { l <- strsplit(readLines(f), "\t"); setNames(lapply(l, `[`, -(1:2)), vapply(l, `[`, "", 1)) }
DB <- list(KEGG = "KEGG", Reactome = "Reactome", WikiPathways = "WikiPathways", GO_BP = "GO_BP")
DB <- lapply(DB, function(k) read_gmt(file.path(gd, paste0("gmt_", k, ".gmt"))))
ora <- function(genes, universe, sets) {
  N <- length(universe); n <- length(genes)
  r <- do.call(rbind, lapply(names(sets), function(s) {
    g <- intersect(sets[[s]], universe); K <- length(g)
    if (K < 10 || K > 500) return(NULL)
    k <- length(intersect(g, genes))
    data.frame(set = s, size = K, overlap = k, expected = n * K / N, p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               genes = paste(head(sort(intersect(g, genes)), 15), collapse = ","))
  }))
  r$padj <- p.adjust(r$p, "BH"); r$fold_enrichment <- r$overlap / r$expected
  r[order(r$p), ]
}
lcpm <- function(x) { y <- normLibSizes(DGEList(x)); list(log = cpm(y, log = TRUE, prior.count = 2), lin = cpm(y)) }

idx <- read.csv(file.path(sd, "pbsub_index.csv"), colClasses = c(cluster = "character"))
idx <- idx[as.character(idx$tested) %in% c("True", "TRUE"), ]  # written by pandas
summ <- list(); top <- list(); otop <- list()
for (i in seq_len(nrow(idx))) {
  t <- idx$cell_type[i]; k <- idx$cluster[i]; nm <- idx$name[i]
  xs <- rd(file.path(sd, sprintf("pbsub_%s_%s.csv.gz", t, k))); xr <- rd(file.path(sd, sprintf("pbrest_%s_%s.csv.gz", t, k)))
  par <- read.csv(gzfile(file.path(gd, paste0("de_", t, ".csv.gz"))))
  # ---- DE within the subcluster ----
  y <- DGEList(xs)
  keep <- filterByExpr(y, group = meta$genotype)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  tt <- topTags(glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = "genotypeI4M"), n = Inf, sort.by = "none")$table
  d <- data.frame(gene = rownames(tt), logCPM = tt$logCPM, logFC = tt$logFC, PValue = tt$PValue, FDR = tt$FDR)
  cs <- rowMeans(t(t(xs) / colSums(xs)) * 1e6)[d$gene]
  oth <- cpm_mean[d$gene, setdiff(colnames(cpm_mean), t), drop = FALSE]; oth[is.na(oth)] <- 0
  d$technical <- is_tech(d$gene)
  d$ambient <- apply(oth, 1, max) >= 5 * cs
  d$shared <- d$gene %in% par$gene[par$shared]
  d$sig <- d$FDR <= FDR & abs(d$logFC) >= LFC
  d$robust <- d$sig & !d$technical & !d$ambient & !d$shared
  m <- match(d$gene, par$gene)
  d$logFC_celltype <- par$logFC[m]; d$FDR_celltype <- par$FDR[m]; d$robust_celltype <- par$robust[m] %in% TRUE
  # ---- specificity: subcluster vs rest of the cell type, within mouse ----
  cs_ <- lcpm(xs); cr_ <- lcpm(xr)
  g2 <- d$gene[d$gene %in% rownames(xr)]
  mr <- rowMeans(cr_$lin)[g2]; ms <- rowMeans(cs_$lin)[g2]
  g2 <- g2[mr >= 1 & mr >= 0.2 * ms]
  dd <- cs_$log[g2, ] - cr_$log[g2, ]
  st <- topTable(eBayes(lmFit(dd, design), trend = TRUE, robust = TRUE), coef = "genotypeI4M", n = Inf, sort.by = "none")
  ms2 <- match(d$gene, rownames(st))
  d$logFC_vs_rest <- st$logFC[ms2]; d$FDR_vs_rest <- st$adj.P.Val[ms2]
  spec <- !is.na(d$FDR_vs_rest) & d$FDR_vs_rest <= FDR & abs(d$logFC_vs_rest) >= LFC & sign(d$logFC_vs_rest) == sign(d$logFC)
  d$specificity <- ifelse(!d$robust, NA, ifelse(spec, "subcluster-specific", ifelse(is.na(d$FDR_vs_rest), "not_testable", "shared in cell type")))
  write.csv(d, gzfile(file.path(sd, sprintf("sde_%s_%s.csv.gz", t, k))), row.names = FALSE)
  # ---- summary ----
  rb <- d$robust
  pc <- d$robust_celltype & !is.na(d$logFC_celltype)
  summ[[i]] <- data.frame(cell_type = t, cluster = k, subcluster = nm, min_nuclei = idx$min_nuclei[i],
    total_nuclei = idx$total_nuclei[i], n_tested = nrow(d), deg = sum(d$sig), robust = sum(rb),
    up = sum(rb & d$logFC > 0), down = sum(rb & d$logFC < 0),
    also_robust_in_celltype = sum(rb & d$robust_celltype), only_in_subcluster = sum(rb & !d$robust_celltype),
    celltype_robust_tested = sum(pc), celltype_robust_same_sign = sum(pc & sign(d$logFC) == sign(d$logFC_celltype)),
    logFC_cor_celltype = round(cor(d$logFC, d$logFC_celltype, use = "complete.obs"), 3),
    logFC_cor_celltype_robust = round(cor(d$logFC[pc], d$logFC_celltype[pc]), 3),
    specific = sum(d$specificity %in% "subcluster-specific"),
    specific_up = sum(d$specificity %in% "subcluster-specific" & d$logFC > 0),
    specific_down = sum(d$specificity %in% "subcluster-specific" & d$logFC < 0))
  print(summ[[i]][, c("cell_type", "subcluster", "robust", "up", "down", "only_in_subcluster", "specific", "logFC_cor_celltype")])
  r <- d[rb, ]; r <- r[order(r$PValue), ]
  if (nrow(r)) top[[i]] <- cbind(cell_type = t, cluster = k, subcluster = nm, head(r[, c("gene", "logFC", "FDR",
                           "logFC_celltype", "FDR_celltype", "robust_celltype", "logFC_vs_rest", "FDR_vs_rest", "specificity")], 50))
  # ---- ORA ----
  for (dir_ in c("Up", "Down")) {
    g <- d$gene[rb & (if (dir_ == "Up") d$logFC > 0 else d$logFC < 0)]
    if (length(g) < 15) next
    for (db in names(DB)) {
      o <- ora(g, d$gene, DB[[db]])
      write.csv(o, gzfile(file.path(sd, sprintf("sora_%s_%s_%s_%s.csv.gz", t, k, db, dir_))), row.names = FALSE)
      tp <- head(o, 5)
      otop[[length(otop) + 1]] <- data.frame(cell_type = t, cluster = k, subcluster = nm, database = db, direction = dir_,
        n_input = length(g), n_sig = sum(o$padj <= 0.05), set = tp$set, size = tp$size, overlap = tp$overlap,
        p = tp$p, padj = tp$padj, genes = tp$genes)
    }
    cat("   ORA", dir_, length(g), "genes; sig sets:",
        paste(sapply(names(DB), function(db) sum(read.csv(gzfile(file.path(sd, sprintf("sora_%s_%s_%s_%s.csv.gz", t, k, db, dir_))))$padj <= 0.05)), collapse = "/"), "\n")
  }
}
write.csv(do.call(rbind, summ), file.path(od, "subcluster_de_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, top), file.path(od, "subcluster_de_top_genes.csv"), row.names = FALSE)
if (length(otop)) write.csv(do.call(rbind, otop), file.path(od, "subcluster_ora_top.csv"), row.names = FALSE)
