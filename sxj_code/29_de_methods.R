# WT vs TNNC1-I4M: compare pseudobulk DE methods; also report the stricter FDR <= 0.05 list of the main method.
#
# Same input, gene filter (filterByExpr), TMM normalisation and model (~ sex + genotype) as 13_pseudobulk_de.R; only
# the test changes:
#   edgeR_QL    edgeR quasi-likelihood F-test (main analysis in 13_pseudobulk_de.R)
#   edgeR_LRT   edgeR likelihood-ratio test (less conservative: ignores uncertainty of the dispersion estimate)
#   voom        limma-voom, robust empirical Bayes
#   voom_qw     limma-voom with sample quality weights (down-weights noisier mice)
# DEG: FDR <= 0.1 and |log2FC| >= 0.25 (main cut-off, as in 13). Flags as in 13 (technical, ambient; shared = DE in the same direction in
# >= 3 cell types for that method). Robust = DEG without flags.
# Consensus: robust in >= 3 of the 4 methods.
# Strict list: robust DEGs of edgeR_QL that also have FDR <= 0.05 (the previous main cut-off).
#
# Outputs: <work_dir>/genotype/dem_<type>.csv.gz   per-gene logFC / P / FDR of every method (local)
#          <work_dir>/genotype/strict_<type>.csv.gz robust DEGs with FDR <= 0.05 (local)
#          <results_dir>/genotype/de_methods_summary.csv, de_methods_overlap.csv (counts only)
# Usage: Rscript 29_de_methods.R <work_dir> <results_dir>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
meta <- read.csv(file.path(wd, "meta.csv"), row.names = 1)
meta$genotype <- factor(meta$genotype, levels = c("WT", "I4M")); meta$sex <- factor(meta$sex)
design <- model.matrix(~ sex + genotype, data = meta)
FDR <- 0.1; LFC <- 0.25
METHODS <- c("edgeR_QL", "edgeR_LRT", "voom", "voom_qw")
TECH <- c("Gm42418", "Cmss1", "Lars2", "Camk1d", "Gphn", "Gm19951", "Cdk8", "AY036118", "Gm26917", "Malat1")
is_tech <- function(g) g %in% TECH | grepl("^mt-", g) | grepl("^Rp[ls][0-9]", g)
types <- sub("^pb_(.*)[.]csv[.]gz$", "\\1", list.files(wd, "^pb_.*[.]csv[.]gz$"))  # all types, as in 13 (needed for the shared flag)
pb <- lapply(setNames(types, types), function(t)
  as.matrix(read.csv(gzfile(file.path(wd, paste0("pb_", t, ".csv.gz"))), row.names = 1, check.names = FALSE))[, rownames(meta)])
allg <- Reduce(union, lapply(pb, rownames))
cpm_mean <- sapply(pb, function(x) { m <- rowMeans(t(t(x) / colSums(x)) * 1e6); m[allg] })
cpm_mean[is.na(cpm_mean)] <- 0; rownames(cpm_mean) <- allg

res <- list()
for (t in types) {
  y <- DGEList(pb[[t]])
  keep <- filterByExpr(y, group = meta$genotype)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  ql <- topTags(glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = "genotypeI4M"), n = Inf, sort.by = "none")$table
  lrt <- topTags(glmLRT(glmFit(y, design), coef = "genotypeI4M"), n = Inf, sort.by = "none")$table
  v <- voom(y, design)
  vo <- topTable(eBayes(lmFit(v, design), robust = TRUE), coef = "genotypeI4M", n = Inf, sort.by = "none")
  vq <- voomWithQualityWeights(y, design)
  vqo <- topTable(eBayes(lmFit(vq, design), robust = TRUE), coef = "genotypeI4M", n = Inf, sort.by = "none")
  d <- data.frame(gene = rownames(y), logCPM = ql$logCPM,
                  edgeR_QL_logFC = ql$logFC, edgeR_QL_P = ql$PValue, edgeR_QL_FDR = ql$FDR,
                  edgeR_LRT_logFC = lrt$logFC, edgeR_LRT_P = lrt$PValue, edgeR_LRT_FDR = lrt$FDR,
                  voom_logFC = vo$logFC, voom_P = vo$P.Value, voom_FDR = vo$adj.P.Val,
                  voom_qw_logFC = vqo$logFC, voom_qw_P = vqo$P.Value, voom_qw_FDR = vqo$adj.P.Val)
  others <- setdiff(colnames(cpm_mean), t)
  d$technical <- is_tech(d$gene)
  d$ambient <- apply(cpm_mean[d$gene, others, drop = FALSE], 1, max) >= 5 * cpm_mean[d$gene, t]
  d$sample_weights <- paste(round(vq$targets$sample.weights, 2), collapse = ";")
  res[[t]] <- d
  cat(t, "weights:", d$sample_weights[1], "\n")
}

# shared flag per method
for (m in METHODS) {
  keys <- unlist(lapply(res, function(d) {
    s <- d[[paste0(m, "_FDR")]] <= FDR & abs(d[[paste0(m, "_logFC")]]) >= LFC
    paste(d$gene[s], sign(d[[paste0(m, "_logFC")]][s]))
  }))
  ns <- table(keys)
  for (t in types) {
    d <- res[[t]]
    s <- d[[paste0(m, "_FDR")]] <= FDR & abs(d[[paste0(m, "_logFC")]]) >= LFC
    k <- as.integer(ns[paste(d$gene, sign(d[[paste0(m, "_logFC")]]))])
    d[[paste0(m, "_sig")]] <- s
    d[[paste0(m, "_robust")]] <- s & !d$technical & !d$ambient & !(!is.na(k) & k >= 3)
    res[[t]] <- d
  }
}

summ <- list(); ov <- list()
for (t in types) {
  d <- res[[t]]
  rb <- sapply(METHODS, function(m) d[[paste0(m, "_robust")]])
  d$n_methods_robust <- rowSums(rb)
  d$consensus <- d$n_methods_robust >= 3
  # exploratory list from the main method
  e <- d$edgeR_QL_robust & d$edgeR_QL_FDR <= 0.05
  ex <- d[e, c("gene", "logCPM", "edgeR_QL_logFC", "edgeR_QL_P", "edgeR_QL_FDR", "edgeR_QL_robust", "n_methods_robust")]
  ex$tier <- rep("strict (FDR<=0.05)", nrow(ex))
  ex <- ex[order(ex$edgeR_QL_P), ]
  write.csv(ex, gzfile(file.path(wd, paste0("strict_", t, ".csv.gz"))), row.names = FALSE)
  write.csv(d, gzfile(file.path(wd, paste0("dem_", t, ".csv.gz"))), row.names = FALSE)
  for (m in METHODS) {
    s <- d[[paste0(m, "_sig")]]; r <- d[[paste0(m, "_robust")]]
    summ[[paste(t, m)]] <- data.frame(cell_type = t, method = m, n_tested = nrow(d), deg = sum(s),
      up_robust = sum(r & d[[paste0(m, "_logFC")]] > 0), down_robust = sum(r & d[[paste0(m, "_logFC")]] < 0),
      robust = sum(r), overlap_with_edgeR_QL = sum(r & d$edgeR_QL_robust),
      logFC_cor_with_edgeR_QL = round(cor(d[[paste0(m, "_logFC")]], d$edgeR_QL_logFC), 3))
  }
  summ[[paste(t, "consensus")]] <- data.frame(cell_type = t, method = "consensus_3of4", n_tested = nrow(d),
    deg = NA, up_robust = sum(d$consensus & d$edgeR_QL_logFC > 0), down_robust = sum(d$consensus & d$edgeR_QL_logFC < 0),
    robust = sum(d$consensus), overlap_with_edgeR_QL = sum(d$consensus & d$edgeR_QL_robust), logFC_cor_with_edgeR_QL = NA)
  summ[[paste(t, "explore")]] <- data.frame(cell_type = t, method = "edgeR_QL_FDR0.05", n_tested = nrow(d),
    deg = sum(d$edgeR_QL_FDR <= 0.05 & abs(d$edgeR_QL_logFC) >= LFC),
    up_robust = sum(e & d$edgeR_QL_logFC > 0), down_robust = sum(e & d$edgeR_QL_logFC < 0),
    robust = sum(e), overlap_with_edgeR_QL = sum(e & d$edgeR_QL_robust), logFC_cor_with_edgeR_QL = NA)
  # overlap pattern of robust sets
  pat <- apply(rb, 1, function(z) paste(METHODS[z], collapse = "+"))
  tb <- table(pat[pat != ""])
  if (length(tb)) ov[[t]] <- data.frame(cell_type = t, methods = names(tb), n = as.integer(tb))
}
S <- do.call(rbind, summ)
write.csv(S, file.path(od, "de_methods_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, ov), file.path(od, "de_methods_overlap.csv"), row.names = FALSE)
print(S[, c("cell_type", "method", "deg", "robust", "overlap_with_edgeR_QL", "logFC_cor_with_edgeR_QL")], row.names = FALSE)
