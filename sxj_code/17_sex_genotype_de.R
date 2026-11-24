# WT vs TNNC1-I4M by sex: is the genotype effect different in males and females?
#
# Same pseudobulk counts as 13_pseudobulk_de.R. Groups WT_M (A, B), WT_F (C, D), I4M_M (E, F), I4M_F (G, H);
# edgeR quasi-likelihood on ~ 0 + group, contrasts:
#   male:        I4M_M - WT_M
#   female:      I4M_F - WT_F
#   interaction: (I4M_F - WT_F) - (I4M_M - WT_M)   (does the genotype effect differ between sexes?)
# Every batch (genotype) contains both sexes, so the interaction is not confounded with the genotype-batch effect.
# With 2 mice per group the power is low: only large, consistent sex differences can reach FDR <= 0.1.
# DEGs flagged as technical / ambient / shared in 13_pseudobulk_de.R are excluded from the "robust" counts.
#
# Outputs: <work_dir>/genotype/sex_<type>.csv.gz (per-gene logFC / FDR for the three contrasts, local)
#          <results_dir>/genotype/sex_summary.csv
#
# Usage: Rscript 17_sex_genotype_de.R <work_dir> <results_dir>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
meta <- read.csv(file.path(wd, "meta.csv"), row.names = 1)
meta$group <- factor(paste(meta$genotype, meta$sex, sep = "_"), levels = c("WT_M", "WT_F", "I4M_M", "I4M_F"))
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
con <- makeContrasts(male = I4M_M - WT_M, female = I4M_F - WT_F,
                     interaction = (I4M_F - WT_F) - (I4M_M - WT_M), levels = design)
FDR <- 0.1  # same DE cut-off as 13_pseudobulk_de.R
LFC <- 0.25

types <- sub("^pb_(.*)[.]csv[.]gz$", "\\1", list.files(wd, "^pb_.*[.]csv[.]gz$"))
summ <- list()
for (t in types) {
  x <- as.matrix(read.csv(gzfile(file.path(wd, paste0("pb_", t, ".csv.gz"))), row.names = 1,
                          check.names = FALSE))[, rownames(meta)]
  y <- DGEList(x)
  keep <- filterByExpr(y, group = meta$group)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  flags <- read.csv(gzfile(file.path(wd, paste0("de_", t, ".csv.gz"))))
  rownames(flags) <- flags$gene
  res <- data.frame(gene = rownames(y))
  for (k in colnames(con)) {
    tt <- topTags(glmQLFTest(fit, contrast = con[, k]), n = Inf, sort.by = "none")$table
    res[[paste0("logFC_", k)]] <- tt$logFC
    res[[paste0("P_", k)]] <- tt$PValue
    res[[paste0("FDR_", k)]] <- tt$FDR
  }
  f <- flags[res$gene, ]
  res$flagged <- (!is.na(f$technical) & f$technical) | (!is.na(f$ambient) & f$ambient) | (!is.na(f$shared) & f$shared)
  res$overall_robust <- !is.na(f$robust) & f$robust
  res$overall_logFC <- f$logFC
  write.csv(res, gzfile(file.path(wd, paste0("sex_", t, ".csv.gz"))), row.names = FALSE)
  ok <- !res$flagged
  r_all <- cor(res$logFC_male, res$logFC_female)
  rb <- res$overall_robust
  summ[[t]] <- data.frame(
    cell_type = t, n_genes = nrow(res),
    male_up = sum(ok & res$FDR_male <= FDR & res$logFC_male >= LFC), male_down = sum(ok & res$FDR_male <= FDR & res$logFC_male <= -LFC),
    female_up = sum(ok & res$FDR_female <= FDR & res$logFC_female >= LFC),
    female_down = sum(ok & res$FDR_female <= FDR & res$logFC_female <= -LFC),
    interaction = sum(ok & res$FDR_interaction <= FDR & abs(res$logFC_interaction) >= LFC),
    n_overall_robust = sum(rb),
    same_direction_frac = if (sum(rb)) mean(sign(res$logFC_male[rb]) == sign(res$logFC_female[rb])) else NA,
    r_logFC_overall_degs = if (sum(rb) > 2) cor(res$logFC_male[rb], res$logFC_female[rb]) else NA,
    r_logFC_all_genes = r_all)
  print(summ[[t]])
}
write.csv(do.call(rbind, summ), file.path(od, "sex_summary.csv"), row.names = FALSE)
