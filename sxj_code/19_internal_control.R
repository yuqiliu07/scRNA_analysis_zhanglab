# Internal-control analysis: is the genotype effect in a target cell type larger than in a reference cell type
# from the same mice?  (cell type x genotype interaction, "difference in differences")
#
# Rationale: genotype is confounded with a processing batch. A batch effect acts on every cell type of a mouse,
# while TNNC1-I4M should act mainly on cardiomyocytes. Pseudobulks of the target and of a reference cell type from
# the same 8 mice are compared within each mouse (target minus reference), so a shift shared by both cell types
# of that mouse (batch, sex, depth) cancels; the genotype coefficient of the difference is
#   genotype effect in target - genotype effect in reference
# References: fibroblast and endothelial (abundant, almost no DEGs on their own).
# A gene can only be controlled if the reference itself expresses it (see fit_pair).
#
# Categories for each gene (target-vs-reference test done against both references):
#   confirmed       interaction FDR < 0.05 against both references, same sign as the target's own change
#   partly          FDR < 0.05 against one reference
#   shared          testable, FDR >= 0.05 against both: the reference changes too (batch-like)
#   not_testable    not expressed by either reference (mean CPM < 1 or < 20% of the target's)
#
# Outputs: <work_dir>/genotype/ic_<target>.csv.gz ; <results_dir>/genotype/ic_summary.csv, ic_top_genes.csv,
#          ic_go_top.csv (GO BP over-representation of confirmed genes, if >= 15 in a direction)
# Usage: Rscript 19_internal_control.R <work_dir> <results_dir>
suppressPackageStartupMessages({ library(edgeR); library(limma); library(AnnotationDbi); library(org.Mm.eg.db) })
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
meta <- read.csv(file.path(wd, "meta.csv"), row.names = 1)
FDR <- 0.05
targets <- c("Cardiomyocyte", "Endocardial")
refs <- c("Fibroblast", "Endothelial")
pb <- function(t) as.matrix(read.csv(gzfile(file.path(wd, paste0("pb_", t, ".csv.gz"))), row.names = 1,
                                     check.names = FALSE))[, rownames(meta)]

# Per-mouse difference of log-CPM (target - reference), tested with limma: d ~ sex + genotype. Each cell type is
# normalised on its own and the gene-wise variance of the difference is estimated directly (empirical Bayes,
# trend), so the different noise levels of the two cell types are handled.
# Only "testable" genes are tested: expressed by the reference itself (mean CPM >= 1 and >= 20% of the target's),
# so that the reference signal is not mainly ambient RNA carried over from the target.
lcpm <- function(x) { y <- normLibSizes(DGEList(x)); list(log = cpm(y, log = TRUE, prior.count = 2), lin = cpm(y)) }
fit_pair <- function(tg, rf) {
  xt <- pb(tg); xr <- pb(rf)
  ct <- lcpm(xt); cr <- lcpm(xr)
  mt <- rowMeans(ct$lin); mr <- rowMeans(cr$lin)
  keep <- filterByExpr(DGEList(xt), group = meta$genotype) & mr >= 1 & mr >= 0.2 * mt
  d <- ct$log[keep, ] - cr$log[keep, ]
  des <- model.matrix(~ sex + genotype, data = transform(meta, genotype = factor(genotype, levels = c("WT", "I4M"))))
  fit <- eBayes(lmFit(d, des), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "genotypeI4M", n = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC_int = tt$logFC, P_int = tt$P.Value, FDR_int = tt$adj.P.Val,
             ref_cpm = mr[rownames(tt)])
}

summ <- list(); top <- list(); gotop <- list()
eg_all <- NULL
for (tg in targets) {
  de <- read.csv(gzfile(file.path(wd, paste0("de_", tg, ".csv.gz"))))
  res <- de[, c("gene", "logFC", "FDR", "robust", "technical", "ambient", "shared")]
  names(res)[2:3] <- c("logFC_target", "FDR_target")
  for (rf in refs) {
    r <- fit_pair(tg, rf)
    rd <- read.csv(gzfile(file.path(wd, paste0("de_", rf, ".csv.gz"))))[, c("gene", "logFC")]
    names(rd)[2] <- paste0("logFC_in_", rf)
    names(r)[-1] <- paste0(names(r)[-1], "_vs_", rf)
    res <- merge(res, r, by = "gene", all.x = TRUE)
    res <- merge(res, rd, by = "gene", all.x = TRUE)
  }
  sigsame <- function(rf) {
    f <- res[[paste0("FDR_int_vs_", rf)]]; l <- res[[paste0("logFC_int_vs_", rf)]]
    !is.na(f) & f <= FDR & abs(l) >= 0.25 & sign(l) == sign(res$logFC_target)
  }
  testable <- function(rf) !is.na(res[[paste0("FDR_int_vs_", rf)]])
  s1 <- sigsame(refs[1]); s2 <- sigsame(refs[2])
  anyt <- testable(refs[1]) | testable(refs[2])
  res$category <- ifelse(s1 & s2, "confirmed", ifelse(s1 | s2, "partly", ifelse(anyt, "shared", "not_testable")))
  # genes passing both references but not significant on their own (gained by the internal control)
  res$new_by_ic <- res$category == "confirmed" & !(res$FDR_target <= FDR)
  write.csv(res, gzfile(file.path(wd, paste0("ic_", tg, ".csv.gz"))), row.names = FALSE)
  rb <- res[res$robust, ]
  tab <- table(factor(rb$category, levels = c("confirmed", "partly", "shared", "not_testable")))
  summ[[tg]] <- data.frame(target = tg, n_robust_degs = nrow(rb), confirmed = tab[["confirmed"]],
                           partly = tab[["partly"]], shared = tab[["shared"]], not_testable = tab[["not_testable"]],
                           confirmed_up = sum(rb$category == "confirmed" & rb$logFC_target > 0),
                           confirmed_down = sum(rb$category == "confirmed" & rb$logFC_target < 0),
                           extra_confirmed_not_in_robust = sum(res$category == "confirmed" & !res$robust &
                                                                 !res$technical & !res$ambient))
  print(summ[[tg]])
  cf <- rb[rb$category == "confirmed", ]
  cf <- cf[order(cf$FDR_target), ]
  if (nrow(cf)) top[[tg]] <- head(cbind(target = tg, cf[, c("gene", "logFC_target", "FDR_target", "logFC_in_Fibroblast",
                                              "logFC_in_Endothelial", "FDR_int_vs_Fibroblast", "FDR_int_vs_Endothelial")]), 40)
  uni_sym <- res$gene[!is.na(res$FDR_int_vs_Fibroblast) | !is.na(res$FDR_int_vs_Endothelial)]
  eg <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, res$gene, "ENTREZID", "SYMBOL", multiVals = "first"))
  uni <- unique(na.omit(eg[uni_sym]))
  for (dd in c("Up", "Down")) {
    g <- cf$gene[if (dd == "Up") cf$logFC_target > 0 else cf$logFC_target < 0]
    if (length(g) < 15) next
    gi <- unique(na.omit(eg[g]))
    ga <- goana(gi, universe = uni, species = "Mm")
    ga <- ga[ga$Ont == "BP" & ga$N >= 10 & ga$N <= 500, ]
    ga$FDR <- p.adjust(ga$P.DE, "BH")
    ga <- head(ga[order(ga$P.DE), ], 15)
    gotop[[paste(tg, dd)]] <- data.frame(target = tg, direction = dd, n_deg = length(g), n_deg_entrez = length(gi),
                                         universe = length(uni), go_id = rownames(ga), term = ga$Term, N = ga$N,
                                         DE = ga$DE, P = ga$P.DE, FDR = ga$FDR)
    cat("  GO", tg, dd, length(g), "genes; FDR<0.05 terms:", sum(ga$FDR < 0.05), "\n")
  }
}
write.csv(do.call(rbind, summ), file.path(od, "ic_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, top), file.path(od, "ic_top_genes.csv"), row.names = FALSE)
if (length(gotop)) write.csv(do.call(rbind, gotop), file.path(od, "ic_go_top.csv"), row.names = FALSE)
