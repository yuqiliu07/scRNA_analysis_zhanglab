# WT vs TNNC1-I4M (+/-): pseudobulk differential expression per cell type, technical flags, GO enrichment.
#
# Input: <work_dir>/genotype/pb_<type>.csv.gz (12_genotype_prep.py: summed CellBender counts per sample), meta.csv.
# Model (edgeR quasi-likelihood): ~ sex + genotype, 4 WT vs 4 I4M (2 M + 2 F each); genes kept by filterByExpr.
# DE: FDR <= 0.1 (Benjamini-Hochberg) and |log2FC| >= 0.25 (about 1.19-fold). Main cut-off set to FDR <= 0.1 at the
#     user's request (2026-11-15; was 0.05): about 10% of the listed genes are expected to be false discoveries.
# Technical flags (genotype is confounded with a library-quality difference: all WT samples have 2-3x more
# ambient / low-quality transcripts and fewer intronic reads than all I4M samples, so part of the DE is technical):
#   technical: mt-, ribosomal protein genes and the ambient / low-quality signature (Gm42418, Cmss1, Lars2, ...)
#   ambient:   gene >= 5x higher (mean CPM) in another cell type - likely carried over in ambient RNA
#   shared:    DE (FDR <= 0.1, |log2FC| >= 0.25) in the same direction in >= 3 cell types - a sample-level, not cell-type, effect
#   A DEG with none of these flags is "robust".
# GO: over-representation (limma::goana, BP, 10-500 genes, universe = tested genes) of robust up / down DEGs, for
#     cell types with >= 15 robust DEGs in that direction.
#
# Outputs: <work_dir>/genotype/de_<type>.csv.gz, go_<type>_<Up|Down>.csv.gz (full tables, kept local)
#          <results_dir>/genotype/de_summary.csv, go_top.csv, de_panel.csv (summary tables)
#
# Usage: Rscript 13_pseudobulk_de.R <work_dir> <results_dir>
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db)
})
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
meta <- read.csv(file.path(wd, "meta.csv"), row.names = 1)
meta$genotype <- factor(meta$genotype, levels = c("WT", "I4M"))
meta$sex <- factor(meta$sex)
design <- model.matrix(~ sex + genotype, data = meta)
FDR <- 0.1  # DE cut-off (user's choice, 2026-11-15; previously 0.05)
LFC <- 0.25  # |log2FC| cut-off (agreed with the user, 2026-11-13)
TECH <- c("Gm42418", "Cmss1", "Lars2", "Camk1d", "Gphn", "Gm19951", "Cdk8", "AY036118", "Gm26917", "Malat1")
is_tech <- function(g) g %in% TECH | grepl("^mt-", g) | grepl("^Rp[ls][0-9]", g)

# curated genes for the biological questions (logFC / FDR reported for each cell type)
panel <- c(
  # thin filament / sarcomere
  "Tnnc1", "Tnni3", "Tnnt2", "Tpm1", "Actc1", "Myh6", "Myh7", "Myl2", "Myl3", "Myl4", "Myl7", "Mybpc3", "Ttn", "Actn2", "Xirp2",
  # calcium handling
  "Atp2a2", "Pln", "Ryr2", "Cacna1c", "Slc8a1", "Casq2", "Camk2d", "Sln", "Trdn",
  # stress / fetal gene programme
  "Nppa", "Nppb", "Acta1", "Ankrd1", "Ccn2", "Rcan1", "Fhl1", "Hspa1a", "Atf3", "Myh7b",
  # metabolism
  "Ppargc1a", "Ppara", "Cpt1b", "Acadm", "Hadha", "Pdk4", "Ucp3", "Slc2a4", "Slc2a1", "Hk2", "Pfkp", "Ldha",
  # fibroblast activation
  "Postn", "Cthrc1", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Acta2", "Lox", "Timp1", "Tnc", "Comp", "Fap", "Meox1", "Ltbp2",
  # endothelial / immune
  "Nos3", "Vegfa", "Kdr", "Apln", "Aplnr", "Klf2", "Klf4", "Icam1", "Vcam1", "Sele",
  "Ccr2", "Lyve1", "Cd163", "Mrc1", "Il1b", "Tnf", "Spp1", "Trem2", "Cd68")

types <- sub("^pb_(.*)[.]csv[.]gz$", "\\1", list.files(wd, "^pb_.*[.]csv[.]gz$"))
pb <- lapply(setNames(types, types), function(t)
  as.matrix(read.csv(gzfile(file.path(wd, paste0("pb_", t, ".csv.gz"))), row.names = 1,
                     check.names = FALSE))[, rownames(meta)])
cpm_mean <- sapply(pb, function(x) rowMeans(t(t(x) / colSums(x)) * 1e6))
eg <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, rownames(cpm_mean), "ENTREZID", "SYMBOL",
                                             multiVals = "first"))

res <- list(); summ <- list(); pan <- list()
for (t in types) {
  y <- DGEList(pb[[t]])
  keep <- filterByExpr(y, group = meta$genotype)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = "genotypeI4M"), n = Inf)$table
  tt$gene <- rownames(tt)
  others <- setdiff(colnames(cpm_mean), t)
  tt$technical <- is_tech(tt$gene)
  tt$ambient <- apply(cpm_mean[tt$gene, others, drop = FALSE], 1, max) >= 5 * cpm_mean[tt$gene, t]
  res[[t]] <- tt
  p <- tt[intersect(panel, tt$gene), c("gene", "logFC", "logCPM", "PValue", "FDR")]
  if (nrow(p)) pan[[t]] <- cbind(cell_type = t, p)
  summ[[t]] <- data.frame(cell_type = t, n_genes_tested = nrow(tt), min_lib = min(y$samples$lib.size),
                          bcv = sqrt(y$common.dispersion))
}

# shared flag: DE in the same direction in >= 3 cell types
sig <- do.call(rbind, lapply(names(res), function(t) {
  d <- res[[t]][res[[t]]$FDR <= FDR & abs(res[[t]]$logFC) >= LFC, ]
  if (nrow(d)) data.frame(key = paste(d$gene, sign(d$logFC)))
}))
nshared <- table(sig$key)

gotop <- list()
for (t in types) {
  tt <- res[[t]]
  ns <- as.integer(nshared[paste(tt$gene, sign(tt$logFC))])
  sig <- tt$FDR <= FDR & abs(tt$logFC) >= LFC
  tt$shared <- !is.na(ns) & ns >= 3 & sig
  tt$robust <- sig & !tt$technical & !tt$ambient & !tt$shared
  write.csv(tt, gzfile(file.path(wd, paste0("de_", t, ".csv.gz"))), row.names = FALSE)
  s <- sig
  summ[[t]] <- cbind(summ[[t]], up = sum(s & tt$logFC > 0), down = sum(s & tt$logFC < 0),
                     up_robust = sum(tt$robust & tt$logFC > 0), down_robust = sum(tt$robust & tt$logFC < 0),
                     flag_technical = sum(s & tt$technical), flag_ambient = sum(s & tt$ambient & !tt$technical),
                     flag_shared = sum(s & tt$shared & !tt$technical & !tt$ambient))
  print(summ[[t]])
  uni <- unique(na.omit(eg[tt$gene]))
  for (d in c("Up", "Down")) {
    g <- tt$gene[tt$robust & (if (d == "Up") tt$logFC > 0 else tt$logFC < 0)]
    if (length(g) < 15) next
    ga <- goana(unique(na.omit(eg[g])), universe = uni, species = "Mm")
    ga <- ga[ga$Ont == "BP" & ga$N >= 10 & ga$N <= 500, ]
    ga$FDR <- p.adjust(ga$P.DE, "BH")
    ga <- ga[order(ga$P.DE), ]
    write.csv(cbind(go_id = rownames(ga), ga), gzfile(file.path(wd, paste0("go_", t, "_", d, ".csv.gz"))),
              row.names = FALSE)
    top <- head(ga, 15)
    gotop[[paste(t, d)]] <- data.frame(cell_type = t, direction = d, n_deg = length(g),
                                       n_deg_entrez = length(unique(na.omit(eg[g]))), universe = length(uni),
                                       go_id = rownames(top),
                                       term = top$Term, N = top$N, DE = top$DE, P = top$P.DE, FDR = top$FDR)
    cat("  GO", t, d, ": n DEG", length(g), " FDR<0.05 terms", sum(ga$FDR < 0.05), "\n")
  }
}
write.csv(do.call(rbind, summ), file.path(od, "de_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, gotop), file.path(od, "go_top.csv"), row.names = FALSE)
write.csv(do.call(rbind, pan), file.path(od, "de_panel.csv"), row.names = FALSE)
