# Pathways by sex: GSEA of the WT vs I4M change within males and within females (17_sex_genotype_de.R contrasts).
#
# Same method as the main GSEA in 23_pathways.R: fgsea (multilevel) on all tested genes ranked by
# sign(logFC) * -log10(P), technical / ambient genes (13_pseudobulk_de.R flags) removed, sets of 15-500 genes, BH within
# each cell type, sex and database. Gene sets: the 5 GMT files written by 27_export_gmt.R.
# With 2 mice per sex and genotype the ranking is noisy; the question is mainly whether male and female pathway
# changes point the same way (NES agreement), not new findings.
#
# Outputs: <work_dir>/genotype/sexgsea_<db>_<type>.csv.gz   NES / padj for males and females (local)
#          <results_dir>/genotype/sex_gsea_summary.csv        per cell type and database: significant sets per sex,
#                                                            NES correlation, % same direction among sets significant
#                                                            in either sex
# Usage: Rscript 43_sex_pathways.R <work_dir> <results_dir>
suppressPackageStartupMessages({ library(fgsea) })
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
set.seed(0)
types <- c("Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial")
dbs <- c("GO_BP", "KEGG", "Reactome", "WikiPathways", "Hallmark")
read_gmt <- function(f) {
  l <- strsplit(readLines(f), "\t")
  setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1))
}
DB <- setNames(lapply(dbs, function(k) read_gmt(file.path(wd, sprintf("gmt_%s.gmt", k)))), dbs)
summ <- list()
for (t in types) {
  sx <- read.csv(file.path(wd, sprintf("sex_%s.csv.gz", t)))
  de <- read.csv(file.path(wd, sprintf("de_%s.csv.gz", t)))
  bad <- de$gene[de$technical | de$ambient]
  sx <- sx[!(sx$gene %in% bad), ]
  rk <- list()
  for (s in c("male", "female")) {
    st <- sign(sx[[paste0("logFC_", s)]]) * -log10(pmax(sx[[paste0("P_", s)]], 1e-300))
    names(st) <- sx$gene
    st <- st[!is.na(st) & !duplicated(names(st))]
    rk[[s]] <- st + rnorm(length(st), 0, 1e-8)
  }
  for (k in dbs) {
    res <- lapply(c("male", "female"), function(s) {
      r <- suppressWarnings(fgseaMultilevel(DB[[k]], rk[[s]], minSize = 15, maxSize = 500, eps = 0))
      data.frame(set = r$pathway, NES = r$NES, padj = r$padj)
    })
    m <- merge(res[[1]], res[[2]], by = "set", suffixes = c("_male", "_female"))
    write.csv(m, gzfile(file.path(wd, sprintf("sexgsea_%s_%s.csv.gz", k, t))), row.names = FALSE)
    sig <- (m$padj_male <= 0.05 | m$padj_female <= 0.05) & !is.na(m$padj_male) & !is.na(m$padj_female)
    summ[[length(summ) + 1]] <- data.frame(cell_type = t, database = k, sets_tested = nrow(m),
      sig_male = sum(m$padj_male <= 0.05, na.rm = TRUE), sig_female = sum(m$padj_female <= 0.05, na.rm = TRUE),
      sig_both_same_direction = sum(m$padj_male <= 0.05 & m$padj_female <= 0.05 & sign(m$NES_male) == sign(m$NES_female), na.rm = TRUE),
      r_NES = round(cor(m$NES_male, m$NES_female, use = "complete.obs"), 3),
      pct_same_direction_sig_either = round(100 * mean(sign(m$NES_male[sig]) == sign(m$NES_female[sig])), 1))
    cat(t, k, "done\n")
  }
}
S <- do.call(rbind, summ)
write.csv(S, file.path(od, "sex_gsea_summary.csv"), row.names = FALSE)
print(S)
