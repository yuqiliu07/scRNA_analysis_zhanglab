# WT vs TNNC1-I4M: pathway analysis across databases, ORA and GSEA, to compare which gives the most useful results.
#
# Databases (mouse gene symbols):
#   GO_BP        org.Mm.eg.db (GO biological process, all descendants)
#   KEGG         KEGG REST via limma::getGeneKEGGLinks / getKEGGPathwayNames (mmu); cached in work_dir/genotype/kegg_mmu.rds
#   Reactome     MSigDB mouse M2 CP:REACTOME (msigdbr)
#   WikiPathways MSigDB mouse M2 CP:WIKIPATHWAYS (msigdbr)
#   Hallmark     MSigDB mouse MH (msigdbr)
# Methods:
#   ORA   hypergeometric test of robust up / down DEGs (13_pseudobulk_de.R: FDR <= 0.05, |log2FC| >= 0.25, not flagged)
#         against the genes tested in that cell type; sets of 10-500 tested genes; cell type x direction with >= 15 DEGs.
#   GSEA  fgsea (multilevel) on all tested genes ranked by sign(logFC) * -log10(P); genes flagged technical / ambient
#         are removed from the ranking so the batch-related ambient signal does not drive the result; sets of 15-500.
# Significance: BH-adjusted p (padj) <= 0.05 within each cell type x database x method (x direction for ORA).
#
# Outputs: <work_dir>/genotype/pw_<method>_<db>_<type>.csv.gz (full tables, local)
#          <results_dir>/genotype/pathway_summary.csv  (number of significant sets)
#          <results_dir>/genotype/pathway_top.csv      (top 10 sets per cell type x database x method x direction)
# Usage: Rscript 23_pathways.R <work_dir> <results_dir>
suppressPackageStartupMessages({
  library(limma); library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db); library(fgsea); library(msigdbr)
})
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
set.seed(0)

# ---------------- gene sets ----------------
sym <- AnnotationDbi::as.list(org.Mm.egSYMBOL)
eg2sym <- function(e) unique(unlist(sym[intersect(e, names(sym))]))
go2eg <- AnnotationDbi::as.list(org.Mm.egGO2ALLEGS)
ont <- AnnotationDbi::Ontology(names(go2eg))
go2eg <- go2eg[!is.na(ont) & ont == "BP"]
GO <- lapply(go2eg, eg2sym)
names(GO) <- paste0(AnnotationDbi::Term(names(GO)), " (", names(GO), ")")

kfile <- file.path(wd, "kegg_mmu.rds")
if (file.exists(kfile)) KEGG <- readRDS(kfile) else {
  kl <- getGeneKEGGLinks(species.KEGG = "mmu")
  kn <- getKEGGPathwayNames(species.KEGG = "mmu", remove.qualifier = TRUE)
  rownames(kn) <- kn$PathwayID
  KEGG <- lapply(split(kl$GeneID, kl$PathwayID), eg2sym)
  names(KEGG) <- kn[names(KEGG), "Description"]
  saveRDS(KEGG, kfile)
}
msig <- function(collection, sub = NULL) {
  d <- tryCatch(msigdbr(db_species = "MM", species = "Mus musculus", collection = collection, subcollection = sub),
                error = function(e) msigdbr(species = "Mus musculus", category = collection, subcategory = sub))
  split(d$gene_symbol, d$gs_name)
}
DB <- list(GO_BP = GO, KEGG = KEGG, Reactome = msig("M2", "CP:REACTOME"), WikiPathways = msig("M2", "CP:WIKIPATHWAYS"),
           Hallmark = msig("MH"))
for (k in names(DB)) cat(k, length(DB[[k]]), "sets\n")

ora <- function(genes, universe, sets) {
  N <- length(universe); n <- length(genes)
  r <- lapply(names(sets), function(s) {
    g <- intersect(sets[[s]], universe); K <- length(g)
    if (K < 10 || K > 500) return(NULL)
    k <- length(intersect(g, genes))
    data.frame(set = s, size = K, overlap = k, expected = n * K / N,
               p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               genes = paste(head(sort(intersect(g, genes)), 15), collapse = ","))
  })
  r <- do.call(rbind, r)
  r$padj <- p.adjust(r$p, "BH")
  r$fold_enrichment <- r$overlap / r$expected
  r[order(r$p), ]
}

types <- sub("^de_(.*)[.]csv[.]gz$", "\\1", list.files(wd, "^de_.*[.]csv[.]gz$"))
summ <- list(); top <- list()
for (t in types) {
  de <- read.csv(gzfile(file.path(wd, paste0("de_", t, ".csv.gz"))))
  uni <- de$gene
  # ---- ORA ----
  for (d in c("Up", "Down")) {
    g <- de$gene[de$robust & (if (d == "Up") de$logFC > 0 else de$logFC < 0)]
    if (length(g) < 15) next
    for (k in names(DB)) {
      r <- ora(g, uni, DB[[k]])
      write.csv(r, gzfile(file.path(wd, sprintf("pw_ORA_%s_%s_%s.csv.gz", k, t, d))), row.names = FALSE)
      summ[[length(summ) + 1]] <- data.frame(cell_type = t, method = "ORA", database = k, direction = d,
                                             n_input = length(g), n_sets_tested = nrow(r), n_sig = sum(r$padj <= 0.05))
      tp <- head(r, 10)
      top[[length(top) + 1]] <- data.frame(cell_type = t, method = "ORA", database = k, direction = d, set = tp$set,
                                           size = tp$size, stat = tp$fold_enrichment, p = tp$p, padj = tp$padj,
                                           genes = tp$genes)
    }
  }
  # ---- GSEA ----
  keep <- !(de$technical | de$ambient)
  st <- sign(de$logFC[keep]) * -log10(pmax(de$PValue[keep], 1e-300))
  names(st) <- de$gene[keep]
  st <- st[!duplicated(names(st))]
  st <- st + rnorm(length(st), 0, 1e-8)  # break ties
  for (k in names(DB)) {
    r <- suppressWarnings(fgseaMultilevel(DB[[k]], st, minSize = 15, maxSize = 500, eps = 0))
    r <- r[order(r$pval), ]
    le <- vapply(r$leadingEdge, function(x) paste(head(x, 15), collapse = ","), "")
    out <- data.frame(set = r$pathway, size = r$size, NES = r$NES, p = r$pval, padj = r$padj, leading_edge = le)
    write.csv(out, gzfile(file.path(wd, sprintf("pw_GSEA_%s_%s.csv.gz", k, t))), row.names = FALSE)
    summ[[length(summ) + 1]] <- data.frame(cell_type = t, method = "GSEA", database = k, direction = "both",
                                           n_input = length(st), n_sets_tested = nrow(out), n_sig = sum(out$padj <= 0.05, na.rm = TRUE))
    for (d in c("Up", "Down")) {
      o <- out[if (d == "Up") out$NES > 0 else out$NES < 0, ]
      tp <- head(o, 10)
      if (nrow(tp)) top[[length(top) + 1]] <- data.frame(cell_type = t, method = "GSEA", database = k, direction = d,
                                                         set = tp$set, size = tp$size, stat = tp$NES, p = tp$p,
                                                         padj = tp$padj, genes = tp$leading_edge)
    }
  }
  cat(t, "done\n")
}
S <- do.call(rbind, summ)
write.csv(S, file.path(od, "pathway_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, top), file.path(od, "pathway_top.csv"), row.names = FALSE)
print(S)
