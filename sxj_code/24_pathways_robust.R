# Robustness layer for 23_pathways.R (genotype is confounded with a processing batch):
#  1) shared-pathway flag: a set significant (GSEA padj <= 0.05) with the same NES sign in >= 3 cell types is
#     "shared" (batch-like, e.g. mRNA splicing / RNA metabolism down in every cell type); the rest is "specific".
#  2) internal-control GSEA: for cardiomyocyte and endocardial, genes are ranked by the target-minus-reference
#     statistic of 19_internal_control.R (sign(logFC_int) * -log10(P_int), vs fibroblast and vs endothelial),
#     so changes shared with the reference cell type cancel. A set is "IC-supported" if padj <= 0.05 with the same
#     sign as in the target's own GSEA against at least one reference.
# Outputs: <results_dir>/genotype/pathway_robust_summary.csv, pathway_specific_top.csv, pathway_ic_top.csv
# Usage: Rscript 24_pathways_robust.R <work_dir> <results_dir>
suppressPackageStartupMessages({
  library(limma); library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db); library(fgsea); library(msigdbr)
})
args <- commandArgs(TRUE)
wd <- file.path(args[1], "genotype"); od <- file.path(args[2], "genotype")
set.seed(0)
dbs <- c("GO_BP", "KEGG", "Reactome", "WikiPathways", "Hallmark")
types <- sub("^de_(.*)[.]csv[.]gz$", "\\1", list.files(wd, "^de_.*[.]csv[.]gz$"))

# ---- 1) shared vs specific ----
G <- do.call(rbind, lapply(dbs, function(k) do.call(rbind, lapply(types, function(t) {
  f <- file.path(wd, sprintf("pw_GSEA_%s_%s.csv.gz", k, t))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(gzfile(f)); d$database <- k; d$cell_type <- t; d
}))))
sig <- G[!is.na(G$padj) & G$padj <= 0.05, ]
key <- paste(sig$database, sig$set, sign(sig$NES))
nct <- table(key)
G$n_celltypes_same_dir <- as.integer(nct[paste(G$database, G$set, sign(G$NES))])
G$n_celltypes_same_dir[is.na(G$n_celltypes_same_dir)] <- 0
G$sig <- !is.na(G$padj) & G$padj <= 0.05
G$class <- ifelse(!G$sig, "ns", ifelse(G$n_celltypes_same_dir >= 3, "shared", "specific"))

# ---- 2) internal-control GSEA ----
sym <- AnnotationDbi::as.list(org.Mm.egSYMBOL)
eg2sym <- function(e) unique(unlist(sym[intersect(e, names(sym))]))
go2eg <- AnnotationDbi::as.list(org.Mm.egGO2ALLEGS)
ont <- AnnotationDbi::Ontology(names(go2eg))
go2eg <- go2eg[!is.na(ont) & ont == "BP"]
GO <- lapply(go2eg, eg2sym)
names(GO) <- paste0(AnnotationDbi::Term(names(GO)), " (", names(GO), ")")
msig <- function(collection, sub = NULL) {
  d <- tryCatch(msigdbr(db_species = "MM", species = "Mus musculus", collection = collection, subcollection = sub),
                error = function(e) msigdbr(species = "Mus musculus", category = collection, subcategory = sub))
  split(d$gene_symbol, d$gs_name)
}
DB <- list(GO_BP = GO, KEGG = readRDS(file.path(wd, "kegg_mmu.rds")), Reactome = msig("M2", "CP:REACTOME"),
           WikiPathways = msig("M2", "CP:WIKIPATHWAYS"), Hallmark = msig("MH"))
IC <- list()
for (tg in c("Cardiomyocyte", "Endocardial")) {
  ic <- read.csv(gzfile(file.path(wd, paste0("ic_", tg, ".csv.gz"))))
  for (rf in c("Fibroblast", "Endothelial")) {
    p <- ic[[paste0("P_int_vs_", rf)]]; l <- ic[[paste0("logFC_int_vs_", rf)]]
    ok <- !is.na(p) & !(ic$technical | ic$ambient)
    st <- sign(l[ok]) * -log10(pmax(p[ok], 1e-300)); names(st) <- ic$gene[ok]
    st <- st[!duplicated(names(st))] + rnorm(sum(!duplicated(ic$gene[ok])), 0, 1e-8)
    for (k in dbs) {
      r <- suppressWarnings(fgseaMultilevel(DB[[k]], st, minSize = 15, maxSize = 500, eps = 0))
      IC[[length(IC) + 1]] <- data.frame(cell_type = tg, reference = rf, database = k, set = r$pathway,
                                         NES_ic = r$NES, padj_ic = r$padj)
    }
    cat(tg, "vs", rf, "done\n")
  }
}
IC <- do.call(rbind, IC)
own <- G[G$cell_type %in% c("Cardiomyocyte", "Endocardial"), c("cell_type", "database", "set", "NES", "padj", "class")]
m <- merge(IC, own, by = c("cell_type", "database", "set"))
m$ic_sig <- !is.na(m$padj_ic) & m$padj_ic <= 0.05 & sign(m$NES_ic) == sign(m$NES) & m$padj <= 0.05
sup <- aggregate(ic_sig ~ cell_type + database + set, data = m, FUN = any)
names(sup)[4] <- "ic_supported"
G <- merge(G, sup, by = c("cell_type", "database", "set"), all.x = TRUE)
G$ic_supported[is.na(G$ic_supported)] <- FALSE

S <- aggregate(cbind(sig = G$sig, shared = G$class == "shared", specific = G$class == "specific",
                     ic_supported = G$ic_supported & G$sig) ~ cell_type + database, data = G, FUN = sum)
write.csv(S, file.path(od, "pathway_robust_summary.csv"), row.names = FALSE)
sp <- G[G$class == "specific", ]
sp <- sp[order(sp$cell_type, sp$database, sp$p), ]
sp <- do.call(rbind, lapply(split(sp, list(sp$cell_type, sp$database), drop = TRUE), head, 8))
write.csv(sp[, c("cell_type", "database", "set", "size", "NES", "padj", "ic_supported", "leading_edge")],
          file.path(od, "pathway_specific_top.csv"), row.names = FALSE)
ics <- G[G$ic_supported & G$sig, ]
write.csv(ics[order(ics$cell_type, ics$database, ics$p), c("cell_type", "database", "set", "size", "NES", "padj", "class",
                                                            "leading_edge")],
          file.path(od, "pathway_ic_top.csv"), row.names = FALSE)
print(S)
