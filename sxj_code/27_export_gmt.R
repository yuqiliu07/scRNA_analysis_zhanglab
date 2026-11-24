# Export the gene sets used in 23_pathways.R as GMT files (work_dir/genotype/gmt_<db>.gmt), for plotting in Python.
# Usage: Rscript 27_export_gmt.R <work_dir>
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Mm.eg.db); library(GO.db); library(msigdbr) })
wd <- file.path(commandArgs(TRUE)[1], "genotype")
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
for (k in names(DB)) {
  s <- DB[[k]]
  s <- s[lengths(s) >= 5 & lengths(s) <= 2000]
  writeLines(vapply(names(s), function(n) paste(c(gsub("\t", " ", n), "na", unique(s[[n]])), collapse = "\t"), ""),
             file.path(wd, paste0("gmt_", k, ".gmt")))
  cat(k, length(s), "\n")
}
