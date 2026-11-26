#!/usr/bin/env Rscript
# Compare our UMAP / clusters (09_harmony_cluster.R) with Novogene's delivered UMAP / clusters.
#
# Inputs
#   <HARMONY_DIR>/harmony_umap_per_cell.csv.gz   ours: cell "<S>_<barcode>", UMAP before/after, cluster
#   <NOVO_DIR>/Combination1.cluster.UMAP.xls     Novogene: barcode "<barcode>_<k>" (k = 1..8 = A..H), UMAP_1/2
#   <NOVO_DIR>/Combination1.clusters_info.xls    Novogene: barcode, orig.ident, seurat_clusters
#   (extracted from 03.Result_*.zip: 3.Identification_subpopulation/{UMAP,Clustering}/)
#
# Comparisons
#   1. Cell overlap per sample (ours = pass-QC singlets; Novogene = their own QC)
#   2. Cluster agreement on shared cells: contingency table, ARI, best-match fraction
#   3. Sample mixing in each 2-D embedding (same footing for all three, since only
#      Novogene's 2-D UMAP is available): fraction of the K_NN nearest neighbours in
#      UMAP space from the same sample, divided by the sample's share of cells
#      (1 = perfectly mixed)
#
# Outputs (<OUT_DIR>/)
#   novogene_vs_ours_umap.png        Novogene UMAP | ours before | ours after, coloured by sample
#   novogene_vs_ours_by_sample_split.png  per-sample panels, 3 rows (Novogene / ours before / after)
#   novogene_vs_ours_clusters.png    ours (harmony UMAP) coloured by Novogene cluster, and
#                                    Novogene UMAP coloured by our cluster (shared cells)
#   novogene_vs_ours_cluster_heatmap.png  row-normalised contingency (Novogene cluster -> ours)
#   compare_cells_per_sample.csv, compare_cluster_contingency.csv,
#   compare_cluster_best_match.csv, compare_umap_mixing.csv, compare_summary.csv
#
# Usage: sbatch 10_compare_novogene_umap.sh  (login node runs out of memory)
#        or: Rscript 10_compare_novogene_umap.R inside a job
# Env (optional): HARMONY_DIR, NOVO_DIR, OUT_DIR, K_NN (30)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

scratch     <- Sys.getenv("SCRATCH", "/u/scratch/y/yuqi")
harmony_dir <- Sys.getenv("HARMONY_DIR", file.path(scratch, "FSU_seq/work_dir/harmony_mm10"))
novo_dir    <- Sys.getenv("NOVO_DIR",    file.path(scratch, "compare/novo_umap"))
out_dir     <- Sys.getenv("OUT_DIR",     file.path(scratch, "FSU_seq/work_dir/compare_novogene"))
k_nn        <- as.integer(Sys.getenv("K_NN", "30"))
SEED        <- 42
samples     <- LETTERS[1:8]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- read ----
ours <- fread(file.path(harmony_dir, "harmony_umap_per_cell.csv.gz"))
ours[, barcode := sub("^[A-H]_", "", cell)]
ours[, key := paste(sample, barcode, sep = "_")]
ours[, cluster := as.character(cluster)]

nu <- fread(file.path(novo_dir, "Combination1.cluster.UMAP.xls"))
nc <- fread(file.path(novo_dir, "Combination1.clusters_info.xls"))
stopifnot(identical(nu$barcode, nc$barcode))
novo <- data.table(nbc = nu$barcode, sample = nc$orig.ident,
                   UMAP_1 = nu$UMAP_1, UMAP_2 = nu$UMAP_2,
                   cluster = as.character(nc$seurat_clusters))
idx <- as.integer(sub(".*_", "", novo$nbc))
stopifnot(all(samples[idx] == novo$sample))              # suffix _k == sample k
novo[, barcode := sub("_[0-9]+$", "", nbc)]
novo[, key := paste(sample, barcode, sep = "_")]
message(sprintf("ours: %d cells, %d clusters | Novogene: %d cells, %d clusters",
                nrow(ours), uniqueN(ours$cluster), nrow(novo), uniqueN(novo$cluster)))

# ---- 1. overlap ----
shared <- intersect(ours$key, novo$key)
cps <- rbindlist(lapply(samples, function(s) {
  o <- ours[sample == s, key]; n <- novo[sample == s, key]
  data.table(sample = s, n_ours = length(o), n_novogene = length(n),
             n_shared = length(intersect(o, n)),
             ours_only = length(setdiff(o, n)), novogene_only = length(setdiff(n, o)))
}))
cps <- rbind(cps, cps[, c(list(sample = "total"), lapply(.SD, sum)), .SDcols = -"sample"])
cps[, pct_ours_in_novogene := round(100 * n_shared / n_ours, 2)]
cps[, pct_novogene_in_ours := round(100 * n_shared / n_novogene, 2)]
fwrite(cps, file.path(out_dir, "compare_cells_per_sample.csv"))

# ---- 2. cluster agreement (shared cells) ----
m <- merge(ours[, .(key, sample, cl_ours = cluster, H1 = UMAP_harmony_1, H2 = UMAP_harmony_2)],
           novo[, .(key, cl_novo = cluster, N1 = UMAP_1, N2 = UMAP_2)], by = "key")
tab <- table(novo = m$cl_novo, ours = m$cl_ours)
ari <- function(t) {
  n <- sum(t); c2 <- function(x) sum(as.numeric(x) * (as.numeric(x) - 1) / 2)
  a <- c2(t); b <- c2(rowSums(t)); c <- c2(colSums(t)); e <- b * c / c2(n)
  (a - e) / ((b + c) / 2 - e)
}
ARI <- ari(tab)
num_sort <- function(x) as.character(sort(as.integer(unique(x))))
ct <- as.data.table(tab)[N > 0]
setnames(ct, "N", "n_cells")
fwrite(ct[order(as.integer(novo), -n_cells)], file.path(out_dir, "compare_cluster_contingency.csv"))

best_novo <- ct[, .(n_novo = sum(n_cells), best_ours = ours[which.max(n_cells)],
                    frac_in_best = round(max(n_cells) / sum(n_cells), 3)), by = novo]
best_ours <- ct[, .(n_ours = sum(n_cells), best_novo = novo[which.max(n_cells)],
                    frac_in_best = round(max(n_cells) / sum(n_cells), 3)), by = ours]
fwrite(rbind(best_novo[, .(direction = "novogene->ours", cluster = novo, n_cells = n_novo,
                           best_match = best_ours, frac_in_best)][order(as.integer(cluster))],
             best_ours[, .(direction = "ours->novogene", cluster = ours, n_cells = n_ours,
                           best_match = best_novo, frac_in_best)][order(as.integer(cluster))]),
       file.path(out_dir, "compare_cluster_best_match.csv"))

# ---- 3. sample mixing in each 2-D embedding ----
mix2d <- function(x, y, smp) {
  nn <- RANN::nn2(cbind(x, y), k = k_nn + 1)$nn.idx[, -1, drop = FALSE]
  same <- rowMeans(matrix(smp[nn], nrow = nrow(nn)) == smp)
  data.table(sample = smp, same = same, ratio = same / (as.numeric(table(smp)[smp]) / length(smp)))
}
mix <- rbind(
  cbind(embedding = "Novogene UMAP (no integration)", mix2d(novo$UMAP_1, novo$UMAP_2, novo$sample)),
  cbind(embedding = "Ours before integration (PCA)", mix2d(ours$UMAP_pca_1, ours$UMAP_pca_2, ours$sample)),
  cbind(embedding = "Ours after Harmony", mix2d(ours$UMAP_harmony_1, ours$UMAP_harmony_2, ours$sample)))
mix_s <- mix[, .(n_cells = .N, same_frac_mean = round(mean(same), 4),
                 ratio_mean = round(mean(ratio), 3), ratio_median = round(median(ratio), 3)),
             by = .(embedding, sample)]
mix_all <- mix[, .(sample = "all", n_cells = .N, same_frac_mean = round(mean(same), 4),
                   ratio_mean = round(mean(ratio), 3), ratio_median = round(median(ratio), 3)),
               by = embedding]
mix_s <- rbind(mix_s, mix_all)
fwrite(mix_s, file.path(out_dir, "compare_umap_mixing.csv"))

summ <- data.table(
  n_cells_ours = nrow(ours), n_cells_novogene = nrow(novo), n_shared = length(shared),
  n_clusters_ours = uniqueN(ours$cluster), n_clusters_novogene = uniqueN(novo$cluster),
  ARI_shared_cells = round(ARI, 4), k_nn_mixing = k_nn,
  mixing_ratio_novogene = mix_all[embedding %like% "Novogene", ratio_mean],
  mixing_ratio_ours_before = mix_all[embedding %like% "before", ratio_mean],
  mixing_ratio_ours_harmony = mix_all[embedding %like% "Harmony", ratio_mean])
fwrite(summ, file.path(out_dir, "compare_summary.csv"))

# ---- plots ----
lab <- c(A = "A (WT, M)", B = "B (WT, M)", C = "C (WT, F)", D = "D (WT, F)",
         E = "E (I4M, M)", F = "F (I4M, M)", G = "G (I4M, F)", H = "H (I4M, F)")
pal <- setNames(c("#257de0", "#f47230", "#00b894", "#f2a600",
                  "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"), samples)
theme_set(theme_classic(base_size = 11) +
          theme(axis.line = element_line(colour = "grey60", linewidth = 0.3),
                axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
                axis.text = element_text(colour = "grey30"),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold", hjust = 0)))
emb_lv <- c("Novogene UMAP (no integration)", "Ours before integration (PCA)", "Ours after Harmony")
long <- rbind(novo[, .(sample, x = UMAP_1, y = UMAP_2, emb = emb_lv[1])],
              ours[, .(sample, x = UMAP_pca_1, y = UMAP_pca_2, emb = emb_lv[2])],
              ours[, .(sample, x = UMAP_harmony_1, y = UMAP_harmony_2, emb = emb_lv[3])])
long[, emb := factor(emb, levels = emb_lv)]
mr <- setNames(mix_all$ratio_mean, mix_all$embedding)
facet_lab <- setNames(sprintf("%s\n%s cells; kNN same-sample ratio %.2f (1 = mixed)", emb_lv,
                              format(long[, .N, by = emb][match(emb_lv, emb), N], big.mark = ","),
                              mr[emb_lv]), emb_lv)
set.seed(SEED)
pd <- long[sample(.N)]
p1 <- ggplot(pd, aes(x, y, colour = sample)) +
  geom_point(size = 0.04, alpha = 0.35, stroke = 0) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~emb, scales = "free", labeller = as_labeller(facet_lab)) +
  labs(x = "UMAP_1", y = "UMAP_2",
       title = "UMAP coloured by sample: Novogene delivery vs ours (random draw order)")
ggsave(file.path(out_dir, "novogene_vs_ours_umap.png"), p1, width = 21, height = 7.5, dpi = 170)

long[, panel := factor(lab[sample], levels = lab)]
bg <- long[, .(x, y, emb)]
p2 <- ggplot(long, aes(x, y)) +
  geom_point(data = bg, colour = "grey85", size = 0.02, stroke = 0) +
  geom_point(aes(colour = sample), size = 0.05, alpha = 0.5, stroke = 0, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  facet_grid(emb ~ panel, scales = "free",
             labeller = labeller(emb = as_labeller(setNames(sub(" \\(", "\n(", emb_lv), emb_lv)))) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Each sample (coloured) on all other cells (grey)")
ggsave(file.path(out_dir, "novogene_vs_ours_by_sample_split.png"), p2, width = 24, height = 10, dpi = 140)

cl_cols <- function(lv) setNames(grDevices::hcl(h = seq(15, 375, length.out = length(lv) + 1)[seq_along(lv)],
                                                c = 70, l = 62), lv)
lv_n <- num_sort(m$cl_novo); lv_o <- num_sort(m$cl_ours)
set.seed(SEED)
md <- m[sample(.N)]
cl_plot <- rbind(md[, .(x = H1, y = H2, cl = cl_novo, lab = cl_novo, panel = "a")],
                 md[, .(x = N1, y = N2, cl = paste0("o", cl_ours), lab = cl_ours, panel = "b")])
cl_plot[, panel := factor(panel, levels = c("a", "b"),
                          labels = c(sprintf("Our Harmony UMAP, coloured by Novogene cluster (%d)", length(lv_n)),
                                     sprintf("Novogene UMAP, coloured by our Leiden cluster (%d)", length(lv_o))))]
cen <- cl_plot[, .(x = median(x), y = median(y)), by = .(panel, cl, lab)]
cols <- c(cl_cols(lv_n), setNames(cl_cols(lv_o), paste0("o", lv_o)))
p3 <- ggplot(cl_plot, aes(x, y, colour = cl)) +
  geom_point(size = 0.04, alpha = 0.4, stroke = 0, show.legend = FALSE) +
  geom_label(data = cen, aes(x, y, label = lab), inherit.aes = FALSE, size = 2.8,
             linewidth = 0, fill = alpha("white", 0.75), label.padding = unit(0.1, "lines")) +
  scale_colour_manual(values = cols) +
  facet_wrap(~panel, scales = "free") +
  labs(x = "UMAP_1", y = "UMAP_2",
       title = sprintf("Cluster cross-check on %s shared cells (ARI = %.3f)",
                       format(nrow(m), big.mark = ","), ARI))
ggsave(file.path(out_dir, "novogene_vs_ours_clusters.png"), p3, width = 18, height = 8.5, dpi = 170)

hm <- ct[, frac := n_cells / sum(n_cells), by = novo]
# order our clusters by their best-matching Novogene cluster so the diagonal is visible
ord_o <- best_ours[order(as.integer(best_novo), -frac_in_best), ours]
p4 <- ggplot(hm, aes(factor(ours, levels = ord_o), factor(novo, levels = rev(lv_n)), fill = frac)) +
  geom_tile() +
  scale_fill_gradient(low = "#eef5fd", high = "#1757ac", limits = c(0, 1),
                      name = "fraction of\nNovogene\ncluster") +
  labs(x = "our Leiden cluster (ordered by best-matching Novogene cluster)", y = "Novogene cluster",
       title = sprintf("Where each Novogene cluster's cells fall in our clustering (shared cells, ARI = %.3f)", ARI)) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7))
ggsave(file.path(out_dir, "novogene_vs_ours_cluster_heatmap.png"), p4, width = 13, height = 7.5, dpi = 160)

cat("\n===== cells =====\n"); print(cps)
cat("\n===== summary =====\n"); print(t(summ))
cat("\n===== mixing (2-D UMAP, kNN) =====\n"); print(dcast(mix_s, sample ~ embedding, value.var = "ratio_mean"))
cat("\n===== Novogene cluster -> best our cluster =====\n"); print(best_novo[order(as.integer(novo))])
message("done")
