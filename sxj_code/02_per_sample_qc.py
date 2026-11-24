"""Step B (workflow docs/mouse_sc_qc_workflow_zh.md sec. 1.1, 1.3, 1.4): per-sample QC metrics, MAD filters,
doublet detection (Scrublet, per 10x lane = per sample) and a quick clustering used by SoupX.

Candidate barcodes = EmptyDrops calls (FDR < 0.001, 01_emptydrops.R). Every candidate is kept in the output
with flags, so the per-step cell counts can be audited; nothing is dropped silently.

Filters (per sample, never pooled):
  n_genes >= 200 (discovery_qc_strategy_zh.md sec. 5);
  log1p(total UMI) and log1p(n_genes) within median +/- 5 MAD;
  pct_mt <= median + 3 MAD and <= MT_CAP (single-nucleus data: hard cap 5%);
  Scrublet doublet removed, capped at the expected multiplet rate (uncapped call and score kept for step D).
Stress genes (Fos/Jun/Hspa1a ...) are scored only, never used to remove cells.

Usage: python 02_per_sample_qc.py <data_dir> <work_dir> <sample> [<sample> ...]
Output: <work_dir>/per_sample/<S>.h5ad (raw integer counts in X, all candidates, flags in obs)
        <work_dir>/per_sample/<S>_clusters.csv (QC-pass cells -> quick cluster, for SoupX)
        <work_dir>/per_sample/<S>_thresholds.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scanpy.preprocessing._scrublet import core as _scrublet_core
from sklearn.utils import check_random_state


def _sample_comb_int64(dims, nsamp, *, random_state=None, method="auto"):
    # scanpy's sample_comb (sklearn sample_without_replacement on prod(dims)) uses a C long, which is 32-bit
    # on Windows and overflows once a sample has > 46,340 nuclei (n_obs^2 > 2^31). Same distribution
    # (uniform grid cells without repetition), drawn in int64 numpy.
    rs = check_random_state(random_state)
    n = int(np.prod(dims, dtype=np.int64))
    seen = np.empty(0, dtype=np.int64)
    while seen.size < nsamp:
        draw = rs.randint(0, n, size=2 * (nsamp - seen.size), dtype=np.int64)
        _, first = np.unique(np.concatenate([seen, draw]), return_index=True)
        seen = np.concatenate([seen, draw])[np.sort(first)]
    return np.vstack(np.unravel_index(seen[:nsamp], dims)).T


_scrublet_core.sample_comb = _sample_comb_int64

MIN_GENES = 200
NMAD_COUNTS = 5
NMAD_MT = 3
MT_CAP = 5.0          # % ; single-nucleus data (median pct_mt 0.5-1.3% in all samples)
STRESS = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2", "Ier3", "Dusp1",
          "Hspa1a", "Hspa1b", "Hspa8", "Hsph1", "Dnajb1", "Nr4a1", "Klf2", "Klf4", "Zfp36", "Socs3"]
SEED = 0


def mad_bounds(x, nmad):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return med - nmad * mad, med + nmad * mad, med, mad


def expected_doublet_rate(n_cells):
    # 10x GEM-X (3' v4) multiplet rate ~0.4% per 1,000 cells recovered; capped for heavily over-called lanes
    return float(min(0.004 * n_cells / 1000, 0.15))


def run(data_dir, work_dir, s):
    out = Path(work_dir) / "per_sample"
    out.mkdir(parents=True, exist_ok=True)
    ed = pd.read_csv(Path(work_dir) / "emptydrops" / f"{s}_emptydrops.csv")
    ed = ed[ed["called"]].set_index("barcode")

    ad = sc.read_10x_h5(Path(data_dir) / s / "raw_feature_bc_matrix.h5")
    ad.var["symbol"] = ad.var_names
    ad.var_names_make_unique()
    ad = ad[ad.obs_names.isin(ed.index)].copy()
    ad.X = sp.csr_matrix(ad.X, dtype=np.float32)
    ad.obs["sample"] = s
    ad.obs["ed_FDR"] = ed.loc[ad.obs_names, "FDR"].values
    ad.obs["in_cellranger"] = ed.loc[ad.obs_names, "in_cellranger"].values

    sym = ad.var["symbol"]
    ad.var["mt"] = sym.str.startswith("mt-")
    ad.var["ribo"] = sym.str.match(r"^Rp[sl]\d")
    ad.var["hb"] = sym.str.match(r"^Hb[ab]-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
    ad.obs.rename(columns={"total_counts": "total", "n_genes_by_counts": "n_genes",
                           "pct_counts_mt": "pct_mt", "pct_counts_ribo": "pct_ribo",
                           "pct_counts_hb": "pct_hb"}, inplace=True)
    ad.obs = ad.obs[["sample", "ed_FDR", "in_cellranger", "total", "n_genes", "pct_mt", "pct_ribo", "pct_hb"]]
    ad.obs["pct_malat1"] = (np.asarray(ad[:, sym == "Malat1"].X.sum(1)).ravel() / ad.obs["total"] * 100)

    # thresholds
    lt, lg = np.log1p(ad.obs["total"]), np.log1p(ad.obs["n_genes"])
    t_lo, t_hi, _, _ = mad_bounds(lt, NMAD_COUNTS)
    g_lo, g_hi, _, _ = mad_bounds(lg, NMAD_COUNTS)
    _, m_hi, m_med, m_mad = mad_bounds(ad.obs["pct_mt"].values, NMAD_MT)
    mt_cut = min(m_hi, MT_CAP)
    o = ad.obs
    o["pass_min_genes"] = o["n_genes"] >= MIN_GENES
    o["pass_umi"] = (lt >= t_lo) & (lt <= t_hi)
    o["pass_genes"] = (lg >= g_lo) & (lg <= g_hi)
    o["pass_mt"] = o["pct_mt"] <= mt_cut
    o["pass_basic"] = o["pass_min_genes"] & o["pass_umi"] & o["pass_genes"] & o["pass_mt"]

    # raw counts kept untouched in X; work on a copy for scores / clustering
    w = ad[o["pass_basic"].values].copy()
    w.layers["counts"] = w.X.copy()

    rate = expected_doublet_rate(w.n_obs)
    sc.pp.scrublet(w, expected_doublet_rate=rate, random_state=SEED, n_prin_comps=30)

    sc.pp.normalize_total(w, target_sum=1e4)
    sc.pp.log1p(w)
    stress = [g for g in STRESS if g in w.var_names]
    sc.tl.score_genes(w, stress, score_name="stress_score", random_state=SEED)
    sc.pp.highly_variable_genes(w, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.pca(w, n_comps=30, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(w, n_neighbors=20, random_state=SEED)
    sc.tl.leiden(w, resolution=0.8, key_added="quick_cluster", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)

    # Scrublet's automatic threshold (threshold_minimum on a bimodal fit) flagged 33-44% of nuclei in five
    # samples here, concentrated in whole high-UMI cardiomyocyte clusters (large polyploid nuclei on top of a
    # heavy ambient CM signal). Calls are therefore capped at the expected multiplet rate: a nucleus is a
    # doublet only if its score exceeds both the automatic threshold and the (1 - expected rate) quantile.
    # The uncapped call is kept as predicted_doublet_auto for the per-cluster review (step D).
    auto_thr = float(w.uns["scrublet"]["threshold"])
    cap_thr = float(np.quantile(w.obs["doublet_score"], 1 - rate))
    w.obs["predicted_doublet_auto"] = w.obs["doublet_score"] > auto_thr
    w.obs["predicted_doublet"] = w.obs["doublet_score"] > max(auto_thr, cap_thr)
    for c in ["doublet_score", "predicted_doublet", "predicted_doublet_auto", "stress_score", "quick_cluster"]:
        o[c] = w.obs[c].reindex(o.index)
    o["predicted_doublet"] = o["predicted_doublet"].fillna(False).astype(bool)
    o["predicted_doublet_auto"] = o["predicted_doublet_auto"].fillna(False).astype(bool)
    o["quick_cluster"] = o["quick_cluster"].astype(str).replace("nan", "")
    o["qc_pass"] = o["pass_basic"] & ~o["predicted_doublet"]
    ad.uns["scrublet_threshold"] = max(auto_thr, cap_thr)
    ad.uns["scrublet_threshold_auto"] = auto_thr

    thr = pd.DataFrame([dict(sample=s, n_candidates=ad.n_obs,
                             umi_lo=np.expm1(t_lo), umi_hi=np.expm1(t_hi),
                             genes_lo=max(np.expm1(g_lo), MIN_GENES), genes_hi=np.expm1(g_hi),
                             mt_median=m_med, mt_mad=m_mad, mt_mad_cut=m_hi, mt_cut=mt_cut,
                             scrublet_expected_rate=rate,
                             scrublet_threshold=ad.uns["scrublet_threshold"],
                             scrublet_threshold_auto=auto_thr,
                             n_pass_min_genes=int(o["pass_min_genes"].sum()),
                             n_pass_basic=int(o["pass_basic"].sum()),
                             n_doublet_auto=int(o["predicted_doublet_auto"].sum()),
                             n_doublet=int(o["predicted_doublet"].sum()),
                             n_qc_pass=int(o["qc_pass"].sum()),
                             n_qc_pass_in_cellranger=int((o["qc_pass"] & o["in_cellranger"]).sum()))])
    thr.to_csv(out / f"{s}_thresholds.csv", index=False)
    o.loc[o["pass_basic"], ["quick_cluster"]].to_csv(out / f"{s}_clusters.csv", index_label="barcode")
    ad.write_h5ad(out / f"{s}.h5ad", compression="gzip")
    print(thr.T.to_string(header=False), flush=True)


if __name__ == "__main__":
    data_dir, work_dir, *samples = sys.argv[1:]
    for s in samples:
        print(f"=== {s} ===", flush=True)
        run(data_dir, work_dir, s)
