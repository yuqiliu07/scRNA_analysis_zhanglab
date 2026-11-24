"""Name the subclusters from 09_refine.py and build the final nucleus set.

Each subcluster gets a state label and a flag:
  ok                 - kept
  low_quality        - kept but marked: markers are ambient / mitochondrial genes (Gm42418, Cmss1, Lars2, mt-*),
                       low UMI; these look like provider cluster 2
  suspected_doublet  - removed: markers of a second lineage plus a high Scrublet rate

The names were assigned by reviewing refine_subclusters.csv / refine_submarkers.csv (seeded run, so the
subcluster numbers are reproducible). To guard against a rerun renumbering them, every subcluster is checked
against marker genes that must be among its top 30 markers; the script stops if a check fails.

Usage: python 10_name_subclusters.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import importlib.util


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s05 = _load("step05", "05_qc_report.py")

# subcluster -> (state, flag, genes that must be in its top-30 Wilcoxon markers,
#                canonical genes of which at least one must be detected in >= MIN_DETECT of its nuclei)
# The canonical-detection check guards against Wilcoxon ranking a gene high only because it is absent elsewhere
# (an earlier run named a low-UMI macrophage subcluster "lymphocyte" from Pax5 / Skap1 detected in 1-2% of nuclei).
MIN_DETECT = 0.2
NAMES = {
    "CM_0": ("vCM_high_UMI", "ok", ["Kcnd2", "Myom1"], ["Ryr2", "Ttn"]),
    "CM_1": ("vCM_ambient", "low_quality", ["Gm42418", "Cmss1"], ["Gm42418"]),
    "CM_2": ("vCM_Rbm20", "ok", ["Rbm20", "Fhl2"], ["Ryr2", "Ttn"]),
    "Fib_0": ("Fib_high_UMI", "ok", ["Gm42418", "Cmss1"], ["Dcn", "Abca8a", "Pdgfra"]),
    "Fib_1": ("Fib_low_UMI", "ok", ["Commd7"], ["Dcn", "Abca8a", "Pdgfra"]),
    "EC_0": ("EC_capillary_high_UMI", "ok", ["Ptprg", "Plcb4"], ["Flt1", "Ptprb", "Cdh5"]),
    "EC_1": ("Endocardial_1", "ok", ["Tmem108", "Npr3"], ["Npr3", "Tmem108"]),
    "EC_2": ("EC_arterial", "ok", ["Arl15", "Eln"], ["Eln", "Mecom"]),
    "EC_3": ("Endocardial+Fib", "suspected_doublet", ["Abca8a", "Dlc1"], ["Abca8a"]),
    "EC_4": ("Endocardial_2", "ok", ["Tmem108", "Npr3"], ["Npr3", "Tmem108"]),
    "EC_5": ("Lymphatic_EC", "ok", ["Reln", "Mmrn1"], ["Mmrn1", "Reln", "Flt4"]),
    "EC_6": ("EC_capillary_low_UMI", "ok", ["Ptprg", "Flt1"], ["Flt1", "Ptprb", "Cdh5"]),
    "Mural_0": ("Pericyte_1", "ok", ["Mdga2", "Nrxn1"], ["Abcc9", "Pdgfrb", "Rgs5"]),
    "Mural_1": ("Pericyte_ambient", "low_quality", ["Gm42418", "Cmss1"], ["Gm42418"]),
    "Mural_2": ("Pericyte_2", "ok", ["Dlc1", "Ebf1"], ["Abcc9", "Pdgfrb", "Rgs5"]),
    "Mural_3": ("SMC", "ok", ["Myh11", "Kcnab1"], ["Myh11"]),
    "Mural_4": ("Pericyte+CM", "suspected_doublet", ["Ttn", "Myh6"], ["Ttn"]),
    "Immune_0": ("Macrophage_1", "ok", ["F13a1", "Mir99ahg"], ["F13a1", "Mrc1"]),
    "Immune_1": ("Macrophage+EC", "suspected_doublet", ["Ptprg", "Ptprm"], ["Ptprg"]),
    "Immune_2": ("Macrophage_ambient", "low_quality", ["Gm42418", "F13a1"], ["Gm42418"]),
    "Immune_3": ("Macrophage_ambient_2", "low_quality", ["Gm42418", "Cmss1"], ["Gm42418"]),
    # Ciita / Ccr2 / Cd74 rank high but are detected in 4-10% of nuclei; median 111 UMI
    "Immune_4": ("Macrophage_low_UMI", "low_quality", ["Malat1", "Ciita"], ["Mrc1"]),
    "Immune_5": ("Macrophage_Mrc1", "ok", ["F13a1", "Mrc1"], ["F13a1", "Mrc1"]),
    "Immune_6": ("Lymphocyte", "ok", ["Skap1", "Bach2"], ["Skap1", "Pax5", "Ptprc"]),
    "Immune_7": ("DC_cDC1", "ok", ["Wdfy4", "Irf8"], ["Clec9a", "Flt3", "Ciita"]),
}
# coarse type for the final composition (overrides cell_type where the subcluster resolves a new population)
NEW_TYPE = {"Immune_6": "Lymphocyte", "Immune_7": "Dendritic"}
# per-nucleus labels from 07 / 09 that are kept regardless of the subcluster flag
KEEP_TYPES = {"Epicardial", "Lymphocyte"}


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    tdir, fdir = res_dir / "tables", res_dir / "figures"
    mk = pd.read_csv(tdir / "refine_submarkers.csv")
    top = mk.groupby("group")["names"].agg(lambda x: set(x.iloc[:30]))
    bad = [(k, g) for k, (_, _, gs, _) in NAMES.items() for g in gs if g not in top.get(k, set())]
    missing = set(top.index) - set(NAMES)
    if bad or missing:
        raise SystemExit(f"subcluster check failed: markers {bad}; unnamed {sorted(missing)} - review names")

    a = ad.read_h5ad(work_dir / "merged" / "refined.h5ad")
    sub = a.obs["subcluster"].astype(str)
    det_rows = []
    for k, (_, _, _, canon) in NAMES.items():
        g = [x for x in canon if x in a.var_names]
        mm = np.where((sub == k).values)[0]
        d = np.asarray((a[mm][:, g].layers["counts"] > 0).mean(0)).ravel() if len(mm) else np.zeros(len(g))
        det_rows.append({"subcluster": k, **{f"det_{x}": v for x, v in zip(g, d)}, "max_detect": d.max()})
    det = pd.DataFrame(det_rows).set_index("subcluster")
    det.to_csv(tdir / "final_name_detection_check.csv")
    low = det.index[det["max_detect"] < MIN_DETECT].tolist()
    if low:
        raise SystemExit(f"canonical-marker detection < {MIN_DETECT} in {low} - review names")
    a.obs["state"] = sub.map(lambda k: NAMES[k][0] if k in NAMES else k)
    a.obs["flag"] = sub.map(lambda k: NAMES[k][1] if k in NAMES else "ok")
    # epicardial / mesothelial nuclei (07: Muc16 cluster) express fibroblast genes and sit inside EC_3;
    # keep them under their own label instead of the subcluster's doublet flag
    # lymphocytes (09, per nucleus) scattered in other subclusters likewise keep their own label
    for t in KEEP_TYPES:
        mm = (a.obs["cell_type"].astype(str) == t).values
        a.obs.loc[mm, "state"] = t
        a.obs.loc[mm, "flag"] = "ok"
    ct = a.obs["cell_type"].astype(str).copy()
    for k, t in NEW_TYPE.items():
        ct[(sub == k).values] = t
    a.obs["cell_type_final"] = ct

    summ = a.obs.astype({"subcluster": str, "state": str}).groupby(["subcluster", "state"], observed=True).agg(flag=("flag", "first"),
                                          cell_type_final=("cell_type_final", "first"), n=("state", "size"))
    summ.to_csv(tdir / "final_subcluster_names.csv")
    print(summ.to_string(), flush=True)
    print(a.obs["flag"].value_counts().to_string(), flush=True)

    f = a[(a.obs["flag"] != "suspected_doublet").values].copy()
    order = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
             "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"]
    order += sorted(set(f.obs["cell_type_final"]) - set(order))
    cnt = pd.crosstab(f.obs["sample"], f.obs["cell_type_final"]).reindex(columns=order, fill_value=0)
    cnt.loc["total"] = cnt.sum()
    cnt.to_csv(tdir / "final_celltype_counts.csv")
    frac = cnt.div(cnt.sum(1), axis=0)
    frac.to_csv(tdir / "final_celltype_fractions.csv")
    st = pd.crosstab(f.obs["state"], f.obs["sample"])
    st["total"] = st.sum(1)
    st.to_csv(tdir / "final_state_counts.csv")
    print(cnt.to_string(), "\n", frac.round(3).to_string(), flush=True)
    print("final nuclei", f.n_obs, flush=True)

    fig, axes = plt.subplots(1, 2, figsize=(17, 7.5))
    s05.scatter(axes[0], f.obsm["X_umap"], f.obs["cell_type_final"].values, "Final cell type", cats=order)
    um = f.obsm["X_umap"]
    for t in ["Lymphocyte", "Dendritic", "Neural_Schwann", "Epicardial"]:
        mm = (f.obs["cell_type_final"] == t).values
        if mm.sum():
            axes[0].annotate(t, np.median(um[mm], 0), fontsize=9, weight="bold")
    q = f.obs["flag"].values
    axes[1].scatter(um[:, 0], um[:, 1], s=0.3, c=np.where(q == "low_quality", "#d62728", "#c8c8c8"),
                    rasterized=True)
    axes[1].set_title("Low-quality (ambient-high) subclusters in red")
    axes[1].set_xticks([]), axes[1].set_yticks([])
    fig.tight_layout()
    fig.savefig(fdir / "final_umap.png", dpi=150)
    plt.close(fig)
    f.write_h5ad(work_dir / "merged" / "final.h5ad", compression="gzip")


if __name__ == "__main__":
    main(*sys.argv[1:3])
