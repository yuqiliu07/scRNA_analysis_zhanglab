"""tSNE of the final nuclei (counterpart of the provider's "Combination cluster tSNE"), split by genotype.

Computed with openTSNE (FFT-accelerated) on the Harmony-corrected PCs used for the UMAP, for the flag == "ok"
nuclei; perplexity 30, seed 0. Coordinates are cached in work_dir/genotype/tsne_ok.csv.gz.
Figure: results_dir/genotype/styled/tsne_split_genotype.png (same style as the split UMAP, 15_genotype_figures_styled.py).

Usage: python 16_tsne.py <work_dir> <results_dir>
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    s15 = _load("step15", "15_genotype_figures_styled.py")
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    obs = a.obs[["sample", "cell_type_final", "flag"]].copy()
    pcs = np.asarray(a.obsm["X_pca_harmony"])
    a.file.close()
    keep = (obs["flag"] == "ok").values
    obs, pcs = obs[keep], pcs[keep]
    cache = work_dir / "genotype" / "tsne_ok.csv.gz"
    if cache.exists():
        ts = pd.read_csv(cache, index_col=0).reindex(obs.index).values
    else:
        from openTSNE import TSNE
        print("tSNE on", pcs.shape, flush=True)
        ts = np.asarray(TSNE(perplexity=30, n_jobs=8, random_state=0, verbose=True).fit(pcs))
        pd.DataFrame(ts, index=obs.index, columns=["tSNE_1", "tSNE_2"]).to_csv(cache)
    b = ad.AnnData(obs=obs.assign(cell_type=obs["cell_type_final"].astype(str),
                                  genotype=obs["sample"].map(meta["genotype"]).astype(str)))
    b.obsm["X_tsne"] = ts
    out = res_dir / "genotype" / "styled"
    out.mkdir(parents=True, exist_ok=True)
    s15.umap_split(b, out, basis="X_tsne", xlab="tSNE_1", ylab="tSNE_2", fname="tsne_split_genotype.png")
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
