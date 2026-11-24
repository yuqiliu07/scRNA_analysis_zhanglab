"""Top-10 marker heatmap of the final cell types, drawn like the provider's Seurat DoHeatmap for comparison.

Rows: top 10 markers per type (08_compare_provider.py, cmp_markers_ours.csv: Wilcoxon on CellBender counts),
duplicates dropped. Columns: nuclei grouped by type; 10% of each type is drawn (at least 100, or all if fewer),
so widths are roughly proportional to type size as in the provider's figure, except that each block is stretched
to at least 3% of the width so the small types stay readable. Values: log-normalized expression
z-scored per gene over all final nuclei and clipped to [-2.5, 2.5] (the DoHeatmap defaults).

Usage: python 11_marker_heatmap.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ORDER = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Endocardial", "Lymphatic_EC", "Pericyte", "SMC",
         "Macrophage", "Lymphocyte", "Dendritic", "Epicardial", "Neural_Schwann"]
ZH = {"Cardiomyocyte": "心肌", "Fibroblast": "成纤维", "Endothelial": "内皮", "Endocardial": "心内膜",
      "Lymphatic_EC": "淋巴管内皮", "Pericyte": "周细胞", "SMC": "平滑肌", "Macrophage": "巨噬",
      "Lymphocyte": "淋巴细胞", "Dendritic": "树突状细胞", "Epicardial": "心外膜", "Neural_Schwann": "神经/施万"}
FRAC, MIN_N, CLIP, SEED = 0.10, 100, 2.5, 0
MIN_W = 0.03  # minimum display width of a type block (fraction of the drawn nuclei), so small types stay readable
# same colours as the numbered final UMAP (tab20: dark shades for types 1-10, then light blue / light orange)
COLORS = [plt.get_cmap("tab20")(i) for i in list(range(0, 20, 2)) + [1, 3]]


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    mk = pd.read_csv(res_dir / "provider_compare" / "cmp_markers_ours.csv").set_index("cell_type")
    genes = []
    for t in ORDER:
        for g in mk.loc[t, "top_markers"].split(", ")[:10]:
            if g not in genes:
                genes.append(g)
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    ct = a.obs["cell_type_final"].astype(str).values
    x = a[:, genes].to_memory().X
    x = np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float32)
    a.file.close()
    z = np.clip((x - x.mean(0)) / x.std(0).clip(1e-6), -CLIP, CLIP)

    rng = np.random.default_rng(SEED)
    cols, bounds = [], [0]
    for t in ORDER:
        idx = np.where(ct == t)[0]
        n = min(len(idx), max(MIN_N, int(round(FRAC * len(idx)))))
        cols.append(rng.choice(idx, n, replace=False))
        bounds.append(bounds[-1] + n)
    cols = np.concatenate(cols)
    m = z[cols].T

    # display positions: each block is stretched to at least MIN_W of the total width
    w = np.diff(bounds).astype(float)
    w = np.maximum(w, MIN_W * w.sum())
    xs = np.concatenate([[0], np.cumsum(w)])
    ng = len(genes)
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_axes([0.08, 0.04, 0.80, 0.80])
    axb = fig.add_axes([0.08, 0.845, 0.80, 0.02])
    for i in range(len(ORDER)):
        im = ax.imshow(m[:, bounds[i]:bounds[i + 1]], aspect="auto", cmap="RdBu_r", vmin=-CLIP, vmax=CLIP,
                       interpolation="nearest", extent=(xs[i], xs[i + 1], ng - 0.5, -0.5))
    for b in xs[1:-1]:
        ax.axvline(b, color="#808080", lw=1)
    ax.set_xlim(0, xs[-1])
    ax.set_ylim(ng - 0.5, -0.5)
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes, fontsize=5)
    ax.set_xticks([])
    for i, t in enumerate(ORDER):
        axb.axvspan(xs[i], xs[i + 1], color=COLORS[i])
        axb.text((xs[i] + xs[i + 1]) / 2, 1.3, f"{i + 1} {ZH[t]}", rotation=45, ha="left", va="bottom",
                 fontsize=9)
    axb.set_xlim(0, xs[-1])
    axb.axis("off")
    cax = fig.add_axes([0.90, 0.25, 0.015, 0.25])
    fig.colorbar(im, cax=cax).set_label("缩放后表达量", fontsize=9)
    fig.savefig(res_dir / "figures" / "final_top10_marker_heatmap.png", dpi=200)
    plt.close(fig)
    print("genes", len(genes), "nuclei drawn", len(cols), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
