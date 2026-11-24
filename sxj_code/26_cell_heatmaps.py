"""Single-cell heatmaps (Seurat DoHeatmap look) of selected genes, WT vs TNNC1-I4M.

Columns: individual nuclei, 250 sampled per mouse (equal weight per mouse), ordered WT (A-D) then I4M (E-H) and by
mouse; colour bars on top for genotype and mouse. Rows: genes, grouped; only genes detected in >= 2% of the nuclei
of that cell type are drawn (many activation genes are sparse in single nuclei). Values: log-normalised expression scaled per gene over the drawn nuclei, clipped to
[-1.5, 2]; magenta (low) - black - yellow (high).

Gene choice (ours, not from an example figure):
  fibroblast    activated-state genes of the project's fibroblast-state method (src/v3/fib_states_v2.py, REF) and
                core-matrisome / ECM-regulator genes (as in 22_fibroblast_states.py)
  cardiomyocyte the curated panel used for the cardiomyocyte questions (sarcomere, calcium, stress, metabolism)
  endocardial / endothelial  the per-mouse panels of 40_celltype_panels.py (same genes, single nuclei)

Outputs: results_dir/genotype/styled/cellheatmap_fibroblast.png, cellheatmap_cardiomyocyte.png,
         cellheatmap_endocardial.png, cellheatmap_endothelial.png
Usage: python 26_cell_heatmaps.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
import numpy as np
import pandas as pd
import importlib.util

_spec = importlib.util.spec_from_file_location("panels", Path(__file__).with_name("40_celltype_panels.py"))
panels = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(panels)

plt.rcParams.update({"font.family": "Arial"})
GCOL = {"WT": "#B8682E", "I4M": "#3F8F6E"}
MCOL = ["#E6A570", "#C77B3F", "#9C5A26", "#6E3C14", "#8FD0B0", "#5FB08A", "#3A8A66", "#1F6045"]
N_PER_MOUSE, MIN_DET, SEED = 250, 0.02, 0
CMAP = LinearSegmentedColormap.from_list("purple_yellow", ["#FF00FF", "#000000", "#FFFF00"])

GENES = {
    "Fibroblast": {
        "activated state": ["Postn", "Cthrc1", "Comp", "Fap", "Tnc", "Lrrc15", "Col11a1", "Adam12", "Aspn"],
        "collagens": ["Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col6a1", "Col14a1", "Col15a1"],
        "ECM glycoproteins": ["Fn1", "Fbn1", "Eln", "Sparc", "Thbs2"],
        "ECM regulators / pro-fibrotic": ["Lox", "Loxl1", "Timp1", "Serpinh1", "Pcolce", "Ccn2", "Tgfb1"],
    },
    "Cardiomyocyte": {
        "thin filament / sarcomere": ["Tnnc1", "Tnni3", "Tnnt2", "Tpm1", "Actc1", "Myh6", "Myh7", "Myl2", "Myl3",
                                      "Mybpc3", "Ttn", "Actn2"],
        "calcium handling": ["Atp2a2", "Pln", "Ryr2", "Cacna1c", "Slc8a1", "Casq2", "Camk2d", "Trdn"],
        "stress / hypertrophy": ["Nppa", "Nppb", "Acta1", "Ankrd1", "Rcan1", "Fhl1", "Myh7b"],
        "metabolism": ["Ppargc1a", "Ppara", "Cpt1b", "Acadm", "Hadha", "Pdk4", "Ucp3", "Slc2a4", "Hk2", "Ldha"],
    },
    "Endocardial": panels.GROUPS,
    "Endothelial": panels.EC_GROUPS,
}


def draw(a, meta, cell_type, groups, out, fname):
    b = a[(a.obs["cell_type_final"].astype(str) == cell_type).values]
    rng = np.random.default_rng(SEED)
    order = list(meta.index[meta.genotype == "WT"]) + list(meta.index[meta.genotype == "I4M"])
    idx = []
    for s in order:
        m = np.where((b.obs["sample"] == s).values)[0]
        idx.append(rng.choice(m, min(N_PER_MOUSE, len(m)), replace=False))
    sizes = [len(i) for i in idx]
    idx = np.concatenate(idx)
    genes, gstart = [], []
    for g, gl in groups.items():
        gl = [x for x in gl if x in b.var_names]
        det = np.asarray((b[:, gl].X > 0).mean(0)).ravel()
        gl = [x for x, d in zip(gl, det) if d >= MIN_DET]
        gstart.append((g, len(genes), len(gl)))
        genes += gl
    X = b[idx][:, genes].X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float).T
    Z = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True).clip(1e-9)
    Z = Z.clip(-1.5, 2)

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(13, 0.24 * len(genes) + 1.4))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[0.32, 0.2, 0.24 * len(genes)], width_ratios=[0.74, 0.2, 0.06],
                  hspace=0.06, wspace=0.02, left=0.11, right=0.97, top=0.93, bottom=0.04)
    axg = fig.add_subplot(gs[0, 0]); axm = fig.add_subplot(gs[1, 0]); ax = fig.add_subplot(gs[2, 0])
    axl = fig.add_subplot(gs[2, 1]); cax_host = fig.add_subplot(gs[2, 2])
    im = ax.imshow(Z, aspect="auto", cmap=CMAP, vmin=-1.5, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(genes)), genes, fontsize=9, style="italic")
    ax.set_xticks([])
    for _, st, n in gstart[1:]:
        ax.axhline(st - 0.5, color="white", lw=1.5)
    axl.set_ylim(len(genes) - 0.5, -0.5); axl.set_xlim(0, 1); axl.axis("off")
    for g, st, n in gstart:
        if n:
            axl.plot([0.03, 0.03], [st - 0.3, st + n - 0.7], color="#555555", lw=1)
            axl.text(0.08, st + n / 2 - 0.5, g, fontsize=9, va="center", ha="left")
    bounds = np.cumsum([0] + sizes)
    gcol = [to_rgb(GCOL[meta.loc[s, "genotype"]]) for s in order for _ in range(sizes[order.index(s)])]
    axg.imshow(np.array([gcol]), aspect="auto")
    nwt = sum(sizes[:4])
    axg.set_title("WT" + " " * 60 + "I4M", fontsize=12)
    mcol = [to_rgb(MCOL[i]) for i, s in enumerate(order) for _ in range(sizes[i])]
    axm.imshow(np.array([mcol]), aspect="auto")
    for i, s in enumerate(order):
        axm.text((bounds[i] + bounds[i + 1]) / 2, 0, f"{s} ({meta.loc[s, 'sex']})", ha="center", va="center",
                 fontsize=7, color="white")
    for a_ in (axg, axm):
        a_.set_xticks([]); a_.set_yticks([])
        for s_ in a_.spines.values():
            s_.set_visible(False)
    ax.axvline(nwt - 0.5, color="white", lw=2)
    cax_host.axis("off")
    cax = cax_host.inset_axes([0.25, 0.55, 0.35, 0.4])
    fig.colorbar(im, cax=cax)
    cax.set_title("Expression", fontsize=9, loc="left")
    fig.text(0.11, 0.005, f"{cell_type} nuclei: {N_PER_MOUSE} per mouse; expression scaled per gene; genes detected in "
             f">= {int(MIN_DET * 100)}% of nuclei", fontsize=8, color="#555555")
    fig.savefig(out / fname, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(cell_type, len(genes), "genes", len(idx), "nuclei", flush=True)


def main(work_dir, res_dir):
    work_dir, res_dir = Path(work_dir), Path(res_dir)
    out = res_dir / "genotype" / "styled"
    meta = pd.read_csv(Path(__file__).with_name("sample_sheet.csv")).set_index("sample")
    a = ad.read_h5ad(work_dir / "merged" / "final.h5ad", backed="r")
    keep = ((a.obs["flag"] == "ok") & a.obs["cell_type_final"].astype(str).isin(list(GENES))).values
    a = a[keep].to_memory()
    draw(a, meta, "Fibroblast", GENES["Fibroblast"], out, "cellheatmap_fibroblast.png")
    draw(a, meta, "Cardiomyocyte", GENES["Cardiomyocyte"], out, "cellheatmap_cardiomyocyte.png")
    draw(a, meta, "Endocardial", GENES["Endocardial"], out, "cellheatmap_endocardial.png")
    draw(a, meta, "Endothelial", GENES["Endothelial"], out, "cellheatmap_endothelial.png")


if __name__ == "__main__":
    main(*sys.argv[1:3])
