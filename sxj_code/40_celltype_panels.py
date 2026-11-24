"""Whole-cell-type gene panels (per-mouse pseudobulk), same style as the cardiomyocyte / fibroblast panels.

  genotype/styled/heatmap_endocardial_panel.png   signals to the myocardium, Notch, flow / shear-responsive genes,
                                                  transcription factors, endocardial markers (reference)
  genotype/styled/heatmap_endothelial_panel.png   basement membrane, VEGF / angiogenesis, flow / vasotone, capillary
                                                  lipid uptake, inflammation, proliferation, EC markers (reference)
  genotype/styled/heatmap_fibroblast_panel.png    resting / activation states, TGF-b / PDGF receptors, growth factors,
                                                  immediate-early genes, proliferation, fibroblast markers (reference;
                                                  ECM genes are in heatmap_fibroblast_matrisome.png)
Gene groups were chosen from known biology, not from the significance ranking.
Usage: python 40_celltype_panels.py <work_dir> <results_dir>
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("sty", HERE / "15_genotype_figures_styled.py")
sty = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sty)

GROUPS = {
    "signals to myocardium": ["Nrg1", "Bmp6", "Pdgfc", "Nrp1", "Epha7", "Edn1", "Efnb2"],
    "Notch": ["Notch1", "Dll4", "Jag1", "Hes1"],
    "flow / shear-responsive": ["Klf2", "Klf4", "Piezo1", "Nos3", "Nfatc1"],
    "transcription factors": ["Gata6", "Meis2", "Sox5", "Runx1t1", "Ikzf2", "Gata4"],
    "endocardial markers": ["Npr3", "Cdh11", "Tmem108", "Hmcn1"],
}

EC_GROUPS = {
    "basement membrane": ["Col4a1", "Col4a2", "Lama4", "Hspg2", "Sparc"],
    "VEGF / angiogenesis": ["Kdr", "Flt1", "Nrp1", "Nrp2", "Apln", "Aplnr", "Angpt2", "Tek"],
    "flow / vasotone": ["Klf2", "Klf4", "Nos3", "Piezo1", "Edn1"],
    "capillary lipid uptake": ["Cd36", "Car4", "Gpihbp1", "Fabp4", "Lpl"],
    "inflammation / adhesion": ["Icam1", "Vcam1", "Selp"],
    "proliferation": ["Mki67", "Top2a"],
    "EC markers": ["Pecam1", "Cdh5", "Vwf"],
}

FIB_GROUPS = {
    "resting / universal": ["Pi16", "Dpt", "Cd248", "Tcf21", "Scara5", "Pdgfra"],
    "activation / myofibroblast": ["Postn", "Acta2", "Fap", "Cilp", "Meox1"],
    "TGF-b / PDGF receptors": ["Tgfbr1", "Tgfbr2", "Smad3", "Smad7", "Pdgfrb"],
    "growth factors / paracrine": ["Igf1", "Igf1r", "Vegfc", "Ptn", "Cxcl12", "Serpina3n"],
    "immediate-early / stress": ["Egr1", "Junb", "Nr4a1", "Hif1a"],
    "proliferation": ["Mki67", "Top2a"],
    "fibroblast markers": ["Dcn", "Col1a2", "Gsn"],
}


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"
    meta = pd.read_csv(HERE / "sample_sheet.csv").set_index("sample")
    pb = pd.read_csv(wd / "pb_Endocardial.csv.gz", index_col=0)[meta.index]
    de = pd.read_csv(wd / "de_Endocardial.csv.gz")
    sty.heatmap(pb, de, meta, GROUPS, out, "heatmap_endocardial_panel.png",
                "Endocardium: signalling, Notch, flow-responsive and transcription-factor genes")
    pb = pd.read_csv(wd / "pb_Endothelial.csv.gz", index_col=0)[meta.index]
    de = pd.read_csv(wd / "de_Endothelial.csv.gz")
    sty.heatmap(pb, de, meta, EC_GROUPS, out, "heatmap_endothelial_panel.png",
                "Endothelial cells: basement membrane, angiogenesis, flow, lipid uptake and inflammation genes")
    pb = pd.read_csv(wd / "pb_Fibroblast.csv.gz", index_col=0)[meta.index]
    de = pd.read_csv(wd / "de_Fibroblast.csv.gz")
    sty.heatmap(pb, de, meta, FIB_GROUPS, out, "heatmap_fibroblast_panel.png",
                "Fibroblasts: state, receptor, growth-factor, immediate-early and proliferation genes")
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
