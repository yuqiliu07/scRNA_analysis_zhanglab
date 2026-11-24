"""Heart-failure (HF) gene panel and pathway list shared by 46_heart_failure.py.

GENES: classic HF / remodelling markers, grouped by process, with the direction they usually take in failing hearts
(+1 up, -1 down) and the cell type in which they are read.
THEMES: HF-related gene sets in the five pathway databases (GO BP, KEGG, Reactome, WikiPathways, MSigDB Hallmark).
EARLY: processes that already change in I4M cardiomyocytes at 4 weeks (electrical, fatty-acid oxidation, sarcomere).
"""

# (gene, expected direction in HF, cell type)
GENES = {
    "Stress / fetal genes": [("Nppa", 1, "Cardiomyocyte"), ("Nppb", 1, "Cardiomyocyte"), ("Myh7", 1, "Cardiomyocyte"),
                             ("Acta1", 1, "Cardiomyocyte"), ("Ankrd1", 1, "Cardiomyocyte"),
                             ("Xirp2", 1, "Cardiomyocyte"), ("Fhl1", 1, "Cardiomyocyte"), ("Myh6", -1, "Cardiomyocyte")],
    "Hypertrophy signalling": [("Rcan1", 1, "Cardiomyocyte"), ("Atf3", 1, "Cardiomyocyte"),
                               ("Ccn2", 1, "Cardiomyocyte")],
    "Ca2+ handling": [("Atp2a2", -1, "Cardiomyocyte"), ("Pln", -1, "Cardiomyocyte"), ("Ryr2", -1, "Cardiomyocyte"),
                      ("Slc8a1", 1, "Cardiomyocyte"), ("Camk2d", 1, "Cardiomyocyte"), ("Cacna1c", -1, "Cardiomyocyte")],
    "Electrical (K+ / Na+ / gap junction)": [("Kcnd2", -1, "Cardiomyocyte"), ("Kcnip2", -1, "Cardiomyocyte"),
                                             ("Kcnj2", -1, "Cardiomyocyte"), ("Kcnh2", -1, "Cardiomyocyte"),
                                             ("Kcnq1", -1, "Cardiomyocyte"), ("Scn5a", -1, "Cardiomyocyte"),
                                             ("Gja1", -1, "Cardiomyocyte")],
    "Energy metabolism": [("Cpt1b", -1, "Cardiomyocyte"), ("Ppargc1a", -1, "Cardiomyocyte"),
                          ("Ppara", -1, "Cardiomyocyte"), ("Acadm", -1, "Cardiomyocyte"),
                          ("Hadha", -1, "Cardiomyocyte"), ("Slc27a1", -1, "Cardiomyocyte"),
                          ("Slc2a1", 1, "Cardiomyocyte"), ("Hk2", 1, "Cardiomyocyte"),
                          ("Bckdhb", -1, "Cardiomyocyte"), ("Ppm1k", -1, "Cardiomyocyte")],
    "Sarcomere": [("Ttn", -1, "Cardiomyocyte"), ("Tnni3", -1, "Cardiomyocyte"), ("Tnnt2", -1, "Cardiomyocyte"),
                  ("Tnnc1", -1, "Cardiomyocyte"), ("Mybpc3", -1, "Cardiomyocyte"), ("Myl2", -1, "Cardiomyocyte"),
                  ("Actc1", -1, "Cardiomyocyte")],
    "Fibrosis (fibroblasts)": [("Postn", 1, "Fibroblast"), ("Col1a1", 1, "Fibroblast"), ("Col3a1", 1, "Fibroblast"),
                               ("Acta2", 1, "Fibroblast"), ("Ccn2", 1, "Fibroblast"), ("Tgfb1", 1, "Fibroblast"),
                               ("Lox", 1, "Fibroblast"), ("Timp1", 1, "Fibroblast"), ("Fn1", 1, "Fibroblast")],
    "Inflammation (macrophages)": [("Ccr2", 1, "Macrophage"), ("Il1b", 1, "Macrophage"), ("Ccl2", 1, "Macrophage"),
                                   ("Tnf", 1, "Macrophage"), ("Spp1", 1, "Macrophage"), ("Lyz2", 1, "Macrophage")],
}

# theme -> [(database, set name)]
THEMES = {
    "Stress / hypertrophy / cardiomyopathy": [
        ("GO_BP", "cardiac muscle hypertrophy (GO:0003300)"),
        ("GO_BP", "regulation of cardiac muscle hypertrophy (GO:0010611)"),
        ("GO_BP", "response to muscle stretch (GO:0035994)"),
        ("KEGG", "Hypertrophic cardiomyopathy"), ("KEGG", "Dilated cardiomyopathy"),
        ("KEGG", "Arrhythmogenic right ventricular cardiomyopathy"), ("KEGG", "Adrenergic signaling in cardiomyocytes"),
        ("KEGG", "cGMP-PKG signaling pathway"), ("KEGG", "Apelin signaling pathway"),
        ("WikiPathways", "WP_HYPERTROPHY_MODEL"), ("WikiPathways", "WP_MICRORNAS_IN_CARDIOMYOCYTE_HYPERTROPHY"),
        ("WikiPathways", "WP_REGULATION_OF_CARDIAC_HYPERTROPHY_BY_MIR208")],
    "Contraction / sarcomere": [
        ("GO_BP", "sarcomere organization (GO:0045214)"), ("GO_BP", "myofibril assembly (GO:0030239)"),
        ("GO_BP", "cardiac muscle contraction (GO:0060048)"), ("GO_BP", "muscle filament sliding (GO:0030049)"),
        ("KEGG", "Cardiac muscle contraction"), ("Reactome", "REACTOME_STRIATED_MUSCLE_CONTRACTION"),
        ("Reactome", "REACTOME_MUSCLE_CONTRACTION"), ("WikiPathways", "WP_STRIATED_MUSCLE_CONTRACTION"),
        ("Hallmark", "HALLMARK_MYOGENESIS")],
    "Ca2+ handling": [
        ("GO_BP", "regulation of cardiac muscle contraction by calcium ion signaling (GO:0010882)"),
        ("GO_BP", "calcium ion transport into cytosol (GO:0060402)"), ("KEGG", "Calcium signaling pathway"),
        ("Reactome", "REACTOME_ION_HOMEOSTASIS"), ("WikiPathways", "WP_CALCIUM_REGULATION_IN_CARDIAC_CELLS")],
    "Electrical remodelling": [
        ("GO_BP", "cardiac conduction (GO:0061337)"), ("GO_BP", "cardiac muscle cell action potential (GO:0086001)"),
        ("GO_BP", "membrane repolarization (GO:0086009)"), ("GO_BP", "potassium ion transmembrane transport (GO:0071805)"),
        ("GO_BP", "regulation of heart rate (GO:0002027)"), ("Reactome", "REACTOME_CARDIAC_CONDUCTION"),
        ("Reactome", "REACTOME_POTASSIUM_CHANNELS"), ("Reactome", "REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS"),
        ("Reactome", "REACTOME_PHASE_0_RAPID_DEPOLARISATION"), ("Reactome", "REACTOME_PHASE_3_RAPID_REPOLARISATION")],
    "Energy metabolism": [
        ("GO_BP", "fatty acid beta-oxidation (GO:0006635)"), ("GO_BP", "fatty acid oxidation (GO:0019395)"),
        ("GO_BP", "tricarboxylic acid cycle (GO:0006099)"), ("GO_BP", "oxidative phosphorylation (GO:0006119)"),
        ("GO_BP", "branched-chain amino acid catabolic process (GO:0009083)"),
        ("KEGG", "Fatty acid degradation"), ("KEGG", "PPAR signaling pathway"), ("KEGG", "Oxidative phosphorylation"),
        ("KEGG", "Citrate cycle (TCA cycle)"), ("KEGG", "Valine, leucine and isoleucine degradation"),
        ("KEGG", "Glycolysis / Gluconeogenesis"),
        ("Reactome", "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"),
        ("Reactome", "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA"),
        ("Reactome", "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"),
        ("Reactome", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"), ("Reactome", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"),
        ("WikiPathways", "WP_FATTY_ACID_BETAOXIDATION"), ("WikiPathways", "WP_PPAR_SIGNALING_PATHWAY"),
        ("WikiPathways", "WP_OXIDATIVE_PHOSPHORYLATION"), ("WikiPathways", "WP_TCA_CYCLE"),
        ("Hallmark", "HALLMARK_FATTY_ACID_METABOLISM"), ("Hallmark", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
        ("Hallmark", "HALLMARK_GLYCOLYSIS"), ("Hallmark", "HALLMARK_ADIPOGENESIS")],
    "Protein stress / hypoxia": [
        ("GO_BP", "response to unfolded protein (GO:0006986)"), ("GO_BP", "cellular response to hypoxia (GO:0071456)"),
        ("Reactome", "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR"), ("Reactome", "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA"),
        ("WikiPathways", "WP_OXIDATIVE_STRESS_RESPONSE"), ("Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
        ("Hallmark", "HALLMARK_HYPOXIA"), ("Hallmark", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")],
    "Fibrosis": [
        ("GO_BP", "extracellular matrix organization (GO:0030198)"), ("GO_BP", "collagen fibril organization (GO:0030199)"),
        ("GO_BP", "transforming growth factor beta receptor signaling pathway (GO:0007179)"),
        ("GO_BP", "wound healing (GO:0042060)"), ("KEGG", "ECM-receptor interaction"),
        ("KEGG", "TGF-beta signaling pathway"), ("Reactome", "REACTOME_COLLAGEN_FORMATION"),
        ("Reactome", "REACTOME_ECM_PROTEOGLYCANS"), ("Reactome", "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX"),
        ("WikiPathways", "WP_TGFBETA_SIGNALING_PATHWAY"), ("Hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
        ("Hallmark", "HALLMARK_TGF_BETA_SIGNALING")],
    "Inflammation": [
        ("GO_BP", "inflammatory response (GO:0006954)"), ("Reactome", "REACTOME_INTERLEUKIN_1_SIGNALING"),
        ("WikiPathways", "WP_INFLAMMATORY_RESPONSE_PATHWAY"), ("WikiPathways", "WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE"),
        ("Hallmark", "HALLMARK_INFLAMMATORY_RESPONSE"), ("Hallmark", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
        ("Hallmark", "HALLMARK_IL6_JAK_STAT3_SIGNALING")],
}
