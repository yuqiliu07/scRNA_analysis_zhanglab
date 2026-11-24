"""Representative plot for each pathway method (23_pathways.R results; gene sets from 27_export_gmt.R).

  styled/pw_ora_dotplot.png     ORA: compareCluster-style dot plot (x = comparison, size = gene ratio, colour = padj)
  styled/pw_ora_cnet.png        ORA: gene-concept network (pathways linked to their DEGs, genes coloured by log2FC)
  styled/pw_gsea_curves.png     GSEA: running enrichment-score plots (ES curve, hit barcode, ranked metric)
  styled/pw_gsea_emap.png       GSEA: enrichment map (pathways linked by shared genes; node colour = NES), CM + Endo
  styled/pw_gsea_emap_fib_ec.png  same for fibroblasts + endothelial cells
  styled/pw_hallmark_heatmap.png  MSigDB Hallmark GSEA: NES across all cell types (* padj <= 0.05; shared sets marked)

Usage: python 28_pathway_method_plots.py <work_dir> <results_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import networkx as nx
import numpy as np
import pandas as pd

plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False, "axes.spines.right": False})
DBL = {"GO_BP": "GO", "KEGG": "KEGG", "Reactome": "Reactome", "WikiPathways": "WikiPW", "Hallmark": "Hallmark"}
CTS = ["Cardiomyocyte", "Endocardial", "Fibroblast", "Endothelial", "Lymphatic_EC", "Pericyte", "SMC", "Macrophage",
       "Lymphocyte", "Neural_Schwann"]
SHORT = {"Cardiomyocyte": "CM", "Endocardial": "Endo", "Fibroblast": "FB"}


def short(name, n=45):
    s = name.replace("HALLMARK_", "").replace("REACTOME_", "").replace("WP_", "").replace("_", " ")
    s = s if not s.isupper() else s.capitalize()
    return s if len(s) <= n else s[:n - 3] + "..."


def read_gmt(wd, db):
    out = {}
    for line in open(wd / f"gmt_{db}.gmt", encoding="utf8"):
        p = line.rstrip("\n").split("\t")
        out[p[0]] = set(p[2:])
    return out


def ranking(wd, ct):
    d = pd.read_csv(wd / f"de_{ct}.csv.gz")
    d = d[~(d.technical | d.ambient)]
    st = np.sign(d.logFC) * -np.log10(d.PValue.clip(1e-300))
    s = pd.Series(st.values, index=d.gene).groupby(level=0).first().sort_values(ascending=False)
    return s


def running_es(rank, genes):
    hit = rank.index.isin(genes)
    w = np.abs(rank.values) * hit
    p_hit = np.cumsum(w) / max(w.sum(), 1e-12)
    p_miss = np.cumsum(~hit) / max((~hit).sum(), 1)
    return p_hit - p_miss, hit


def main(work_dir, res_dir):
    wd = Path(work_dir) / "genotype"
    out = Path(res_dir) / "genotype" / "styled"

    # ---------------- ORA dot plot ----------------
    groups = [("Cardiomyocyte", "Up"), ("Cardiomyocyte", "Down"), ("Endocardial", "Up"), ("Endocardial", "Down"),
              ("Fibroblast", "Up")]
    rows = []
    for ct, d in groups:
        tabs = []
        for db in DBL:
            f = wd / f"pw_ORA_{db}_{ct}_{d}.csv.gz"
            if f.exists():
                t = pd.read_csv(f)
                t["db"] = db
                tabs.append(t)
        t = pd.concat(tabs).sort_values("p")
        n_in = None
        de = pd.read_csv(wd / f"de_{ct}.csv.gz")
        n_in = int((de.robust & ((de.logFC > 0) if d == "Up" else (de.logFC < 0))).sum())
        for r in t.head(4).itertuples():
            rows.append({"group": f"{SHORT[ct]} {d.lower()}\n(n={n_in})", "term": f"{short(r.set, 48)} [{DBL[r.db]}]",
                         "ratio": r.overlap / n_in, "padj": r.padj, "count": r.overlap})
    o = pd.DataFrame(rows)
    terms = list(dict.fromkeys(o.term))
    grps = list(dict.fromkeys(o.group))
    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    cm = LinearSegmentedColormap.from_list("ora", ["#D7263D", "#F2A3A8", "#3B5BA9"])
    sc_ = ax.scatter([grps.index(g) for g in o.group], [terms.index(t) for t in o.term], s=o.ratio * 2500,
                     c=np.log10(o.padj.clip(1e-6)), cmap=cm, vmin=-4, vmax=0, edgecolors="#333333", lw=0.4)
    ax.set_xticks(range(len(grps)), grps)
    ax.set_yticks(range(len(terms)), terms, fontsize=8)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    ax.set_xlim(-0.6, len(grps) - 0.4)
    ax.grid(color="#EEEEEE", lw=0.6); ax.set_axisbelow(True)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.03, pad=0.02)
    cb.set_ticks([-4, -3, -2, -1, 0], labels=["1e-4", "1e-3", "0.01", "0.1", "1"])
    cb.set_label("p.adjust")
    hs = [plt.scatter([], [], s=v * 2500, c="#999999") for v in [0.02, 0.05, 0.1]]
    ax.legend(hs, ["0.02", "0.05", "0.10"], title="GeneRatio", loc="upper left", bbox_to_anchor=(1.18, 1),
              frameon=False, labelspacing=1.3)
    ax.set_title("ORA (over-representation analysis): top 4 pathways per comparison, all databases", fontsize=11)
    fig.savefig(out / "pw_ora_dotplot.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---------------- ORA cnet plot (cardiomyocyte, down) ----------------
    de = pd.read_csv(wd / "de_Cardiomyocyte.csv.gz").set_index("gene")
    picks = [("GO_BP", "homophilic cell-cell adhesion (GO:0007156)"), ("KEGG", "Protein processing in endoplasmic reticulum"),
             ("KEGG", "Antigen processing and presentation"), ("Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")]
    G = nx.Graph()
    for db, s in picks:
        t = pd.read_csv(wd / f"pw_ORA_{db}_Cardiomyocyte_Down.csv.gz").set_index("set")
        if s not in t.index:
            continue
        lab = f"{short(s, 44)}\n[{DBL[db]}, padj {t.loc[s, 'padj']:.1g}]"
        G.add_node(lab, kind="term", size=int(t.loc[s, "overlap"]))
        for g in str(t.loc[s, "genes"]).split(","):
            G.add_node(g, kind="gene")
            G.add_edge(lab, g)
    # radial layout: pathways on a row; genes found in only one pathway on a circle around it (labels pushed outward
    # along the radius); shared genes on two rows below, at least 1 unit apart so that labels do not overlap
    terms_n = [n for n, a in G.nodes(data=True) if a["kind"] == "term"]
    genes_n = [n for n, a in G.nodes(data=True) if a["kind"] == "gene"]
    tpos = {t: np.array([i * 4.3, 0.0]) for i, t in enumerate(terms_n)}
    pos, lab_pos = dict(tpos), {}
    own = {t: [g for g in G.neighbors(t) if G.degree(g) == 1] for t in terms_n}
    R = 1.45
    for t, gl in own.items():
        for k, g in enumerate(sorted(gl)):
            ang = 2 * np.pi * k / max(len(gl), 1) + np.pi / 2
            u = np.array([np.cos(ang), np.sin(ang)])
            pos[g] = tpos[t] + R * u
            ha = "center" if abs(u[0]) < 0.3 else ("left" if u[0] > 0 else "right")
            va = "center" if abs(u[1]) < 0.5 else ("bottom" if u[1] > 0 else "top")
            lab_pos[g] = (tpos[t] + (R + 0.22) * u, ha, va)
    shared_g = sorted([g for g in genes_n if G.degree(g) > 1],
                      key=lambda g: np.mean([tpos[t][0] for t in G.neighbors(g)]))
    last = {0: -1e9, 1: -1e9}
    for k, g in enumerate(shared_g):
        row = k % 2
        cx = np.mean([tpos[t][0] for t in G.neighbors(g)])
        x = max(cx, last[row] + 1.15)
        last[row] = x
        pos[g] = np.array([x, -2.4 - 0.95 * row])
        lab_pos[g] = (pos[g] + np.array([0, -0.2]), "center", "top")
    fig, ax = plt.subplots(figsize=(4.3 * len(terms_n) + 1.5, 6.4))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#BBBBBB", width=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=terms_n, node_color="#E8C27A", ax=ax, node_size=900, edgecolors="#555555")
    lf = [de.loc[g, "logFC"] if g in de.index else 0 for g in genes_n]
    nn = nx.draw_networkx_nodes(G, pos, nodelist=genes_n, node_color=lf, cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                                node_size=170, ax=ax, edgecolors="#555555", linewidths=0.4)
    for n in terms_n:
        ax.text(pos[n][0], pos[n][1] + R + 0.75, n, fontsize=8.5, ha="center", va="bottom", weight="bold")
    texts = []
    for n in genes_n:
        (x, y), ha, va = lab_pos[n]
        texts.append(ax.text(x, y, n, fontsize=7, ha=ha, va=va, style="italic"))
    ax.margins(0.08, 0.15)
    cb = fig.colorbar(nn, ax=ax, fraction=0.03)
    cb.set_label("log2FC in cardiomyocytes (I4M vs WT)")
    ax.set_title("ORA gene-concept network (cnet plot): cardiomyocyte, genes down in I4M", fontsize=11, pad=36)
    ax.axis("off")
    fig.canvas.draw()
    bb = [t.get_window_extent() for t in texts]
    clash = [(texts[i].get_text(), texts[j].get_text()) for i in range(len(bb)) for j in range(i + 1, len(bb))
             if bb[i].overlaps(bb[j])]
    print("cnet label overlaps:", clash or "none", flush=True)
    fig.savefig(out / "pw_ora_cnet.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---------------- GSEA running enrichment curves ----------------
    gmt = {db: read_gmt(wd, db) for db in DBL}
    curves = [("Cardiomyocyte", "KEGG", "Valine, leucine and isoleucine degradation"),
              ("Cardiomyocyte", "Hallmark", "HALLMARK_INFLAMMATORY_RESPONSE"),
              ("Endocardial", "GO_BP", "angiogenesis (GO:0001525)"),
              ("Cardiomyocyte", "Reactome", "REACTOME_MRNA_SPLICING")]
    fig = plt.figure(figsize=(17, 4.8))
    gs = fig.add_gridspec(3, len(curves), height_ratios=[3, 0.5, 1.4], hspace=0.05, wspace=0.28)
    for j, (ct, db, s) in enumerate(curves):
        rank = ranking(wd, ct)
        es, hit = running_es(rank, gmt[db][s])
        g = pd.read_csv(wd / f"pw_GSEA_{db}_{ct}.csv.gz").set_index("set")
        a0, a1, a2 = fig.add_subplot(gs[0, j]), fig.add_subplot(gs[1, j]), fig.add_subplot(gs[2, j])
        x = np.arange(len(rank))
        a0.plot(x, es, color="#2CA02C", lw=1.8)
        a0.axhline(0, color="#999999", lw=0.6)
        a0.set_xlim(0, len(rank)); a0.set_xticks([])
        a0.set_ylabel("Enrichment score" if j == 0 else "")
        tag = " (shared, batch-like)" if "SPLICING" in s else ""
        a0.set_title(f"{ct}: {short(s, 40)}\n[{DBL[db]}]{tag}  NES {g.loc[s, 'NES']:+.2f}, padj {g.loc[s, 'padj']:.2g}",
                     fontsize=9.5)
        a1.vlines(x[hit], 0, 1, color="black", lw=0.5)
        a1.set_xlim(0, len(rank)); a1.set_yticks([]); a1.set_xticks([])
        for sp in a1.spines.values():
            sp.set_visible(False)
        a2.fill_between(x, rank.values, color="#9E9E9E", lw=0)
        a2.set_xlim(0, len(rank))
        a2.set_xlabel("Genes ranked: up in I4M  ->  down in I4M", fontsize=8.5)
        a2.set_ylabel("Ranked metric" if j == 0 else "", fontsize=8.5)
    fig.suptitle("GSEA (gene set enrichment analysis): running enrichment score", fontsize=12, y=1.04)
    fig.savefig(out / "pw_gsea_curves.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ---------------- GSEA enrichment map ----------------
    allg = []
    for ct in CTS:
        for db in DBL:
            f = wd / f"pw_GSEA_{db}_{ct}.csv.gz"
            if f.exists():
                t = pd.read_csv(f); t["db"] = db; t["cell_type"] = ct
                allg.append(t)
    A = pd.concat(allg)
    sig = A[A.padj <= 0.05].copy()
    sig["key"] = sig.db + "|" + sig.set + "|" + np.sign(sig.NES).astype(str)
    nshare = sig.groupby("key").cell_type.nunique()
    A["key"] = A.db + "|" + A.set + "|" + np.sign(A.NES).astype(str)
    A["n_same"] = A.key.map(nshare).fillna(0)
    # each module (connected component) gets a number next to its strongest pathway; the names are listed beside the
    # network instead of being written on it, so no text overlaps
    def emap(cts, fname):
        fig = plt.figure(figsize=(21, 9.5))
        gsp = fig.add_gridspec(1, 5, width_ratios=[1.25, 1.0, 1.25, 1.0, 0.06], wspace=0.04)
        for k_, ct in enumerate(cts):
            ax, axl = fig.add_subplot(gsp[0, 2 * k_]), fig.add_subplot(gsp[0, 2 * k_ + 1])
            t = A[(A.cell_type == ct) & (A.padj <= 0.05) & (A.n_same < 3) & (A.db != "WikiPathways")].sort_values("p").head(45)
            sets = {f"{r.db}|{r.set}": gmt[r.db].get(r.set, set()) for r in t.itertuples()}
            G = nx.Graph()
            for r in t.itertuples():
                G.add_node(f"{r.db}|{r.set}", NES=r.NES, size=r.size, label=f"{short(r.set, 48)} [{DBL[r.db]}]")
            keys = list(sets)
            for i in range(len(keys)):
                for k in range(i + 1, len(keys)):
                    a_, b_ = sets[keys[i]], sets[keys[k]]
                    if a_ and b_:
                        oc = len(a_ & b_) / min(len(a_), len(b_))  # overlap coefficient (EnrichmentMap default 0.5)
                        if oc >= 0.5:
                            G.add_edge(keys[i], keys[k], w=oc)
            comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), -max(abs(G.nodes[n]["NES"]) for n in c)))
            ncol = int(np.ceil(np.sqrt(len(comps) * 1.2)))
            pos = {}
            for ci, c in enumerate(comps):
                sub = G.subgraph(c)
                p_ = nx.spring_layout(sub, seed=1, k=0.6, iterations=200) if len(c) > 1 else {next(iter(c)): np.zeros(2)}
                arr = np.array(list(p_.values())); ctr = arr.mean(0)
                sc = 0.34 / max(np.abs(arr - ctr).max(), 1e-9) if len(c) > 1 else 0
                r, q = divmod(ci, ncol)
                for n_, v in p_.items():
                    pos[n_] = (np.asarray(v) - ctr) * sc + np.array([q, -r])
            nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#BBBBBB", width=[1.5 * G.edges[e]["w"] for e in G.edges])
            nodes = list(G.nodes)
            nn = nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax, node_color=[G.nodes[n]["NES"] for n in nodes],
                                        cmap="RdBu_r", vmin=-2.5, vmax=2.5, edgecolors="#555555", linewidths=0.4,
                                        node_size=[40 + G.nodes[n]["size"] * 1.0 for n in nodes])
            lines = []
            for ci, c in enumerate(comps, 1):
                top = max(c, key=lambda n: abs(G.nodes[n]["NES"]))
                cx = np.mean([pos[n][0] for n in c]); ymax = max(pos[n][1] for n in c)
                ax.text(cx, ymax + 0.2, str(ci), fontsize=9, weight="bold", ha="center", va="bottom",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#888888", lw=0.5))
                more = f"  (+{len(c) - 1})" if len(c) > 1 else ""
                arrow = "↑" if G.nodes[top]["NES"] > 0 else "↓"
                lines.append((f"{ci:>2}  {arrow} {G.nodes[top]['label']}{more}", "#B2182B" if arrow == "↑" else "#2166AC"))
            ax.set_xlim(-0.6, ncol - 0.4); ax.set_ylim(-(len(comps) - 1) // ncol - 0.6, 0.75)
            ax.set_title(f"{ct}: cell-type-specific pathways", fontsize=12, loc="left")
            ax.axis("off")
            axl.axis("off")
            fs_ = 8.5 if len(lines) <= 26 else 7.5
            for i, (txt_, col) in enumerate(lines):
                axl.text(0.0, 1 - i / max(len(lines), 26) * 0.97, txt_, fontsize=fs_, color=col, va="top", transform=axl.transAxes)
        cax = fig.add_subplot(gsp[0, 4])
        cb = fig.colorbar(nn, cax=cax)
        cb.set_label("NES (red = up in I4M, blue = down)")
        fig.suptitle("GSEA enrichment map: nodes = pathways (size = genes in set), edges = shared genes (overlap >= 0.5); "
                     "numbers = modules, listed with their strongest pathway (+n = further pathways in the module)", fontsize=11.5)
        fig.savefig(out / fname, dpi=220, bbox_inches="tight")
        plt.close(fig)

    emap(["Cardiomyocyte", "Endocardial"], "pw_gsea_emap.png")
    emap(["Fibroblast", "Endothelial"], "pw_gsea_emap_fib_ec.png")

    # ---------------- Hallmark NES heatmap ----------------
    H = A[A.db == "Hallmark"]
    keep = H[H.padj <= 0.05].set.unique()
    M = H[H.set.isin(keep)].pivot_table(index="set", columns="cell_type", values="NES").reindex(columns=CTS)
    P = H[H.set.isin(keep)].pivot_table(index="set", columns="cell_type", values="padj").reindex(columns=CTS)
    shared = H[H.set.isin(keep)].groupby("set").n_same.max()
    M = M.loc[M.mean(1).sort_values().index]
    fig, ax = plt.subplots(figsize=(9.5, 0.33 * len(M) + 2))
    im = ax.imshow(M.values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if P.loc[M.index[i], M.columns[j]] <= 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=10)
    ax.set_xticks(range(len(CTS)), [c.replace("_", " ") for c in CTS], rotation=40, ha="right")
    ax.set_yticks(range(len(M)), [short(s, 40) + ("  (shared)" if shared[s] >= 3 else "") for s in M.index], fontsize=8.5)
    for lab_, s in zip(ax.get_yticklabels(), M.index):
        if shared[s] >= 3:
            lab_.set_color("#888888")
    cb = fig.colorbar(im, ax=ax, fraction=0.03)
    cb.set_label("NES (I4M vs WT)")
    ax.set_title("MSigDB Hallmark, GSEA: NES in every cell type\n(* padj <= 0.05; grey = significant in >= 3 cell types, "
                 "batch-like)", fontsize=10.5)
    fig.savefig(out / "pw_hallmark_heatmap.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("done", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
