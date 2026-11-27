#!/usr/bin/env python3
"""Part 4: pathways and biology.

Main result : GSEA on the full ranked gene list (DESeq2 Wald stat).
Supporting  : g:Profiler over-representation on genes with raw p < 0.01 (exploratory:
              only 6 genes pass FDR 0.05, so an FDR-based list would be useless).

Run from Part4: python3 scripts/part4_pathways_v2.py
Each run writes to a fresh directory; the submitted results are not overwritten.
"""
import argparse, io, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
from pathway_utils import LIBRARIES, get_gmt, rank_genes, new_dir, lead_genes

# hypoxia / stress genes that Part 2 reports as driven by P7_T2D and P8_T2D.
STRESS = ["AK4", "LOX", "IGFBP3", "LUCAT1", "ANGPTL4", "SLC16A3", "GDF15",
          "THBS2", "PLOD2", "G0S2", "TMEM158", "TREM1", "HIC1"]
# beta-cell / known T2D islet genes that Part 2 reports as trending down.
BETA = ["SLC2A2", "RASGRP1", "PPP1R1A", "GLRA1", "CAPN13", "CHL1", "PDX1",
        "IAPP", "ABCC8", "SLC30A8", "INS", "GCG", "TCF7L2", "MAFA", "KCNJ11"]

ap = argparse.ArgumentParser()
ap.add_argument("--de", default="../Part2_statistics/results/de_results.csv")
ap.add_argument("--output-dir", default="reruns/primary")
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--skip-ora", action="store_true")
ap.add_argument("--nperm", type=int, default=10000)
ap.add_argument("--p-cut", type=float, default=0.01)
a = ap.parse_args()

root = new_dir(a.output_dir)
RES, FIG = str(root / "results"), str(root / "figures")
for d in (RES, FIG):
    os.makedirs(d, exist_ok=True)

log_lines = []
def say(*x):
    s = " ".join(str(i) for i in x)
    print(s, flush=True)
    log_lines.append(s)

# ---------------------------------------------------------------- ranked list
de = pd.read_csv(a.de)
rk = rank_genes(de)
rk.to_csv(f"{RES}/ranked_list.csv", index=False)
say(f"ranked genes: {len(rk)} of {len(de)} tested "
    f"({de.symbol.isna().sum()} have no symbol, {de.symbol.dropna().duplicated().sum()} duplicate symbols)")

# ---------------------------------------------------------------- GSEA
libs = list(LIBRARIES)
allres = []
for lib in libs:
    name = os.path.basename(lib).replace(".gmt", "")
    out_csv = f"{RES}/gsea_{name}.csv"
    path = get_gmt(lib)
    pre = gp.prerank(rnk=rk, gene_sets=path, min_size=15, max_size=500,
                     permutation_num=a.nperm, seed=7, threads=a.threads,
                     outdir=None, verbose=False)
    r = pre.res2d.copy()
    r.insert(0, "library", name)
    r.to_csv(out_csv, index=False)
    allres.append(r)
    fdr = r["FDR q-val"].astype(float)
    say(f"{name}: {len(r)} sets tested, FDR<0.05: {int((fdr<0.05).sum())}")
    for i, t in enumerate(r.assign(f=fdr).sort_values("f").head(3)["Term"].tolist()):
        try:
            ax = pre.plot(terms=t, show_ranking=True, figsize=(5, 6))
            ax.savefig(f"{FIG}/gsea_{name}_top{i+1}.png", dpi=200, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            raise RuntimeError(f"Enrichment plot failed for {t}") from e

G = pd.concat(allres, ignore_index=True)
for c in ["NES", "NOM p-val", "FDR q-val"]:
    G[c] = pd.to_numeric(G[c], errors="coerce")
G["n_lead"] = G["Lead_genes"].map(lambda s: len(lead_genes(s)))
G["stress_frac"] = G["Lead_genes"].map(
    lambda s: round(sum(g in STRESS for g in lead_genes(s)) / max(len(lead_genes(s)), 1), 2))
G.sort_values("FDR q-val").to_csv(f"{RES}/gsea_all.csv", index=False)

# ---------------------------------------------------------------- dot plots
def dotplot(d, title, path, n=10):
    if d.empty:
        return False
    d = d.copy()
    d["label"] = d["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True).str.slice(0, 55)
    d = pd.concat([d[d.NES > 0].nlargest(n, "NES"), d[d.NES < 0].nsmallest(n, "NES")]).drop_duplicates("Term").sort_values("NES")
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 1.8))
    sc = ax.scatter(d["NES"], d["label"], s=20 + 12 * d["n_lead"],
                    c=d["FDR q-val"], cmap="viridis_r", edgecolor="k", lw=.4)
    ax.axvline(0, c="grey", lw=.8)
    lo, hi = d.NES.min(), d.NES.max()
    pad = max(.35, .12*(hi-lo))
    ax.set_xlim(min(0,lo)-pad,max(0,hi)+pad)
    ax.margins(y=.06)
    ax.set_xlabel("NES (positive = higher in T2D)")
    ax.set_title(title, fontsize=10)
    fig.colorbar(sc, label="FDR q-value", shrink=.6)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
    return True

for lib, d in G.groupby("library"):
    sig = d[d["FDR q-val"] < 0.05]
    note = " (FDR < 0.05)"
    if sig.empty:  # fall back so slide 10 always gets a figure
        sig = d.reindex(d["NES"].abs().sort_values(ascending=False).index).head(20)
        note = " (no set reached FDR 0.05; exploratory NES)"
    if dotplot(sig, f"GSEA {lib}: T2D vs ND{note}", f"{FIG}/dotplot_gsea_{lib}.png"):
        say(f"dot plot written for {lib}{note}")

# ---------------------------------------------------------------- g:Profiler ORA
if not a.skip_ora:
    try:
        from gprofiler import GProfiler
        gpf = GProfiler(return_dataframe=True)
        bg = de["ensembl"].dropna().unique().tolist()
        sets = {"up":   de[(de.pvalue < a.p_cut) & (de.log2FoldChange > 0)]["ensembl"],
                "down": de[(de.pvalue < a.p_cut) & (de.log2FoldChange < 0)]["ensembl"]}
        ora = []
        for k, ids in sets.items():
            ids = ids.dropna().tolist()
            say(f"ORA {k}: {len(ids)} genes at raw p<{a.p_cut}")
            if not ids:
                say("  empty query, skipped")
                continue
            r = gpf.profile(organism="hsapiens", query=ids, background=bg,
                            domain_scope="custom", sources=["GO:BP", "GO:MF", "GO:CC", "REAC"],
                            significance_threshold_method="g_SCS")
            r.to_csv(f"{RES}/ora_gprofiler_{k}_query.csv", index=False)
            r.insert(0, "direction", k)
            ora.append(r)
            say(f"  {len(r)} significant terms")
        if ora:
            O = pd.concat(ora, ignore_index=True)
            O["parents"] = O["parents"].astype(str)
            O.to_csv(f"{RES}/ora_gprofiler.csv", index=False)
            for col in ("p_value", "intersection_size", "query_size", "term_size"):
                O[col] = pd.to_numeric(O[col], errors="raise")
            O["gene_ratio"] = O["intersection_size"] / O["query_size"]
            O.nsmallest(20, "p_value")[["direction", "source", "name", "term_size", "intersection_size", "p_value"]].to_csv(f"{RES}/ora_top20.csv", index=False)
            for k, d in O.groupby("direction"):
                d = d.nsmallest(15, "p_value").sort_values("gene_ratio")
                fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 1.8))
                sc = ax.scatter(d["gene_ratio"], d["name"].str.slice(0, 55),
                                s=20 + 12 * d["intersection_size"],
                                c=-np.log10(d["p_value"].clip(lower=np.finfo(float).tiny)), cmap="viridis", edgecolor="k", lw=.4)
                ax.margins(x=.15, y=.10)
                ax.set_xlabel("gene ratio (overlap / recognized query genes)")
                ax.set_title(f"g:Profiler ORA, genes {k} in T2D (exploratory, raw p<{a.p_cut})\nDot area reflects overlap count", fontsize=10)
                fig.colorbar(sc, label="-log10 adjusted p", shrink=.6)
                fig.tight_layout(); fig.savefig(f"{FIG}/dotplot_ora_{k}.png", dpi=200); plt.close(fig)
    except Exception as e:
        raise RuntimeError("g:Profiler did not complete; rerun in a fresh directory") from e

# ---------------------------------------------------------------- checks for interpretation
fad = ["TSPAN33", "NT5E", "TMED6", "PAK5"]   # PAK7 was renamed PAK5 by HGNC; Ensembl 110 uses PAK5
chk = de.loc[de.symbol.isin(fad), ["symbol", "ensembl", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"]].copy()
chk.insert(0, "fadista_name", chk.symbol.replace({"PAK5": "PAK7 (now PAK5)"}))
chk = chk.sort_values("pvalue")
chk.to_csv(f"{RES}/fadista_silenced_genes_check.csv", index=False)
missing = [g for g in fad if g not in set(chk.symbol)]

pos = {s: i + 1 for i, s in enumerate(rk.symbol)}
marker = pd.DataFrame({"symbol": STRESS + BETA,
                       "set": ["stress"] * len(STRESS) + ["beta_cell"] * len(BETA)})
marker["rank"] = marker.symbol.map(pos)
marker["of_total"] = len(rk)
marker = marker.merge(de[["symbol", "log2FoldChange", "pvalue", "padj"]], on="symbol", how="left")
marker.to_csv(f"{RES}/marker_gene_ranks.csv", index=False)

# ---------------------------------------------------------------- compact summary
with io.open(f"{RES}/part4_summary.txt", "w", encoding="utf-8") as f:
    w = lambda s="": f.write(s + "\n")
    w("PART 4 SUMMARY"); w("=" * 70)
    w(f"GSEA permutations: {a.nperm}; seed: 7; set size: 15-500")
    for s in log_lines:
        w(s)
    for lib, d in G.groupby("library"):
        w(); w("=" * 70); w(f"LIBRARY: {lib}")
        fdr = d["FDR q-val"]
        w(f"sets tested {len(d)} | FDR<0.05 {int((fdr<0.05).sum())}")
        for lab, dd in (("TOP UP in T2D (ranked by NES; see FDR)", d[d.NES > 0].nlargest(12, "NES")),
                        ("TOP DOWN in T2D (ranked by NES; see FDR)", d[d.NES < 0].nsmallest(12, "NES"))):
            w(); w(lab + "   (stress_frac = share of leading-edge genes that are Part 2 stress genes)")
            w(f"{'NES':>6} {'FDRq':>9} {'nLead':>6} {'strFr':>6}  Term")
            for _, r in dd.iterrows():
                w(f"{r['NES']:6.2f} {r['FDR q-val']:9.4f} {r['n_lead']:6d} {r['stress_frac']:6.2f}  {str(r['Term'])[:58]}")
    w(); w("=" * 70); w("FADISTA SILENCED GENES (TSPAN33, NT5E, TMED6, PAK5 = former PAK7)")
    w(chk.to_string(index=False))
    if missing:
        w(f"not present in de_results.csv: {', '.join(missing)}")
    w(); w("=" * 70); w("MARKER GENE RANKS  (rank 1 = most up in T2D, "
                        f"{len(rk)} = most down)")
    w(marker.to_string(index=False))
    w(); w("=" * 70)
    w("files written:")
    for d_ in (RES, FIG):
        for fn in sorted(os.listdir(d_)):
            w(f"  {d_}/{fn}")

say("")
say(f"Done. Summary saved to {RES}/part4_summary.txt.")
