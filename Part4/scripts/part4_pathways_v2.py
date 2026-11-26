#!/usr/bin/env python3
"""Part 4: pathways and biology.

Main result : GSEA on the full ranked gene list (DESeq2 Wald stat).
Supporting  : g:Profiler over-representation on genes with raw p < 0.01 (exploratory:
              only 6 genes pass FDR 0.05, so an FDR-based list would be useless).

Run from the Part4 folder:   python3 scripts/part4_pathways.py

Resumable: a library whose results/gsea_<name>.csv already exists is skipped.
Use --force to redo everything.
Offline test: --gmt toy.gmt --skip-ora
"""
import argparse, io, os, sys, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp

ENRICHR_GMT = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"

# hypoxia / stress genes that Part 2 reports as driven by P7_T2D and P8_T2D.
STRESS = ["AK4", "LOX", "IGFBP3", "LUCAT1", "ANGPTL4", "SLC16A3", "GDF15",
          "THBS2", "PLOD2", "G0S2", "TMEM158", "TREM1", "HIC1"]
# beta-cell / known T2D islet genes that Part 2 reports as trending down.
BETA = ["SLC2A2", "RASGRP1", "PPP1R1A", "GLRA1", "CAPN13", "CHL1", "PDX1",
        "IAPP", "ABCC8", "SLC30A8", "INS", "GCG", "TCF7L2", "MAFA", "KCNJ11"]

ap = argparse.ArgumentParser()
ap.add_argument("--de", default="../Part2_statistics/results/de_results.csv")
ap.add_argument("--gmt", nargs="*", help="local .gmt files instead of Enrichr libraries")
ap.add_argument("--skip-ora", action="store_true")
ap.add_argument("--force", action="store_true", help="redo libraries already finished")
ap.add_argument("--nperm", type=int, default=10000)
ap.add_argument("--p-cut", type=float, default=0.01)
a = ap.parse_args()

RES, FIG, GMT = "results", "figures", "genesets"
for d in (RES, FIG, GMT):
    os.makedirs(d, exist_ok=True)

log_lines = []
def say(*x):
    s = " ".join(str(i) for i in x)
    print(s, flush=True)
    log_lines.append(s)

# ---------------------------------------------------------------- ranked list
de = pd.read_csv(a.de)
rk = de.dropna(subset=["symbol", "stat"]).copy()
rk["abs"] = rk["stat"].abs()
rk = rk.sort_values("abs", ascending=False).drop_duplicates("symbol")
rk = rk[["symbol", "stat"]].sort_values("stat", ascending=False).reset_index(drop=True)
rk.to_csv(f"{RES}/ranked_list.csv", index=False)
say(f"ranked genes: {len(rk)} of {len(de)} tested "
    f"({de.symbol.isna().sum()} have no symbol, {de.symbol.dropna().duplicated().sum()} duplicate symbols)")

# ---------------------------------------------------------------- gene sets
def get_gmt(name):
    """Return a local .gmt path. Downloads from Enrichr once and caches as UTF-8.

    Doing the download here instead of inside gseapy avoids the Windows cp1252
    crash on Greek letters in pathway names.
    """
    if name.endswith(".gmt"):
        return name
    path = os.path.join(GMT, name + ".gmt")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    import requests
    for attempt in (1, 2, 3):
        try:
            r = requests.get(ENRICHR_GMT.format(name), timeout=180)
            r.raise_for_status()
            text = r.content.decode("utf-8", errors="replace")
            if len(text) < 1000 or "\t" not in text:
                raise ValueError("library came back empty, check the name")
            with io.open(path, "w", encoding="utf-8", newline="\n") as f:
                f.write(text)
            say(f"  downloaded {name}: {text.count(chr(10))} sets")
            return path
        except Exception as e:
            say(f"  download attempt {attempt} for {name} failed: {type(e).__name__}: {e}")
            if attempt < 3:
                time.sleep(5)
    return None

libs = a.gmt or ["MSigDB_Hallmark_2020", "Reactome_Pathways_2024", "GO_Biological_Process_2025"]

# ---------------------------------------------------------------- GSEA
def lead_genes(s):
    return [g for g in str(s).split(";") if g]

allres, failed = [], []
for lib in libs:
    name = os.path.basename(lib).replace(".gmt", "")
    out_csv = f"{RES}/gsea_{name}.csv"
    if os.path.exists(out_csv) and not a.force:
        r = pd.read_csv(out_csv)
        allres.append(r)
        say(f"{name}: already done, {len(r)} sets (use --force to redo)")
        continue
    path = get_gmt(lib)
    if path is None:
        failed.append(name)
        say(f"{name}: SKIPPED, could not get the gene sets")
        continue
    try:
        pre = gp.prerank(rnk=rk, gene_sets=path, min_size=15, max_size=500,
                         permutation_num=a.nperm, seed=7, threads=4,
                         outdir=None, verbose=False)
    except Exception as e:
        failed.append(name)
        say(f"{name}: FAILED during GSEA: {type(e).__name__}: {e}")
        continue
    r = pre.res2d.copy()
    r.insert(0, "library", name)
    r.to_csv(out_csv, index=False)
    allres.append(r)
    fdr = r["FDR q-val"].astype(float)
    say(f"{name}: {len(r)} sets tested, FDR<0.05: {int((fdr<0.05).sum())}, FDR<0.25: {int((fdr<0.25).sum())}")
    for i, t in enumerate(r.assign(f=fdr).sort_values("f").head(3)["Term"].tolist()):
        try:
            ax = pre.plot(terms=t, show_ranking=True, figsize=(5, 6))
            ax.savefig(f"{FIG}/gsea_{name}_top{i+1}.png", dpi=200, bbox_inches="tight")
            plt.close("all")
        except Exception as e:
            say(f"  enrichment plot failed for {t}: {type(e).__name__}")

if not allres:
    say("No GSEA results at all. Nothing further to do.")
    with io.open(f"{RES}/part4_summary.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines) + "\n")
    sys.exit(1)

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
    d = pd.concat([d.nlargest(n, "NES"), d.nsmallest(n, "NES")]).drop_duplicates("label").sort_values("NES")
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 1.8))
    sc = ax.scatter(d["NES"], d["label"], s=20 + 12 * d["n_lead"],
                    c=d["FDR q-val"], cmap="viridis_r", edgecolor="k", lw=.4)
    ax.axvline(0, c="grey", lw=.8)
    ax.set_xlabel("NES (positive = higher in T2D)")
    ax.set_title(title, fontsize=10)
    fig.colorbar(sc, label="FDR q-value", shrink=.6)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
    return True

for lib, d in G.groupby("library"):
    sig = d[d["FDR q-val"] < 0.25]
    note = ""
    if sig.empty:  # fall back so slide 10 always gets a figure
        sig = d.reindex(d["NES"].abs().sort_values(ascending=False).index).head(20)
        note = " (no set reached FDR 0.25; showing the largest NES)"
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
            if len(ids) < 5:
                say(f"  too few genes, skipped")
                continue
            r = gpf.profile(organism="hsapiens", query=ids, background=bg,
                            domain_scope="custom", sources=["GO:BP", "GO:MF", "GO:CC", "REAC"],
                            significance_threshold_method="g_SCS")
            r.insert(0, "direction", k)
            ora.append(r)
            say(f"  {len(r)} significant terms")
        if ora:
            O = pd.concat(ora, ignore_index=True)
            O["parents"] = O["parents"].astype(str)
            O.to_csv(f"{RES}/ora_gprofiler.csv", index=False)
            O["gene_ratio"] = O["intersection_size"] / O["query_size"]
            for k, d in O.groupby("direction"):
                d = d.nsmallest(15, "p_value").sort_values("gene_ratio")
                fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 1.8))
                sc = ax.scatter(d["gene_ratio"], d["name"].str.slice(0, 55),
                                s=20 + 12 * d["intersection_size"],
                                c=-np.log10(d["p_value"]), cmap="viridis", edgecolor="k", lw=.4)
                ax.set_xlabel("gene ratio")
                ax.set_title(f"g:Profiler ORA, genes {k} in T2D (exploratory, raw p<{a.p_cut})", fontsize=10)
                fig.colorbar(sc, label="-log10 adjusted p", shrink=.6)
                fig.tight_layout(); fig.savefig(f"{FIG}/dotplot_ora_{k}.png", dpi=200); plt.close(fig)
    except Exception as e:
        say(f"g:Profiler step FAILED: {type(e).__name__}: {e}")

# ---------------------------------------------------------------- checks for interpretation
fad = ["TSPAN33", "NT5E", "TMED6", "PAK5"]   # PAK7 was renamed PAK5 by HGNC; Ensembl 110 uses PAK5
chk = de[de.symbol.isin(fad)][["symbol", "baseMean", "log2FoldChange", "pvalue", "padj"]]
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
    for s in log_lines:
        w(s)
    if failed:
        w(); w("LIBRARIES THAT FAILED: " + ", ".join(failed))
    for lib, d in G.groupby("library"):
        w(); w("=" * 70); w(f"LIBRARY: {lib}")
        fdr = d["FDR q-val"]
        w(f"sets tested {len(d)} | FDR<0.05 {int((fdr<0.05).sum())} | FDR<0.25 {int((fdr<0.25).sum())}")
        for lab, dd in (("TOP UP in T2D", d.nlargest(12, "NES")),
                        ("TOP DOWN in T2D", d.nsmallest(12, "NES"))):
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
say("done. Summary saved to results/part4_summary.txt.")
