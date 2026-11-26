#!/usr/bin/env python3
"""Part 4 robustness check: does the pathway signal depend on P7_T2D and P8_T2D?

Part 2 found that two T2D donors carry a strong hypoxia/stress signature and that
dropping them removes all six significant genes. This script repeats the Part 4
pathway analysis on a DESeq2 ranking computed without those two donors, and prints
the two runs side by side.

Inputs : results/de_results_noP7P8.csv  (DESeq2, ~ batch + sex + group, 14 donors)
         results/gsea_all.csv           (the main 16-donor run, for comparison)
Outputs: results_noP7P8/ and figures_noP7P8/   (the main run is not touched)

Run from the Part4 folder:   python3 scripts/part4_robustness.py
"""
import os, io, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
from scipy.stats import mannwhitneyu

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"
LIBS = ["MSigDB_Hallmark_2020", "Reactome_Pathways_2024", "GO_Biological_Process_2025"]
RES, FIG, GMT = "results_noP7P8", "figures_noP7P8", "genesets"
for d in (RES, FIG, GMT): os.makedirs(d, exist_ok=True)

# a priori beta-cell identity set, from the literature. Deliberately NOT the genes
# Part 2 flagged as down, which would make the test circular.
BETA = ["INS","IAPP","PDX1","NKX6-1","MAFA","NEUROD1","PAX6","NKX2-2","FOXA2","SLC2A2",
        "GCK","ABCC8","KCNJ11","PCSK1","PCSK2","CHGA","CHGB","SCG2","G6PC2","UCN3",
        "SLC30A8","ISL1","PAX4","RFX6","GLIS3","INSM1","SYT4","ERO1B","PPP1R1A","DLK1"]
STRESS = ["AK4","LOX","IGFBP3","LUCAT1","ANGPTL4","SLC16A3","GDF15","THBS2","PLOD2",
          "G0S2","TMEM158","TREM1","HIC1"]

def rank_of(df):
    d = df.dropna(subset=["symbol","stat"]).copy()
    d["a"] = d["stat"].abs()
    d = d.sort_values("a", ascending=False).drop_duplicates("symbol")
    return d[["symbol","stat"]].sort_values("stat", ascending=False).reset_index(drop=True)

def set_test(rk, genes, name):
    pos = {s: i+1 for i, s in enumerate(rk.symbol)}
    pres = [g for g in genes if g in pos]
    inset = rk[rk.symbol.isin(pres)].stat.values
    rest  = rk[~rk.symbol.isin(pres)].stat.values
    _, p = mannwhitneyu(inset, rest, alternative="less")
    med = int(np.median([pos[g] for g in pres]))
    return dict(set=name, n=len(pres), median_rank=med, of=len(rk),
                pct_down=round(100*med/len(rk), 1), p_shifted_down=p)

# ---------------------------------------------------------------- rank-level test
de_d = pd.read_csv("results/de_results_noP7P8.csv")
de_m = pd.read_csv("../Part2_statistics/results/de_results.csv")
rk_d, rk_m = rank_of(de_d), rank_of(de_m)
rows = []
for tag, rk in [("main 16 donors", rk_m), ("without P7/P8 (14)", rk_d)]:
    for nm, gs in [("beta-cell identity", BETA), ("Part 2 stress genes", STRESS)]:
        r = set_test(rk, gs, nm); r["run"] = tag; rows.append(r)
T = pd.DataFrame(rows)[["run","set","n","median_rank","of","pct_down","p_shifted_down"]]
T.to_csv(f"{RES}/rank_test.csv", index=False)
print("\nRANK-LEVEL TEST (Mann-Whitney, set shifted DOWN against all other genes)")
print(T.to_string(index=False))

# ---------------------------------------------------------------- GSEA
def get_gmt(name):
    p = os.path.join(GMT, name + ".gmt")
    if os.path.exists(p) and os.path.getsize(p) > 1000: return p
    import requests
    for i in (1, 2, 3):
        try:
            r = requests.get(ENRICHR.format(name), timeout=180); r.raise_for_status()
            t = r.content.decode("utf-8", errors="replace")
            if len(t) < 1000: raise ValueError("empty library")
            io.open(p, "w", encoding="utf-8", newline="\n").write(t)
            return p
        except Exception as e:
            print(f"  download {i} failed for {name}: {e}"); time.sleep(5)
    return None

allr = []
for lib in LIBS:
    out = f"{RES}/gsea_{lib}.csv"
    if os.path.exists(out):
        allr.append(pd.read_csv(out)); print(f"{lib}: already done"); continue
    g = get_gmt(lib)
    if g is None: print(f"{lib}: SKIPPED"); continue
    pre = gp.prerank(rnk=rk_d, gene_sets=g, min_size=15, max_size=500,
                     permutation_num=10000, seed=7, threads=4, outdir=None, verbose=False)
    r = pre.res2d.copy(); r.insert(0, "library", lib); r.to_csv(out, index=False)
    allr.append(r)
    f = pd.to_numeric(r["FDR q-val"], errors="coerce")
    print(f"{lib}: {len(r)} sets, FDR<0.05: {int((f<0.05).sum())}")

if allr:
    G = pd.concat(allr, ignore_index=True)
    for c in ["NES","FDR q-val"]: G[c] = pd.to_numeric(G[c], errors="coerce")
    G.sort_values("FDR q-val").to_csv(f"{RES}/gsea_all.csv", index=False)

    # side-by-side with the main run
    if os.path.exists("results/gsea_all.csv"):
        M = pd.read_csv("results/gsea_all.csv")
        for c in ["NES","FDR q-val"]: M[c] = pd.to_numeric(M[c], errors="coerce")
        watch = ["Hypoxia","TNF-alpha Signaling via NF-kB","Epithelial Mesenchymal Transition",
                 "Glycolysis","Inflammatory Response","Pancreas Beta Cells","Protein Secretion",
                 "Collagen Formation","tRNA Processing in the Mitochondrion"]
        cmp = []
        for t in watch:
            a = M[M.Term == t]; b = G[G.Term == t]
            if len(a) and len(b):
                cmp.append(dict(Term=t, NES_main=round(a.NES.iloc[0],2), FDR_main=round(a["FDR q-val"].iloc[0],4),
                                NES_noP7P8=round(b.NES.iloc[0],2), FDR_noP7P8=round(b["FDR q-val"].iloc[0],4)))
        if cmp:
            C = pd.DataFrame(cmp); C.to_csv(f"{RES}/comparison_key_sets.csv", index=False)
            print("\nKEY SETS, MAIN RUN vs WITHOUT P7/P8")
            print(C.to_string(index=False))

    # dot plot of the new run
    for lib, d in G.groupby("library"):
        sig = d[d["FDR q-val"] < 0.25]
        note = ""
        if sig.empty:
            sig = d.reindex(d.NES.abs().sort_values(ascending=False).index).head(20)
            note = " (no set reached FDR 0.25)"
        s = sig.copy()
        s["label"] = s.Term.astype(str).str.replace(r"\s*\(GO:\d+\)", "", regex=True).str.slice(0,55)
        s["n_lead"] = s["Lead_genes"].astype(str).str.count(";") + 1
        sel = pd.concat([s.nlargest(10,"NES"), s.nsmallest(10,"NES")]).drop_duplicates("label").sort_values("NES")
        if sel.empty: continue
        fig, ax = plt.subplots(figsize=(8.5, 0.34*len(sel)+1.8))
        sc = ax.scatter(sel.NES, sel.label, s=20+12*sel.n_lead, c=sel["FDR q-val"],
                        cmap="viridis_r", edgecolor="k", lw=.4)
        ax.axvline(0, c="grey", lw=.8)
        lo, hi = sel.NES.min(), sel.NES.max(); pad = max(0.35, 0.12*(hi-lo))
        ax.set_xlim(lo-pad, hi+pad)
        ax.set_xlabel("NES (positive = higher in T2D)")
        ax.set_title(f"GSEA {lib}, without P7_T2D and P8_T2D{note}", fontsize=10)
        fig.colorbar(sc, label="FDR q-value", shrink=.6); fig.tight_layout()
        fig.savefig(f"{FIG}/dotplot_gsea_{lib}.png", dpi=200); plt.close(fig)
        print("wrote", f"{FIG}/dotplot_gsea_{lib}.png")
print("\ndone")
