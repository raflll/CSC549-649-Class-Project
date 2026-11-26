#!/usr/bin/env python3
"""Compare complete GSEA libraries using the same settings in both cohorts.

Uses the archived full-cohort and 14-donor DESeq2 tables. The 14-donor fit is an
input, not re-estimated here. The separate donor-influence script refits counts.
Run from Part4; outputs go to a fresh directory rather than overwriting results.
"""
import argparse
import numpy as np
import pandas as pd
from pathway_utils import LIBRARIES, get_gmt, rank_genes, new_dir

LABELS = {"MSigDB_Hallmark_2020": "MSigDB Hallmark",
          "Reactome_Pathways_2024": "Reactome", "GO_Biological_Process_2025": "GO:BP"}
WATCH = [('raised in T2D', 'GO:BP', 'Cellular Response to Hypoxia (GO:0071456)'), ('raised in T2D', 'GO:BP', 'Cellular Response to Decreased Oxygen Levels (GO:0036294)'), ('raised in T2D', 'GO:BP', 'Extracellular Matrix Organization (GO:0030198)'), ('raised in T2D', 'GO:BP', 'Collagen Fibril Organization (GO:0030199)'), ('raised in T2D', 'GO:BP', 'Inflammatory Response (GO:0006954)'), ('raised in T2D', 'Reactome', 'Collagen Formation'), ('raised in T2D', 'Reactome', 'Extracellular Matrix Organization'), ('raised in T2D', 'Reactome', 'Collagen Degradation'), ('reduced in T2D', 'Reactome', 'tRNA Processing in the Mitochondrion'), ('reduced in T2D', 'Reactome', 'rRNA Processing in the Mitochondrion'), ('reduced in T2D', 'Reactome', 'Mitochondrial RNA Degradation'), ('exploratory', 'Reactome', 'Regulation of Insulin Secretion'), ('exploratory', 'GO:BP', 'Regulation of Insulin Secretion (GO:0050796)'), ('exploratory', 'GO:BP', 'Insulin Secretion (GO:0030073)'), ('exploratory', 'GO:BP', 'Mitochondrial Translation (GO:0032543)'), ('reduced in T2D', 'MSigDB Hallmark', 'Pancreas Beta Cells')]

def compare(main, reduced):
    cols = ["Term", "NES", "NOM p-val", "FDR q-val"]
    a, b = main[cols].copy(), reduced[cols].copy()
    for frame in (a, b):
        for col in cols[1:]:
            frame[col] = pd.to_numeric(frame[col], errors="raise")
    m = a.merge(b, on="Term", how="outer", suffixes=("_16", "_14"), validate="one_to_one", indicator=True)
    if not m["_merge"].eq("both").all():
        raise ValueError("Gene-set membership differs between runs; inspect the rankings and GMT before comparing.")
    m = m.drop(columns="_merge")
    x, y = m["FDR q-val_16"] < .05, m["FDR q-val_14"] < .05
    m["verdict"] = np.select([x & y, x & ~y, ~x & y], ["stays sig", "loses sig", "gains sig"], default="ns both")
    stats = dict(sets=len(m), sig16=int(x.sum()), sig14=int(y.sum()), lose=int((x & ~y).sum()),
                 stay=int((x & y).sum()), gain=int((~x & y).sum()),
                 sign_flips=int((np.sign(m.NES_16) != np.sign(m.NES_14)).sum()),
                 median_abs_NES_16=round(m.NES_16.abs().median(),3), median_abs_NES_14=round(m.NES_14.abs().median(),3))
    return m.sort_values(["FDR q-val_16", "NES_16"], ascending=[True, False]), stats

def main():
    import gseapy as gp
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="reruns/sensitivity")
    ap.add_argument("--de14", default="results/de_results_noP7P8.csv")
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    out = new_dir(args.output_dir)
    r16 = rank_genes(pd.read_csv("../Part2_statistics/results/de_results.csv"))
    r14 = rank_genes(pd.read_csv(args.de14))
    if set(r16.symbol) != set(r14.symbol):
        raise ValueError("Ranked gene universes differ; do not compare silently.")
    comparisons, stats = {}, []
    for name in LIBRARIES:
        gmt = get_gmt(name)
        runs = []
        for n, rank in [(16,r16),(14,r14)]:
            pre = gp.prerank(rnk=rank, gene_sets=gmt, min_size=15, max_size=500,
                             permutation_num=1000, seed=7, threads=args.threads,
                             outdir=None, verbose=False)
            r = pre.res2d.copy()
            r.to_csv(out / f"gsea_{name}_{n}donors.csv", index=False)
            runs.append(r)
        comp, counts = compare(*runs)
        label = LABELS[name]
        comparisons[label] = comp
        stats.append(dict(library=label, **counts))
        comp.to_csv(out / f"comparison_{name}.csv", index=False)
        print(label, counts, flush=True)
    comparisons["MSigDB Hallmark"].to_csv(out / "hallmark_sensitivity.csv", index=False)
    pd.DataFrame(stats).to_csv(out / "library_wide_sensitivity.csv", index=False)
    long_verdict = {"stays sig":"stays significant", "loses sig":"loses significance",
                    "gains sig":"gains significance", "ns both":"not significant in either"}
    rows = []
    for theme, lib, term in WATCH:
        r = comparisons[lib].set_index("Term").loc[term]
        rows.append(dict(theme=theme, library=lib, gene_set=term,
                         NES_16donors=round(r.NES_16,3), FDR_16donors=round(r["FDR q-val_16"],4),
                         NES_14donors=round(r.NES_14,3), FDR_14donors=round(r["FDR q-val_14"],4),
                         verdict=long_verdict[r.verdict]))
    pd.DataFrame(rows).to_csv(out / "robustness_sensitivity.csv", index=False)
    print(f"Saved matched 1,000-permutation sensitivity analyses to {out}")

if __name__ == "__main__":
    main()
