#!/usr/bin/env python3
"""Draw the g:Profiler dot plot from results/ora_gprofiler.csv.

Recreate the figure and top-20 table from the archived ORA results without a live query.

Run from the Part4 folder:   python3 scripts/fix_ora_plot.py
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "results/ora_gprofiler.csv"
if not os.path.exists(SRC):
    raise SystemExit(f"{SRC} not found. Run from the Part4 folder.")

O = pd.read_csv(SRC)
for c in ["p_value", "intersection_size", "query_size", "term_size"]:
    O[c] = pd.to_numeric(O[c], errors="coerce")
O = O.dropna(subset=["p_value", "intersection_size", "query_size"])
O["gene_ratio"] = O["intersection_size"] / O["query_size"]
print(f"{len(O)} terms: " + ", ".join(f"{k} {v}" for k, v in O.direction.value_counts().items()))

for k, d in O.groupby("direction"):
    d = d.nsmallest(15, "p_value").sort_values("gene_ratio")
    if d.empty:
        print(f"{k}: nothing to plot")
        continue
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 1.8))
    sc = ax.scatter(d["gene_ratio"], d["name"].astype(str).str.slice(0, 55),
                    s=20 + 12 * d["intersection_size"],
                    c=-np.log10(d["p_value"]), cmap="viridis", edgecolor="k", lw=.4)
    ax.margins(x=.15, y=.10)
    ax.set_xlabel("gene ratio (overlap / recognized query genes)")
    ax.set_title(f"g:Profiler over-representation, genes {k} in T2D\n"
                 f"exploratory: raw p<0.01 list, {int(d['query_size'].iloc[0])} genes; area reflects overlap", fontsize=10)
    fig.colorbar(sc, label="-log10 adjusted p", shrink=.6)
    fig.tight_layout()
    out = f"figures/dotplot_ora_{k}.png"
    fig.savefig(out, dpi=200); plt.close(fig)
    print("wrote", out)

top = O.nsmallest(20, "p_value")[["direction", "source", "name", "term_size",
                                  "intersection_size", "p_value"]]
top.to_csv("results/ora_top20.csv", index=False)
print("\nTOP 20 g:Profiler TERMS")
print(top.to_string(index=False))
