#!/usr/bin/env python3
"""Redraw the GSEA dot plots from results/gsea_all.csv with room for the largest dots.

The first version clipped the biggest circles against the right-hand axis.
Nothing is recomputed; this only redraws the figures.

Run from the Part4 folder:   python3 scripts/replot_dotplots.py
"""
import os
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "results/gsea_all.csv"
if not os.path.exists(SRC):
    raise SystemExit(f"{SRC} not found. Run this from the Part4 folder.")
os.makedirs("figures", exist_ok=True)

G = pd.read_csv(SRC)
for c in ["NES", "FDR q-val"]:
    G[c] = pd.to_numeric(G[c], errors="coerce")
if "n_lead" not in G.columns:
    G["n_lead"] = G["Lead_genes"].astype(str).str.count(";") + 1
G = G.dropna(subset=["NES", "FDR q-val"])

for lib, d in G.groupby("library"):
    sig = d[d["FDR q-val"] < 0.25]
    note = ""
    if sig.empty:
        sig = d.reindex(d["NES"].abs().sort_values(ascending=False).index).head(20)
        note = " (no set reached FDR 0.25; showing the largest NES)"
    sig = sig.copy()
    sig["label"] = sig["Term"].astype(str).str.replace(r"\s*\(GO:\d+\)", "", regex=True).str.slice(0, 55)
    sel = pd.concat([sig.nlargest(10, "NES"), sig.nsmallest(10, "NES")]
                    ).drop_duplicates("label").sort_values("NES")
    if sel.empty:
        continue
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(sel) + 1.8))
    sc = ax.scatter(sel["NES"], sel["label"], s=20 + 12 * sel["n_lead"],
                    c=sel["FDR q-val"], cmap="viridis_r", edgecolor="k", lw=.4)
    ax.axvline(0, c="grey", lw=.8)
    lo, hi = sel["NES"].min(), sel["NES"].max()
    pad = max(0.35, 0.12 * (hi - lo))          # keep the largest circles off the spines
    ax.set_xlim(lo - pad, hi + pad)
    ax.margins(y=0.04)
    ax.set_xlabel("NES (positive = higher in T2D)")
    ax.set_title(f"GSEA {lib}: T2D vs ND{note}", fontsize=10)
    fig.colorbar(sc, label="FDR q-value", shrink=.6)
    fig.tight_layout()
    out = f"figures/dotplot_gsea_{lib}.png"
    fig.savefig(out, dpi=200); plt.close(fig)
    print("wrote", out, f"({len(sel)} sets)")
