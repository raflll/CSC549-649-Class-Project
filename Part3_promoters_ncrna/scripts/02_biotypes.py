"""Step 2: gene lists and gene types (biotypes) of the changed genes.

What: split the tested genes into up (T2D higher), down (T2D lower) and background; count gene types in each;
      test whether any type is over- or under-represented among changed genes.
Why:  the Part 3 question asks if non-coding RNAs are over-represented among the changed genes.
How:  changed gene = padj < 0.05. If fewer than 30 genes per direction pass, use nominal pvalue < 0.05 instead
      (rule set after seeing that DESeq2 gave 1 gene at padj < 0.05, before any motif or biotype result).
      Background = all other tested genes (not all 62k genes, because the expression filter alone removes many lncRNAs).
      Fisher exact test per type (changed vs background), Benjamini-Hochberg correction.
      Low-expression genes (lncRNAs) have less power to be called changed, so the test is repeated inside five
      expression (baseMean) bins and combined (Mantel-Haenszel). Read the expression-matched columns first.
Usage: python3 scripts/02_biotypes.py results/genes.tsv ../Part2_statistics/results/de_results.csv results/
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

genes_f, de_f, out = sys.argv[1:4]
out = out.rstrip("/") + "/"
genes = pd.read_csv(genes_f, sep="\t").set_index("gene_id")
de = pd.read_csv(de_f).rename(columns={"ensembl": "gene_id"}).set_index("gene_id")[["baseMean", "log2FoldChange", "pvalue", "padj"]]
de = de.join(genes[["symbol", "biotype", "group"]], how="left")
print("tested genes without annotation:", de["group"].isna().sum())
de = de.dropna(subset=["group"])

ok = ((de.padj < 0.05) & (de.log2FoldChange > 0)).sum() >= 30 and ((de.padj < 0.05) & (de.log2FoldChange < 0)).sum() >= 30
col = "padj" if ok else "pvalue"
sig = de[col] < 0.05
de["status"] = "background"
de.loc[sig & (de.log2FoldChange > 0), "status"] = "up"
de.loc[sig & (de.log2FoldChange < 0), "status"] = "down"
print("rule:", col, "< 0.05 ->", de.status.value_counts().to_dict())
for s in ["up", "down"]:
    de[de.status == s].index.to_series().to_csv(out + s + "_genes.txt", index=False, header=False)
de[de.status == "background"].index.to_series().to_csv(out + "background_genes.txt", index=False, header=False)

# biotype table: counts and percent for each status
tab = pd.crosstab(de["group"], de["status"])[["up", "down", "background"]]
pct = (tab / tab.sum() * 100).round(1).add_suffix("_pct")
tab = pd.concat([tab, pct], axis=1)

# Fisher test: changed (up+down) vs background, one row per type
changed = de.status != "background"
de["bin"] = pd.qcut(de.baseMean, 5, labels=False)
rows = []
for g in tab.index:
    a = int(((de.group == g) & changed).sum());  b = int(((de.group != g) & changed).sum())
    c = int(((de.group == g) & ~changed).sum()); d = int(((de.group != g) & ~changed).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    tabs = []                                   # one 2x2 table per expression bin
    for _, x in de.groupby("bin"):
        ing, ch = x.group == g, x.status != "background"
        t2 = [[(ing & ch).sum(), (~ing & ch).sum()], [(ing & ~ch).sum(), (~ing & ~ch).sum()]]
        if min(sum(t2[0]), sum(t2[1])) > 0 and min(t2[0][0] + t2[1][0], t2[0][1] + t2[1][1]) > 0:
            tabs.append(np.array(t2, float) + 0.0)
    st = StratifiedTable(tabs)
    rows.append((g, odds, p, st.oddsratio_pooled, st.test_null_odds().pvalue))
t = pd.DataFrame(rows, columns=["group", "odds_ratio", "p", "odds_ratio_expr_matched", "p_expr_matched"]).set_index("group")
t["q"] = multipletests(t.p, method="fdr_bh")[1]
t["q_expr_matched"] = multipletests(t.p_expr_matched, method="fdr_bh")[1]
tab = tab.join(t).round(4)
tab.to_csv(out + "biotype_table.csv")
print(tab.to_string())

# figure
fig, ax = plt.subplots(figsize=(6, 3.5))
w = 0.27
x = range(len(tab))
for k, (s, colr) in enumerate([("up", "#b2182b"), ("down", "#2166ac"), ("background", "#999999")]):
    ax.bar([i + (k - 1) * w for i in x], tab[s + "_pct"], w, label=f"{s} (n={int(tab[s].sum())})", color=colr)
ax.set_xticks(list(x)); ax.set_xticklabels(tab.index, rotation=20)
ax.set_ylabel("% of genes"); ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/biotype_bar.png", dpi=200)

# figure 2: share of each gene type that was called changed, against the 5% expected by chance
rate = (tab["up"] + tab["down"]) / (tab["up"] + tab["down"] + tab["background"]) * 100
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.bar(rate.index, rate.values, color="#555555")
ax.axhline(5, color="#b2182b", ls="--", lw=1); ax.text(len(rate) - 0.5, 5.3, "5% expected by chance", ha="right", color="#b2182b", fontsize=8)
ax.set_ylabel("% of genes called changed (p < 0.05)")
ax.tick_params(axis="x", rotation=20); ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/biotype_changed_rate.png", dpi=200)
print((rate.round(1)).to_dict())
