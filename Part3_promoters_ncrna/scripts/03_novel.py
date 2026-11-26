"""Step 3: new (not in Ensembl) transcripts in Part 1's merged.gtf.

What: sort the StringTie transcripts that Ensembl does not have into classes, and flag lncRNA-like candidates.
Why:  the outline asks for novel transcripts and non-coding RNAs. gffcompare is not installed, so this is a small
      Python version of the same idea.
How:  a transcript without ref_gene_id is new. If its merged locus (gene_id MSTRG.n) also holds Ensembl transcripts, it is a
      new isoform of a known gene. Otherwise it sits in a new locus: antisense if it overlaps an Ensembl gene on the
      opposite strand, same strand (own class), or none (intergenic). lncRNA-like candidate = new locus, >= 2 exons, spliced length >= 200 nt.
      Coding potential is NOT tested here (needs CPAT or CPC2), so the word is 'candidate'.
Usage: python3 scripts/03_novel.py ../Part1_data_pipeline/merged.gtf.gz results/genes.tsv results/
"""
import gzip, re, sys
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

merged, genes_f, out = sys.argv[1:4]
tx = {}                                    # transcript_id -> dict
with gzip.open(merged, "rt") as f:
    for line in f:
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] not in ("transcript", "exon"):
            continue
        a = c[8]
        tid = re.search(r'transcript_id "([^"]+)"', a).group(1)
        d = tx.setdefault(tid, {"chrom": c[0], "strand": c[6], "locus": re.search(r'gene_id "([^"]+)"', a).group(1),
                                "known": "ref_gene_id" in a, "exons": 0, "len": 0, "start": 10**12, "end": 0})
        if c[2] == "exon":
            d["exons"] += 1
            d["len"] += int(c[4]) - int(c[3]) + 1
        d["start"] = min(d["start"], int(c[3])); d["end"] = max(d["end"], int(c[4]))

locus_known = defaultdict(bool)
for d in tx.values():
    locus_known[d["locus"]] |= d["known"]

g = pd.read_csv(genes_f, sep="\t")
by_chr = {k: v[["start", "end", "strand"]].values for k, v in g.groupby("chrom")}

def overlap(d):
    """'same', 'opposite' or None: overlap with an Ensembl gene on the same or the opposite strand."""
    found = None
    for s, e, st in by_chr.get(d["chrom"], []):
        if s <= d["end"] and e >= d["start"]:
            if st == d["strand"]:
                return "same"
            found = "opposite"
    return found

rows = []
for tid, d in tx.items():
    if d["known"]:
        continue
    if locus_known[d["locus"]]:
        cls = "new isoform of known gene"
    elif d["strand"] == ".":
        cls = "new locus, no strand"
    else:
        ov = overlap(d)
        cls = {"same": "new locus, overlaps a gene on the same strand", "opposite": "new locus, antisense",
               None: "new locus, intergenic"}[ov]
    cand = cls in ("new locus, antisense", "new locus, intergenic") and d["exons"] >= 2 and d["len"] >= 200
    rows.append((tid, d["locus"], d["chrom"], d["strand"], d["exons"], d["len"], cls, cand))
nv = pd.DataFrame(rows, columns=["transcript_id", "locus", "chrom", "strand", "exons", "length", "class", "lncRNA_candidate"])
nv.to_csv(out + "novel_transcripts.csv", index=False)

summ = nv["class"].value_counts()
with open(out + "novel_summary.txt", "w") as f:
    f.write(f"transcripts in merged.gtf: {len(tx)}\nnot in Ensembl: {len(nv)}\n")
    for k, v in summ.items():
        f.write(f"  {k}: {v} ({100*v/len(nv):.1f}%)\n")
    f.write(f"lncRNA candidates: {nv.lncRNA_candidate.sum()} transcripts in {nv[nv.lncRNA_candidate].locus.nunique()} loci\n")
print(open(out + "novel_summary.txt").read())

fig, ax = plt.subplots(figsize=(6, 3))
summ[::-1].plot.barh(ax=ax, color="#555555")
ax.set_xscale("log"); ax.set_xlabel("transcripts (log scale)")
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/novel_classes.png", dpi=200)
