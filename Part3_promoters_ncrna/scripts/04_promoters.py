"""Step 4: promoter DNA of the up, down and background genes.

What: cut the promoter window out of the genome for each gene and write one FASTA file per gene list.
Why:  motifs are searched in promoter DNA.
How:  window = 1000 bp before to 200 bp after the TSS (1201 bp), read on the gene strand (minus-strand genes are
      reverse-complemented). The genome FASTA is read one chromosome at a time. Windows that run off a chromosome end or
      have > 10% N are skipped and counted. Only chromosomes 1-22, X, Y, MT are used.
Usage: python3 scripts/04_promoters.py ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz results/genes.tsv results/
"""
import gzip, sys
import pandas as pd

fa, genes_f, out = sys.argv[1:4]
out = out.rstrip("/") + "/"
UP, DOWN = 1000, 200
comp = str.maketrans("ACGTN", "TGCAN")

genes = pd.read_csv(genes_f, sep="\t").set_index("gene_id")
lists = {n: [l.strip() for l in open(f"{out}{n}_genes.txt")] for n in ["up", "down", "background"]}
which = {g: n for n, ids in lists.items() for g in ids}
wanted = genes.loc[list(which)].groupby("chrom")

files = {n: open(f"{out}promoters_{n}.fa", "w") for n in lists}
count = {n: 0 for n in lists}
skipped = 0

def write_chrom(name, seq):
    global skipped
    if name not in wanted.groups:
        return
    for gid, r in wanted.get_group(name).iterrows():
        a, b = (r.tss - UP, r.tss + DOWN) if r.strand == "+" else (r.tss - DOWN, r.tss + UP)
        if a < 1 or b > len(seq):
            skipped += 1; continue
        s = seq[a - 1:b]
        if s.count("N") > 0.1 * len(s):
            skipped += 1; continue
        if r.strand == "-":
            s = s.translate(comp)[::-1]
        n = which[gid]
        files[n].write(f">{gid}|{r.symbol}|{r.chrom}:{a}-{b}({r.strand})\n{s}\n")
        count[n] += 1

name, parts = None, []
with gzip.open(fa, "rt") as f:
    for line in f:
        if line[0] == ">":
            if name:
                write_chrom(name, "".join(parts).upper())
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
write_chrom(name, "".join(parts).upper())
for f in files.values():
    f.close()
print("promoters written:", count, " skipped (edge, N, or other contigs):", skipped)
