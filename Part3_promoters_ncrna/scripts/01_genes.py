"""Step 1: gene table from the Ensembl release 110 GTF.

What: one row per gene with Ensembl ID, symbol, biotype (gene type), chromosome, strand, start, end, TSS.
Why:  biotypes tell us which genes are non-coding. The TSS (transcription start site) is where the promoter window sits.
How:  read the 'gene' lines of the GTF. TSS = start on the + strand, end on the - strand.
Usage: python3 scripts/01_genes.py ref/Homo_sapiens.GRCh38.110.gtf.gz results/genes.tsv
"""
import gzip, re, sys
import pandas as pd

gtf, out = sys.argv[1], sys.argv[2]

# Biotype groups used in the report.
GROUPS = {"protein_coding": "protein_coding", "lncRNA": "lncRNA", "miRNA": "small_RNA", "snRNA": "small_RNA",
          "snoRNA": "small_RNA", "misc_RNA": "small_RNA", "rRNA": "small_RNA", "scaRNA": "small_RNA",
          "Mt_rRNA": "small_RNA", "Mt_tRNA": "small_RNA", "ribozyme": "small_RNA", "sRNA": "small_RNA", "scRNA": "small_RNA"}

def group(bt):
    if bt in GROUPS:
        return GROUPS[bt]
    if "pseudogene" in bt:
        return "pseudogene"
    return "other"

rows = []
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "gene":
            continue
        a = c[8]
        gid = re.search(r'gene_id "([^"]+)"', a).group(1)
        name = re.search(r'gene_name "([^"]+)"', a)
        bt = re.search(r'gene_biotype "([^"]+)"', a).group(1)
        start, end, strand = int(c[3]), int(c[4]), c[6]
        rows.append((gid, name.group(1) if name else "", bt, group(bt), c[0], strand, start, end,
                     start if strand == "+" else end))

df = pd.DataFrame(rows, columns=["gene_id", "symbol", "biotype", "group", "chrom", "strand", "start", "end", "tss"])
df.to_csv(out, sep="\t", index=False)
print(len(df), "genes written to", out)
print(df["group"].value_counts().to_string())
