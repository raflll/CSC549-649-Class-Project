# Part 3: Promoters and non-coding RNAs

Do the genes that change in T2D share promoter motifs, and are non-coding RNAs over-represented among them?
Methods, results and caveats are in `note.txt`.

## Files

| file | content |
|------|---------|
| `figures/biotype_bar.png` | gene types of up, down and background genes |
| `figures/biotype_changed_rate.png` | share of each gene type called changed |
| `figures/motif_top.png` | ten best motifs per promoter set (uncorrected p value) |
| `figures/novel_classes.png` | classes of transcripts that are not in Ensembl |
| `results/biotype_table.csv` | gene type counts, odds ratios, p and q values |
| `results/motif_results.csv` | all 879 motifs for the up, down and changed sets |
| `results/promoters_up.fa`, `promoters_down.fa`, `promoters_background.fa.gz` | promoter sequences (TSS -1000 to +200) |
| `results/up_genes.txt`, `down_genes.txt`, `background_genes.txt` | changed genes (p < 0.05) and the other tested genes |
| `results/novel_transcripts.csv`, `novel_summary.txt` | transcripts not in Ensembl and lncRNA candidates |
| `scripts/` | `00_reference.sh` downloads the references; `01` to `05` make the files above; `run_all.sh` runs them |

## Run

From this folder (reads `../Part1_data_pipeline/` and `../Part2_statistics/`):

```
bash scripts/00_reference.sh     # Ensembl GTF, GRCh38 genome, JASPAR motifs into ref/ (about 1 GB)
bash scripts/run_all.sh          # about 15 minutes
```

## Software

Python 3 with pandas, numpy, scipy, statsmodels, matplotlib.
