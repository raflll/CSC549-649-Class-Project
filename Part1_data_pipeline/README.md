# Part 1: Data and pipeline

Islet RNA-seq, type 2 diabetes (T2D) vs non-diabetic (ND).
Methods, results and caveats are in `note.txt`.

## Data

Gene Expression Omnibus (GEO) GSE50244 / Sequence Read Archive (SRA) SRP029262 (Fadista et al., PNAS 2014, doi:10.1073/pnas.1402665111).
16 donors, grouped by glycated hemoglobin (HbA1c): 8 T2D (HbA1c >= 6.5%) and 8 ND (HbA1c < 6.0%). 4 women and 4 men in each group, paired on sex, age and body mass index (BMI).

| group | n | women | age, mean (range) | BMI, mean (range) | HbA1c %, mean (range) |
|-------|---|-------|-------------------|-------------------|-----------------------|
| T2D   | 8 | 4     | 59.0 (44-77)      | 28.3 (24.7-30.9)  | 7.3 (6.7-8.6)         |
| ND    | 8 | 4     | 61.8 (44-74)      | 27.8 (24.7-30.5)  | 5.5 (4.6-5.8)         |

## Files

| file | content |
|------|---------|
| `Table1_donors.csv` | the 16 donors and their runs |
| `SRR_list.txt` | the 16 run accessions, input of the scripts |
| `mapping_table.csv` | per run: read pairs, % mapped and flow cell (batch) |
| `gene_count_matrix.csv` | raw counts, 62,754 genes x 16 runs |
| `merged.gtf.gz` | transcripts assembled from the reads, including new ones |
| `fastqc_summary.txt` | FastQC result (PASS / WARN / FAIL) of all 32 FASTQ files |
| `fastqc/` | FastQC reports of SRR957961 (batch A) and SRR957957 (batch B) |
| `figures/` | workflow, and mapping rate per donor (y-axis from 0% and from 80%) |
| `scripts/` | the scripts |

## Notes

- Batch: the donors were sequenced on two flow cells (column `flowcell` in `mapping_table.csv`). Batch should be in the DESeq2 model (`design = ~ batch + group`).
- Column order: in `gene_count_matrix.csv` the runs are in alphabetical order, not in the order of `Table1_donors.csv`.
- Gene names: the first column is `ENSG...|symbol`, or only the Ensembl ID when there is no symbol. Some symbols (for example Y_RNA, U6) belong to more than one Ensembl ID.
- New transcripts: in `merged.gtf`, transcripts named `MSTRG...` are not in Ensembl (31,934 of 284,533).

## Run

```
./00_reference.sh
nohup ./01_download.sh > download.log 2>&1 &
nohup ./run_all.sh > run_all.log 2>&1 &
python 05_mapping_table.py Table1_donors.csv mapping_table.csv
python 06_figure_mapping.py mapping_table.csv fig_mapping_from0.png 0
python 06_figure_mapping.py mapping_table.csv fig_mapping_from80.png 80
python 07_figure_workflow.py fig_workflow.png
```

`run_all.sh` runs `02_fastqc.sh`, `03_hisat2.sh` and `04_stringtie.sh` after the download.
All commands were run from the project folder, next to `fastq/`, `bam/` and `ref/`.

## Software

SRA Toolkit 3.4.1, FastQC 0.13.0, HISAT2 2.2.3, SAMtools 1.25, StringTie 2.1.7, Python 3 with matplotlib.
Reference: Ensembl GRCh38, HISAT2 index `grch38_genome`, Ensembl annotation release 110.
All steps ran in Windows Subsystem for Linux (WSL, Ubuntu) on one laptop, in a conda environment.
