# CSC549/649 Class Project: islet RNA-seq, type 2 diabetes vs non-diabetic

Pancreatic islet RNA-seq from GEO GSE50244 (Fadista et al., PNAS 2014): 8 T2D and 8 matched non-diabetic (ND) donors.

| folder | part | content |
|--------|------|---------|
| [`Part1_data_pipeline/`](Part1_data_pipeline) | 1. Data and pipeline | download, FastQC, HISAT2 mapping, StringTie counts, `gene_count_matrix.csv` |
| [`Part2_statistics/`](Part2_statistics) | 2. Statistics | islet QC, PCA, DESeq2 differential expression, volcano plot |
| [`Part3/`](Part3) | 3. | not started |
| [`Part4/`](Part4) | 4. | not started |

Each part has a `README.md` (files and how to run) and a `note.txt` (what, why, how, results, caveats).

Raw FASTQ, BAM and reference files (~323 GB) are not in the repository; see `.gitignore`.
