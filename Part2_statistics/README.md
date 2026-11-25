# Part 2: Statistics

Which genes differ between T2D and non-diabetic (ND) islets?
Methods, results and caveats are in `note.txt`.

## Files

| file | content |
|------|---------|
| `figures/pca.png` | PCA of the 16 donors (VST, 500 most variable genes) |
| `figures/pca_batch_sex_removed.png` | the same after removing batch and sex (display only) |
| `figures/volcano.png` | volcano plot, T2D vs ND |
| `figures/ma_plot.png` | MA plot, apeglm-shrunken fold changes |
| `figures/qc_purity.png` | islet purity: endocrine vs acinar vs ductal reads |
| `figures/qc_sex.png` | sex check: XIST vs chrY genes |
| `results/de_results.csv` | DESeq2 results for all 19,072 tested genes, sorted by p value |
| `results/up_genes.txt`, `down_genes.txt` | genes higher / lower in T2D at FDR < 0.05 (also as `.csv`) |
| `results/expressed_genes.txt` | all tested genes (background for enrichment) |
| `results/qc_islet_purity_sex.csv` | per-donor QC |
| `results/pca_associations.csv` | PC1-5 vs group, batch, sex and purity |
| `results/sensitivity_models.csv`, `sensitivity_top_genes.csv` | other models and leaving out two donors |
| `results/known_t2d_genes.csv` | genes from earlier T2D islet studies |
| `results/summary.txt` | printed summary and R session info |
| `scripts/part2_statistics.R` | the script that makes all of the above |

## Run

From this folder (reads `../Part1_data_pipeline/`):

```
Rscript scripts/part2_statistics.R
```

## Software

R 4.3.2, DESeq2 1.42.1, apeglm 1.24.0, limma 3.58.1, ggplot2 3.5.2, ggrepel 0.9.6.
