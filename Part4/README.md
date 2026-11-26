# Part 4: Pathways and biology

What biology do the genes that differ between T2D and non-diabetic (ND) islets point to, and
does it fit type 2 diabetes? Methods, results and caveats are in `note.txt`. The written
interpretation and references are in `interpretation_notes.md`. The audit and validation
report, including what was checked and what was corrected, is in `AUDIT.md`.

## Files

| file | content |
|------|---------|
| `results/ranked_list.csv` | all 17,358 symbol-annotated genes ranked by the DESeq2 Wald statistic |
| `results/gsea_MSigDB_Hallmark_2020.csv` | GSEA, 50 Hallmark sets |
| `results/gsea_Reactome_Pathways_2024.csv` | GSEA, 1,257 Reactome sets |
| `results/gsea_GO_Biological_Process_2025.csv` | GSEA, 2,305 GO:BP sets |
| `results/gsea_all.csv` | the three above combined, with leading-edge size and stress fraction |
| `results/ora_gprofiler.csv` | g:Profiler over-representation, genes at raw p < 0.01 (exploratory) |
| `results/ora_top20.csv` | the 20 most significant g:Profiler terms |
| `results/marker_gene_ranks.csv` | rank of each stress and beta-cell marker in the ranked list |
| `results/fadista_silenced_genes_check.csv` | TSPAN33, NT5E, TMED6 and PAK5 in our results |
| `results/de_results_noP7P8.csv` | DESeq2 re-fitted without P7_T2D and P8_T2D (14 donors) |
| `results/donor_influence.csv` | 21 DESeq2 refits; NES and mean log2FC of four fixed gene sets per refit |
| `results/robustness_sensitivity.csv` | selected gene sets, NES and FDR with 16 donors against 14 |
| `results/hallmark_sensitivity.csv` | all 50 Hallmark sets, 16 donors against 14 |
| `results/library_wide_sensitivity.csv` | per-library counts of sets lost, kept, gained and sign-flipped |
| `results/ora_gprofiler_rerun.csv` | repeat query of g:Profiler, 274 up genes, for verification |
| `results/ora_gprofiler_down_rerun.csv` | repeat query of g:Profiler, 32 down genes; header only, 0 significant terms |
| `results/part4_summary.txt` | printed summary of the main run |
| `figures/dotplot_gsea_*.png` | GSEA dot plot per library (NES, FDR, leading-edge size) |
| `figures/gsea_*_top1-3.png` | enrichment plots for the three most significant sets per library |
| `figures/dotplot_ora_up.png` | g:Profiler dot plot, genes higher in T2D |
| `figures/robustness_sensitivity.png` | selected gene sets, 16 donors against 14 |
| `figures/donor_influence.png` | donor influence: leave-one-out and matched two-donor exclusions |
| `figures/fig1_hallmark_report.png` | Hallmark enrichment dot plot for the report and presentation |
| `figures/fig2_donor_report.png` | four fixed gene sets across 21 refits |
| `scripts/part4_pathways_v2.py` | main analysis: ranking, GSEA, g:Profiler, marker checks |
| `scripts/fix_ora_plot.py` | redraws the g:Profiler figure from the saved table |
| `scripts/replot_dotplots.py` | redraws the GSEA dot plots with corrected axis margins |
| `scripts/part4_robustness.py` | 14-donor sensitivity analysis, writes to `results_noP7P8/` |
| `scripts/part4_donor_influence.py` | 21 DESeq2 refits for the donor influence analysis |
| `scripts/make_report_figures.py` | regenerate the two report figures from saved results |
| `requirements.txt` | pinned dependency versions |
| `GENESETS.md` | gene-set source URLs, versions and checksums |
| `AUDIT.md` | what was verified, what was corrected, what is a known limitation |

## Run

From this folder (reads `../Part2_statistics/results/de_results.csv`):

    python3 -m pip install -r requirements.txt
    python3 scripts/part4_pathways_v2.py        # main analysis, 16 donors
    python3 scripts/fix_ora_plot.py             # g:Profiler figure
    python3 scripts/replot_dotplots.py          # GSEA dot plots
    python3 scripts/part4_robustness.py         # 14-donor sensitivity
    python3 scripts/part4_donor_influence.py    # 21 refits, ~20 min, needs pydeseq2
    python3 scripts/make_report_figures.py      # report figures from saved results

Needs internet: gene sets are downloaded from Enrichr and cached in `genesets/` (not committed),
and g:Profiler is queried over its API. A library whose `results/gsea_<name>.csv` already exists
is skipped; use `--force` to redo it.

## Software

Python 3.14.0, gseapy 1.3.2, gprofiler-official 1.0.0, pandas 3.1.1, numpy 2.5.5,
scipy 1.18.1, matplotlib 3.11.2. The audit reproduction additionally used pydeseq2 0.7.0.

Gene sets: MSigDB Hallmark 2020, Reactome Pathways 2024 and GO Biological Process 2025 via
Enrichr; GO:BP, GO:MF, GO:CC and Reactome via the g:Profiler API.

GSEA settings: 10,000 permutations for the main analysis, 1,000 for the sensitivity comparisons
(identical for both rankings being compared), set size 15 to 500, seed 7. GSEApy estimates FDR
the standard GSEA way, from the observed and permutation NES distributions, not by
Benjamini-Hochberg correction. Preranked GSEA permutes gene labels rather than sample labels; the GSEA User Guide
notes this mode is less conservative than phenotype permutation and that FDR < 0.05, not 0.25,
is the appropriate threshold for it, which is what we use.

g:Profiler: custom background of the tested genes, g:SCS correction (g:Profiler's own correction
for GO and Reactome graph structure, not Benjamini-Hochberg).
