# Part 4: Pathways and biology

What biology do the genes that differ between T2D and non-diabetic (ND) islets point to, and
does it fit type 2 diabetes? Methods, results and caveats are in `note.txt`. The written
interpretation and references are in `interpretation_notes.md`. The audit and validation
report, including what was checked and what was corrected, is in `AUDIT.md`.

## Files

| file | content |
|------|---------|
| `results/ranked_list.csv` | 17,358 unique gene symbols ranked by the DESeq2 Wald statistic |
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
| `figures/donor_influence.png` | donor influence: archived overview of leave-one-out and same-size two-donor exclusions |
| `figures/fig1_hallmark_report.png` | Hallmark enrichment dot plot for the report and presentation |
| `figures/fig2_donor_report.png` | four fixed gene sets across 21 refits |
| `scripts/part4_pathways_v2.py` | main analysis: ranking, GSEA, g:Profiler, marker checks |
| `scripts/fix_ora_plot.py` | redraws the g:Profiler figure and top-20 table from archived results |
| `scripts/replot_dotplots.py` | redraws the GSEA dot plots with corrected axis margins |
| `scripts/part4_robustness.py` | full-library sensitivity analysis, both rankings at 1,000 permutations |
| `scripts/part4_donor_influence.py` | 21 DESeq2 refits for the donor influence analysis |
| `scripts/make_report_figures.py` | regenerate the report figures and selected sensitivity plot from saved results |
| `scripts/pathway_utils.py` | shared ranking, GMT validation and fresh output-directory helpers |
| `scripts/check_results.py` | offline checks of archived inputs and outputs |
| `requirements.txt` | pinned dependency versions |
| `GENESETS.md` | gene-set source URLs, versions and checksums |
| `AUDIT.md` | what was verified, what was corrected, what is a known limitation |

## Run

From this folder, install dependencies and verify the archived outputs:

    python3 -m pip install -r requirements.txt
    python3 scripts/check_results.py

To redraw figures from the archived results (no internet needed):

    python3 scripts/fix_ora_plot.py
    python3 scripts/replot_dotplots.py
    python3 scripts/make_report_figures.py

To repeat the analysis, run these separately:

    python3 scripts/part4_pathways_v2.py
    python3 scripts/part4_robustness.py
    python3 scripts/part4_donor_influence.py

The primary run writes to `reruns/primary/`; the 1,000-permutation sensitivity run to
`reruns/sensitivity/`; and the 21-refit donor analysis to `reruns/donor_influence/`.
The sensitivity script uses the archived 14-donor DE table by default; use `--de14
reruns/donor_influence/de_results_noP7P8.csv` to use a new refit. All scripts require
fresh output directories (use `--output-dir` to choose another). None silently reuses
previous enrichment results. The donor refits can be expensive. The stored result files
are the original analysis; compare reruns before replacing them. PyDESeq2 refits are an
approximate cross-check of the R DESeq2 model, not an identical implementation.

The gene-set files must match the documented version and checksum prefix in `GENESETS.md`.
The script downloads them if absent and stops if a version differs. g:Profiler queries the
live API, so future versions may return different results. The 2,000-permutation fixed-set
NES calculation is distinct from full-library GSEA. Its corrected nominal p-value
calculation uses the same-sign null; the archived nominal p columns are historical and are
not used to support the conclusions.

## Software

Python 3.14.0, gseapy 1.3.2, gprofiler-official 1.0.0, pandas 3.1.1, numpy 2.5.5,
scipy 1.18.1, matplotlib 3.11.2. The original audit reproduction reported pydeseq2 0.7.0 (now pinned).

Gene sets: MSigDB Hallmark 2020, Reactome Pathways 2024 and GO Biological Process 2025 via
Enrichr; GO:BP, GO:MF, GO:CC and Reactome via the g:Profiler API.

GSEA settings: 10,000 permutations for the main analysis, 1,000 for the sensitivity comparisons
(identical for both rankings being compared), set size 15 to 500, seed 7. GSEApy estimates FDR
the standard GSEA way, from the observed and permutation NES distributions, not by
Benjamini-Hochberg correction. Preranked GSEA permutes gene labels rather than sample labels; the GSEA User Guide
notes this mode is less conservative than phenotype permutation and that FDR < 0.05, not 0.25,
is the appropriate threshold for it, which is what we use.

g:Profiler: custom background of the tested genes, g:SCS correction (g:Profiler's own correction
for the annotated gene-set structure, not Benjamini-Hochberg). The GSEA FDR is estimated separately within each library; g:Profiler queries use their own correction. Finite Wald statistics were retained even when Part 2 suppressed a nominal p-value (73 ranked symbols have such flags). The archived ranking is preserved; see `note.txt` for the limitation.
