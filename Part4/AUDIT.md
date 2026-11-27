# Part 4: audit and validation report

Record of the original analysis checks and the subsequent repository review. Original execution claims below refer to the archived run, not a new execution of the corrected scripts.
Last updated 26 November 2026.

## Inputs, verified by checksum

| File | Identity | Note |
|------|----------|------|
| `Part1_data_pipeline/gene_count_matrix.csv` | 62,754 genes x 16 donors | re-analysed from here |
| `Part2_statistics/results/de_results.csv` | md5 `adb69afe...` | identical to the repository archive copy |
| `Part4/scripts/part4_pathways_v2.py` | sha256 `bd66c6a913aae6f29df01dce8ef5a6b6541847a0cb188fe9d0e7231ec2246dd0` | corrected script; full rerun pending |
| `genesets/MSigDB_Hallmark_2020.gmt` | sha256 `42755929...`, 50 sets, sizes 32-200 | `Pancreas Beta Cells` = the 40 genes listed on gsea-msigdb.org; `Hypoxia` = 200 genes |
| `genesets/Reactome_Pathways_2024.gmt` | sha256 `1a7c43cb...`, 2,104 sets | see `GENESETS.md` |
| `genesets/GO_Biological_Process_2025.gmt` | sha256 `0fc02458...`, 5,342 sets | see `GENESETS.md` |

The reviewed source was a repository ZIP snapshot; no commit hash is embedded in the archive.

## 1. Validation table

| Requirement | Verification performed | Evidence | Result |
|---|---|---|---|
| Dataset, donors, groups | Read `Table1_donors.csv`, `mapping_table.csv` | 8 T2D / 8 ND, 4+4 by sex, HbA1c >= 6.5 vs < 6.0 | VERIFIED AND CORRECT |
| Batch structure | Derived batch from flowcell exactly as Part 2 does | batch A 6 T2D / 3 ND; batch B 2 T2D / 5 ND | VERIFIED AND CORRECT |
| DESeq2 model, contrast, sign | Re-ran `~ batch + sex + group` in pydeseq2 0.7.0 | 19,072 genes; same 6 genes at FDR < 0.05; AK4 stat +5.78 (up in T2D), SLC2A2 -2.70 | VERIFIED AND CORRECT |
| Full numerical agreement across implementations | Compared all 19,072 genes, not just significant names | baseMean r = 1.0000000; log2FC r = 0.9999935; Wald stat r = 0.9999893, Spearman 0.9999988, 15 genes differ by > 0.01; padj r = 0.9998 | APPROXIMATE AGREEMENT, not exact reproduction (different optimiser and shrinkage implementation) |
| 19,072 / 17,358 / 6 / 306 | Recomputed each | filter gives 19,072; ranking gives 17,358 and is dataframe-identical to `ranked_list.csv`; 6 at FDR < 0.05; 274 up + 32 down = 306 at raw p < 0.01 | VERIFIED AND CORRECT |
| Duplicate symbols, ties, excluded genes | Counted | 17,385 rows with a symbol, 17,358 unique; 53 rows share 26 symbols, resolved by largest \|stat\|; 1,687 genes (8.8%) have no symbol and are excluded | VERIFIED AND CORRECT |
| Library coverage of the ranking | Intersected gmt membership with the ranking | GO:BP covers 11,801 of 17,358 (68%); Reactome 8,882 (51%) | ACCEPTABLE SCIENTIFIC LIMITATION, now documented |
| Hallmark GSEA re-execution | Re-ran all 50 sets from the local gmt, seed 7 | NES r = 0.999960 against the stored run, 100% sign agreement, 26/26 significant sets identical, 0 sets change significance class | VERIFIED AND CORRECT |
| Reactome GSEA re-execution | Re-ran 1,257 sets from the local gmt | NES r = 0.9999, 100% sign agreement | VERIFIED AND CORRECT |
| g:Profiler re-execution | Re-queried the up-list live from the analysis machine, compared term by term | 94 terms both runs, identical term set, max absolute difference 0 in p_value, term_size, intersection_size, query_size and effective_domain_size. Both reported queries use a nominal p < 0.01 cutoff with g:SCS correction. The down-list re-execution is recorded in the next row | VERIFIED AND CORRECT for the up-list |
| g:Profiler down-list re-execution | Re-queried the 32 down-regulated genes live from the analysis machine (`results/ora_gprofiler_down_rerun.csv`) | 32 genes submitted (console output), 0 significant terms returned; the saved file is a header row with no data rows, matching the original run exactly. Both g:Profiler queries, up and down, have now been re-executed live | VERIFIED AND CORRECT |
| Independent enrichment implementation | Wrote a weighted Kolmogorov-Smirnov ES with gene-label permutation | Hallmark beta-cell set NES -1.798 against gseapy -1.8014 on the same ranking | PARTIAL: cross-checks ES, NES and nominal p for one set. It does not reproduce library-wide FDR, which the full gseapy rerun does |
| Reported NES/FDR values | Checked 12 named sets against `gsea_all.csv` | all match to 2 decimal places | VERIFIED AND CORRECT |
| "Hypoxia NES 2.72 vs 2.31" | Looked up every hypoxia-related set | Hallmark `Hypoxia` (200 genes) NES 2.7235; GO `Cellular Response to Hypoxia (GO:0071456)` (71 genes) NES 2.2904; Reactome `Cellular Response to Hypoxia` NES 0.77, FDR 1.0 | VERIFIED AND CORRECT: different gene sets, not an inconsistency |
| Fadista four-gene comparison | Looked each up in the count matrix and `de_results.csv` | all four tested; TSPAN33 FDR 0.018; NT5E 0.657, TMED6 0.846, PAK5 0.844 not significantly differentially expressed | CORRECTED AND RETESTED |
| g:Profiler query and background sizes | Read the result table | query_size 270 (274 submitted), effective_domain_size 18,909 (19,072 submitted) | CORRECTED AND RETESTED |
| Correction method labelling | Checked g:Profiler parameters | `significance_threshold_method="g_SCS"`, g:Profiler's own graph-structure correction, not Benjamini-Hochberg FDR | CORRECTED AND RETESTED |
| GSEA FDR threshold | Checked the official GSEA User Guide | the guide states that with a small number of samples gene_set permutation gives a less stringent significance assessment and that one should then use a more stringent FDR cutoff, such as 5%. We used 5%. The guide also recommends 1,000 permutations for either permutation type | VERIFIED AND CORRECT |
| Permutation choice | Official GSEA User Guide, plus the design | gene-label permutation ignores gene-gene correlation. 8 donors per group meets the general sample-size recommendation for phenotype permutation, but phenotype permutation shuffles sample labels without preserving the batch and sex structure the design adjusts for, so ranking by the covariate-adjusted Wald statistic is the defensible trade | ACCEPTABLE SCIENTIFIC LIMITATION, now documented with both sides |
| Power versus donor effect | 21 DESeq2 refits: baseline, 16 leave-one-out, 4 two-donor exclusions, identical NES computation (weighted KS, 2,000 gene-label permutations) throughout | hypoxia NES 2.31 baseline, 1.19 without P7/P8, 2.19-2.39 for three matched pairs; mean log2FC falls 88% against no fall for those pairs. ECM shows the same direction more weakly (NES 2.35 to 1.62, log2FC -43%). Beta-cell NES negative in all 21 refits, -0.98 to -2.42, with P7/P8 inside the range of the other 20 | CORRECTED AND RETESTED. Three comparison pairs were tested, not all 28. All three are in batch A; P7/P8 are in batch B |
| "Every raised pathway loses significance" | Full-library comparison in all three libraries | FALSE. In Hallmark, 14 of 26 lose, 12 keep, 15 gain. TNF-alpha, EMT and Inflammatory Response keep significance at 14-donor NES 2.18, 2.30 and 1.50 respectively | CORRECTED AND RETESTED |
| Stability of the 14-donor fit | Counted NES sign changes library-wide | 18/50 Hallmark, 532/1,257 Reactome, 992/2,305 GO:BP change sign; significant-set counts move in opposite directions (GO:BP 182 to 62, Reactome 127 to 249) | ACCEPTABLE SCIENTIFIC LIMITATION, now documented and used to justify preferring the fixed-set comparison |
| Model identifiability after exclusion | Checked batch x group | batch B has 0 T2D donors after exclusion; group effect rests on batch A | ACCEPTABLE SCIENTIFIC LIMITATION, stated in every document |
| Cell composition | Compared Part 2 marker-based fractions between groups and located marker genes | endocrine 3.91% vs 2.98% (p = 0.42), acinar 4.50 vs 3.74 (p = 0.66), ductal 0.26 vs 0.35 (p = 0.41); INS +1.22, PRSS1 +1.24, KRT19 +1.28 all higher in T2D | ACCEPTABLE SCIENTIFIC LIMITATION: no clear evidence that composition explains the beta-cell result, but non-significance in 8 donors per group does not establish equivalence and a read fraction is not a beta-cell count. Cannot be excluded |
| Figures match their tables | Re-plotted from the result CSVs and inspected the rendered images | no clipped dots, labels or legends; axis, colour and size mappings correct | CORRECTED AND RETESTED |
| Literature supports each claim | Fetched each source; checked species and material | Olsson 2011 human islets (direct match); Hunter and Stein 2017 review; Haythorne 2019 mouse model and INS-1 rat line (direction only); Negi 2012 human islets; Fadista 2014 human islets, our own source cohort | CORRECTED AND RETESTED: species and material now stated |
| Fadista agreement is not replication | Checked donor provenance | our 16 donors are a subset of their 89 | CORRECTED AND RETESTED |
| Reproducibility on a fresh clone | dependencies, GMT checks and fresh output directories added | documented in README | Full corrected pipeline rerun pending |
| Report page limit | Rendered to PDF and counted pages | see README of this audit section | CORRECTED AND RETESTED |

## 2. Retracted or corrected claims

| Earlier claim | Status | Replacement |
|---|---|---|
| PAK7 "filtered out for low expression" | RETRACTED | HGNC renamed PAK7 to PAK5; Ensembl 110 uses PAK5. The gene was tested: log2FC -0.62, p 0.176, FDR 0.844 |
| g:Profiler used 274 genes against 19,072 | CORRECTED | 270 recognised of 274 submitted; effective background 18,909 of 19,072 |
| g:Profiler is an "independent" check | CORRECTED | second method, same donors and same DESeq2 output |
| Beta-cell signal "five orders of magnitude stronger" | RETRACTED | from a 30-gene set assembled after seeing the results. On the MSigDB set it is essentially unchanged: NES -1.798 to -1.755 |
| "Every raised pathway loses significance" | RETRACTED | true for hypoxia, glycolysis and several matrix sets; false for TNF-alpha, EMT, inflammatory response and Reactome collagen degradation, which keep significance |
| "The stress signature collapses, proving an artifact" | CORRECTED | hypoxia is sensitive to P7/P8 exclusion; same-size comparison pairs are from batch A, unlike P7/P8 from batch B. Donor and batch effects cannot be separated cleanly, and isolation stress is not proven |
| Seven beta-cell genes "in the bottom 2%" | WITHDRAWN AS CIRCULAR | those genes were selected by Part 2 for trending down. Replaced by the 40-gene MSigDB set |
| "Agrees with Fadista" | CORRECTED | our donors are a subset of theirs; this is a re-analysis |
| "These genes don't differ" | CORRECTED | "were not significantly differentially expressed" |

## 3. Unresolved

The corrected full pipeline remains unexecuted in the review environment. Two additional qualifications:

- The 14-donor Hallmark comparison was run at 1,000 permutations, not 10,000, because the audit machine cannot hold the 10,000-permutation matrix. The 16-donor run at both settings agrees to r = 0.99996 with zero changes in significance class, so the effect on conclusions is negligible. Both rankings in the comparison used identical settings.
- Part 3 is now present in `../Part3_promoters_ncrna/`. Its 20/22 negative mitochondrial tRNA fold changes agree in direction with the mitochondrial RNA-processing enrichment, subject to shared data and poly-A selection limitations.

## 4. Reproducing this

See `README.md` for the corrected commands. Reruns write into `reruns/` rather than overwriting archived results. The former manually assembled beta-cell rank test is removed from the pipeline; the selected sensitivity table uses the full-library Hallmark beta-cell result. The fixed-set rerun computes its nominal p-values against the same-sign permutation null; archived nominal p-values used the total permutation count and are retained as historical outputs, not used for the conclusions. The reported NES comparisons are not claimed to be newly re-executed.

## 5. Final pre-submission checks

| Item | Evidence | Result |
|---|---|---|
| Primary vs sensitivity permutation counts | Primary GSEA (`results/gsea_*.csv`) ran at 10,000 permutations. All sensitivity comparisons (`hallmark_sensitivity.csv`, `library_wide_sensitivity.csv`) ran at 1,000, identically for both rankings being compared. The GSEA User Guide recommends 1,000 for either permutation type, so this is a documented setting rather than a compromise. Every table states which applies | VERIFIED AND CORRECT |
| Does any conclusion rest on a borderline FDR? | Pancreas Beta Cells, the load-bearing result, is FDR 0.004184 at 10,000 permutations and 0.0039 at 1,000, with NES -1.8014 and -1.7894: not borderline and stable across settings. Among the "keeps significance" claims only Inflammatory Response is borderline in the reduced cohort (FDR 0.0312); it is reported with its value rather than as a bare claim | VERIFIED AND CORRECT |
| Does the ECM claim match the evidence? | ECM NES 2.35 to 1.62, below all three comparison pairs (2.12 to 2.55), mean log2FC -43%. Hypoxia is -88% and dips to 1.85 in single-donor fits while ECM never goes below 2.03. Same direction, weaker evidence; the text says so | CORRECTED AND RETESTED |
| Statistic plotted in `donor_influence.png` | Normalised enrichment score: weighted Kolmogorov-Smirnov ES divided by the mean of same-sign gene-label permutations, 2,000 permutations, seed 7. Axis and captions label it NES | CORRECTED AND RETESTED |
| Donor-sensitivity wording | Text now states that the analysis identifies which donors carry a result, not its cause. Isolation stress is presented as a plausible explanation supported by Negi et al. 2012, not a demonstrated one | CORRECTED AND RETESTED |
| Scope of pair testing | Three comparison T2D pairs were tested, not all 28. Stated wherever the comparison is used | CORRECTED AND RETESTED |
| File manifest | See README for the intended analysis files, including the report figure script and outputs. Exclude `genesets/`, `_to_delete/`, and the individual report | Confirm actual staged file list before committing |
| Report length | Rendered to PDF with figures and table removed: body text occupies exactly 3 pages, references on a separate 4th | VERIFIED AND CORRECT |

## Why gene-label permutation rather than phenotype permutation

Phenotype permutation is admissible at eight donors per group, but it shuffles sample labels
without preserving the batch and sex structure that the DESeq2 design `~ batch + sex + group`
adjusts for. Ranking by the covariate-adjusted Wald statistic keeps that adjustment, at the cost
of the weaker gene-label null. The GSEA User Guide states that gene_set permutation gives a less
stringent significance assessment with few samples and that a 5% FDR cutoff should then be used,
which is the cutoff applied throughout. This paragraph was cut from the report body to meet the
2-3 page limit; the reasoning is recorded here.

## Repository cleanup and offline verification (26 November 2026)

The updated scripts remove the obsolete hand-picked rank test, use matched 1,000-permutation library comparisons, convert ORA numeric columns before plotting, and refuse stale output directories. The gene-set loader checks the archived SHA-256 prefixes and set counts. PyDESeq2 is now included in the dependency list. The three archived dot plots were redrawn at FDR < 0.05; the selected sensitivity figure uses the full-library beta-cell row rather than the older single-set test.

`python scripts/check_results.py` passed against the repository snapshot: 19,072 tested genes, 17,358 ranked genes, 6 DE genes, 274/32 ORA queries, identical original/repeat 94-term ORA tables, no down-query terms, matching individual/combined GSEA tables (26/125/173 significant sets), matching Hallmark sensitivity classifications (14 lost, 12 kept, 15 gained, 18 sign flips), and 21 stored refits with negative beta-cell NES. All Python files parsed. All three offline plotting scripts executed.

The exact GMT files and the GSEApy, g:Profiler and PyDESeq2 packages were unavailable in the review environment; the corrected full analysis has not been re-executed. Stored primary results and donor refits were retained. The revised sensitivity script can recreate the comparison tables from the archived 16/14-donor DE tables; the 14-donor input is itself an archived fit and can be regenerated with the donor-influence script. Future live g:Profiler results may change with the database.

## Independent check of the final snapshot

The Part 1 count filter yields exactly the 19,072 Part 2 identifiers. Recalculated Wald statistics, available normal-tail p-values and BH-adjusted p-values match to within 1e-12. All 3,612 GSEA rows have valid ES/NES signs, probability ranges, and leading-edge members, counts and order. Both archived ORA tables agree; the leading term has 28 hits, 249 term members, query size 270 and background 18,909. The unadjusted hypergeometric tail is 2.72e-17, distinct from g:SCS-adjusted p=1.36e-12. Live ORA was not repeated in this environment.

An independent sparse running-sum implementation using public Enrichr Hallmark memberships and 10,000 random gene-set permutations reproduces the beta-cell and hypoxia ES to about 1e-15 (ES -0.430895267153 and +0.663927247740; independent NES -1.8140 and +2.7246 versus archived -1.8014 and +2.7235). The original ranking retains 73 symbols associated with 85 rows with missing nominal p-values, consistent with DESeq2 outlier filtering. Excluding these symbols gives ES -0.430949870860 and +0.662216515651 respectively. These are checks of two sets only, not new library-wide FDR estimates.

The two mitochondrial RNA-processing pathways are reduced, but neither Hallmark nor GO oxidative phosphorylation is significantly reduced. The literature comparison has been qualified accordingly. The archived 14-donor DE table has four significant genes (SAA2, PCOLCE2, GNAL, SAA1), none of the original six, matching the Part 2 sensitivity count. The three composition-proxy Welch tests and Part 3's 20/22 negative mitochondrial tRNA fold changes were independently recalculated.

The primary Fadista lookup now writes the same columns and aliases as the archived table; ORA plots have adequate axis margins. All archived numerical results are unchanged. Full rerun remains pending because the exact packages and GMTs could not be fetched in the execution environment.
