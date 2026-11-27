# Part 4: interpretation

## What the pathways say

In the full 16-donor analysis two groups of pathways move in opposite directions.

**Reduced in T2D.** Beta-cell identity (Hallmark Pancreas Beta Cells, NES -1.80, FDR 0.004), protein secretion (NES -2.13) and mitochondrial RNA processing (Reactome tRNA Processing in the Mitochondrion, NES -3.12; rRNA Processing in the Mitochondrion, NES -2.95). These are coordinated transcript-level shifts, consistent with altered beta-cell-associated biology but not direct measurements of cell identity, secretion or mitochondrial activity.

**Raised in T2D.** Hypoxia is the strongest (Hallmark Hypoxia, NES 2.72), with inflammation, epithelial-mesenchymal transition, glycolysis, and collagen and extracellular matrix alongside it. g:Profiler, run on a fixed gene list rather than the ranking, puts response to hypoxia first (g:SCS-adjusted p = 1.4e-12). That is a second method, not independent evidence: it uses the same donors and the same DESeq2 output.

## Which of this survives a change of donors

Part 2 found that two T2D donors, P7_T2D and P8_T2D, carry hypoxia and stress genes ten to fifteen times higher than every other donor, and that none of the original six differentially expressed genes stays significant after removing them; four other genes reach FDR < 0.05. The obvious test, refitting without them, has a trap: any result can lose significance in a smaller, more imbalanced cohort simply through lost power. So we refitted DESeq2 twenty-one times, once for all 16 donors, once for each single-donor exclusion, and once for each of four T2D pairs, and scored four fixed gene sets identically every time.

The comparison helps assess sensitivity beyond loss of sample size (`figures/donor_influence.png`).

- **Hypoxia is sensitive to exclusion of P7/P8.** Cellular Response to Hypoxia scores NES 2.31 in the full cohort, 1.19 without P7/P8, and 2.19 to 2.39 for the three other T2D pairs we tested, which remove the same number of samples. The gene-level effect moves the same way: mean log2 fold change of the set falls 88%, from 0.54 to 0.06, while the comparison pairs show no fall at all (0.55 to 0.64). This pattern is consistent with influence beyond sample size alone. The three comparison pairs are all in batch A, whereas P7/P8 are in batch B. Thus donor influence cannot be separated from changed batch representation; only three alternative pairs were examined (four of 28 total pairs).
- **Extracellular matrix moves in the same direction but less strongly.** NES 2.35 to 1.62, which is below all three comparison pairs (2.12 to 2.55), with mean log2 fold change falling 43% against no fall for those pairs. Across the sixteen single-donor exclusions, hypoxia drops as far as 1.85 whereas ECM never goes below 2.03. So the matrix signature is donor-influenced in the same direction, but on weaker evidence than hypoxia, and we report it that way. Reactome mitochondrial tRNA processing behaves similarly, -3.13 to -1.95 against -2.89 to -3.15.
- **The beta-cell set remains negative after excluding these donors.** Its NES is negative in all 21 refits, between -0.98 and -2.42. The P7/P8 value of -1.72 sits inside the range the other twenty refits produce, and the four two-donor exclusions span only -1.56 to -1.72. It is not insensitive to every donor: removing P4_ND weakens it to -0.98. The point is narrower than "unmoved" - it is that P7_T2D and P8_T2D are not what carries it.

At full-library scale the picture is messier and does not support a blanket claim. In the Hallmark library 14 of 26 significant sets lose significance without P7/P8, 12 keep it and 15 others gain it. Hypoxia and glycolysis lose it, dropping to NES 0.89 and -1.22. TNF-alpha signalling via NF-kB keeps it at 2.18, epithelial-mesenchymal transition at 2.30, and inflammatory response at 1.50, though that last is borderline at FDR 0.031. So it is **not** true that every raised pathway depends on those donors; TNF-alpha, EMT and inflammatory response retain significance in this comparison. Protein secretion reverses sign (-2.08 to +0.98) and is not directionally stable; mitochondrial tRNA and rRNA processing remain negative and significant (FDR 0.038 and 0.0459), while mitochondrial RNA degradation loses significance. These comparison runs used 1,000 permutations for both cohorts; the primary results above used 10,000.

The 14-donor fit should also be read with care. It is unstable rather than merely weaker: 18 of 50 Hallmark sets, 532 of 1,257 Reactome sets and 992 of 2,305 GO:BP sets change the sign of their enrichment score, and the count of significant sets moves in opposite directions between libraries (GO:BP 182 to 62, Reactome 127 to 249). Excluding those donors also leaves batch B with no T2D donors, so the group effect then rests on batch A alone. The fixed-set comparison complements the full-library FDR table, but same-size exclusions do not control the batch imbalance.

## Does it fit type 2 diabetes?

The beta-cell-associated and mitochondrial RNA-processing signatures fit some aspects of published T2D biology, but the agreement is incomplete.

- Loss of mature beta-cell markers has been implicated in T2D (Hunter and Stein 2017, a review covering human and rodent work). Here PDX1 is negative but not individually significant (FDR 0.765), whereas MAFA and INS have positive, nonsignificant fold changes. The set-level finding is not proof of dedifferentiation.
- Reduced oxidative phosphorylation gene expression has been reported directly in **human** T2D islets against non-diabetic donors, with 21 downregulated OXPHOS genes against about 7 expected by chance (Olsson et al. 2011). However, our oxidative phosphorylation sets are not significantly reduced (Hallmark NES +0.83, FDR 0.873; GO NES +1.19, FDR 0.445). The specific reduction in mitochondrial RNA-processing transcripts is related but does not reproduce that oxidative phosphorylation finding.
- Haythorne et al. (2019) report increased glycolytic expression and impaired mitochondrial metabolism in a **mouse** model (betaV59M) and the INS-1 rat beta-cell line. Our Hallmark glycolysis signal is positive but loses significance after P7/P8 exclusion. Our mitochondrial RNA-processing results do not measure metabolism; the agreement is therefore partial.

The raised direction has an alternative explanation that fits just as well. Negi et al. (2012) used laser capture on human islets before isolation, after isolation and after three days in culture, and found a strong inflammatory response induced by isolation and continuing in culture, together with loss of the mature endocrine phenotype. Islet isolation stress is not by itself a confound here, because all 16 donors went through the same procedure; it becomes one only where it differs between groups, which remains a plausible explanation for the pattern in these two donors.

One more observation points the same way. The differential expression result is globally asymmetric: 56.1% of all genes have a positive fold change, 90% of genes at raw p < 0.01 are higher in T2D, and the median Wald statistic is +0.17 rather than near zero. This asymmetry alone cannot distinguish disease biology from a technical or donor-related shift.

## Could the beta-cell result just be cell composition?

This is bulk islet RNA-seq, so a lower beta-cell score could mean fewer beta cells in the sample rather than beta cells expressing less. Part 2's marker-based fractions give no clear evidence for that: endocrine 3.91% in T2D against 2.98% in ND (Welch p = 0.42), acinar 4.50 against 3.74 (p = 0.66), ductal 0.26 against 0.35 (p = 0.41).

Those numbers have to be read carefully. A p-value of 0.42 in 8 donors per group does not establish that composition is equivalent, only that we cannot detect a difference; and an endocrine read fraction computed from INS, GCG, SST, PPY and GHRL is a marker proxy, not a count of beta cells. The endocrine fraction is numerically higher in the diabetic samples, but cannot identify beta-cell abundance. Insulin message itself is higher in T2D (log2FC +1.22), as are the exocrine markers PRSS1 (+1.24) and KRT19 (+1.28). The reduced set score therefore comes from the regulatory and secretory-machinery members of the set, the leading edge in `results/gsea_MSigDB_Hallmark_2020.csv`, rather than from less hormone message.

That is a refinement, not a resolution. **The marker analysis does not give clear evidence that composition explains the beta-cell result, but it cannot exclude it.** We have not measured insulin secretion, mitochondrial respiration, or anything about an individual beta cell. Deconvolution against a single-cell reference could help distinguish the explanations.

## Comparison with Fadista et al. 2014

This is a comparison, not a replication. Our 16 donors are a subset of their 89, so agreement is not independent evidence.

Fadista et al. identified 35 eQTL genes whose expression differed between normoglycemic and hyperglycemic donors, and showed that silencing TSPAN33, NT5E, TMED6 and PAK7 in INS1 cells reduced glucose-stimulated insulin secretion. All four were tested in our data. TSPAN33 is lower in our T2D donors (log2FC -1.01, FDR 0.018). NT5E, TMED6 and PAK5 were not significantly differentially expressed (FDR 0.66, 0.85, 0.84). PAK7 appears in our results under its current HGNC symbol PAK5, which is what the Ensembl 110 annotation uses.

## Connection to Part 3

Part 3 reports negative fold changes for 20 of 22 mitochondrial tRNAs, consistent with the reduced mitochondrial RNA-processing signal. No promoter motif reached corrected significance. Both analyses use the same cohort, so this is complementary evidence, not replication; poly-A selection also complicates small-RNA interpretation.

## Limitations

- 8 donors per group, so the gene-level analysis is underpowered.
- Batch is partly confounded with group: 6 of the 8 T2D donors are in batch A.
- Hypoxia and ECM are sensitive to P7/P8 exclusion, which leaves batch B with no T2D donors. All three comparison pairs are from batch A.
- Islet purity varies about 15-fold between donors, so differences in cell composition can look like differences in expression.
- T2D was defined by HbA1c alone; treatment, culture time and cause of death are unknown.
- The reduced direction rests on GSEA only, since g:Profiler returned no significant terms for the 32 genes at raw p < 0.01. This does not rule out enrichment.
- 1,687 tested genes (8.8%) have no gene symbol and are excluded from GSEA; the libraries themselves cover only 68% (GO:BP) and 51% (Reactome) of the ranked genes.
- The g:Profiler list uses raw p < 0.01 rather than FDR, so that analysis is exploratory. The ranking retains finite Wald statistics for 73 symbols associated with missing Part 2 nominal p-values (consistent with Cook-distance filtering); outlier-sensitive statistics are a limitation. Excluding these symbols barely changes the independent beta-cell and hypoxia ES checks, but other pathways have not been re-evaluated under this restriction. GSEA FDR is estimated separately per library.
- Enrichment is correlation, not mechanism.

## Next steps, for the synthesis slide

1. Deconvolve cell composition against a published single-cell islet reference. The marker-based check cannot exclude composition effects; deconvolution could help distinguish these explanations.
2. Consider full GSEA with phenotype permutation. Our 8 donors per group meet the general sample-size recommendation, but it is not a free upgrade: phenotype permutation shuffles sample labels without preserving the batch and sex structure our design adjusts for, so it would test a null that ignores the batch confounding. Ranking by the covariate-adjusted Wald statistic, as we did, keeps that adjustment at the cost of the less stringent gene-label null. A design that balanced batch against group would let both be had at once.
3. Replicate in a cohort balanced across sequencing batches, with islet handling and culture time recorded, so the donor-level stress effect can be modelled rather than excluded.
4. Take a candidate back to the bench. TSPAN33 is the obvious one: it is the most robust single gene in Part 2, it is lower in our T2D donors, and Fadista et al. showed that silencing it reduces glucose-stimulated insulin secretion in INS1 cells.

## References

1. Fadista J, Vikman P, Ottosson Laakso E, et al. Global genomic and transcriptomic analysis of human pancreatic islets reveals novel genes influencing glucose metabolism. Proc Natl Acad Sci USA. 2014;111(38):13924-13929. doi:10.1073/pnas.1402665111
2. Negi S, Jetha A, Aikin R, Hasilo C, Sladek R, Paraskevas S. Analysis of beta-cell gene expression reveals inflammatory signaling and evidence of dedifferentiation following human islet isolation and culture. PLoS One. 2012;7(1):e30415. doi:10.1371/journal.pone.0030415
3. Olsson AH, Yang BT, Hall E, Taneera J, Salehi A, Dekker Nitert M, Ling C. Decreased expression of genes involved in oxidative phosphorylation in human pancreatic islets from patients with type 2 diabetes. Eur J Endocrinol. 2011;165(4):589-595. doi:10.1530/EJE-11-0282
4. Haythorne E, Rohm M, van de Bunt M, et al. Diabetes causes marked inhibition of mitochondrial metabolism in pancreatic beta-cells. Nat Commun. 2019;10(1):2474. doi:10.1038/s41467-019-10189-x
5. Hunter CS, Stein RW. Evidence for loss in identity, de-differentiation, and trans-differentiation of islet beta-cells in type 2 diabetes. Front Genet. 2017;8:35. doi:10.3389/fgene.2017.00035
