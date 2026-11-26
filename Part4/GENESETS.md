# Gene set provenance

The `.gmt` files are not committed, because Enrichr redistributes libraries derived from
MSigDB and Reactome under their own licences. They are re-downloaded by the analysis script
and cached in `Part4/genesets/`. If the download fails or the libraries change, fetch them by
hand and place them in that folder under exactly these names.

| File | Source URL | Sets in file | sha256 (first 32 hex) | Downloaded |
|------|-----------|--------------|------------------------|------------|
| `MSigDB_Hallmark_2020.gmt` | `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020` | 50 | 4275592957a1587652092bb398cf7721 | 26 Nov 2026 |
| `Reactome_Pathways_2024.gmt` | `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_Pathways_2024` | 2104 | 1a7c43cb503d9db23494605169614a7b | 25 Nov 2026 |
| `GO_Biological_Process_2025.gmt` | `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2025` | 5342 | 0fc02458765cfce6c9348dce9f7a9397 | 25 Nov 2026 |

Save each file with UTF-8 encoding and Unix line endings. Writing them with the Windows
default cp1252 encoding fails, because some pathway names contain Greek letters.

`HALLMARK_PANCREAS_BETA_CELLS` membership used for the single-set robustness test was taken
from https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_PANCREAS_BETA_CELLS.html
(40 genes).

g:Profiler is queried live over its API (https://biit.cs.ut.ee/gprofiler/api/), so its
database version is whatever is current at run time. The stored result in
`results/ora_gprofiler.csv` is from 25 November 2026.
