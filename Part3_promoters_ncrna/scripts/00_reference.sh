#!/bin/bash
# Download the reference files used by Part 3 into ref/ (run from the Part3 folder).
# Same Ensembl release (110) and genome build (GRCh38) as Part 1.
mkdir -p ref
cd ref
wget -nc https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
wget -nc https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -nc https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
