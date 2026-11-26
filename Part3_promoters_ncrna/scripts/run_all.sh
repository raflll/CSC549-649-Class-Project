#!/bin/bash
# Run from the Part3 folder: bash scripts/run_all.sh
# Needs ref/ from scripts/00_reference.sh and Part 1 and Part 2 outputs.
set -e
mkdir -p results figures
python3 scripts/01_genes.py ref/Homo_sapiens.GRCh38.110.gtf.gz results/genes.tsv
python3 scripts/02_biotypes.py results/genes.tsv ../Part2_statistics/results/de_results.csv results/
python3 scripts/03_novel.py ../Part1_data_pipeline/merged.gtf.gz results/genes.tsv results/
python3 scripts/04_promoters.py ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz results/genes.tsv results/
python3 scripts/05_motifs.py ref/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt results/
