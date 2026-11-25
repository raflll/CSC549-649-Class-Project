#!/bin/bash

mkdir -p ref
cd ref

wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzvf grch38_genome.tar.gz

wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
