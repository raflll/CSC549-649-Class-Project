#!/bin/bash

list=$1
gtf=ref/Homo_sapiens.GRCh38.110.gtf

mkdir -p quant assembly

for run in `cat $list`
do
    echo "start $run"
    date
    stringtie -e -p 16 -G $gtf -o quant/$run/$run.gtf bam/$run.bam
    stringtie -p 16 -G $gtf -l $run -o assembly/$run.gtf bam/$run.bam
done

prepDE.py -i quant -l 101 -g gene_count_matrix.csv

ls assembly/SRR*.gtf > assembly/mergelist.txt
stringtie --merge -p 16 -G $gtf -o assembly/merged.gtf assembly/mergelist.txt

echo "all done"
date
