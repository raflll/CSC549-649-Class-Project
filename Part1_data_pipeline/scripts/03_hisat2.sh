#!/bin/bash

list=$1
fqdir=$2

mkdir -p bam

for run in `cat $list`
do
    echo "start $run"
    date
    hisat2 -p 16 --dta -x ref/grch38/genome -1 $fqdir/${run}_1.fastq -2 $fqdir/${run}_2.fastq 2> bam/$run.log | samtools sort -@ 4 -o bam/$run.bam -
    samtools index bam/$run.bam
done

echo "all done"
date
