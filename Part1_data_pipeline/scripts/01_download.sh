#!/bin/bash

mkdir -p fastq

for run in `cat SRR_list.txt`
do
    echo "start $run"
    date
    fastq-dump --split-files -O fastq $run
done

echo "all done"
date
