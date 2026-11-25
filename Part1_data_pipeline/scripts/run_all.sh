#!/bin/bash

while ! grep -q "all done" download.log
do
    sleep 300
done

./02_fastqc.sh fastq > fastqc.log 2>&1
./03_hisat2.sh SRR_list.txt fastq > hisat2.log 2>&1
./04_stringtie.sh SRR_list.txt > stringtie.log 2>&1

echo "run_all done"
date
