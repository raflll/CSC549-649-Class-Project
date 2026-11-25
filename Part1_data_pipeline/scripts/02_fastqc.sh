#!/bin/bash

fqdir=$1

mkdir -p fastqc

fastqc -t 8 --outdir fastqc $fqdir/*.fastq
