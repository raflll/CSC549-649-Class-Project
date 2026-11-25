#!/usr/bin/python

import sys

runs = []
pair = {}
group = {}
table = open(sys.argv[1], 'r')
table.readline()
for line in table:
    f = line.strip().split(',')
    runs.append(f[4])
    pair[f[4]] = f[0]
    group[f[4]] = f[1]
table.close()

out = open(sys.argv[2], 'w')
out.write('run,pair,group,flowcell,lane,read_pairs,unique_pct,multi_pct,overall_pct\n')

for run in runs:
    fq = open('fastq/' + run + '_1.fastq', 'r')
    header = fq.readline().split()[1].split(':')
    fq.close()
    flowcell = header[2]
    lane = header[3]

    log = open('bam/' + run + '.log', 'r')
    for l in log:
        l = l.strip()
        if 'reads; of these:' in l:
            pairs = l.split()[0]
        if 'aligned concordantly exactly 1 time' in l:
            unique = l.split()[1].strip('()%')
        if 'aligned concordantly >1 times' in l:
            multi = l.split()[1].strip('()%')
        if 'overall alignment rate' in l:
            overall = l.split()[0].strip('%')
    log.close()

    out.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (run, pair[run], group[run], flowcell, lane, pairs, unique, multi, overall))

out.close()
