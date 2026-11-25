#!/usr/bin/python

import sys
import matplotlib.pyplot as plt

plt.figure(figsize=(11, 4))
plt.xlim(0, 11)
plt.ylim(0.4, 3.7)
plt.axis('off')

def box(x, y, text, color):
    plt.text(x, y, text, ha='center', va='center', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.5', facecolor=color, edgecolor='black'))

def arrow(x1, y1, x2, y2, tool):
    plt.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='->', lw=1.5))
    plt.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.25, tool, ha='center', va='bottom', fontsize=8, style='italic')

box(1.0, 3.0, 'SRA\nGSE50244 / SRP029262\n16 donors (8 T2D, 8 ND)', 'lightgrey')
box(3.9, 3.0, 'FASTQ\n32 files\n2 x 101 bases', 'lightyellow')
box(6.8, 3.0, 'BAM\n16 sorted files', 'lightyellow')
box(9.7, 3.0, 'Count matrix\n62,754 genes x 16 donors', 'lightgreen')

arrow(2.15, 3.0, 3.25, 3.0, 'fastq-dump')
arrow(4.55, 3.0, 6.05, 3.0, 'HISAT2\n+ GRCh38 index')
arrow(7.55, 3.0, 8.5, 3.0, 'StringTie -e\n+ Ensembl GTF')

box(3.9, 1.0, 'Quality reports\n(FastQC)', 'white')
box(6.8, 1.0, 'mapping_table.csv\n(mapping rate, batch)', 'white')
box(9.7, 1.0, 'merged.gtf\n(new transcripts)', 'white')

arrow(3.9, 2.45, 3.9, 1.45, '')
arrow(6.8, 2.6, 6.8, 1.45, '')
arrow(7.3, 2.65, 9.0, 1.45, '')
plt.text(4.0, 1.95, 'FastQC', fontsize=8, style='italic')
plt.text(6.9, 1.95, 'HISAT2 logs', fontsize=8, style='italic')
plt.text(8.35, 2.05, 'StringTie (no -e)\n+ --merge', fontsize=8, style='italic')

plt.tight_layout()
plt.savefig(sys.argv[1], dpi=300)
