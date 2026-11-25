#!/usr/bin/python

import sys
import matplotlib.pyplot as plt

labels = []
rates = []
colors = []

table = open(sys.argv[1], 'r')
table.readline()
for line in table:
    f = line.strip().split(',')
    labels.append(f[1] + '\n' + f[2])
    rates.append(float(f[8]))
    if f[3] == 'C11YEACXX':
        colors.append('tab:orange')
    else:
        colors.append('tab:blue')
table.close()

ymin = float(sys.argv[3])

plt.figure(figsize=(9, 4.5))
plt.bar(range(len(rates)), rates, color=colors)
plt.xticks(range(len(rates)), labels, fontsize=8)
plt.ylim(ymin, 100 + (100 - ymin) * 0.06)
plt.ylabel('Reads mapped to GRCh38 (%)')

if ymin == 0:
    step = 20
else:
    step = 5
plt.yticks(range(int(ymin), 101, step))

for i in range(len(rates)):
    plt.text(i, rates[i] + (100 - ymin) * 0.01, '%.1f' % rates[i], ha='center', va='bottom', fontsize=7)

plt.bar(0, 0, color='tab:orange', label='Batch A (flow cell C11YEACXX)')
plt.bar(0, 0, color='tab:blue', label='Batch B (flow cell D088TACXX)')
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=8, frameon=False)

plt.tight_layout()
plt.savefig(sys.argv[2], dpi=300)
