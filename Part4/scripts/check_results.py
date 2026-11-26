"""Offline integrity checks for the archived Part 4 results (no enrichment rerun)."""
from pathlib import Path
import ast
import numpy as np
import pandas as pd
from pathway_utils import rank_genes, LIBRARIES
from part4_robustness import compare

R = Path('results')
de = pd.read_csv('../Part2_statistics/results/de_results.csv')
rank = rank_genes(de)
archived = pd.read_csv(R / 'ranked_list.csv')
assert len(de) == 19072 and len(rank) == 17358
pd.testing.assert_frame_equal(rank, archived, check_dtype=False, rtol=1e-12, atol=1e-12)
assert int((de.padj < .05).sum()) == 6
up = (de.pvalue < .01) & (de.log2FoldChange > 0)
down = (de.pvalue < .01) & (de.log2FoldChange < 0)
assert (int(up.sum()), int(down.sum())) == (274,32)
print('PASS: Part 2 counts, ranking, 6 significant genes, 274/32 ORA queries')

G = pd.read_csv(R / 'gsea_all.csv')
for name, expected in zip(LIBRARIES, [(50,26),(1257,125),(2305,173)]):
    d = pd.read_csv(R / f'gsea_{name}.csv')
    assert (len(d),int((d['FDR q-val'] < .05).sum())) == expected
    a = d.sort_values('Term').reset_index(drop=True)
    b = G.loc[G.library.eq(name), d.columns].sort_values('Term').reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-12,atol=1e-12)
    print(f'PASS: {name}, {expected[0]} sets, {expected[1]} significant; combined table matches')

ora = pd.read_csv(R / 'ora_gprofiler.csv')
repeat = pd.read_csv(R / 'ora_gprofiler_rerun.csv')
cols = ['source','native','p_value','term_size','intersection_size','query_size','effective_domain_size']
assert len(ora)==len(repeat)==94 and set(ora.direction)=={'up'}
pd.testing.assert_frame_equal(ora[cols].sort_values(['source','native']).reset_index(drop=True),repeat[cols].sort_values(['source','native']).reset_index(drop=True),check_dtype=False)
assert pd.read_csv(R / 'ora_gprofiler_down_rerun.csv').empty
print('PASS: 94 ORA terms agree; down-query CSV has zero data rows')

h = pd.read_csv(R / 'hallmark_sensitivity.csv')
fields = ['NES','NOM p-val','FDR q-val']
def side(n):
    return h[['Term']+[f'{c}_{n}' for c in fields]].rename(columns={f'{c}_{n}':c for c in fields})
calc, stats = compare(side(16),side(14))
assert (stats['sets'],stats['sig16'],stats['sig14'],stats['lose'],stats['stay'],stats['gain'],stats['sign_flips']) == (50,26,27,14,12,15,18)
assert calc.set_index('Term').verdict.sort_index().equals(h.set_index('Term').verdict.sort_index())
print('PASS: sensitivity comparison logic matches archived Hallmark results')

b = h.set_index('Term').loc['Pancreas Beta Cells']
sel = pd.read_csv(R / 'robustness_sensitivity.csv').set_index('gene_set').loc['Pancreas Beta Cells']
assert np.isclose(sel.NES_16donors, b.NES_16, atol=.0005) and np.isclose(sel.NES_14donors,b.NES_14,atol=.0005)
assert np.isclose(sel.FDR_16donors,b['FDR q-val_16'],atol=.00005) and np.isclose(sel.FDR_14donors,b['FDR q-val_14'],atol=.00005)
D = pd.read_csv(R / 'donor_influence.csv')
assert len(D)==21 and D.n.value_counts().to_dict()=={15:16,14:4,16:1}
assert (D.HALLMARK_PANCREAS_BETA_CELLS_NES < 0).all()
print('PASS: selected beta-cell comparison matches full-library result; 21 refits, all beta-cell NES negative')

for path in sorted(Path('scripts').glob('*.py')):
    ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
print('PASS: all Python files parse')
print('Offline checks complete. This does not repeat GSEA, g:Profiler or DESeq2 refits.')
