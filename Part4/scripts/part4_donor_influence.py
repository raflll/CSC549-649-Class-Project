#!/usr/bin/env python3
"""Part 4: donor influence analysis.

Refits DESeq2 21 times (all 16 donors, each single-donor exclusion, and four two-donor
T2D exclusions) and scores four fixed gene sets identically each time, so that a donor
effect can be told apart from simple loss of power. Also records the mean log2 fold
change of each set, so effect size is compared and not only whether FDR crosses 0.05.

Run from the Part4 folder:   python3 scripts/part4_donor_influence.py
Needs pydeseq2 (see requirements.txt) and the cached gmt files in genesets/.
Takes roughly 20 minutes.
"""
import warnings, sys, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
B="../Part1_data_pipeline/"
U=".."
raw=pd.read_csv(B+"gene_count_matrix.csv",index_col=0)
don=pd.read_csv(B+"Table1_donors.csv")
mapt=pd.read_csv(B+"mapping_table.csv")
meta=don.merge(mapt[["run","flowcell"]],on="run")
meta["batch"]=np.where(meta.flowcell=="C11YEACXX","A","B")
meta["label"]=meta.pair+"_"+meta.group
meta=meta.set_index("run")
ids=raw.index.to_series().str.split("|",expand=True)
gmap=pd.DataFrame({"symbol":ids[1].values},index=ids[0].values); gmap=gmap[~gmap.index.duplicated()]
raw.index=ids[0].values; raw=raw[~raw.index.duplicated()]; raw=raw[meta.index]
keep=(raw>=10).sum(axis=1)>=8
cnt=raw[keep].T

def gmt_set(path,name):
    for line in open(path,encoding="utf-8"):
        p=line.rstrip("\n").split("\t")
        if p[0]==name: return [x for x in p[2:] if x]
    return []
SETS={
 "HALLMARK_PANCREAS_BETA_CELLS":["ABCC8","AKT3","CHGA","DCX","DPP4","ELP4","FOXA2","FOXO1","G6PC2","GCG",
   "GCK","HNF1A","IAPP","INS","INSM1","ISL1","LMO2","MAFB","NEUROD1","NEUROG3","NKX2-2","NKX6-1","PAK3",
   "PAX4","PAX6","PCSK1","PCSK2","PDX1","PKLR","SCGN","SEC11A","SLC2A2","SPCS1","SRP14","SRP9","SRPRB",
   "SST","STXBP1","SYT13","VDR"],
 "GO_Cellular_Response_to_Hypoxia":gmt_set(f"genesets/GO_Biological_Process_2025.gmt",
   "Cellular Response to Hypoxia (GO:0071456)"),
 "GO_Extracellular_Matrix_Organization":gmt_set(f"genesets/GO_Biological_Process_2025.gmt",
   "Extracellular Matrix Organization (GO:0030198)"),
 "Reactome_tRNA_Processing_Mitochondrion":gmt_set(f"genesets/Reactome_Pathways_2024.gmt",
   "tRNA Processing in the Mitochondrion"),
}
for k,v in SETS.items(): print(f"  set {k}: {len(v)} genes", flush=True)

rng=np.random.default_rng(7)
def es_of(idx,w,N,n):
    hit=np.zeros(N); hit[idx]=w[idx]; Nr=hit.sum()
    if Nr==0: return 0.0
    run=np.cumsum(hit/Nr-(1-(hit>0))/(N-n))
    return run[np.argmax(np.abs(run))]
def score(rk,genes,nperm=2000):
    sym=rk.symbol.values; w=np.abs(rk.stat.values.astype(float)); N=len(sym)
    idx=np.where(np.isin(sym,genes))[0]; n=len(idx)
    if n<10: return np.nan,np.nan,np.nan
    es=es_of(idx,w,N,n)
    null=np.array([es_of(rng.choice(N,n,replace=False),w,N,n) for _ in range(nperm)])
    same=null[np.sign(null)==np.sign(es)]
    nes=es/np.abs(same).mean() if len(same) else np.nan
    p=((null<=es).sum()+1)/(nperm+1) if es<0 else ((null>=es).sum()+1)/(nperm+1)
    # gene-level effect size: mean |log2FC| of set members, to separate power from signal
    return es,nes,p

def run(samples,tag):
    m=meta.loc[samples].copy()
    for c in ["batch","sex","group"]: m[c]=pd.Categorical(m[c])
    design="~batch + sex + group" if m.batch.nunique()>1 else "~sex + group"
    dds=DeseqDataSet(counts=cnt.loc[samples],metadata=m,design=design,refit_cooks=True,quiet=True)
    dds.deseq2(); st=DeseqStats(dds,contrast=["group","T2D","ND"],quiet=True); st.summary()
    r=st.results_df.copy(); r["symbol"]=gmap.reindex(r.index).symbol.values
    d=r.dropna(subset=["symbol","stat"]).copy(); d["a"]=d.stat.abs()
    d=d.sort_values("a",ascending=False).drop_duplicates("symbol")
    rk=d[["symbol","stat","log2FoldChange"]].sort_values("stat",ascending=False).reset_index(drop=True)
    out={"excluded":tag,"n":len(samples),"n_T2D":int((m.group=="T2D").sum()),"design":design}
    for k,g in SETS.items():
        es,nes,p=score(rk,g)
        out[k+"_NES"]=round(nes,3); out[k+"_p"]=round(p,4)
        out[k+"_meanLFC"]=round(rk[rk.symbol.isin(g)].log2FoldChange.mean(),3)
    return out

allS=meta.index.tolist()
rows=[run(allS,"none (all 16)")]
print("baseline done",flush=True)
for s in allS:
    rows.append(run([x for x in allS if x!=s], meta.loc[s,"label"]))
    print("  loo",meta.loc[s,"label"],flush=True)
p7p8=meta[meta.label.isin(["P7_T2D","P8_T2D"])].index.tolist()
rows.append(run([x for x in allS if x not in p7p8],"P7_T2D + P8_T2D"))
# three other T2D pairs for comparison
t2d=meta[meta.group=="T2D"].index.tolist()
for pair in [("P1_T2D","P2_T2D"),("P3_T2D","P4_T2D"),("P5_T2D","P6_T2D")]:
    dr=meta[meta.label.isin(pair)].index.tolist()
    rows.append(run([x for x in allS if x not in dr]," + ".join(pair)))
    print("  pair",pair,flush=True)
pd.DataFrame(rows).to_csv("results/donor_influence.csv",index=False)
print("written")
