import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.size":13,"font.family":"DejaVu Sans",
    "axes.labelsize":14,"axes.titlesize":14,"xtick.labelsize":13,"ytick.labelsize":13})

# ---------- Figure 1: Hallmark dotplot ----------
import os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","results")
H=os.path.join(R,"gsea_MSigDB_Hallmark_2020.csv")
rows=[r for r in csv.DictReader(open(H, encoding="utf-8", newline="")) if float(r["FDR q-val"])<0.05]
for r in rows:
    r["NES"]=float(r["NES"]); r["FDR"]=float(r["FDR q-val"])
    r["lead"]=len(r["Lead_genes"].split(";"))
rows.sort(key=lambda r:r["NES"])
print("Fig1 sets at FDR<0.05:",len(rows))

fig,ax=plt.subplots(figsize=(9.6,0.40*len(rows)+2.0))
sc=ax.scatter([r["NES"] for r in rows],range(len(rows)),
    s=[28+r["lead"]*1.5 for r in rows],c=[r["FDR"] for r in rows],
    cmap="viridis_r",edgecolors="black",linewidths=0.7,zorder=3)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r["Term"] for r in rows])
ax.axvline(0,color="0.45",lw=1.2,zorder=1)
ax.grid(axis="x",color="0.90",lw=0.8,zorder=0)
lo=min(r["NES"] for r in rows); hi=max(r["NES"] for r in rows); pad=0.45*max(1,(hi-lo)/6)
ax.set_xlim(lo-pad,hi+pad); ax.set_ylim(-0.9,len(rows)-0.1)
ax.set_xlabel("Normalised enrichment score   (positive = higher in T2D)",labelpad=9)
ax.set_title("MSigDB Hallmark gene sets at FDR < 0.05\nT2D vs non-diabetic islets, 16 donors, 10,000 permutations",pad=12)
cb=fig.colorbar(sc,ax=ax,pad=0.02,fraction=0.040); cb.set_label("GSEA FDR q-value",fontsize=13)
for sz,lab in [(28+40*1.5,"40"),(28+110*1.5,"110")]:
    ax.scatter([],[],s=sz,c="white",edgecolors="black",linewidths=0.7,label=lab)
ax.legend(title="leading-edge genes",loc="lower right",frameon=True,fontsize=11,title_fontsize=11,labelspacing=1.1,borderpad=0.8)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","figures","fig1_hallmark_report.png"),dpi=200,bbox_inches="tight"); plt.close(fig)

# ---------- Figure 2: donor influence ----------
D=list(csv.DictReader(open(os.path.join(R,"donor_influence.csv"), encoding="utf-8", newline="")))
sets=[("HALLMARK_PANCREAS_BETA_CELLS","Pancreas beta cells\n(Hallmark)"),
      ("Reactome_tRNA_Processing_Mitochondrion","tRNA processing in the\nmitochondrion (Reactome)"),
      ("GO_Cellular_Response_to_Hypoxia","Cellular response to\nhypoxia (GO:0071456)"),
      ("GO_Extracellular_Matrix_Organization","Extracellular matrix\norganization (GO:0030198)")]
full={k:float([r for r in D if r["n"]=="16"][0][k+"_NES"]) for k,_ in sets}
single=[r for r in D if r["n"]=="15"]
pairs=[r for r in D if r["n"]=="14" and "P7_T2D" not in r["excluded"]]
p78=[r for r in D if r["n"]=="14" and "P7_T2D" in r["excluded"]][0]
print("Fig2 refits:",len(D),"| single:",len(single),"| control pairs:",len(pairs))

fig,axes=plt.subplots(1,4,figsize=(14.5,5.3))
for ax,(k,lab) in zip(axes,sets):
    ax.axhline(full[k],color="0.35",ls="--",lw=1.5,zorder=1)
    ax.scatter([1]*len(single),[float(r[k+"_NES"]) for r in single],s=46,c="0.62",
               edgecolors="0.25",linewidths=0.6,zorder=3)
    ax.scatter([2]*len(pairs),[float(r[k+"_NES"]) for r in pairs],s=74,c="#1f77b4",
               edgecolors="black",linewidths=0.7,zorder=3)
    ax.scatter([2],[float(p78[k+"_NES"])],s=180,marker="D",c="#d62728",
               edgecolors="black",linewidths=0.9,zorder=4)
    ax.axhline(0,color="0.80",lw=1.0,zorder=0)
    ax.set_xlim(0.5,2.5); ax.set_xticks([1,2]); ax.set_xticklabels(["1 donor\nremoved","2 donors\nremoved"])
    ax.set_title(lab,fontsize=13,pad=10)
    ax.grid(axis="y",color="0.92",lw=0.8,zorder=0)
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("Normalised enrichment score",fontsize=14)
handles=[Line2D([],[],marker="o",ls="",ms=8,mfc="0.62",mec="0.25",label="each donor removed (16 refits)"),
         Line2D([],[],marker="o",ls="",ms=9,mfc="#1f77b4",mec="black",label="3 other T2D pairs removed"),
         Line2D([],[],marker="D",ls="",ms=11,mfc="#d62728",mec="black",label="P7_T2D + P8_T2D removed"),
         Line2D([],[],ls="--",color="0.35",lw=1.5,label="all 16 donors")]
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,fontsize=12,bbox_to_anchor=(0.5,-0.015))
fig.suptitle("Donor influence on four fixed gene sets across 21 DESeq2 refits",fontsize=15,y=1.00)
fig.tight_layout(rect=[0,0.07,1,0.97])
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","figures","fig2_donor_report.png"),dpi=200,bbox_inches="tight"); plt.close(fig)
print("done")

# Selected matched-cohort comparison (both rankings: 1,000 permutations).
S = list(csv.DictReader(open(os.path.join(R, "robustness_sensitivity.csv"), encoding="utf-8", newline="")))
fig, ax = plt.subplots(figsize=(10, .38 * len(S) + 2))
labels = [r["gene_set"] + " [" + r["library"] + "]" for r in S]
for i, r in enumerate(S):
    a, b = float(r["NES_16donors"]), float(r["NES_14donors"])
    ax.plot([a,b], [i,i], color="0.65", lw=1)
    ax.scatter(a,i,c="#1f77b4",s=40,zorder=3)
    ax.scatter(b,i,c="#d62728",marker="D",s=40,zorder=3)
ax.set_yticks(range(len(S)), labels, fontsize=10)
ax.invert_yaxis(); ax.axvline(0,color="0.5",lw=.8)
ax.set_xlabel("Normalised enrichment score (NES)")
ax.set_title("Selected gene sets: 16 versus 14 donors (1,000 permutations each)", fontsize=12)
ax.legend(handles=[Line2D([],[],marker="o",ls="",color="#1f77b4",label="16 donors"), Line2D([],[],marker="D",ls="",color="#d62728",label="Without P7/P8")], loc="lower right", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(R,"..","figures","robustness_sensitivity.png"),dpi=200,bbox_inches="tight"); plt.close(fig)
