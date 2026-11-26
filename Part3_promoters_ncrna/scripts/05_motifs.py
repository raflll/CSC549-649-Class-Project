"""Step 5: motif enrichment in promoters of changed genes.

What: for each JASPAR motif, test if promoters of up (or down) genes contain it more often than matched background promoters.
Why:  genes that change together may share transcription factor (TF) binding sites.
How:  1. Each motif is a position weight matrix (PWM, log-odds of each base per position; background 0.25 each, pseudocount 0.8).
      2. Slide it over every promoter on both strands. Relative score = (score - min) / (max - min) of the PWM.
         A promoter is a 'hit' if its best relative score >= 0.85.
      3. Background: for every changed promoter pick the 5 unchanged tested promoters with the closest GC content
         and CpG observed/expected ratio (from all of them, without reuse). This matters: GC-rich promoters hit
         GC-rich motifs by chance, and CpG-island promoters hit NRF1, E2F and SP/KLF motifs, so unmatched groups
         give false motifs. Only matched promoters are scanned.
      4. Fisher exact test (hit / no hit, changed vs matched background), Benjamini-Hochberg correction.
Usage: python3 scripts/05_motifs.py ref/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt results/
Slow step (about 10 minutes): the score matrix is cached in ref/best_scores.npy.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

jaspar, out = sys.argv[1:3]
out = out.rstrip("/") + "/"
THRESH, RATIO = 0.85, 5

def read_fasta(path):
    ids, seqs = [], []
    for line in open(path):
        if line[0] == ">":
            ids.append(line[1:].strip()); seqs.append("")
        else:
            seqs[-1] += line.strip()
    return ids, seqs

def read_jaspar(path):
    motifs, name, rows = [], None, {}
    for line in open(path):
        if line[0] == ">":
            if name:
                motifs.append((name, np.array([rows[b] for b in "ACGT"], float).T))
            name, rows = line[1:].strip().replace("\t", " "), {}
        elif line.strip():
            rows[line[0]] = [float(x) for x in line.split("[")[1].split("]")[0].split()]
    motifs.append((name, np.array([rows[b] for b in "ACGT"], float).T))
    return motifs                                   # list of (name, L x 4 count matrix)

def to_pwm(counts):
    freq = (counts + 0.8 * 0.25) / (counts.sum(axis=1, keepdims=True) + 0.8)
    return np.log2(freq / 0.25)

# sequences: up, down, background
ids, seqs, label = [], [], []
for n in ["up", "down", "background"]:
    i, s = read_fasta(f"{out}promoters_{n}.fa")
    ids += i; seqs += s; label += [n] * len(i)
label = np.array(label)
P = len(seqs[0])
print("promoters:", {n: int((label == n).sum()) for n in ["up", "down", "background"]})
gc = np.array([(s.count("G") + s.count("C")) / len(s) for s in seqs])
cpg = np.array([x.count("CG") * len(x) / max(1, x.count("C") * x.count("G")) for x in seqs])   # CpG observed/expected
feat = np.column_stack([gc / gc.std(), cpg / cpg.std()])                                        # both on the same scale
bg_idx = np.where(label == "background")[0]
tree = cKDTree(feat[bg_idx])

def matched_bg(fg_idx):
    """For each changed promoter, RATIO unused background promoters closest in (GC, CpG o/e)."""
    used, chosen = set(), []
    for i in np.random.default_rng(1).permutation(fg_idx):
        _, nb = tree.query(feat[i], k=RATIO * 40)
        free = [bg_idx[k] for k in nb if bg_idx[k] not in used][:RATIO]
        used.update(free); chosen += free
    return np.array(chosen, int)

up, down = np.where(label == "up")[0], np.where(label == "down")[0]
sets = {"up": up, "down": down, "changed": np.concatenate([up, down])}
bgs = {n: matched_bg(fg) for n, fg in sets.items()}
need = np.unique(np.concatenate([up, down] + list(bgs.values())))      # only these promoters are scanned
row = {int(g): k for k, g in enumerate(need)}
print("promoters scanned:", len(need))

onehot = np.zeros((len(need), P, 4), np.float32)
lut = {"A": 0, "C": 1, "G": 2, "T": 3}
for k, g in enumerate(need):
    idx = np.array([lut.get(ch, -1) for ch in seqs[g]])
    ok = idx >= 0
    onehot[k, np.where(ok)[0], idx[ok]] = 1
motifs = read_jaspar(jaspar)

def best_score(m):
    pwm = to_pwm(m[1])
    lo, hi = pwm.min(axis=1).sum(), pwm.max(axis=1).sum()
    best = np.full(len(need), -np.inf, np.float32)
    for w in (pwm, pwm[::-1, ::-1]):                # forward and reverse-complement (A,C,G,T -> T,G,C,A)
        L = len(w); S = np.zeros((len(need), P - L + 1), np.float32)
        for j in range(L):
            S += onehot[:, j:j + P - L + 1, :] @ w[j].astype(np.float32)
        best = np.maximum(best, S.max(axis=1))
    return (best - lo) / (hi - lo)

cache = "ref/best_scores.npy"
if os.path.exists(cache):
    B = np.load(cache)
else:
    B = np.array([best_score(m) for m in motifs])   # motifs x scanned promoters
    np.save(cache, B)
hit = B >= THRESH

res = []
for sname, fg in sets.items():
    bg = bgs[sname]
    fr, br = [row[int(g)] for g in fg], [row[int(g)] for g in bg]
    print(f"{sname}: {len(fg)} changed promoters, {len(bg)} matched background, GC {gc[fg].mean():.3f} vs {gc[bg].mean():.3f}, CpG o/e {cpg[fg].mean():.3f} vs {cpg[bg].mean():.3f}")
    for m, (name, _) in enumerate(motifs):
        a, b = int(hit[m, fr].sum()), int(hit[m, br].sum())
        odds, p = fisher_exact([[a, len(fg) - a], [b, len(bg) - b]])
        res.append((sname, name, a / len(fg), b / len(bg), odds, p))
r = pd.DataFrame(res, columns=["set", "motif", "frac_changed", "frac_background", "odds_ratio", "p"])
r["q"] = r.groupby("set")["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
r = r.sort_values(["set", "p"]).round(5)
r.to_csv(out + "motif_results.csv", index=False)
print(r.groupby("set").head(8).to_string(index=False))
print("motifs with q < 0.05:", r[r.q < 0.05].groupby("set").size().to_dict())

# figure: top 10 motifs per set by p value
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, sname in zip(axes, sets):
    t = r[r.set == sname].head(10)[::-1]
    ax.barh(t.motif, -np.log10(t.p), color="#555555")
    ax.set_title(f"{sname} (no motif has q < 0.05)" if not (r[r.set == sname].q < 0.05).any() else sname, fontsize=9); ax.set_xlabel("-log10 p, uncorrected")
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/motif_top.png", dpi=200)
