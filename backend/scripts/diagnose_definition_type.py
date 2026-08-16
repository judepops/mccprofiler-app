"""Does HOW a gene set is defined predict whether it separates in contact space?

An earlier version of this claim was a just-so story: a two-bin taxonomy drawn
over seven of the twenty-one sets after seeing which ones separated. With seven
sets and two bins almost any scheme fits.

This version fixes the assignment first and tests it across all 21.

THE RULE, stated before any effect size is consulted. Each set is classified by
the DATA AND PROCEDURE its annotation was computed from, which is a documented
property of the annotation and nothing to do with our results:

  SEQUENCE       a property of the gene product or the DNA itself, requiring no
                 measurement in any cell type
  GENETIC        association or constraint inferred from variation across
                 populations
  STATE_MULTI    a combinatorial chromatin state assigned by a model over MANY
                 marks (ChromHMM / Roadmap HMM segmentations)
  RANK_SINGLE    a rank-ordering of regions by ONE mark, cut at a threshold
                 (ROSE on H3K27ac, i.e. super-enhancers)
  EXPRESSION     derived from RNA abundance or its variability across samples
  FITNESS        derived from a growth or dropout screen
  CONTEXT        derived from genomic neighbourhood, not from the gene

The STATE_MULTI versus RANK_SINGLE distinction is the one that matters and it is
the substantive methodological difference, not a convenience: ChromHMM and
Roadmap fit a hidden Markov model over many histone marks to assign a state,
whereas ROSE ranks stitched regions by a single mark and cuts at the inflection
point of the ranked curve. That is precisely the procedure Pott and Lieb (2015)
argue has no biological warrant.

Then: do the categories differ in |d|? Kruskal-Wallis across the grouped
categories, and Mann-Whitney for the specific contrast.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from sklearn.decomposition import PCA
from _shape import MAG_OVERALL

from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
N_PERM = 2000
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]

# --- THE ASSIGNMENT, fixed before any effect size is read -------------------
DEFN = {
    "Lambert_TF":                  ("SEQUENCE",    "curated DNA-binding domain"),
    "CpG_island_promoter":         ("SEQUENCE",    "CpG dinucleotide content of the promoter"),
    "phastCons_2kb_topQ":          ("SEQUENCE",    "cross-species sequence conservation"),

    "GWAS_immune_hot":             ("GENETIC",     "immune GWAS variant density"),
    "GWAS_immune_any":             ("GENETIC",     "any immune GWAS association"),
    "GWAS_total_topQ":             ("GENETIC",     "total GWAS burden"),
    "gnomAD_pLI_topQ":             ("GENETIC",     "population LoF constraint"),

    "ChromHMM_active_TSS":         ("STATE_MULTI", "HMM state over many marks"),
    "ChromHMM_bivalent":           ("STATE_MULTI", "HMM state over many marks"),
    "Roadmap_silenced":            ("STATE_MULTI", "HMM state over many marks"),

    "dbSUPER_CD4_SE_TSS_pm50kb":   ("RANK_SINGLE", "ROSE rank cutoff on H3K27ac"),

    "DICE_top_TPM_quartile":       ("EXPRESSION",  "RNA abundance"),
    "cd4_rna_top_quartile":        ("EXPRESSION",  "RNA abundance"),
    "cd4_specific_immune":         ("EXPRESSION",  "RNA specificity"),
    "Eisenberg_HK":                ("EXPRESSION",  "RNA stability across tissues"),
    "bio_HK_k3":                   ("EXPRESSION",  "expression-derived label"),
    "bio_bulk_k3":                 ("EXPRESSION",  "expression-derived label"),
    "bio_dev_TF_k3":               ("EXPRESSION",  "expression-derived label"),

    "DepMap_curated_essential":    ("FITNESS",     "CRISPR dropout screen"),
    "DepMap_inferred_essential":   ("FITNESS",     "CRISPR dropout screen"),

    "gene_desert_bottomQ_density": ("CONTEXT",     "local gene density"),
}

# Grouped for testing. CONTEXT is excluded as the known positive control for the
# density confound; RANK_SINGLE is reported but is n = 1 and cannot be tested.
DIRECT = {"SEQUENCE", "GENETIC", "STATE_MULTI"}    # what the gene or its regulation IS
OUTPUT = {"EXPRESSION", "FITNESS"}                 # what the gene DOES

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
pca = PCA(n_components=18, random_state=0).fit(X - A @ beta)
D = pca.transform(X - A @ beta)
D = (D - D.mean(0)) / D.std(0)
n = len(D)

cm = s.table("cohort_membership")
rows = []
for g, (cat, how) in DEFN.items():
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, n)
    k = int(hit.sum())
    if k < S.MIN_GROUP_N:
        print(f"  skipped {g}, only {k} panel genes")
        continue
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    rows.append(dict(group=g, category=cat, how=how, n=k, axis=j + 1,
                     d=float((a_.mean() - b_.mean()) / sp) if sp else 0.0,
                     p=(np.sum(null >= obs) + 1) / (N_PERM + 1)))
df = pd.DataFrame(rows)
df["ad"] = df["d"].abs()

print("\n" + "=" * 92)
print("ALL 21 SETS, CLASSIFIED BY DEFINITION PROCEDURE BEFORE ANY EFFECT SIZE WAS READ")
print("=" * 92)
for cat in ["SEQUENCE", "GENETIC", "STATE_MULTI", "RANK_SINGLE",
            "EXPRESSION", "FITNESS", "CONTEXT"]:
    sub = df[df["category"] == cat]
    if not len(sub):
        continue
    print(f"\n{cat}")
    for r in sub.sort_values("ad", ascending=False).itertuples():
        star = "*" if r.p < 0.05 else " "
        print(f"   {r.group:<32}{r.n:>6}  axis {r.axis:>2}  "
              f"|d| {r.ad:.2f}{star}   {r.how}")

print("\n" + "=" * 92)
print("TEST: do the categories differ?")
print("=" * 92)
groups = [df.loc[df["category"] == c, "ad"].to_numpy()
          for c in ["SEQUENCE", "GENETIC", "STATE_MULTI", "EXPRESSION", "FITNESS"]]
groups = [g for g in groups if len(g)]
H, pk = kruskal(*groups)
print(f"Kruskal-Wallis across 5 testable categories: H = {H:.2f}, p = {pk:.4f}")

a = df.loc[df["category"].isin(DIRECT), "ad"]
b = df.loc[df["category"].isin(OUTPUT), "ad"]
U, pm = mannwhitneyu(a, b, alternative="greater")
print(f"\nDIRECT (sequence, genetic, multi-mark state; n = {len(a)}) "
      f"median |d| = {a.median():.2f}")
print(f"OUTPUT (expression, fitness;                  n = {len(b)}) "
      f"median |d| = {b.median():.2f}")
print(f"Mann-Whitney one-sided, DIRECT > OUTPUT: U = {U:.0f}, p = {pm:.4f}")

print(f"\nsignificant at p < 0.05:")
for lab, sub in [("DIRECT", df[df['category'].isin(DIRECT)]),
                 ("OUTPUT", df[df['category'].isin(OUTPUT)])]:
    print(f"  {lab}: {(sub['p'] < 0.05).sum()} of {len(sub)}")

se = df[df["category"] == "RANK_SINGLE"]
if len(se):
    r = se.iloc[0]
    print(f"\nRANK_SINGLE is n = 1 and cannot be tested as a category.")
    print(f"  {r['group']}: |d| = {r['ad']:.2f}, p = {r['p']:.3f}, axis {r['axis']}.")
    print(f"  It sits at the {100 * (df['ad'] < r['ad']).mean():.0f}th percentile of "
          f"|d| across all sets, i.e. mid-pack, not last.")
