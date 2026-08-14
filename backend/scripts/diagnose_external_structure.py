"""Is there ANY positional structure to the external sets, or are they everywhere?

The enrichment grid says 14 of 21 sets clear no cell on PC2 x PC3. Two ways to
under-call structure from that, both fixed here.

It only tests the two axes on screen. A set strongly positioned on PC7 looks
flat on PC2 x PC3. So test all 19 components above the noise ceiling.

It only tests patches. A shifted centre of mass with heavy overlap produces no
patch anywhere while still being a real displacement, which is exactly what a
continuum would look like. So test the centroid, not the cells.

Two tests per set:

  1. UNIVARIATE. Standardised mean difference on each of the 19 dimensions,
     with a label-permutation null and Benjamini-Hochberg across all
     21 x 19 tests.
  2. MULTIVARIATE. Distance of the set's centroid from the panel centroid in
     the full 19-dimensional space, against the same permutation null. This
     asks "is this set displaced at all, in any direction" and is the honest
     single answer to whether the sets are positioned.

Displacement and separation are different claims and both are reported. A set
can be reliably displaced and still overlap almost completely, which is what a
continuum with weak biological correlates looks like.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from app.store import get_store
from app import store_schema as S

N_PERM = 2000
np.random.seed(0)

s = get_store()
emb = s.table("embeddings")
scree = s.table("pc_scree")
dims = [f"pc{int(r.pc)}" for r in scree.itertuples() if r.above_noise]
dims = [d for d in dims if d in emb.columns]

genes = s.genes[["gene_id", "symbol_key"]].merge(emb[["gene_id"] + dims], on="gene_id")
X = genes[dims].to_numpy(float)
X = (X - X.mean(0)) / X.std(0)          # per-dimension z, so all dims weigh equally
sym = genes["symbol_key"].to_numpy()
n = len(genes)
print(f"{n} genes x {len(dims)} dimensions above the noise ceiling\n")

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
sets = sorted(sizes[sizes >= S.MIN_GROUP_N].index)

rng = np.random.default_rng(0)
rows, multi = [], []

for g in sets:
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, n)
    k = int(hit.sum())
    if k < S.MIN_GROUP_N:
        continue

    obs_mean = X[hit].mean(0)
    # Displacement of the centroid, in SD units of the null. Under permutation
    # a centroid of k genes wanders by ~1/sqrt(k) per dimension, so the raw
    # distance is not comparable between sets of different size; the z below is.
    obs_dist = float(np.linalg.norm(obs_mean))

    perm_d = np.empty(N_PERM)
    perm_max = np.empty((N_PERM, len(dims)))
    for b in range(N_PERM):
        idx = rng.choice(n, size=k, replace=False)
        m = X[idx].mean(0)
        perm_d[b] = np.linalg.norm(m)
        perm_max[b] = m
    z_multi = (obs_dist - perm_d.mean()) / perm_d.std()
    p_multi = float((np.sum(perm_d >= obs_dist) + 1) / (N_PERM + 1))
    multi.append(dict(group=g, n=k, dist=obs_dist, z=z_multi, p=p_multi))

    for j, d in enumerate(dims):
        a, b_ = X[hit, j], X[~hit, j]
        # Cohen's d, pooled.
        sp = np.sqrt(((k - 1) * a.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
        dd = float((a.mean() - b_.mean()) / sp) if sp else 0.0
        p = float((np.sum(np.abs(perm_max[:, j]) >= abs(obs_mean[j])) + 1) / (N_PERM + 1))
        rows.append(dict(group=g, dim=d, n=k, cohens_d=dd, p=p))

df = pd.DataFrame(rows)
# Benjamini-Hochberg over every set x dimension test.
df = df.sort_values("p").reset_index(drop=True)
m = len(df)
df["q"] = (df["p"] * m / (df.index + 1)).cummin().clip(upper=1.0)

print("=" * 78)
print("1. MULTIVARIATE: is the set displaced at all, in any direction?")
print("=" * 78)
md = pd.DataFrame(multi).sort_values("z", ascending=False)
print(f"{'set':<32}{'n':>6}{'z':>9}{'p':>10}")
for r in md.itertuples():
    star = "  ***" if r.p < 0.001 else ("  **" if r.p < 0.01 else ("  *" if r.p < 0.05 else ""))
    print(f"{r.group:<32}{r.n:>6}{r.z:>9.1f}{r.p:>10.4f}{star}")
sig = int((md["p"] < 0.05).sum())
print(f"\n{sig} of {len(md)} sets are displaced from the panel centroid at p < 0.05.")

print("\n" + "=" * 78)
print("2. UNIVARIATE: the strongest set-by-dimension effects (BH q < 0.05)")
print("=" * 78)
hits = df[df["q"] < 0.05].reindex(df[df["q"] < 0.05]["cohens_d"].abs()
                                  .sort_values(ascending=False).index)
print(f"{'set':<32}{'dim':>6}{'n':>6}{'cohen d':>10}{'q':>10}")
for r in hits.head(18).itertuples():
    print(f"{r.group:<32}{r.dim:>6}{r.n:>6}{r.cohens_d:>10.2f}{r.q:>10.2g}")
print(f"\n{len(hits)} of {m} set-by-dimension tests survive BH correction.")
print(f"sets with at least one surviving dimension: {hits['group'].nunique()} of {len(sets)}")

print("\n" + "=" * 78)
print("3. HOW BIG? displacement vs separation")
print("=" * 78)
print("Cohen's d -> overlap of the two distributions (normal approximation):")
big = hits["cohens_d"].abs().max() if len(hits) else 0.0
from math import erf, sqrt
def overlap(d):     # overlapping coefficient for two equal-variance normals
    return 2 * (1 - 0.5 * (1 + erf(abs(d) / (2 * sqrt(2)))))
for d in [0.2, 0.4, 0.6, 0.8, 1.0, big]:
    print(f"   d = {d:4.2f}  ->  {overlap(d):5.1%} of the two distributions overlap")
print(f"\nLargest surviving effect anywhere: d = {big:.2f}")
