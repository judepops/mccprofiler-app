"""Can the axes be made nameable? An experiment, not a pipeline change.

The problem, measured in diagnose_pc_names.py: each PC name covers only 22 to
35 percent of its own axis, and PC3 is 62 percent correlated with total signal,
so a shape claim read off it is partly an amount claim.

Neither is a naming failure. PCA maximises variance, not interpretability, and
nothing in it asks a component to be about one thing. So the question is whether
a different basis over the same subspace is more nameable.

Four bases, all spanning the same 19-dimensional space so gene geometry is
preserved:

  A  PCA                    what we ship now
  B  PCA + varimax          rotate for simple structure
  C  amount-removed + PCA   regress total_mcc out of every feature first
  D  amount-removed + varimax

Scored on:
  concentration  share of a factor's squared loading mass in its top 5 features.
                 Higher means fewer features carry the factor, so a short name
                 can be faithful to it.
  n_salient      how many features load above 0.20. Lower is nameable.
  |r| amount     worst correlation with total_mcc across factors. This is the
                 contamination that makes PC3 hard to quote.
  coverage       variance retained, so a gain in nameability can be priced.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from app.store import get_store

N_COMP = 19  # the components above the parallel-analysis noise ceiling
np.random.seed(0)


def varimax(L, gamma=1.0, q=200, tol=1e-6):
    """Kaiser varimax. Rotates loadings toward simple structure.

    Orthogonal, so the subspace and all gene-to-gene distances are preserved
    exactly. Only the choice of axes within it changes.
    """
    L = L.copy()
    p, k = L.shape
    R = np.eye(k)
    d = 0.0
    for _ in range(q):
        d_old = d
        Lam = L @ R
        u, s, vt = np.linalg.svd(
            L.T @ (Lam ** 3 - (gamma / p) * Lam @ np.diag(np.diag(Lam.T @ Lam)))
        )
        R = u @ vt
        d = s.sum()
        if d_old != 0 and d / d_old < 1 + tol:
            break
    return L @ R


def score(loadings, scores, amount, label):
    """loadings: features x factors. scores: genes x factors."""
    k = loadings.shape[1]
    conc, salient, ramt = [], [], []
    for j in range(k):
        l2 = loadings[:, j] ** 2
        l2 = l2 / l2.sum()
        conc.append(np.sort(l2)[::-1][:5].sum())
        salient.append(int((np.abs(loadings[:, j]) > 0.20).sum()))
        ramt.append(abs(np.corrcoef(scores[:, j], amount)[0, 1]))
    return dict(basis=label,
                conc_mean=float(np.mean(conc)),
                conc_best=float(np.max(conc)),
                salient_mean=float(np.mean(salient)),
                amt_max=float(np.max(ramt)),
                amt_mean=float(np.mean(ramt)))


s = get_store()
emb = s.table("embeddings").set_index("gene_id")
X = (s.table("features")
       .pivot(index="gene_id", columns="feature", values="z")
       .loc[emb.index])
feat_names = list(X.columns)
Xv = X.values
amount = X["total_mcc"].values
print(f"{Xv.shape[0]} genes x {Xv.shape[1]} features, {N_COMP} components\n")

results, bases = [], {}

# ---- A / B: as-is ---------------------------------------------------------
p1 = PCA(n_components=N_COMP, random_state=0).fit(Xv)
S1 = p1.transform(Xv)
L1 = p1.components_.T * np.sqrt(p1.explained_variance_)
results.append(score(L1, S1, amount, "A  PCA"))
bases["A"] = (L1, S1, p1.explained_variance_ratio_.sum())

L2 = varimax(L1)
# Rotated scores: project the data onto the rotated loading directions.
S2 = Xv @ np.linalg.pinv(L2).T
results.append(score(L2, S2, amount, "B  PCA + varimax"))
bases["B"] = (L2, S2, p1.explained_variance_ratio_.sum())

# ---- C / D: amount removed first -----------------------------------------
# Regress total_mcc out of every feature. Shape is then whatever is left when
# how much signal a gene has is already accounted for.
a = (amount - amount.mean()) / amount.std()
Xr = np.empty_like(Xv)
for j in range(Xv.shape[1]):
    b = np.dot(a, Xv[:, j]) / np.dot(a, a)
    Xr[:, j] = Xv[:, j] - b * a

p3 = PCA(n_components=N_COMP, random_state=0).fit(Xr)
S3 = p3.transform(Xr)
L3 = p3.components_.T * np.sqrt(p3.explained_variance_)
results.append(score(L3, S3, amount, "C  amount-removed PCA"))
bases["C"] = (L3, S3, p3.explained_variance_ratio_.sum())

L4 = varimax(L3)
S4 = Xr @ np.linalg.pinv(L4).T
results.append(score(L4, S4, amount, "D  amount-removed + varimax"))
bases["D"] = (L4, S4, p3.explained_variance_ratio_.sum())

print("=" * 84)
print(f"{'basis':<30}{'top5 mass':>11}{'best':>8}{'salient':>9}"
      f"{'max|r|amt':>11}{'mean|r|amt':>12}")
print("=" * 84)
for r in results:
    print(f"{r['basis']:<30}{r['conc_mean']:>10.1%}{r['conc_best']:>8.0%}"
          f"{r['salient_mean']:>9.1f}{r['amt_max']:>11.3f}{r['amt_mean']:>12.3f}")

print("\nvariance retained by the 19-dim subspace:")
for k in "ABCD":
    print(f"  {k}: {bases[k][2]:.1%}")

# ---- what the rotated factors are actually about --------------------------
for key, title in [("B", "PCA + VARIMAX"), ("D", "AMOUNT-REMOVED + VARIMAX")]:
    L, S, _ = bases[key]
    print("\n" + "=" * 84)
    print(f"{title}: the six most concentrated factors, and what loads on them")
    print("=" * 84)
    l2n = (L ** 2) / (L ** 2).sum(axis=0, keepdims=True)
    conc = np.sort(l2n, axis=0)[::-1][:5].sum(axis=0)
    for j in np.argsort(conc)[::-1][:6]:
        ramt = np.corrcoef(S[:, j], amount)[0, 1]
        print(f"\n  factor {j+1}   top-5 mass {conc[j]:.0%}   r(amount) {ramt:+.2f}")
        for i in np.argsort(np.abs(L[:, j]))[::-1][:6]:
            print(f"      {L[i, j]:+.3f}  {feat_names[i]}")
