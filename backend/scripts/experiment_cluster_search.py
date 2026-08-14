"""Does clustering structure appear once the substrate is cleaned or conditioned?

Every clustering test in this project so far has been run on all 91 features, on
the whole panel, at once. Three things measured in REVIEW.md say that is the
least favourable setting available, and none of them has been tried:

  1. Eight of the nineteen retained components are dominated by feature families
     that reproduce at rho 0.27 to 0.42, so a good fraction of the space is
     structured measurement noise. Clustering in noise dimensions cannot help
     and can actively hide structure by inflating distances.
  2. Amount dominates the geometry and carries most of the biological signal
     (nested baselines, shape-corrected column). Shape clusters could exist
     inside an amount stratum and be washed out globally.
  3. Nobody has clustered the amount-free space on its own.

So: four substrates x the whole panel, plus the same substrates within amount
tertiles. If clusters exist anywhere reachable from the current features, this
finds them. If they do not, the negative result is far stronger than the current
one because it survives the conditions most favourable to finding them.

Tests per substrate, chosen so no single method's failure mode decides it:

  silhouette      KMeans k=2..8, against a permuted-feature null so the number
                  is read against what structureless data of this shape gives
  gap statistic   vs a uniform null in PCA space, the standard k selector
  HDBSCAN         density clustering, which can return "no clusters" honestly
  dip test        Hartigan unimodality on each retained component

Run in cd4env (needs sklearn, hdbscan, diptest).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import hdbscan
import diptest

from app.store import get_store

warnings.filterwarnings("ignore")
SEED = 0
KS = range(2, 9)
N_PC = 19
TRUST = 0.7


def prep(X, n_pc=N_PC, seed=SEED):
    n_pc = min(n_pc, X.shape[1])
    p = PCA(n_components=n_pc, random_state=seed).fit(X)
    return p.transform(X), p


def permuted(X, rng):
    """Per-feature permutation: marginals kept, joint structure destroyed."""
    Z = X.copy()
    for j in range(Z.shape[1]):
        rng.shuffle(Z[:, j])
    return Z


def gap_statistic(S, ks=KS, n_ref=10, seed=SEED):
    """Tibshirani gap against a uniform box null in the PCA space."""
    rng = np.random.default_rng(seed)
    lo, hi = S.min(0), S.max(0)
    out = {}
    for k in ks:
        wk = np.log(KMeans(k, n_init=10, random_state=seed).fit(S).inertia_)
        refs = []
        for _ in range(n_ref):
            R = rng.uniform(lo, hi, size=S.shape)
            refs.append(np.log(KMeans(k, n_init=5, random_state=seed).fit(R).inertia_))
        out[k] = float(np.mean(refs) - wk)
    return out


def evaluate(X, label, rng):
    S, _ = prep(X)
    Sp, _ = prep(permuted(X, rng))

    sil = {k: silhouette_score(S, KMeans(k, n_init=10, random_state=SEED).fit_predict(S))
           for k in KS}
    sil_null = {k: silhouette_score(Sp, KMeans(k, n_init=10, random_state=SEED).fit_predict(Sp))
                for k in KS}
    gaps = gap_statistic(S)
    best_gap = max(gaps, key=gaps.get)
    gap_rising = gaps[2] < gaps[max(KS)]

    lab = hdbscan.HDBSCAN(min_cluster_size=25).fit_predict(S)
    n_hdb = int(len(set(lab)) - (1 if -1 in lab else 0))
    noise = float((lab == -1).mean())

    dips = [diptest.diptest(S[:, i])[1] for i in range(min(5, S.shape[1]))]

    best_k = max(sil, key=sil.get)
    print(f"\n{'=' * 74}\n{label}   n={X.shape[0]}  features={X.shape[1]}\n{'=' * 74}")
    print("  silhouette  " + "  ".join(f"k{k}:{sil[k]:.3f}" for k in KS))
    print("  null        " + "  ".join(f"k{k}:{sil_null[k]:.3f}" for k in KS))
    print(f"  best k by silhouette {best_k} ({sil[best_k]:.3f}), "
          f"null gives {sil_null[best_k]:.3f}, "
          f"excess {sil[best_k] - sil_null[best_k]:+.3f}")
    print(f"  gap statistic best k = {best_gap}"
          f"{'  (gap still rising at k=8, so k=1 is not excluded)' if gap_rising else ''}")
    print(f"  HDBSCAN clusters = {n_hdb}, {noise:.0%} unassigned")
    print(f"  dip test p on PC1-5: " + ", ".join(f"{p:.2f}" for p in dips)
          + ("   unimodal" if min(dips) > 0.05 else "   NON-UNIMODAL"))
    return dict(label=label, best_k=best_k, sil=sil[best_k], sil_null=sil_null[best_k],
                excess=sil[best_k] - sil_null[best_k], gap_k=best_gap,
                hdbscan=n_hdb, dip_min=min(dips), n=X.shape[0], p=X.shape[1])


# ---------------------------------------------------------------------------
s = get_store()
feat = s.table("features")
W = feat.pivot(index="gene_id", columns="feature", values="z")
gene_ids = W.index.to_numpy()
names = list(W.columns)
X_all = W.to_numpy(float)

rep = s.table("reproducibility_pairs")
rho = {}
for f, sub in rep.groupby("feature"):
    if len(sub) > 10 and sub["gw"].std() > 0 and sub["immune"].std() > 0:
        rho[f] = float(np.corrcoef(sub["gw"].rank(), sub["immune"].rank())[0, 1])
trusted = [f for f in names if rho.get(f, -1) > TRUST]
print(f"{len(names)} features, {len(rho)} measured twice, "
      f"{len(trusted)} reproduce above rho {TRUST}")

amount = W["total_mcc"].to_numpy(float)
a = (amount - amount.mean()) / amount.std()


def deamount(X):
    Z = X.copy()
    for j in range(Z.shape[1]):
        Z[:, j] -= (np.dot(a, Z[:, j]) / np.dot(a, a)) * a
    return Z


X_tr = W[trusted].to_numpy(float)
SUBSTRATES = [
    ("A  all 91 features", X_all),
    ("B  trusted only (rho > 0.7)", X_tr),
    ("C  all 91, amount removed", deamount(X_all)),
    ("D  trusted, amount removed", deamount(X_tr)),
]

rng = np.random.default_rng(SEED)
rows = [evaluate(X, lab, rng) for lab, X in SUBSTRATES]

# ---- conditional: within amount strata -----------------------------------
print(f"\n\n{'#' * 74}\nCONDITIONAL: clustering WITHIN amount tertiles\n"
      f"If amount dominates the geometry, shape clusters could exist inside a\n"
      f"stratum and be invisible globally. This is the test that has never run.\n{'#' * 74}")
q = np.quantile(amount, [1 / 3, 2 / 3])
strata = {"low amount": amount <= q[0],
          "mid amount": (amount > q[0]) & (amount <= q[1]),
          "high amount": amount > q[1]}
for sname, mask in strata.items():
    for lab, X in [("all 91", X_all), ("trusted", X_tr)]:
        rows.append(evaluate(X[mask], f"{sname}, {lab}", rng))

print(f"\n\n{'=' * 74}\nSUMMARY\n{'=' * 74}")
df = pd.DataFrame(rows)
print(df[["label", "n", "p", "best_k", "sil", "sil_null", "excess",
          "gap_k", "hdbscan", "dip_min"]].round(3).to_string(index=False))
print("\nA real cluster structure would show: silhouette well above its null,")
print("gap statistic peaking at k > 1, HDBSCAN returning >1 cluster, and a dip")
print("test rejecting unimodality. Read the four together, not one at a time.")
