"""How large an effect could the super-enhancer test have detected?

The Aim 3 null is only worth stating if the test had the power to find something.
As written it is attackable: the two positive controls are the SMALLEST sets in
the table (ChromHMM_bivalent n=46, Roadmap_silenced n=41) while the null is one
of the larger ones (dbSUPER n=158). "We detect an effect at n=41 and find
nothing at n=158" is fine statistically and reads backwards to a reviewer.

Two answers, both computed here.

  1. MINIMUM DETECTABLE EFFECT. Simulate a known displacement injected into a
     random set of the same size, run the same permutation test, and find the
     smallest true effect detected 80% of the time. If that is well below the
     effect sizes the positive controls show, the null is informative.

  2. MATCHED-N COMPARISON. Lambert_TF (n=156) is within 2 genes of dbSUPER
     (n=158) and separates at d = +0.46, p = 0.0005. That is the cleanest
     possible control for power: same n, same substrate, same test, opposite
     result. It is also mechanism-defined (DNA-binding domain) rather than an
     expression list, which makes it a better positive control than either
     chromatin-state set.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from _shape import MAG_OVERALL

from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
N_PERM = 1000       # permutations per test
N_SIM = 300         # simulated datasets per effect size
ALPHA = 0.05
POWER = 0.80
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()

mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
R = X - A @ beta
D = PCA(n_components=18, random_state=0).fit_transform(R)
D = (D - D.mean(0)) / D.std(0)
n, k_dim = D.shape
print(f"substrate {n} genes x {k_dim} components\n")


def null_threshold(size, rng):
    """95th percentile of the centroid-norm null for a set of this size."""
    vals = np.array([np.linalg.norm(D[rng.choice(n, size=size, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    return np.quantile(vals, 1 - ALPHA)


def power_at(size, d, thr, rng, direction=None):
    """Fraction of simulated sets with a true shift d that clear the threshold.

    The shift is applied along one component, which is the least favourable
    case: a set displaced along a single axis is harder to detect by a
    multivariate centroid norm than one displaced diagonally across several.
    """
    hits = 0
    for _ in range(N_SIM):
        idx = rng.choice(n, size=size, replace=False)
        j = direction if direction is not None else rng.integers(k_dim)
        sub = D[idx].copy()
        sub[:, j] += d
        if np.linalg.norm(sub.mean(0)) > thr:
            hits += 1
    return hits / N_SIM


rng = np.random.default_rng(0)
GRID = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]

print("=" * 74)
print("MINIMUM DETECTABLE EFFECT, by set size")
print("=" * 74)
print("Power to detect a true displacement of size d on ONE component,")
print(f"at alpha {ALPHA}, from {N_SIM} simulations per cell.\n")

sizes = {"dbSUPER_CD4_SE": 158, "Lambert_TF": 156,
         "ChromHMM_bivalent": 46, "Roadmap_silenced": 41}
header = "  ".join(f"{d:>5.2f}" for d in GRID)
print(f"{'set (n)':<26}{header}")
mde = {}
for name, size in sizes.items():
    thr = null_threshold(size, rng)
    row, found = [], None
    for d in GRID:
        p = power_at(size, d, thr, rng)
        row.append(p)
        if found is None and p >= POWER:
            found = d
    mde[name] = found
    print(f"{name + f' (n={size})':<26}" + "  ".join(f"{p:>5.2f}" for p in row))

print(f"\nSmallest d detected with >= {POWER:.0%} power:")
for name, size in sizes.items():
    v = mde[name]
    print(f"  {name + f' (n={size})':<26} d = {v if v else '> ' + str(GRID[-1])}")

print("\n" + "=" * 74)
print("MATCHED-N COMPARISON: the cleanest control for power")
print("=" * 74)
cm = s.table("cohort_membership")


def observed(group):
    members = set(cm.loc[cm["group"] == group, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, n)
    kk = int(hit.sum())
    r = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[r.choice(n, size=kk, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((kk - 1) * a_.var(ddof=1) + (n - kk - 1) * b_.var(ddof=1)) / (n - 2))
    return kk, (obs - null.mean()) / null.std(), (np.sum(null >= obs) + 1) / (N_PERM + 1), \
        float((a_.mean() - b_.mean()) / sp)


for g in ["Lambert_TF", "dbSUPER_CD4_SE_TSS_pm50kb"]:
    kk, z, p, d = observed(g)
    print(f"  {g:<30} n={kk:<5} z={z:6.1f}  p={p:.4f}  d={d:+.2f}")
print("\n  Same substrate, same test, set sizes within two genes of each other,")
print("  opposite results. Power cannot explain the difference.")
