"""Is feature-set fragility a WEAKNESS of the partition or a STRENGTH of the features?

Jude's objection to 6g, 2026-08-17, and it is a fair one: if dropping any 12 of
77 features disturbs the partition, that may simply mean every feature carries
non-redundant information. Non-redundancy is a property people work hard to get.

The two readings are not symmetric, though, and three measurements separate them.

  1. DROP-SIZE SWEEP. If information is spread evenly across features, damage
     should scale smoothly with how many are removed. If even one feature
     collapses the partition, the cut is knife-edge and the features are not
     what is being measured.

  2. AGREEMENT, NOT JUST ARI. ARI is harsh on a k=2 split of a continuum. A
     percentage of genes keeping their label is the number a reader can weigh,
     and section 6g quoted only ARI, which overstates the damage.

  3. THE DECISIVE ONE: IS THE FULL FEATURE SET PRIVILEGED? Reproducibility
     across the 116 twice-captured genes is the only evidence licensing the
     regions. If the full substrate reproduces at 0.741 and random subsets
     reproduce far worse, the features ARE good and the chosen partition is the
     right one, so Jude is right. If random subsets reproduce just as well, then
     many equally reproducible partitions exist and there is no reason to prefer
     this one, which is what a continuum looks like.

Only test 3 can adjudicate. Tests 1 and 2 bound how bad the situation is.

Run in cd4env from the app root.
"""
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np, pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score as ari
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED, DRAWS = 0, 20
K_PER_HALF = {"contained": 2, "extended": 3}

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")


def build(cols, seed=SEED):
    Wc = W[cols]
    mag = [m for m in MAG_OVERALL if m in Wc.columns]
    keep = [c for c in Wc.columns if c not in mag]
    A = np.column_stack([np.ones(len(Wc)), Wc[mag].to_numpy(float)])
    X = Wc[keep].to_numpy(float)
    X = X - A @ np.linalg.lstsq(A, X, rcond=None)[0]
    X = (X - X.mean(0)) / X.std(0)
    S = PCA(20, random_state=seed).fit_transform(X)
    top = KMeans(2, n_init=25, random_state=seed).fit_predict(X)
    if S[top == 0, 0].mean() > S[top == 1, 0].mean():
        top = 1 - top
    reach = np.where(top == 0, "contained", "extended")
    comp = np.empty(len(X), dtype=object)
    for half in ("contained", "extended"):
        m = reach == half
        comp[m] = [f"{half}-{c}" for c in
                   KMeans(K_PER_HALF[half], n_init=50, random_state=seed).fit_predict(X[m])]
    return reach, comp


def agree(a, b):
    """Best achievable label agreement, maximised over label permutations.

    Raw equality is meaningless here: cluster ids are arbitrary, so an identical
    partition can score 0%. Without this the sweep would report nonsense.
    """
    from itertools import permutations
    ua, ub = list(pd.unique(a)), list(pd.unique(b))
    if len(ub) < len(ua):
        ua, ub = ub, ua
        a, b = b, a
    best = 0.0
    for p in permutations(ub, len(ua)):
        m = dict(zip(ua, p))
        best = max(best, float(np.mean([m[x] == y for x, y in zip(a, b)])))
    return best


full_reach, full_comp = build(list(W.columns))
rng = np.random.default_rng(SEED)

print("=" * 78)
print("TEST 1 + 2: does damage scale with how many features are dropped?")
print("=" * 78)
print(f"  {'dropped':>8}{'reach ARI':>22}{'reach agreement':>22}")
print(f"  {'':>8}{'median':>10}{'p10':>7}{'p90':>7}{'median':>12}{'p10':>7}{'p90':>7}")
sweep = {}
for k in (1, 3, 6, 12, 24, 38):
    A_, G_ = [], []
    for _ in range(DRAWS):
        drop = set(rng.choice(list(W.columns), size=k, replace=False))
        r, c = build([x for x in W.columns if x not in drop])
        A_.append(ari(full_reach, r))
        G_.append(agree(full_reach, r))
    A_, G_ = np.array(A_), np.array(G_)
    sweep[k] = (A_, G_)
    print(f"  {k:>8}{np.median(A_):>10.3f}{np.quantile(A_,.1):>7.3f}"
          f"{np.quantile(A_,.9):>7.3f}"
          f"{np.median(G_)*100:>11.1f}%{np.quantile(G_,.1)*100:>6.1f}%"
          f"{np.quantile(G_,.9)*100:>6.1f}%")
print("\n  Smooth decline => information is spread across features (Jude's reading).")
print("  A cliff at k=1 => the cut is knife-edge and the features are not the")
print("  thing being measured.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("TEST 3 (DECISIVE): is the FULL feature set privileged for reproducibility?")
print("=" * 78)
rep = s.table("reproducibility_pairs")
gw = rep.pivot(index="symbol_key", columns="feature", values="gw")
im = rep.pivot(index="symbol_key", columns="feature", values="immune")
sh = [c for c in gw.columns if c in im.columns]
gw, im = gw[sh].dropna(), im[sh].dropna()
cm = gw.index.intersection(im.index)
G = ((gw.loc[cm] - gw.loc[cm].mean()) / gw.loc[cm].std()).to_numpy(float)
I = ((im.loc[cm] - im.loc[cm].mean()) / im.loc[cm].std()).to_numpy(float)
cols = list(gw.columns)
print(f"  {len(cm)} twice-captured genes, {len(cols)} shared features")


def xcap(idx, k, seed=SEED):
    a = KMeans(k, n_init=25, random_state=seed).fit_predict(G[:, idx])
    b = KMeans(k, n_init=25, random_state=seed).fit_predict(I[:, idx])
    return ari(a, b)


allidx = np.arange(len(cols))
for lvl, k in (("reach", 2), ("composition", 3)):
    obs = xcap(allidx, k)
    null = []
    for _ in range(60):
        # Same SIZE as the ablation removed proportionally, so the comparison is
        # about WHICH features, not how many.
        idx = rng.choice(allidx, size=len(cols) - 10, replace=False)
        null.append(xcap(idx, k))
    null = np.array(null)
    p = float((null >= obs).mean())
    print(f"\n  {lvl} (k={k})")
    print(f"     full {len(cols)} features        ARI {obs:.3f}")
    print(f"     random {len(cols)-10}-feature subsets  median {np.median(null):.3f}"
          f"   p10 {np.quantile(null,.1):.3f}   p90 {np.quantile(null,.9):.3f}")
    print(f"     fraction of random subsets reproducing AT LEAST as well: {p:.2f}")
    print(f"     -> {'FULL SET IS NOT PRIVILEGED (continuum reading)' if p > 0.2 else 'full set is privileged (features reading)'}")
