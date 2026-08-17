"""Item 5: the robustness sweep, done properly.

The 2026-08-17 first pass used 20 draws at ONE seed each, which cannot separate
two different things: sensitivity to the SUBSTRATE (which features are present)
from residual sensitivity to the OPTIMISER (which KMeans init won). Reported
figures were therefore provisional and must not be quoted.

This version varies the seed WITHIN each draw, so the two are separable:

  within-draw seed spread   = optimiser sensitivity at a fixed feature set
  between-draw spread       = substrate sensitivity, net of the above

It also runs a full drop-size sweep and a seed-only control (drop nothing, vary
the seed), which is the floor any substrate effect has to clear.

Run in cd4env from the app root. Writes audit/scripts/output/robustness_sweep.tsv.
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
BASE_SEED = 0
DRAWS = 25            # feature subsets per drop size
SEEDS = 5             # KMeans seeds per subset
DROPS = (0, 1, 3, 6, 12, 24, 38)
K_PER_HALF = {"contained": 2, "extended": 3}
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")


def build(cols, seed):
    Wc = W[cols]
    mag = [m for m in MAG_OVERALL if m in Wc.columns]
    keep = [c for c in Wc.columns if c not in mag]
    A = np.column_stack([np.ones(len(Wc)), Wc[mag].to_numpy(float)])
    X = Wc[keep].to_numpy(float)
    X = X - A @ np.linalg.lstsq(A, X, rcond=None)[0]
    X = (X - X.mean(0)) / X.std(0)
    S = PCA(20, random_state=BASE_SEED).fit_transform(X)
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


ref_reach, ref_comp = build(list(W.columns), BASE_SEED)
rng = np.random.default_rng(BASE_SEED)
rows = []
for k in DROPS:
    for draw in range(1 if k == 0 else DRAWS):
        drop = set() if k == 0 else set(rng.choice(list(W.columns), size=k, replace=False))
        cols = [c for c in W.columns if c not in drop]
        for sd in range(SEEDS):
            r, c = build(cols, BASE_SEED + sd)
            rows.append({"n_dropped": k, "draw": draw, "seed": sd,
                         "reach_ari": ari(ref_reach, r),
                         "comp_ari": ari(ref_comp, c)})
    print(f"  drop {k:>2}: done", flush=True)

d = pd.DataFrame(rows)
d.to_csv(OUT / "robustness_sweep.tsv", sep="\t", index=False)

print("\n" + "=" * 78)
print("SEED-ONLY FLOOR: vary the optimiser, hold the substrate fixed")
print("=" * 78)
z = d[d.n_dropped == 0]
print(f"  reach ARI       median {z.reach_ari.median():.3f}"
      f"   min {z.reach_ari.min():.3f}   max {z.reach_ari.max():.3f}")
print(f"  composition ARI median {z.comp_ari.median():.3f}"
      f"   min {z.comp_ari.min():.3f}   max {z.comp_ari.max():.3f}")
print("  Any substrate effect must clear THIS to be real.")

print("\n" + "=" * 78)
print("SUBSTRATE vs OPTIMISER, separated")
print("=" * 78)
print(f"  {'drop':>5}{'reach ARI':>12}{'within-draw':>14}{'between-draw':>14}"
      f"{'comp ARI':>11}")
print(f"  {'':>5}{'median':>12}{'(seed) sd':>14}{'(feature) sd':>14}{'median':>11}")
for k in DROPS:
    g = d[d.n_dropped == k]
    if g.empty:
        continue
    within = g.groupby("draw")["reach_ari"].std().mean()
    between = g.groupby("draw")["reach_ari"].mean().std()
    print(f"  {k:>5}{g.reach_ari.median():>12.3f}"
          f"{(0.0 if np.isnan(within) else within):>14.3f}"
          f"{(0.0 if np.isnan(between) else between):>14.3f}"
          f"{g.comp_ari.median():>11.3f}")
print("\n  If between-draw sd exceeds within-draw sd, the sensitivity is to the")
print("  FEATURE SET and not to the optimiser. If they are comparable, the")
print("  original 'unstable to substrate' claim has no support at all.")
print(f"\nwrote {OUT / 'robustness_sweep.tsv'}  ({len(d)} fits)")
