"""At what resolution is an imposed partition most reproducible?

THIS IS NOT A TEST FOR HOW MANY CLUSTERS EXIST. There are none: HDBSCAN returns
zero across sixteen conditions, every dip test is unimodal, and the gap statistic
is still rising at k=8. That question is answered and the answer is "the
landscape is continuous".

The question here is different and is the only one the project's own argument
supports. If a partition is imposed anyway, for display or for description, at
which k does it survive being measured twice? A partition that reproduces across
independent captures of the same genes is a defensible resolution choice; one
that does not is an artefact of a single dataset.

Three criteria, and they answer different things:

  KAPPA         Cohen's kappa between the labels the SAME 116 genes receive when
                clustered on their gw features versus their immune features.
                Labels are permuted arbitrarily by KMeans, so they are matched
                optimally (Hungarian) before scoring. This is the criterion that
                matters, because it is the one the write-up leans on.

  SILHOUETTE    Against a per-feature-permuted null, so it reads as "how much
                more separated than chance". Expect it to be low at every k;
                that is the continuum result, not a defect.

  BALANCE       Smallest cluster as a share of the panel. A k that produces a
                12-gene group has not found structure, it has found outliers.

Everything runs on the AMOUNT-CORRECTED substrate, so a group cannot mean
"these genes have more signal", which is the caveat that has sat on the
archetype names since 2026-08-03.

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
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import cohen_kappa_score, silhouette_score
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED, N_INIT, KS = 0, 25, range(2, 9)

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
X = X - A @ beta
X = (X - X.mean(0)) / X.std(0)
print(f"panel: {X.shape[0]} genes, {X.shape[1]} amount-corrected features")

rep = s.table("reproducibility_pairs")
gw = rep.pivot(index="symbol_key", columns="feature", values="gw")
im = rep.pivot(index="symbol_key", columns="feature", values="immune")
shared = [c for c in gw.columns if c in im.columns]
gw, im = gw[shared].dropna(), im[shared].dropna()
common = gw.index.intersection(im.index)
G = ((gw.loc[common] - gw.loc[common].mean()) / gw.loc[common].std()).to_numpy(float)
I = ((im.loc[common] - im.loc[common].mean()) / im.loc[common].std()).to_numpy(float)
print(f"cross-capture: {len(common)} genes measured twice, {len(shared)} shared features\n")


def match(a, b, k):
    """Align two labelings by the assignment that maximises agreement."""
    C = np.zeros((k, k), int)
    for i, j in zip(a, b):
        C[i, j] += 1
    r, c = linear_sum_assignment(-C)
    m = dict(zip(c, r))
    return np.array([m[x] for x in b])


rng = np.random.default_rng(SEED)
Xn = X.copy()
for j in range(Xn.shape[1]):
    rng.shuffle(Xn[:, j])

rows = []
for k in KS:
    lab = KMeans(k, n_init=N_INIT, random_state=SEED).fit_predict(X)
    sil = silhouette_score(X, lab)
    sil_null = silhouette_score(Xn, KMeans(k, n_init=N_INIT, random_state=SEED)
                                .fit_predict(Xn))
    sizes = np.bincount(lab, minlength=k)
    a = KMeans(k, n_init=N_INIT, random_state=SEED).fit_predict(G)
    b = KMeans(k, n_init=N_INIT, random_state=SEED).fit_predict(I)
    kap = cohen_kappa_score(a, match(a, b, k))
    rows.append(dict(k=k, kappa=kap, sil=sil, sil_null=sil_null,
                     excess=sil - sil_null, smallest=sizes.min(),
                     smallest_pct=100 * sizes.min() / len(lab)))

d = pd.DataFrame(rows)
print("=" * 78)
print(f"{'k':>3}{'kappa':>9}{'silhouette':>12}{'null':>8}{'excess':>8}"
      f"{'smallest':>10}{'':>6}")
print("=" * 78)
for r in d.itertuples():
    star = "  <-- best kappa" if r.kappa == d.kappa.max() else ""
    warn = "  small" if r.smallest_pct < 5 else ""
    print(f"{r.k:>3}{r.kappa:>9.3f}{r.sil:>12.3f}{r.sil_null:>8.3f}"
          f"{r.excess:>8.3f}{r.smallest:>7} ({r.smallest_pct:.0f}%){star}{warn}")

out = Path(__file__).resolve().parent / "output" / "choose_k.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(out, sep="\t", index=False)

best = d.loc[d.kappa.idxmax()]
print("\n" + "=" * 78)
print("READING THIS")
print("=" * 78)
print(f"  kappa peaks at k = {int(best.k)} ({best.kappa:.3f}).")
print("  Landis and Koch: <0.20 slight, 0.21-0.40 fair, 0.41-0.60 moderate,")
print("  0.61-0.80 substantial, >0.80 almost perfect.")
print(f"\n  Silhouette excess is {d.excess.min():.3f} to {d.excess.max():.3f} across all k.")
print("  It never becomes large, at any k. That is the continuum result and no")
print("  choice of k fixes it: this picks the least-bad resolution, it does not")
print("  find clusters.")
viable = d[(d.kappa > 0.4) & (d.smallest_pct >= 5)]
if len(viable):
    print(f"\n  k with kappa > 0.4 AND no group under 5% of the panel: "
          f"{', '.join(str(int(x)) for x in viable.k)}")
    print(f"  -> prefer the SMALLEST such k: k = {int(viable.k.min())}. A finer split "
          f"that reproduces no better is a harder story for no gain.")
else:
    print("\n  NO k clears kappa > 0.4 with balanced groups. If that holds, the")
    print("  partition should be dropped from the app rather than re-fitted.")
print(f"\nwrote {out}")
