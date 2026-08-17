"""Leiden community detection on the contact-architecture substrate.

Replaces the KMeans sweep in `choose_k.py`, for three reasons that matter on a
continuum:

  * KMeans assumes spherical, equal-variance clusters. Nothing here is spherical.
  * KMeans forces every gene into a group, so at k >= 4 it spent its clusters
    isolating 3-gene outlier pockets rather than describing the panel.
  * KMeans takes k as given. Leiden takes a RESOLUTION, and lets the number of
    communities fall out of the graph, which is the more honest parameterisation
    when the number is not a real property of the data.

WHAT THE PREVIOUS SWEEP GOT WRONG, recorded so it is not repeated. It picked k
by maximising cross-capture kappa, which selected k=2 at kappa 0.859. That
partition puts 97% of its between-centroid variance on sPC1: it is one axis with
a threshold on it, and it reproduces beautifully precisely BECAUSE it is
trivial. Kappa measures stability, not information. A criterion that rewards
thresholding the dominant axis will always do so.

So a partition is judged on three things here, and must pass all of them:

  REPRODUCIBLE   ARI between the communities the same 116 twice-captured genes
                 receive from their gw features and from their immune features.
                 ARI rather than kappa because Leiden need not return the same
                 number of communities on both panels, and ARI handles that.
  MULTIVARIATE   share of between-centroid variance carried by sPC1 alone. Above
                 ~0.75 the partition is a dichotomised axis and adds nothing over
                 simply reporting the axis.
  BALANCED       smallest community as a share of the panel. A 3-gene community
                 is an outlier pocket, not a group.

Also runs the HIERARCHICAL alternative: split on the dominant reach axis first,
then find communities within each half separately. If the panel's structure is
"reach, and then something else within each regime", that will show up as
sub-communities that are NOT themselves sPC1 splits.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED, KNN = 0, 15
RES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
genes = s.genes.set_index("gene_id").loc[W.index]
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
X = X - A @ beta
X = (X - X.mean(0)) / X.std(0)
S = PCA(20, random_state=SEED).fit_transform(X)
print(f"{X.shape[0]} genes, {X.shape[1]} amount-corrected features\n")


def knn_graph(M, k=KNN):
    """Shared-neighbour weighted kNN graph, the standard Leiden substrate."""
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(M))).fit(M)
    _, idx = nn.kneighbors(M)
    edges, seen = [], set()
    for i, row in enumerate(idx):
        for j in row[1:]:
            e = (min(i, j), max(i, j))
            if e not in seen:
                seen.add(e)
                edges.append(e)
    g = ig.Graph(n=len(M), edges=edges)
    return g


def leiden(M, res, seed=SEED):
    g = knn_graph(M)
    part = la.find_partition(g, la.RBConfigurationVertexPartition,
                             resolution_parameter=res, seed=seed, n_iterations=-1)
    return np.asarray(part.membership)


def spc1_share(lab, scores):
    ks = np.unique(lab)
    C = np.array([scores[lab == c].mean(0) for c in ks])
    w = np.array([(lab == c).sum() for c in ks])[:, None]
    v = (w * (C - scores.mean(0)) ** 2).sum(0)
    return float(v[0] / v.sum()), int(np.searchsorted(np.cumsum(v / v.sum()), 0.90) + 1)


rep = s.table("reproducibility_pairs")
gw = rep.pivot(index="symbol_key", columns="feature", values="gw")
im = rep.pivot(index="symbol_key", columns="feature", values="immune")
sh = [c for c in gw.columns if c in im.columns]
gw, im = gw[sh].dropna(), im[sh].dropna()
cm = gw.index.intersection(im.index)
G = ((gw.loc[cm] - gw.loc[cm].mean()) / gw.loc[cm].std()).to_numpy(float)
I = ((im.loc[cm] - im.loc[cm].mean()) / im.loc[cm].std()).to_numpy(float)

print("=" * 92)
print("LEIDEN RESOLUTION SWEEP")
print("=" * 92)
print(f"{'res':>6}{'n comm':>8}{'ARI':>8}{'sPC1':>7}{'ax90':>6}{'smallest':>10}"
      f"{'sil':>7}   verdict")
rows = []
for r in RES:
    lab = leiden(X, r)
    n = len(np.unique(lab))
    if n < 2:
        print(f"{r:>6}{n:>8}   single community")
        continue
    share, ax90 = spc1_share(lab, S)
    sm = int(np.bincount(lab).min())
    ari = adjusted_rand_score(leiden(G, r), leiden(I, r))
    sil = silhouette_score(X, lab) if n > 1 else np.nan
    v = ("1-D, uninformative" if share > 0.75 else
         "outlier pockets" if sm < 0.05 * len(X) else "multivariate + balanced")
    rows.append(dict(res=r, n=n, ari=ari, spc1=share, ax90=ax90, smallest=sm,
                     sil=sil, verdict=v))
    print(f"{r:>6}{n:>8}{ari:>8.3f}{share:>7.0%}{ax90:>6}{sm:>10}{sil:>7.3f}   {v}")

d = pd.DataFrame(rows)
out = Path(__file__).resolve().parent / "output" / "choose_k_leiden.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(out, sep="\t", index=False)

print("\n" + "=" * 92)
print("HIERARCHICAL: split on the reach axis first, then find communities inside each half")
print("=" * 92)
print("If the structure is 'reach, and then something else within each regime',")
print("the sub-communities should NOT themselves be sPC1 splits.\n")
top = KMeans(2, n_init=25, random_state=SEED).fit_predict(X)
for half in (0, 1):
    m = top == half
    lo, hi = S[m, 0].min(), S[m, 0].max()
    side = "far-reaching" if S[m, 0].mean() > S[:, 0].mean() else "mid-range"
    print(f"  half {half} ({side}), n = {m.sum()}")
    for r in (0.2, 0.4, 0.6):
        sub = leiden(X[m], r)
        n = len(np.unique(sub))
        if n < 2:
            print(f"     res {r}: 1 community")
            continue
        sh_, ax_ = spc1_share(sub, S[m])
        sm = int(np.bincount(sub).min())
        print(f"     res {r}: {n} communities, sPC1 share {sh_:.0%}, "
              f"axes to 90% {ax_}, smallest {sm}"
              + ("   <- still a reach split" if sh_ > 0.6 else "   <- new structure"))
    print()

if len(d):
    ok = d[(d.spc1 <= 0.75) & (d.smallest >= 0.05 * len(X))]
    print("=" * 92)
    if len(ok):
        b = ok.loc[ok.ari.idxmax()]
        print(f"BEST: resolution {b.res}, {int(b.n)} communities, ARI {b.ari:.3f}, "
              f"sPC1 {b.spc1:.0%}, smallest {int(b.smallest)}")
        print("Passes all three: reproducible, multivariate, balanced.")
    else:
        print("NO resolution is simultaneously multivariate and balanced.")
        print("That is itself the result: the graph has no community structure")
        print("beyond the dominant axis, which is what a continuum looks like.")
print(f"\nwrote {out}")
