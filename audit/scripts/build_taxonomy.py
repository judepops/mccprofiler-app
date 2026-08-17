"""Build the reach x composition taxonomy, as MEMBERSHIPS rather than labels.

METHOD NOTE. Leiden was tried first and rejected on stability: 0.703 mean
agreement over 50 seeds, so roughly a third of genes changed composition group
with the random seed. KMeans within each reach half reaches 0.868 on the same
test. Leiden remains the better tool for ASKING whether communities exist (it
returns a single community below resolution 0.4, which is the continuum result
arriving by a third independent method) but it is the worse tool for assigning
a stable label once the decision to impose one has been taken.

Phase 1 of audit/PLAN_TAXONOMY.md. Two levels, both named by architecture only:

    reach        RELATIVE extension           contained vs extended
                 NOT absolute bands. Measured 2026-08-17: the clustered reach
                 label matches a gene's dominant distance band for only 20% of
                 genes, BELOW chance for three bands, because "extended" means
                 "more extended than the panel average", not "the bulk of
                 contacts sits beyond 250 kb". Both halves are in fact
                 proximal-dominant (56.8% and 54.5% of signal inside 50 kb);
                 what separates them is signal beyond 250 kb, 13.8% against
                 6.0%. The old names mid-range/far-reaching asserted the
                 absolute reading and were wrong.
    composition  which element class dominates ctcf / enhancer / promoter

The levels are orthogonal: the composition split carries 1-4% of its variance on
sPC1, which is why no flat clustering found it. Reach dominates the distance
metric and swamps composition entirely.

WHY MEMBERSHIPS, NOT LABELS. Leiden returns a SINGLE community at every
resolution below 0.4, HDBSCAN returns zero clusters across sixteen conditions,
and every dip test is unimodal. These six regions are a resolution choice on a
continuum, so hard-assigning a gene to one of them asserts exactly the thing the
project spent two months disproving. Each gene therefore gets a membership
vector over the six regions, and its hard label is only ever the argmax, carried
with its own mixedness.

Membership is a softmax over negative squared distance to each region centroid
in the amount-corrected space, with the temperature set to the median
within-region spread so the scale is a property of the data rather than a knob.
A gene sitting between two regions gets a genuine blend; a gene at a centroid
gets a near-pure vector.

THE GO/NO-GO. Section 3 of the plan: the kappa 0.72 figure belongs to the old
KMeans k=3-4 on 91 features and does NOT transfer. This script re-measures
agreement across the 116 twice-captured genes for each level separately. If the
composition level does not reproduce, it must not ship, because reproducibility
of an imposed partition is the entire justification for imposing one.

Run in cd4env. Writes to audit/scripts/output/ only; nothing here touches the
store or the pipeline.
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
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED, KNN, RES = 0, 15, 0.6
N_SEEDS = 50
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")

# ---------------------------------------------------------------------------
# EXCLUDE UNDER-EVIDENCED GENES BEFORE FITTING, added 2026-08-17 (Jude's call).
#
# Nine genes have too few called peaks for the taxonomy to describe them, and
# leaving them in did active harm rather than merely adding noise.
#
# ASB9, PADI4 and TTPAL have ZERO peaks. Every per-class and element-class
# feature therefore sits at a default, so they land far from every centroid, and
# because the membership softmax squares the distance they came out with the
# THREE HIGHEST top weights in the whole panel (48.7%, 49.2%, 44.9% against a
# panel median of 28.3%). Six more genes have one or two peaks and behave the
# same way in miniature. The eight most confident-looking genes in the panel were
# all at the 96th to 100th percentile for centroid distance, and two of them had
# no peaks at all: confidence was reading as evidence when it measured the
# absence of it.
#
# This is not a quiet-locus problem. Their MCC signal is normal and their
# enrichment sits inside the range of genes that WERE called: TTPAL reaches O/E
# max 59.4 with 73 bins above O/E 5, against 20.8 to 46.5 for matched called
# genes. They are also absent from unfiltered.tsv, which rules out every
# downstream threshold. It is an annotation gap between the viewpoint list and
# the peak output. See audit/scripts/audit_missing_peaks.py.
#
# The threshold is config.ACTIVE_NPEAKS_MIN = 3, the pipeline's own existing
# floor for calling a gene active, rather than a number chosen here.
#
# They are dropped from the FIT and receive no region, rather than being
# projected onto the centroids afterwards. A projected label would be exactly the
# unsupported claim this removes: there is no architecture to compare.
MIN_PEAKS = 3
_pk = s.table("peaks")
_npk = _pk.groupby("symbol_key").size()
_sym = s.genes.set_index("gene_id")["symbol_key"]
_counts = W.index.to_series().map(_sym).map(_npk).fillna(0)
EXCLUDED = list(W.index[_counts.to_numpy() < MIN_PEAKS])
if EXCLUDED:
    _names = sorted(s.genes.set_index("gene_id").loc[EXCLUDED, "gene_symbol"])
    print(f"EXCLUDED {len(EXCLUDED)} under-evidenced genes "
          f"(<{MIN_PEAKS} called peaks): {_names}")
    W = W.drop(index=EXCLUDED)
print(f"fitting on {len(W)} genes")

genes = s.genes.set_index("gene_id").loc[W.index]
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
X = X - A @ beta
X = (X - X.mean(0)) / X.std(0)
S = PCA(20, random_state=SEED).fit_transform(X)
print(f"{X.shape[0]} genes, {X.shape[1]} amount-corrected features")


def leiden(M, res=RES, seed=SEED, knn=KNN):
    nn = NearestNeighbors(n_neighbors=min(knn + 1, len(M))).fit(M)
    _, idx = nn.kneighbors(M)
    e = set()
    for i, row in enumerate(idx):
        for j in row[1:]:
            e.add((min(i, j), max(i, j)))
    g = ig.Graph(n=len(M), edges=list(e))
    return np.asarray(la.find_partition(
        g, la.RBConfigurationVertexPartition,
        resolution_parameter=res, seed=seed, n_iterations=-1).membership)


ELEMENTS = ["ctcf", "enhancer", "promoter"]


def name_composition(D, sub):
    """Name each sub-community by the element class it is richest in.

    UNIQUE assignment, via the score matrix. Naming each cluster independently
    by its own argmax let two clusters claim the same element and silently
    merge: the first run returned 3 mid-range communities but only 2 names
    (promoter 546, enhancer 528), collapsing a real split into a label
    collision. Hungarian assignment forces a one-to-one map instead.

    Architecture only. The old pipeline named clusters by best Eisenberg Fisher
    p, which makes any later "this group is enriched for housekeeping genes"
    circular; external sets stay strictly post-hoc here.
    """
    from scipy.optimize import linear_sum_assignment

    ks = list(np.unique(sub))
    Sc = np.zeros((len(ks), len(ELEMENTS)))
    for i, c in enumerate(ks):
        m = sub == c
        for j, el in enumerate(ELEMENTS):
            cols = [f for f in D.columns
                    if f.endswith(f"_{el}") or f.startswith(f"{el}_")]
            if cols:
                Sc[i, j] = float((D.loc[m, cols].mean() - D.loc[~m, cols].mean()).mean())
    r, c = linear_sum_assignment(-Sc)
    return {ks[i]: (ELEMENTS[j], Sc[i, j]) for i, j in zip(r, c)}


# --- level 1: reach ---------------------------------------------------------
top = KMeans(2, n_init=25, random_state=SEED).fit_predict(X)
if S[top == 0, 0].mean() > S[top == 1, 0].mean():
    top = 1 - top                                    # 0 = mid-range, 1 = far
# contained / extended, not mid / far: these are RELATIVE positions on the
# reach axis, and naming them after absolute distance bands produced a label
# that matched a gene's dominant band only 20% of the time.
reach = np.where(top == 1, "extended", "contained")
print(f"\nreach: contained {np.sum(reach == 'contained')}, "
      f"extended {np.sum(reach == 'extended')}")

# --- level 2: composition, within each reach half ---------------------------
D = pd.DataFrame(X, columns=keep, index=W.index)
composition = np.empty(len(X), dtype=object)
# k CHOSEN PER HALF, not forced symmetric. Measured on this substrate:
#   mid-range   k=2 -> smallest 519, seed stability 0.998
#               k=3 -> smallest 4 (a pocket), stability 0.796
#   far-reaching k=3 -> smallest 238, stability 0.950
# So the taxonomy is asymmetric and that asymmetry is the finding: a
# CTCF-dominated contact profile only becomes a distinct group among
# far-reaching genes, which is what CTCF biology predicts, since CTCF loops are
# long-range structural contacts. Mid-range genes divide on promoter vs
# enhancer alone. Forcing a symmetric 3x2 grid manufactures a 4-gene group.
K_PER_HALF = {"contained": 2, "extended": 3}
for half in ("contained", "extended"):
    m = reach == half
    sub = KMeans(K_PER_HALF[half], n_init=50, random_state=SEED).fit_predict(X[m])
    named = name_composition(D[m], sub)
    lab = np.array([named[c][0] for c in sub])
    composition[m] = lab
    counts = pd.Series(lab).value_counts().to_dict()
    print(f"  {half:<4} -> {len(np.unique(sub))} communities: {counts}")

region = np.array([f"{r}-{c}" for r, c in zip(reach, composition)])
print("\nsix regions:")
print(pd.Series(region).value_counts().to_string())

# --- memberships ------------------------------------------------------------
regions = sorted(set(region))
C = np.array([X[region == r].mean(0) for r in regions])
d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
# Temperature from the NEAREST-centroid distance, not the within-region mean.
# The first version used the latter, which in 77 dimensions is large enough to
# flatten the softmax to near-uniform (median top weight 0.22 against 0.167 for
# six equal regions) and destroyed the very information the mixture is for. The
# distance-ratio test says the structure is really there: median
# nearest/second-nearest is 0.854 and 65% of genes sit below 0.9.
tau = float(np.median(np.sort(d2, 1)[:, 0]))
M = np.exp(-d2 / tau)
M = M / M.sum(1, keepdims=True)
top1 = M.max(1)
ent = -(M * np.log(M + 1e-12)).sum(1) / np.log(len(regions))
# TWO DECIMALS, NOT ZERO. `.0%` rounded 99.78% up to "100%" and printed that
# 2026-08-17, which is the same class of error as everything else found that day:
# a display choice turning "almost all" into "all". The claim licensing this
# taxonomy is that genes are blends rather than labels, so the four exceptions are
# named rather than rounded away.
_n_hi = int((top1 >= 0.5).sum())
print(f"\nmembership: median top weight {np.median(top1):.2f}, "
      f"{np.mean(top1 < 0.5):.2%} of genes below 0.5 (i.e. genuinely mixed)")
if _n_hi:
    _hi_names = genes["gene_symbol"].to_numpy()[np.argsort(-top1)[:_n_hi]]
    print(f"  {_n_hi} gene(s) at or above 0.5: "
          f"{', '.join(f'{n} {w:.3f}' for n, w in zip(_hi_names, np.sort(top1)[::-1][:_n_hi]))}")

# Distance to every centroid, kept alongside the weights so a gene's membership
# can be INTERROGATED rather than taken on trust. Showing 0.38 without showing
# what produced it is the same failure as showing an archetype name without its
# loadings, and that one went unnoticed for six weeks.
dist = pd.DataFrame(np.sqrt(d2), columns=[f"d_{r}" for r in regions], index=W.index)

tax = pd.DataFrame({
    "gene_id": W.index, "symbol_key": genes["symbol_key"].to_numpy(),
    "reach": reach, "composition": composition, "region": region,
    "top_weight": np.round(top1, 4), "mixedness": np.round(ent, 4),
})
for i, r in enumerate(regions):
    tax[f"w_{r}"] = np.round(M[:, i], 4)
for c in dist.columns:
    tax[c] = np.round(dist[c].to_numpy(), 3)
tax["tau"] = round(float(tau), 3)
tax["under_evidenced"] = False

# The excluded genes are written back as EXPLICIT rows with a null region rather
# than left absent. A missing row is indistinguishable from a join failure, and
# this store has already lost 11 genes once to a silent case-mismatch join (the
# Cxorfnn bug). An explicit row with region = NA and under_evidenced = True says
# "deliberately not described" in a way a downstream merge cannot misread, and it
# keeps the table at full panel coverage so `_merge_checked` still demands 100%.
if EXCLUDED:
    extra = pd.DataFrame({"gene_id": EXCLUDED})
    extra["symbol_key"] = extra["gene_id"].map(
        s.genes.set_index("gene_id")["symbol_key"])
    for c in tax.columns:
        if c not in extra.columns:
            extra[c] = np.nan
    extra["under_evidenced"] = True
    extra["tau"] = round(float(tau), 3)
    tax = pd.concat([tax, extra[tax.columns]], ignore_index=True)
    print(f"\n{len(EXCLUDED)} under-evidenced genes written with region = NA")

tax.to_csv(OUT / "taxonomy_labels.tsv", sep="\t", index=False)

# --- GO/NO-GO: does each level reproduce across captures? -------------------
rep = s.table("reproducibility_pairs")
gw = rep.pivot(index="symbol_key", columns="feature", values="gw")
im = rep.pivot(index="symbol_key", columns="feature", values="immune")
sh = [c for c in gw.columns if c in im.columns]
gw, im = gw[sh].dropna(), im[sh].dropna()
cm = gw.index.intersection(im.index)
G = ((gw.loc[cm] - gw.loc[cm].mean()) / gw.loc[cm].std()).to_numpy(float)
I = ((im.loc[cm] - im.loc[cm].mean()) / im.loc[cm].std()).to_numpy(float)

print("\n" + "=" * 74)
print(f"GO / NO-GO: reproducibility across {len(cm)} twice-captured genes")
print("=" * 74)
res = {}
for lvl, k in (("reach", 2), ("composition", 3)):
    a = KMeans(k, n_init=25, random_state=SEED).fit_predict(G)
    b = KMeans(k, n_init=25, random_state=SEED).fit_predict(I)
    res[lvl] = adjusted_rand_score(a, b)
    print(f"  {lvl:<12} ARI {res[lvl]:.3f}")
print("  (ARI: 0 = chance, 1 = identical. Above ~0.3 is meaningful agreement")
print("   for a partition of a continuum; below that it is not reproducible.)")

# --- stability across Leiden seeds -----------------------------------------
print("\nstability of the composition level across Leiden seeds:")
runs = []
for sd in range(N_SEEDS):
    comp = np.empty(len(X), dtype=object)
    for half in ("contained", "extended"):
        m = reach == half
        sub = KMeans(K_PER_HALF[half], n_init=10, random_state=sd).fit_predict(X[m])
        named = name_composition(D[m], sub)
        comp[m] = [named[c][0] for c in sub]
    runs.append(comp)
agree = np.mean([[np.mean(runs[i] == runs[j]) for j in range(N_SEEDS)]
                 for i in range(N_SEEDS)])
print(f"  mean pairwise agreement over {N_SEEDS} seeds: {agree:.3f}")

print("\n" + "=" * 74)
ok = res["reach"] > 0.3 and res["composition"] > 0.3 and agree > 0.8
print("VERDICT:", "SHIP" if ok else "DO NOT SHIP")
if not ok:
    if res["composition"] <= 0.3:
        print("  composition does not reproduce across captures; ship reach alone")
    if agree <= 0.8:
        print("  composition is seed-unstable; the communities are not robust")

# --- radar profiles: the "FIFA stats" comparison across regions -------------
# Six interpretable axes rather than 77 features, because the point of this
# plot is to be read at a glance. Values are panel percentiles of the
# amount-corrected feature, so 50 is the panel median and no axis can mean
# "this region has more signal".
RADAR = {
    "reach": "mean_distance_to_viewpoint_all",
    "promoter share": "promoter_signal_fraction",
    "enhancer share": "enhancer_signal_fraction_raw",
    "CTCF share": "ctcf_signal_fraction",
    "dispersion": "signal_entropy",
    "concentration": "frac_signal_in_top_peak",
}
rows = []
for axis, feat in RADAR.items():
    if feat not in D.columns:
        print(f"  radar: {feat} not in substrate, skipped")
        continue
    pct = D[feat].rank(pct=True) * 100
    for r in regions:
        m = region == r
        rows.append(dict(region=r, axis=axis, feature=feat,
                         median=round(float(pct[m].median()), 1),
                         q25=round(float(pct[m].quantile(0.25)), 1),
                         q75=round(float(pct[m].quantile(0.75)), 1),
                         n=int(m.sum())))
radar = pd.DataFrame(rows)
radar.to_csv(OUT / "taxonomy_radar.tsv", sep="\t", index=False)
print("\nradar profiles (panel percentile, 50 = median):")
print(radar.pivot(index="region", columns="axis", values="median").to_string())

print(f"\nwrote {OUT / 'taxonomy_labels.tsv'}")
print(f"wrote {OUT / 'taxonomy_radar.tsv'}")
