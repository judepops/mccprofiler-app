"""Does the composition axis rest on features that carry contact decay?

Asked 2026-08-17. The app's display grid weights peaks by `peak_max`, which
falls 2.9-fold from the proximal to the far band (rho -0.555 with distance), and
element class is confounded with distance in this panel (median offset: enhancer
57.8 kb, promoter 92.8 kb, CTCF 112.5 kb). So a "70% enhancer" gene is partly
reporting that its peaks are close.

That is a display problem. The question this script answers is whether the same
confound is IN THE CLUSTERING: the substrate contains raw per-class signal
fractions alongside O/E ones, and if the composition split rests on the raw
family then "enhancer-dominated" partly means "contacts are close", which would
undercut the claim that reach and composition are orthogonal levels.

Method: rebuild the taxonomy exactly as build_taxonomy.py does, once with the
full substrate and once with the decay-carrying class features dropped, and
compare the resulting partitions by ARI. High ARI = the split does not depend on
them. Reach is checked too, as a control: reach SHOULD depend on distance.

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
SEED = 0
K_PER_HALF = {"contained": 2, "extended": 3}

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")

# The decay-carrying class features: per-class shares of RAW signal. Their O/E
# counterparts are distance-corrected and stay in.
DECAY = [c for c in W.columns
         if ("signal_fraction" in c and c.endswith("_raw"))
         or c.startswith("raw_log2_enrichment")
         or (c.startswith("raw_peak_max"))]
print(f"{W.shape[1]} features in the substrate")
print(f"{len(DECAY)} decay-carrying class features dropped in the ablation:")
for c in sorted(DECAY):
    print(f"    {c}")


def build(cols):
    """build_taxonomy.py's pipeline, verbatim, on a chosen column set."""
    Wc = W[cols]
    mag = [m for m in MAG_OVERALL if m in Wc.columns]
    keep = [c for c in Wc.columns if c not in mag]
    A = np.column_stack([np.ones(len(Wc)), Wc[mag].to_numpy(float)])
    X = Wc[keep].to_numpy(float)
    X = X - A @ np.linalg.lstsq(A, X, rcond=None)[0]
    X = (X - X.mean(0)) / X.std(0)
    S = PCA(20, random_state=SEED).fit_transform(X)
    top = KMeans(2, n_init=25, random_state=SEED).fit_predict(X)
    if S[top == 0, 0].mean() > S[top == 1, 0].mean():
        top = 1 - top
    reach = np.where(top == 0, "contained", "extended")
    comp = np.empty(len(X), dtype=object)
    for half in ("contained", "extended"):
        m = reach == half
        comp[m] = [f"{half}-{c}" for c in
                   KMeans(K_PER_HALF[half], n_init=50,
                          random_state=SEED).fit_predict(X[m])]
    return reach, comp, X.shape[1]


full_reach, full_comp, nf = build(list(W.columns))
abl_reach, abl_comp, na = build([c for c in W.columns if c not in DECAY])
print(f"\nfull substrate: {nf} amount-corrected features")
print(f"ablated:        {na}")

print("\n" + "=" * 72)
print("DOES DROPPING THE DECAY-CARRYING CLASS FEATURES CHANGE THE PARTITION?")
print("=" * 72)
print(f"  reach        ARI = {ari(full_reach, abl_reach):.3f}"
      "   (SHOULD be high: reach is meant to be about distance,")
print("                                and these features are not the reach signal)")
print(f"  composition  ARI = {ari(full_comp, abl_comp):.3f}"
      "   (if LOW, composition rests on decay-carrying features)")

# Within-half composition, which is the claim actually at stake.
for half in ("contained", "extended"):
    m = (full_reach == half) & (abl_reach == half)
    if m.sum() > 10:
        print(f"    {half:<10} (n={m.sum():>4}, genes in this half both ways) "
              f"ARI = {ari(full_comp[m], abl_comp[m]):.3f}")

print("\n  reach agreement, gene by gene: "
      f"{(full_reach == abl_reach).mean() * 100:.1f}%")

# ---------------------------------------------------------------------------
# THE CONTROL, without which the ablation above means nothing.
#
# Dropping 12 of 77 features changes the substrate by 16%, and KMeans on 66
# columns can disagree with KMeans on 77 for reasons that have nothing to do
# with decay. So: drop 12 RANDOM non-decay features, many times, and see what
# ARI that produces. If the random null sits at the same place as the decay
# ablation, the finding is "composition is fragile to any feature removal", not
# "composition rests on decay-carrying features". Those call for completely
# different fixes.
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("CONTROL: dropping 12 RANDOM non-decay features instead")
print("=" * 72)
pool = [c for c in W.columns if c not in DECAY]
rng = np.random.default_rng(SEED)
rows = []
for i in range(20):
    drop = set(rng.choice(pool, size=len(DECAY), replace=False))
    r, c, _ = build([x for x in W.columns if x not in drop])
    m_c = (full_reach == "contained") & (r == "contained")
    m_e = (full_reach == "extended") & (r == "extended")
    rows.append((ari(full_reach, r), ari(full_comp, c),
                 ari(full_comp[m_c], c[m_c]) if m_c.sum() > 10 else np.nan,
                 ari(full_comp[m_e], c[m_e]) if m_e.sum() > 10 else np.nan))
R = np.array(rows, float)
names = ["reach", "composition", "  contained", "  extended"]
print(f"  {'':<14}{'decay ablation':>16}{'random null':>26}")
print(f"  {'':<14}{'':>16}{'median':>10}{'p10':>8}{'p90':>8}")
obs = [ari(full_reach, abl_reach), ari(full_comp, abl_comp),
       ari(full_comp[(full_reach == 'contained') & (abl_reach == 'contained')],
           abl_comp[(full_reach == 'contained') & (abl_reach == 'contained')]),
       ari(full_comp[(full_reach == 'extended') & (abl_reach == 'extended')],
           abl_comp[(full_reach == 'extended') & (abl_reach == 'extended')])]
for k, nm in enumerate(names):
    col = R[:, k]
    col = col[np.isfinite(col)]
    p = (col <= obs[k]).mean()
    print(f"  {nm:<14}{obs[k]:>16.3f}{np.median(col):>10.3f}"
          f"{np.quantile(col,0.1):>8.3f}{np.quantile(col,0.9):>8.3f}"
          f"   {'DECAY IS SPECIAL' if p <= 0.1 else 'within the null'} (p={p:.2f})")
print("\n  Read the p column: it is the fraction of random 12-feature drops that")
print("  disturb the partition AS MUCH AS dropping the decay features. A high p")
print("  means the decay features are not special and the axis is simply")
print("  fragile; a low p means the axis genuinely depends on them.")
