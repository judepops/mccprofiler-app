"""Item 4: does composition reproduce as well WITHOUT the decay-carrying features?

Established 2026-08-17: dropping the 12 raw-signal class features destroys the
composition split within the contained half (ARI 0.000 against a random-drop
null at p=0.05) while leaving reach intact (0.812). So "enhancer-dominated" is
partly "contacts are close", because element class is confounded with distance
in this panel (median offset: enhancer 57.8 kb, promoter 92.8 kb, CTCF 112.5 kb)
and raw peak height decays 2.9-fold across the window.

That tells us the current composition axis USES those features. It does not tell
us whether a decay-corrected composition axis would be BETTER, and that is the
only question that matters for whether to change anything.

The test that decides it is cross-capture reproducibility on the 116
twice-captured genes, because that is the sole evidence licensing the partition.
If O/E-only composition reproduces as well or better, the current axis is buying
its structure with a confound and should be replaced. If it reproduces worse, the
raw features carry real signal that the O/E versions do not, and the confound is
a limitation to document rather than a defect to fix.

NOTHING IS REPLACED HERE. This builds the alternative and compares it.

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
SEED, SEEDS = 0, 10
K_PER_HALF = {"contained": 2, "extended": 3}
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")


def is_decay(c: str) -> bool:
    """Class or intensity features computed from RAW peak height.

    Raw height falls 2.9-fold across the window (rho -0.555 with distance), so
    any per-class quantity built on it inherits the distance confound. The O/E
    counterparts divide the decay out and stay in.
    """
    return (("signal_fraction" in c and c.endswith("_raw"))
            or c.startswith("raw_log2_enrichment")
            or c.startswith("raw_peak_max"))


DECAY = [c for c in W.columns if is_decay(c)]
print(f"substrate {W.shape[1]} features, {len(DECAY)} decay-carrying:")
for c in sorted(DECAY):
    print(f"    {c}")


def fit(cols, seed=SEED):
    Wc = W[cols]
    mag = [m for m in MAG_OVERALL if m in Wc.columns]
    keep = [c for c in Wc.columns if c not in mag]
    A = np.column_stack([np.ones(len(Wc)), Wc[mag].to_numpy(float)])
    X = Wc[keep].to_numpy(float)
    X = X - A @ np.linalg.lstsq(A, X, rcond=None)[0]
    X = (X - X.mean(0)) / X.std(0)
    S = PCA(20, random_state=SEED).fit_transform(X)
    top = KMeans(2, n_init=25, random_state=seed).fit_predict(X)
    if S[top == 0, 0].mean() > S[top == 1, 0].mean():
        top = 1 - top
    reach = np.where(top == 0, "contained", "extended")
    comp = np.empty(len(X), dtype=object)
    for half in ("contained", "extended"):
        m = reach == half
        comp[m] = [f"{half}-{c}" for c in
                   KMeans(K_PER_HALF[half], n_init=50, random_state=seed).fit_predict(X[m])]
    return reach, comp, X.shape[1]


ALL = list(W.columns)
OE = [c for c in W.columns if not is_decay(c)]
r_all, c_all, n_all = fit(ALL)
r_oe, c_oe, n_oe = fit(OE)
print(f"\nfull      {n_all} amount-corrected features")
print(f"O/E only  {n_oe}")

# --- the decisive comparison: cross-capture reproducibility ----------------
rep = s.table("reproducibility_pairs")
gw = rep.pivot(index="symbol_key", columns="feature", values="gw")
im = rep.pivot(index="symbol_key", columns="feature", values="immune")
sh = [c for c in gw.columns if c in im.columns]
cm = gw.index.intersection(im.index)
G_all = gw.loc[cm, sh]
I_all = im.loc[cm, sh]
sh_oe = [c for c in sh if not is_decay(c)]
print(f"\n{len(cm)} twice-captured genes; {len(sh)} shared features, "
      f"{len(sh) - len(sh_oe)} of them decay-carrying")


def xcap(Gd, Id, k, seeds=SEEDS):
    G = ((Gd - Gd.mean()) / Gd.std()).to_numpy(float)
    I = ((Id - Id.mean()) / Id.std()).to_numpy(float)
    out = []
    for sd in range(seeds):
        a = KMeans(k, n_init=25, random_state=sd).fit_predict(G)
        b = KMeans(k, n_init=25, random_state=sd).fit_predict(I)
        out.append(ari(a, b))
    return np.array(out)


print("\n" + "=" * 76)
print("CROSS-CAPTURE REPRODUCIBILITY, full shared set vs O/E-only")
print("=" * 76)
print(f"  {'level':<14}{'full':>22}{'O/E only':>22}")
print(f"  {'':<14}{'median (min-max)':>22}{'median (min-max)':>22}")
res = {}
for lvl, k in (("reach", 2), ("composition", 3)):
    a = xcap(G_all, I_all, k)
    b = xcap(G_all[sh_oe], I_all[sh_oe], k)
    res[lvl] = (a, b)
    print(f"  {lvl:<14}{np.median(a):>10.3f} ({a.min():.3f}-{a.max():.3f})"
          f"{np.median(b):>10.3f} ({b.min():.3f}-{b.max():.3f})")

d_comp = np.median(res["composition"][1]) - np.median(res["composition"][0])
print(f"\n  composition change: {d_comp:+.3f}")
print("  Positive => the decay-corrected axis reproduces BETTER and the current")
print("  one is buying structure with a confound; replace it.")
print("  Negative => the raw features carry signal the O/E versions do not;")
print("  document the confound as a limitation instead.")

# --- what the alternative partition looks like -----------------------------
print("\n" + "=" * 76)
print("WHAT WOULD CHANGE ON THE PANEL")
print("=" * 76)
print(f"  reach agreement       {(r_all == r_oe).mean() * 100:5.1f}%"
      f"   ARI {ari(r_all, r_oe):.3f}")
print(f"  composition           ARI {ari(c_all, c_oe):.3f}")
print(f"\n  current region sizes: {pd.Series(c_all).value_counts().to_dict()}")
print(f"  O/E-only sizes:       {pd.Series(c_oe).value_counts().to_dict()}")

# Seed stability of the alternative, which has to hold for it to be usable.
runs = [fit(OE, seed=sd)[1] for sd in range(SEEDS)]
ag = np.mean([[np.mean(runs[i] == runs[j]) for j in range(SEEDS)] for i in range(SEEDS)])
print(f"\n  O/E-only seed stability over {SEEDS} seeds: {ag:.3f}"
      "   (current substrate reports 0.977)")

pd.DataFrame({"gene_id": W.index, "region_current": c_all, "region_oe_only": c_oe,
              "reach_current": r_all, "reach_oe_only": r_oe}
             ).to_csv(OUT / "oe_only_composition.tsv", sep="\t", index=False)
print(f"\nwrote {OUT / 'oe_only_composition.tsv'}")
