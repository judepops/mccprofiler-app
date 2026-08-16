"""How much of each PC name is actually true?

A name like "local vs long-range" is a claim that one interpretable contrast
explains the axis. Three ways to test it, each answering a different objection:

1. CONCENTRATION. What fraction of the axis's squared loading mass sits in the
   features the name refers to? PCA loadings are a unit vector, so squared
   loadings partition to 1 and the fraction is directly readable as "how much of
   this axis is about the named thing".

2. BEST SINGLE PROXY. Correlate the axis's gene scores against every feature.
   If the top correlate is the named feature at r 0.9, the name is a good
   summary. If the top correlate is r 0.5 and belongs to a different concept,
   the name is a label of convenience.

3. AMOUNT CONTAMINATION. Correlate the axis with total_mcc. Any axis that is
   substantially amount is not measuring the shape contrast it claims, which is
   exactly what was already found for PC3.
"""
import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "backend"))

import numpy as np
import pandas as pd
from app.store import get_store
from app import store_schema as S

s = get_store()
emb = s.table("embeddings").set_index("gene_id")
feat = s.table("features")
load = s.table("pc_loadings")

X = feat.pivot(index="gene_id", columns="feature", values="z")
X = X.loc[emb.index]
print(f"genes {X.shape[0]}, features {X.shape[1]}\n")

# The features each name actually refers to. Chosen from the concept, not from
# the loadings, otherwise the test is circular: it would be scored against the
# very features that produced the name.
CLAIMS = {
    1: ("amount, richness and reach vs emptiness",
        ["total_mcc", "n_peaks_all", "n_peaks_within_100kb", "n_peaks_within_250kb",
         "n_high_consensus_peaks_075", "distal_signal_density", "raw_signal_sum",
         "n_peaks_enhancer", "n_peaks_promoter", "n_peaks_ctcf", "mean_degree",
         "mean_degree_raw", "empty_band_fraction"]),
    2: ("local vs long-range",
        ["frac_local", "frac_distal", "frac_far_distal", "local_to_distal_ratio",
         "frac_promoter_proximal", "distal_signal_density", "signal_entropy_distal",
         "mean_distance_to_viewpoint", "max_distance_to_viewpoint_all",
         "max_distance_to_viewpoint_enhancer", "max_distance_to_viewpoint_promoter",
         "max_distance_to_viewpoint_ctcf", "dominant_peak_distance_bp",
         "n_peaks_within_100kb", "n_peaks_within_250kb", "oe_distal_mean"]),
    3: ("promoter-driven vs enhancer/CTCF",
        ["promoter_signal_fraction", "promoter_signal_fraction_raw",
         "enhancer_signal_fraction", "enhancer_signal_fraction_raw",
         "ctcf_signal_fraction", "ctcf_signal_fraction_raw",
         "n_peaks_promoter", "n_peaks_enhancer", "n_peaks_ctcf",
         "frac_promoter", "frac_enhancer", "frac_ctcf"]),
    4: ("enhancer vs CTCF composition",
        ["enhancer_signal_fraction", "enhancer_signal_fraction_raw",
         "ctcf_signal_fraction", "ctcf_signal_fraction_raw",
         "n_peaks_enhancer", "n_peaks_ctcf", "frac_enhancer", "frac_ctcf",
         "raw_peak_max_max_enhancer", "raw_peak_max_max_ctcf"]),
    5: ("concentrated vs dispersed",
        ["signal_entropy", "signal_entropy_distal", "frac_signal_in_top_peak",
         "dominant_peak_isolation_bp", "empty_band_fraction", "mean_peak_gap_bp",
         "consensus_fraction_max_all", "gini_coefficient"]),
}

avail = set(X.columns)
rows = []
for pc, (name, claimed) in CLAIMS.items():
    col = f"pc{pc}"
    sc = emb[col]
    L = load[load["pc"] == pc].set_index("feature")["loading"]
    L = L.reindex(X.columns).fillna(0.0)
    L2 = (L ** 2)
    L2 = L2 / L2.sum()

    present = [c for c in claimed if c in avail]
    missing = [c for c in claimed if c not in avail]
    conc = float(L2[present].sum())
    # A useful yardstick: if loading mass were spread evenly across all 91
    # features, this block would hold this much by size alone.
    chance = len(present) / len(X.columns)

    r = X.corrwith(sc)
    r_abs = r.abs().sort_values(ascending=False)
    top = r_abs.index[0]
    r_top = r[top]
    best_claimed = r[present].abs().sort_values(ascending=False)
    bc, bcr = best_claimed.index[0], r[best_claimed.index[0]]
    r_amount = r.get("total_mcc", np.nan)

    rows.append(dict(pc=pc, name=name, conc=conc, chance=chance,
                     lift=conc / chance, top=top, r_top=r_top,
                     best_claimed=bc, r_best_claimed=bcr, r_amount=r_amount,
                     n_missing=len(missing)))

    print("=" * 78)
    print(f"PC{pc}  {name}")
    print(f"  loading mass in the named concept : {conc:5.1%}   "
          f"(even spread would give {chance:5.1%}, lift {conc/chance:.2f}x)")
    print(f"  strongest correlate overall       : {top:<34s} r = {r_top:+.3f}")
    print(f"  strongest NAMED correlate         : {bc:<34s} r = {bcr:+.3f}")
    print(f"  correlation with total_mcc        : {r_amount:+.3f}"
          f"    <-- amount contamination")
    if missing:
        print(f"  (not in the 91: {', '.join(missing)})")
    print("  top 6 correlates, whatever they are:")
    for f in r_abs.index[:6]:
        mark = "  [named]" if f in present else ""
        print(f"      {f:<38s} r = {r[f]:+.3f}{mark}")

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
sm = pd.DataFrame(rows)
print(f"{'PC':<4}{'named-concept mass':>20}{'lift':>7}{'best named r':>15}"
      f"{'top r':>9}{'r(amount)':>11}")
for _, x in sm.iterrows():
    print(f"{x.pc:<4}{x.conc:>19.1%}{x.lift:>7.2f}{x.r_best_claimed:>15.3f}"
          f"{x.r_top:>9.3f}{x.r_amount:>11.3f}")

# How much of the panel does the top-5 story even cover?
sc_tbl = s.table("pc_scree")
top5 = sc_tbl[sc_tbl["pc"] <= 5]["variance_pct"].sum()
above = int(sc_tbl["above_noise"].sum())
tot_above = sc_tbl[sc_tbl["above_noise"]]["variance_pct"].sum()
print(f"\nPC1-5 carry {top5:.1f}% of variance. {above} components sit above the "
      f"noise ceiling, together {tot_above:.1f}%.")
print(f"So naming 5 axes describes {top5:.1f}% of the variation and leaves "
      f"{tot_above - top5:.1f}% real-but-unnamed.")
