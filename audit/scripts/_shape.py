"""Shared, correct amount-correction for every shape-only / amount-adjusted analysis.

WHY THIS EXISTS
---------------
The old convention `AMOUNT = [total_mcc, max_mcc, n_peaks_all]` (drop 3, or adjust for 3)
UNDER-corrects. Overall signal magnitude leaks into most "shape" features: q90_mcc, peak
heights, degree, signal_entropy, and even `promoter_signal_fraction` (~70% predictable from
magnitude). Adjusting for only 3 features let magnitude-confounded biology look real.

THE CORRECT BASIS (verified empirically, 2026-07-27)
----------------------------------------------------
`MAG_OVERALL` = OVERALL magnitude only: signal levels, TOTAL peak counts, peak heights,
connectivity, density. It deliberately EXCLUDES:
  - per-class peak counts (n_peaks_enhancer/ctcf/promoter) and per-class heights, and
  - all fractions / ratios / asymmetries / distances / entropies,
because the enhancer-vs-promoter balance IS the architecture we want to keep. Putting
per-class counts into the amount basis is circular over-correction: it wrongly annihilates
the real ME effector/constitutive split (its biology drops to NONE under the 17-feature
basis but SURVIVES under this 11-feature overall-magnitude basis: GWAS_immune OR 2.46 p=0.02,
gene_desert 4.56 p=4e-8, essential/HK depletion p<1e-3).

USAGE (replaces `AMOUNT = [...]; Xshape = drop(AMOUNT)`)
-------------------------------------------------------
    from _shape import MAG_OVERALL, magnitude_confounders, corrected_shape
    Xamt, amt_names = magnitude_confounders(X, feats)   # confounder matrix for biology tests
    Xshape, shape_feats = corrected_shape(X, feats)     # residualised architecture substrate

Every feature is regressed on overall magnitude and replaced by its residual, so per-class
counts survive as per-class *bias* (architecture), not amount. The magnitude basis features
themselves are dropped. Residualisation is fit on the full panel by default.
"""
import numpy as np
from sklearn.preprocessing import StandardScaler

MAG_OVERALL = [
    "total_mcc", "max_mcc", "q90_mcc",              # MCC signal levels
    "n_peaks_all",                                    # total peak count
    "n_peaks_within_100kb", "n_high_consensus_peaks_075",  # overall peak density / count
    "raw_peak_max_max_all", "raw_peak_max_mean_all",  # overall peak heights
    "mean_degree", "mean_degree_raw",                 # overall connectivity, SEE BELOW
    "distal_signal_density",                          # overall density
]

# `mean_degree` and `mean_degree_raw` no longer exist in the feature substrate.
# They were dropped as DEGENERATE on 2026-08-16 (config.py::DEGENERATE_FEATURES)
# because topology.py builds the complete graph on active peaks, making both
# algebraic functions of two peak counts rather than measurements.
#
# DO NOT DELETE THEM FROM THIS LIST. `_basis_idx` filters to what a matrix
# actually contains, so they cost nothing on the current 73-feature substrate
# (which resolves to a 9-feature basis). But `features_raw.pkl` and the held-back
# 791-gene immune panel still carry all 91 features. On those, removing these two
# entries would stop the correction absorbing them as amount and start treating
# them as SHAPE, putting the degenerate artefact back into the corrected
# substrate through the back door.
#
# Consequence: `len(MAG_OVERALL)` is 11 and OVERSTATES the resolved basis on the
# current substrate. Report `len(magnitude_confounders(X, feats)[1])` instead;
# several call sites still print the literal 11 and are wrong to.

# kept for reference / provenance: the old (insufficient) definition
LEGACY_AMOUNT = ["total_mcc", "max_mcc", "n_peaks_all"]


def _basis_idx(feats):
    return [feats.index(f) for f in MAG_OVERALL if f in feats]


def magnitude_confounders(X, feats):
    """Overall-magnitude confounder matrix (n, k) + names, for amount-adjusted biology tests."""
    idx = _basis_idx(feats)
    return np.asarray(X)[:, idx], [feats[i] for i in idx]


def corrected_shape(X, feats, fit_mask=None, standardize=True):
    """Amount-corrected shape substrate.

    Residualise every non-magnitude feature against overall magnitude (linear), drop the
    magnitude basis, z-score. Per-class architecture is preserved as residual bias.

    fit_mask: optional boolean mask to fit the residualisation coefficients on a subset
              (e.g. active genes); the transform is applied to all rows. Default: fit on all.
    Returns (Xshape_corrected, kept_feature_names).
    """
    X = np.asarray(X, float)
    idx = _basis_idx(feats)
    rest = [i for i in range(len(feats)) if i not in idx]
    m = np.ones(len(X), bool) if fit_mask is None else np.asarray(fit_mask, bool)
    Bfit = np.column_stack([np.ones(int(m.sum())), X[np.ix_(m, idx)]])
    beta, _, _, _ = np.linalg.lstsq(Bfit, X[np.ix_(m, rest)], rcond=None)
    Ball = np.column_stack([np.ones(len(X)), X[:, idx]])
    R = X[:, rest] - Ball @ beta
    if standardize:
        R = StandardScaler().fit_transform(R)
    return R, [feats[i] for i in rest]
