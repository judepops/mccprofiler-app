"""Store layout — the contract between the science side and the server.

`build_store.py` (cd4env) writes it. The server reads it and nothing else: it
never imports mccprofiler and never opens the 4.4 GB pickle. See PLAN.md §2.

Rule: **every figure in the app is data-backed.** The store holds numbers, the
frontend renders them. No pre-rendered PNG is ever served — a static image
cannot be brushed, zoomed, filtered or recoloured, and the app's whole point is
that those interactions are live.
"""

from __future__ import annotations

STORE_VERSION = 1

# ---------------------------------------------------------------------------
# channels
# ---------------------------------------------------------------------------

# Pickle key -> short channel name used throughout the API.
# combined_matrix is the MCC contact channel and the ONLY one clustering reads.
# The other seven gate membership or serve as validation; they are displayed but
# must never be described as clustering inputs (PLAN.md trap #5).
CHANNELS: dict[str, str] = {
    "combined_matrix": "mcc",
    "atac_matrix": "atac",
    "ctcf_matrix": "ctcf",
    "rna_plus_matrix": "rna_plus",
    "rna_minus_matrix": "rna_minus",
    "h3k4me1_matrix": "h3k4me1",
    "h3k4me3_matrix": "h3k4me3",
    "h3k27ac_matrix": "h3k27ac",
}

PRIMARY_CHANNEL = "mcc"

CHANNEL_ROLE: dict[str, str] = {
    "mcc": "clustering substrate",
    "atac": "gates gene and peak membership — never a clustering feature",
    "ctcf": "validation",
    "rna_plus": "validation",
    "rna_minus": "validation",
    "h3k4me1": "validation",
    "h3k4me3": "validation",
    "h3k27ac": "validation",
}

# ---------------------------------------------------------------------------
# multiscale pyramid
# ---------------------------------------------------------------------------

# L0 is the native 50 bp grid. Each level averages a fixed factor of the one
# below, so bin i at level n covers bins [i*f, (i+1)*f) at L0.
N_BINS_L0 = 40_000
BIN_SIZE_L0 = 50

LEVELS: dict[int, int] = {0: 1, 1: 5, 2: 20}  # level -> downsample factor

def level_bins(level: int) -> int:
    return N_BINS_L0 // LEVELS[level]

def level_bp(level: int) -> int:
    return BIN_SIZE_L0 * LEVELS[level]

DEFAULT_LEVEL = 2  # 2,000 bins at 1 kb — the default whole-window view

# ---------------------------------------------------------------------------
# HDF5 layout
# ---------------------------------------------------------------------------
#
#   /profiles/{channel}/L{level}   float32 (n_genes, level_bins)
#   /expected/{channel}            float32 (N_BINS_L0,)   panel mean, for O/E
#   /gene_ids                      utf-8   (n_genes,)
#
# Chunked one gene per chunk so a single-gene read touches one chunk. gzip with
# shuffle: MCC profiles are zero-heavy and compress hard.

H5_PROFILES = "profiles"
H5_EXPECTED = "expected"
H5_GENE_IDS = "gene_ids"

COMPRESSION = "gzip"
COMPRESSION_OPTS = 4

# Numerical floor for the expected curve when forming O/E. Matches
# mccprofiler.config.OE_EPS so the app's O/E agrees with the pipeline's.
OE_EPS = 1e-3

# ---------------------------------------------------------------------------
# tables (parquet, under store/tables/)
# ---------------------------------------------------------------------------

TABLES: dict[str, str] = {
    "genes": "one row per gene: id, symbol, source, viewpoint chrom/pos, label",
    "features": "91 MCCProfiler features, z-scored, plus panel percentile",
    "dimensions": "named axes: variance %, loadings, per-anchor rho (91-feat space)",
    "dimension_reproducibility": "per-axis reproducibility (63-shared-feat space)",
    "dimension_scores": "per-gene score on each named axis, plus percentile",
    "posteriors": "per-gene archetype posterior, confidence, entropy",
    "coords": "gcPCA and PC coordinates per gene",
    "peaks": "per-peak element class, signed distance, mapped gene",
    "ps_fits": "per-gene P(s) exponents and fit quality",
    "ps_curves": "log-spaced binned P(s) per gene",
    "anchors": "gene-level biology anchors (gene_table superset)",
    "external_groups": "reference-set stratification vs random floor",
    "reproducibility": "per-feature rho across the 116 twice-captured genes",
    "confounds_locus": "probe-level confound associations per axis",
    "confounds_window": "window-level confound associations per axis, incl. radial",
    "nested_baselines": "91 features vs n_peaks vs magnitude basis",
    "radar": "archetype median and IQR per radar axis",
}

# Reference sets smaller than this are not offered in the cohort view. Of 23
# sets, Roadmap_CTCF_bound has 1 GW gene and ChromHMM_polycomb has 7 —
# enrichments over those are noise (PLAN.md §5 view 5).
MIN_GROUP_N = 25

# Posterior below which a gene renders as "mixture" rather than a named
# archetype. 44% of active genes fall here and that is the point, not a defect.
CORE_POSTERIOR_MIN = 0.8

# ---------------------------------------------------------------------------
# display names — architecture, not borrowed biology (decided 2026-08-03)
# ---------------------------------------------------------------------------
#
# The pipeline's labels name groups after biological categories they do not
# actually track. Measured on the 1,846-gene panel:
#
#   * arch-HK is NOT distinctively housekeeping. Eisenberg-HK fraction 38.2%
#     (OR 1.32 vs rest of panel) against arch-ME-constitutive's 40.2% (OR 1.33) —
#     the name does not discriminate the thing it is named after. arch-HK also
#     has the LOWEST median blood expression of the three active groups
#     (4.5 TPM vs 11.8 and 9.9) and is dead average on promoter-drivenness
#     (z -0.04). Its one real distinction is essentiality, 17.9% vs 13.8%.
#   * arch-ME-constitutive is the promoter-driven group (promoter_signal_fraction
#     z +0.72, the highest of any group), which "multi-enhancer" actively
#     contradicts.
#   * arch-ME-effector is dominated by ONE tall enhancer contact
#     (raw_peak_max_max_enhancer +1.03) with LOW signal entropy (-0.97), so
#     "multi"-enhancer misdescribes it too.
#
# Names below are derived from the top discriminating features (group mean minus
# rest-of-panel mean, in z units) and make no functional claim. Canonical ids
# are preserved for provenance and for joining back to the pipeline; nothing
# upstream is renamed.
ARCHETYPE_DISPLAY: dict[str, dict] = {
    "arch-HK": {
        "display": "dispersed",
        "architecture": "Contacts spread broadly across the window: high signal "
                        "entropy, few empty bands, strong far-distal and distal "
                        "signal, high connectivity.",
        "top_features": ["signal_entropy +1.25", "empty_band_fraction -1.17",
                         "frac_far_distal +1.09", "mean_degree +1.01"],
    },
    "arch-ME-constitutive": {
        "display": "promoter-local",
        "architecture": "Signal anchored on promoter-class peaks and held close: "
                        "highest promoter signal fraction of any group, more "
                        "promoter peaks, local-band enriched, far-distal depleted.",
        "top_features": ["promoter_signal_fraction_raw +1.13", "n_peaks_promoter +0.93",
                         "frac_local +0.80", "frac_far_distal -0.77"],
    },
    "arch-ME-effector": {
        "display": "enhancer-focal",
        "architecture": "One dominant, tall enhancer-class contact rather than many: "
                        "high enhancer peak maximum and enrichment, but LOW signal "
                        "entropy and depleted far-distal signal.",
        "top_features": ["raw_peak_max_max_enhancer +1.03", "signal_entropy -0.97",
                         "enhancer_signal_fraction_raw +0.86", "frac_far_distal -0.91"],
    },
    "arch-sparse": {
        "display": "sparse",
        "architecture": "Few peaks and low connectivity, with what signal exists "
                        "concentrated into a single peak.",
        "top_features": ["n_high_consensus_peaks_075 -1.62", "n_peaks_all -1.57",
                         "mean_degree -1.34", "frac_signal_in_top_peak +1.29"],
    },
    "arch-off": {
        "display": "empty (QC)",
        "architecture": "Near-empty signal — every CTCF-class descriptor sits at its "
                        "floor (z -5 to -8). A quality-control class, not biology.",
        "top_features": ["max_distance_to_viewpoint_ctcf -8.37",
                         "oe_fwhm_bp_max_ctcf -6.98", "consensus_fraction_max_ctcf -6.92"],
        "is_qc": True,
    },
}

# Shown wherever groups are displayed. The groups are regions of a continuum,
# not discovered clusters, and the UI must not let that fall away.
CONTINUUM_CAVEAT = (
    "These are named regions of a continuum, not discovered clusters. The gap "
    "statistic returns k=1 with the gap declining monotonically, HDBSCAN returns "
    "a single cluster, and the dip test is unimodal on every leading component. "
    "The divisions are an imposed resolution choice — defensible because they "
    "reproduce across independent captures (Cohen's kappa 0.72 at k=3-4), not "
    "because the data separates into them. 44% of active genes are mixtures."
)
