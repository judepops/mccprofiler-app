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

# ---------------------------------------------------------------------------
# element classes and feature blocks
# ---------------------------------------------------------------------------

# Three classes from the `RE` column of annotated.tsv. Counts on the GW panel:
# ctcf 11,487 / enhancer 15,035 / promoter 14,465.
ELEMENT_COLOR = {
    "enhancer": "#c2703d",
    "ctcf": "#4a7c59",
    "promoter": "#2b5070",
}

# Feature blocks, per mccprofiler/FEATURE_SPEC.md, which gives the block sizes
# as G1 5, G2 7, G3 6, G4 5, G5 3, P1 12, P2-P5 64, P6 11.
#
# The classifier is exact-name-first, because prefix matching gets this badly
# wrong: an `oe_` prefix rule swept the 60-odd per-peak O/E aggregations into
# G5, which the spec says holds exactly three profile-level features. The
# distinguishing property of the peak arm is the element-class suffix
# (_all / _enhancer / _ctcf / _promoter), since the peak block computes every
# aggregation once per class.
ELEMENT_SUFFIXES = ("_all", "_enhancer", "_ctcf", "_promoter")

BLOCK_DESC: dict[str, str] = {
    "G1": "Intensity — how much signal is there",
    "G2": "Spatial allocation — where in the window it sits",
    "G3": "Distance moments and bimodality",
    "G4": "Shape and inequality",
    "G5": "O/E on the profile",
    "P1": "Peak composition",
    "P2": "Per-peak aggregations, by element class",
    "P6": "Topology — peak-peak graph",
}

EXACT_BLOCK: dict[str, str] = {
    # G1 — intensity
    **{n: "G1" for n in ("total_mcc", "max_mcc", "mean_mcc", "q90_mcc", "q99_mcc")},
    # G2 — spatial allocation
    **{n: "G2" for n in ("frac_promoter_proximal", "frac_local", "frac_distal",
                         "frac_far_distal", "bait_pileup_fraction",
                         "distal_signal_density", "local_to_distal_ratio")},
    # G3 — distance moments
    **{n: "G3" for n in ("mean_contact_distance", "median_contact_distance",
                         "std_contact_distance", "contact_distance_skew",
                         "contact_distance_kurtosis", "bimodality_score",
                         "oe_distance_kurtosis", "corr_oe_distance",
                         "corr_raw_distance")},
    # G4 — shape and inequality
    **{n: "G4" for n in ("gini_mcc", "contact_asymmetry", "signal_entropy",
                         "signal_entropy_distal", "signal_profile_kurtosis",
                         "empty_band_fraction", "frac_signal_in_top_peak",
                         "peak_dominance_index", "spacing_regularity",
                         "mean_peak_gap_bp", "dominant_peak_isolation_bp", "dominant_peak_distance_bp",
                         "distance_to_nearest_peak_bp")},
    # G5 — O/E computed on the profile itself, not on peaks
    **{n: "G5" for n in ("oe_distal_max", "oe_distal_mean", "oe_proximal_max",
                         "raw_proximal_max")},
    # P6 — topology
    **{n: "P6" for n in ("mean_degree", "max_degree", "n_active_pairs",
                         "frac_active_pairs", "n_isolates", "n_active_peaks",
                         "mean_degree_raw", "max_degree_raw", "n_isolates_raw",
                         "frac_active_pairs_raw", "n_active_peaks_raw",
                         "oe_max_cv", "raw_peak_max_cv")},
}

# P1 composition: counts and per-class signal shares.
P1_PREFIXES = ("n_peaks", "n_high_consensus")
P1_SUFFIXES = ("_signal_fraction", "_signal_fraction_raw")


def feature_block(name: str) -> tuple[str, str]:
    """Return (block_code, description) for a feature name."""
    if name in EXACT_BLOCK:
        code = EXACT_BLOCK[name]
        return code, BLOCK_DESC[code]

    if name.startswith(P1_PREFIXES) or name.endswith(P1_SUFFIXES) or name.startswith("frac_"):
        return "P1", BLOCK_DESC["P1"]

    # Anything carrying an element-class suffix is a per-peak aggregation.
    if name.endswith(ELEMENT_SUFFIXES):
        return "P2", BLOCK_DESC["P2"]

    return "other", "Unclassified"


# Kept for the API's block ordering.
FEATURE_BLOCKS = [(c, BLOCK_DESC[c], ()) for c in ("G1", "G2", "G3", "G4", "G5", "P1", "P2", "P6")]


# ---------------------------------------------------------------------------
# embeddings
# ---------------------------------------------------------------------------

# PC labels from audit/Archetype_Tests/14_name_dimensions.py. Verified aligned:
# the store's PCA reproduces dimension_names.tsv's variance shares exactly
# (15.27 / 9.47 / 8.21 / 6.33 / 5.88), so these are the same components.
PC_LABELS: dict[str, str] = {
    "pc1": "amount — richness and reach vs emptiness",
    "pc2": "local vs long-range",
    "pc3": "promoter-driven vs enhancer-driven",
    "pc4": "enhancer vs CTCF composition",
    "pc5": "concentrated vs dispersed",
    "pc6": "PC6",
    "pc7": "PC7",
    "pc8": "PC8",
    "pc9": "PC9",
    "pc10": "PC10",
    "gcpca_b1_1": "gcPC1 (density-free, B1)",
    "gcpca_b1_2": "gcPC2 (density-free, B1)",
    "gcpca_b2_1": "gcPC1 (density-free, B2)",
    "gcpca_b2_2": "gcPC2 (density-free, B2)",
    "umap_1": "UMAP 1",
    "umap_2": "UMAP 2",
    "umap_null_1": "UMAP 1 (permuted null)",
    "umap_null_2": "UMAP 2 (permuted null)",
}

# PC1 is the amount axis, so PC1 x PC2 mostly re-sorts genes by signal depth.
# PC2 x PC3 is the interpretable plane and is the default.
DEFAULT_EMBEDDING_AXES = ("pc2", "pc3")

UMAP_CAVEAT = (
    "UMAP optimises a local-neighbour objective and is well documented to render "
    "continuous data as apparent clusters. Compare against the permuted null, "
    "which preserves every feature's marginal distribution exactly while "
    "destroying all joint structure — there is provably nothing to find in it. "
    "If the null looks similarly grouped, the grouping is the method, not the "
    "biology. PCA is the default here because it is linear and preserves "
    "distances, so a continuum renders as one."
)


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
