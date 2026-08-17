"""Store layout, the contract between the science side and the server.

`build_store.py` (cd4env) writes it. The server reads it and nothing else: it
never imports mccprofiler and never opens the 4.4 GB pickle. See PLAN.md §2.

Rule: **every figure in the app is data-backed.** The store holds numbers, the
frontend renders them. No pre-rendered PNG is ever served, a static image
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
    "atac": "gates gene and peak membership, never a clustering feature",
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

DEFAULT_LEVEL = 2  # 2,000 bins at 1 kb, the default whole-window view

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
# sets, Roadmap_CTCF_bound has 1 GW gene and ChromHMM_polycomb has 7,
# enrichments over those are noise (PLAN.md §5 view 5).
MIN_GROUP_N = 25

# The symmetric rule, added 2026-08-16. A set covering most of the panel cannot
# produce an interpretable "set versus rest" contrast at any effect size,
# because the comparison group is whatever is left over and is defined only by
# exclusion. Three sets exceed this on gw_cd4_1: ChromHMM_active_TSS at 88%
# (rest = 227 genes), DICE_top_TPM_quartile and CpG_island_promoter at 83%.
#
# Flagged rather than dropped. Removing them would invite the question of why a
# set covering most of the panel is missing, and the honest answer is more
# useful than a silent omission: the row exists, it is simply not evidence in
# either direction. It also explains the otherwise surprising near-null for
# CpG-island promoters, which are architecturally distinct in the literature but
# describe 83% of this panel.
#
# The headline results are unaffected: dbSUPER super-enhancers cover 8.6% and
# Lambert_TF 8.4%, so the matched-size comparison sits well inside the bound,
# and the definition-type result is unchanged by dropping all three
# (p = 0.0019 against 0.0018).
MAX_GROUP_COVERAGE = 0.70


def coverage_ok(n_in_set: int, n_panel: int) -> bool:
    """Whether a set/rest contrast is interpretable at all for this set."""
    return MIN_GROUP_N <= n_in_set <= MAX_GROUP_COVERAGE * n_panel

# Posterior below which a gene renders as "mixture" rather than a named
# archetype. 44% of active genes fall here and that is the point, not a defect.
CORE_POSTERIOR_MIN = 0.8

# ---------------------------------------------------------------------------
# display names, architecture, not borrowed biology (decided 2026-08-03)
# ---------------------------------------------------------------------------
#
# The pipeline's labels name groups after biological categories they do not
# actually track. Measured on the 1,846-gene panel:
#
#   * arch-HK is NOT distinctively housekeeping. Eisenberg-HK fraction 38.2%
#     (OR 1.32 vs rest of panel) against arch-ME-constitutive's 40.2% (OR 1.33),
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
        "architecture": "Near-empty signal, every CTCF-class descriptor sits at its "
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
    "G1": "Intensity, how much signal is there",
    "G2": "Spatial allocation, where in the window it sits",
    "G3": "Distance moments and bimodality",
    "G4": "Shape and inequality",
    "G5": "O/E on the profile",
    "P1": "Peak composition",
    "P2": "Per-peak aggregations, by element class",
    "P6": "Topology, peak-peak graph",
}

EXACT_BLOCK: dict[str, str] = {
    # G1, intensity
    **{n: "G1" for n in ("total_mcc", "max_mcc", "mean_mcc", "q90_mcc", "q99_mcc")},
    # G2, spatial allocation
    **{n: "G2" for n in ("frac_viewpoint_proximal", "frac_local", "frac_distal",
                         "frac_far_distal", "bait_pileup_fraction",
                         "distal_signal_density", "local_to_distal_ratio")},
    # G3, distance moments
    **{n: "G3" for n in ("mean_contact_distance", "median_contact_distance",
                         "std_contact_distance", "contact_distance_skew",
                         "contact_distance_kurtosis", "bimodality_score",
                         "oe_distance_kurtosis", "corr_oe_distance",
                         "corr_raw_distance")},
    # G4, shape and inequality
    **{n: "G4" for n in ("gini_mcc", "contact_asymmetry", "signal_entropy",
                         "signal_entropy_distal", "signal_profile_kurtosis",
                         "empty_band_fraction", "frac_signal_in_top_peak",
                         "peak_dominance_index", "spacing_regularity",
                         "mean_peak_gap_bp", "dominant_peak_isolation_bp", "dominant_peak_distance_bp",
                         "distance_to_nearest_peak_bp")},
    # G5, O/E computed on the profile itself, not on peaks
    **{n: "G5" for n in ("oe_distal_max", "oe_distal_mean", "oe_proximal_max",
                         "raw_proximal_max")},
    # P6, topology
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

# =============================================================================
# THE DISPLAY COORDINATE SYSTEM IS THE AMOUNT-CORRECTED ("SHAPE") SPACE.
# =============================================================================
# Switched 2026-08-16. The app previously plotted raw PCs while every Aim 3
# result was computed on the corrected substrate, so the project carried two
# coordinate systems at once. That is the confusion behind the earlier sPC/PC
# mix-ups, and it got worse after the feature removals: on the 73-feature
# substrate raw PC1 correlates 0.505 with the magnitude basis and PC3 correlates
# 0.519, so amount is smeared across two components and neither can be named
# honestly. The default PC2xPC3 plane had been chosen *because* PC1 was the
# amount axis, and that stopped being true.
#
# In the corrected space magnitude is projected out first, so every component
# correlates 0.000 with it by construction, and the axes name cleanly. Note that
# sPC2 and sPC3 are far crisper than the raw PC3 they replace, and that raw PC4's
# wrong name ("enhancer vs CTCF") disappears along with the component.
SPC_LABELS: dict[str, str] = {
    # RE-DERIVED 2026-08-16 against the 85-feature substrate. The previous set
    # was read off the 73-feature loadings and FOUR OF FIVE were wrong after the
    # NaN fix moved the components: sPC2's poles were inverted (it reads
    # promoter-positive, and the label said enhancer-positive), sPC3's negative
    # pole is enhancer not promoter, sPC5 was a different axis entirely, and
    # sPC1's negative pole is the mid-range band, not "short reach".
    #
    # An inverted pole is worse than a vague one: it hands the reader the
    # opposite biology with full confidence. Re-derive these after ANY substrate
    # change and check the sign, not just the words.
    "spc1": "sPC1 (mid-range contacts vs far-reaching and widely spaced)",
    "spc2": "sPC2 (promoter-rich vs enhancer-rich)",
    "spc3": "sPC3 (CTCF-rich vs enhancer-rich)",
    "spc4": "sPC4 (locally enriched at the viewpoint vs dispersed)",
    "spc5": "sPC5 (bait pile-up vs local enrichment further out)",
    **{f"spc{i}": f"sPC{i}" for i in range(6, 31)},
}

SPC_POLES: dict[str, dict[str, str]] = {
    # Keys are "neg"/"pos", matching PC_POLES and the frontend's {neg, pos}.
    # SIGNS VERIFIED against shape_loadings on 2026-08-16; see SPC_LABELS for
    # what went wrong last time.
    "spc1": {"neg": "mid-range contacts (50-250 kb)",
             "pos": "far-reaching, widely spaced peaks"},
    "spc2": {"neg": "enhancer-rich", "pos": "promoter-rich"},
    "spc3": {"neg": "enhancer-rich", "pos": "CTCF-rich"},
    "spc4": {"neg": "locally enriched, piled near the viewpoint",
             "pos": "dispersed, high entropy"},
    "spc5": {"neg": "bait pile-up at the viewpoint",
             "pos": "local enrichment further out"},
}

# The default display plane. sPC1 x sPC2 rather than the old PC2 x PC3: with
# magnitude removed there is no longer any reason to skip the first component,
# and the old plane's justification ("PC1 is amount") no longer holds.
DEFAULT_PLANE = ("spc1", "spc2")

# -----------------------------------------------------------------------------
# RAW PC labels below. RETAINED FOR PROVENANCE AND THE TOGGLE ONLY.
# Do not use these as the primary coordinate system; see SPC_LABELS above.
# The pc1/pc3 comments are pre-removal (91-feature) readings and are stale in
# their specifics: pc1's top loading was mean_degree, which no longer exists.
# -----------------------------------------------------------------------------
# PC labels from audit/Archetype_Tests/14_name_dimensions.py. Verified aligned:
# the store's PCA reproduces dimension_names.tsv's variance shares exactly
# (15.27 / 9.47 / 8.21 / 6.33 / 5.88), so these are the same components.
# STALE as of the 2026-08-16 rebuild: variance is now 17.23 / 10.92 / 9.34 /
# 7.36 / 7.06 and dimension_names.tsv no longer aligns. Do not join on index.
PC_LABELS: dict[str, str] = {
    # Format is "PC<n> (interpretation)" so the component number is never hidden
    # behind the prose. The interpretation is exactly that, a reading of the
    # loadings, and the app shows those loadings so it can be checked.
    # Measured 2026-08-03, and NOT "amount" in the sense of signal. PC1
    # correlates with total_mcc at +0.033, essentially zero. It tracks peak
    # counts and reach: n_high_consensus_peaks_075 +0.787, mean_degree +0.763,
    # n_peaks_all +0.717, max_distance_to_viewpoint_all +0.715. Total signal
    # sits on PC3 (-0.621) and PC6 (-0.433). Calling this axis "amount" implied
    # amount was quarantined here, which is the opposite of the truth.
    "pc1": "PC1 (peak richness and reach vs emptiness)",
    "pc2": "PC2 (local vs long-range)",
    # CAVEAT, measured 2026-08-03 from the stored loadings: the positive pole is
    # cleanly promoter (promoter_signal_fraction +0.278/+0.254), but the negative
    # pole is NOT purely enhancer, its strongest loading is total_mcc (-0.227),
    # i.e. raw amount, with CTCF terms (-0.193, -0.176) as prominent as
    # n_peaks_enhancer (-0.216). The inherited name over-simplifies. Qualify
    # before quoting PC3 as "housekeeping vs developmental from contact alone".
    "pc3": "PC3 (promoter-driven vs enhancer/CTCF + amount)",
    "pc4": "PC4 (enhancer vs CTCF composition)",
    "pc5": "PC5 (concentrated vs dispersed)",
    **{f"pc{i}": f"PC{i}" for i in range(6, 31)},
    "gcpca_b1_1": "gcPC1 (density-free, B1)",
    "gcpca_b1_2": "gcPC2 (density-free, B1)",
    "gcpca_b2_1": "gcPC1 (density-free, B2)",
    "gcpca_b2_2": "gcPC2 (density-free, B2)",
    "umap_1": "UMAP 1",
    "umap_2": "UMAP 2",
    "umap_null_1": "UMAP 1 (permuted null)",
    "umap_null_2": "UMAP 2 (permuted null)",
}

# Which END of each axis is which.
#
# A label like "local vs long-range" names the contrast but not the direction,
# and PCA sign is arbitrary, so without this a reader cannot tell which side of
# the plot is local. Poles below are read off the stored loadings, not assumed:
#
#   PC1  + n_high_consensus_peaks +0.21, mean_degree +0.20, n_peaks_all +0.19
#        - empty_band_fraction -0.16, frac_viewpoint_proximal -0.14
#   PC2  + n_peaks_within_100kb +0.25, frac_local +0.18
#        - frac_far_distal -0.23, oe_distal_mean -0.21
#   PC3  + promoter_signal_fraction +0.28
#        - total_mcc -0.23, n_peaks_enhancer -0.22, raw_peak_max_max_ctcf -0.19
#   PC4  + enhancer_signal_fraction_raw +0.28
#        - frac_ctcf -0.21, raw_peak_max_max_ctcf -0.20
#   PC5  + raw_peak_max_mean_all +0.24, mean_peak_gap_bp +0.24,
#          frac_signal_in_top_peak +0.23
#        - signal_entropy -0.21
PC_POLES: dict[str, dict[str, str]] = {
    "pc1": {"neg": "few peaks, short reach", "pos": "many peaks, long reach"},
    "pc2": {"neg": "long-range", "pos": "local"},
    "pc3": {"neg": "enhancer / CTCF, high amount", "pos": "promoter-driven"},
    "pc4": {"neg": "CTCF", "pos": "enhancer"},
    "pc5": {"neg": "dispersed", "pos": "concentrated"},
}

# PC1 is the amount axis, so PC1 x PC2 mostly re-sorts genes by signal depth.
# PC2 x PC3 is the interpretable plane and is the default.
# Changed 2026-08-16 from ("pc2", "pc3") to the amount-corrected space. The old
# default existed because PC1 was the amount axis and was worth skipping; after
# the feature removals amount sits on PC1 (r 0.505) AND PC3 (r 0.519), so that
# plane is now half-amount on one axis. See SPC_LABELS.
DEFAULT_EMBEDDING_AXES = DEFAULT_PLANE

# Every axis label, both spaces, so a lookup never silently falls back to the
# bare key. Shape axes first: they are the primary system.
ALL_AXIS_LABELS = {**PC_LABELS, **SPC_LABELS}
ALL_AXIS_POLES = {**PC_POLES, **SPC_POLES}

SHAPE_SPACE_NOTE = (
    "Axes are the amount-corrected (shape) space: the magnitude basis is "
    "projected out before PCA, so every component is uncorrelated with overall "
    "signal by construction and position here cannot mean 'this gene has more "
    "signal'. This is the same space every Aim 3 result is computed in. Raw PCs "
    "remain selectable for provenance, but on the current substrate raw PC1 and "
    "PC3 both carry amount (r 0.505 and 0.519) and are not cleanly nameable."
)

UMAP_CAVEAT = (
    "UMAP optimises a local-neighbour objective and is well documented to render "
    "continuous data as apparent clusters. Compare against the permuted null, "
    "which preserves every feature's marginal distribution exactly while "
    "destroying all joint structure, there is provably nothing to find in it. "
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
    "The divisions are an imposed resolution choice, defensible because they "
    "reproduce across independent captures (Cohen's kappa 0.72 at k=3-4), not "
    "because the data separates into them. 44% of active genes are mixtures. "
    "Reproducibility of an imposed partition is not evidence that the partition "
    "is real (Hennig 2015)."
)

# Added 2026-08-16. The labels are STALE as well as imposed, and the app should
# say so where it renders them rather than only in prose. archetype_labels.tsv
# is dated 2026-07-21: KMeans was fit on the 91-feature substrate, 26 days
# before the degenerate and unreproducible features were removed. `arch-HK` was
# named "dispersed" partly on mean_degree +1.01, a feature that no longer exists
# because it was an algebraic function of two peak counts.
#
# The decision (2026-08-16) was to DEMOTE rather than re-derive: re-running
# KMeans at an arbitrary k, on a substrate whose own tests say no k exists,
# would polish the weakest part of the analysis. So the labels stay, marked
# provisional, with the posterior and the mixture fraction beside them.
ARCHETYPE_PROVISIONAL = (
    "Provisional. These labels were fit on the earlier 91-feature substrate "
    "(2026-07-21) and have not been re-derived since 18 features were removed as "
    "degenerate or unreproducible; one of the features that named this grouping "
    "no longer exists. Treat the position and the mixture as the result and the "
    "label as a convenience."
)


# ---------------------------------------------------------------------------
# What this panel is, and is not. Render at the top of every page.
# ---------------------------------------------------------------------------
#
# Added 2026-08-17. `gw_cd4_1` is routinely called the "genome-wide panel" and
# that is wrong in a way that changes how every result below should be read. It
# is genome-DISTRIBUTED (all chromosomes) but not genome-REPRESENTATIVE: it was
# built by TSS extraction followed by ATAC-accessibility selection in CD4, and
# accessible human promoters are overwhelmingly CpG-island and TATA-less. So the
# same gating that decides gene membership also decided which PROMOTER CLASS is
# present.
PANEL_BIAS = {
    "headline": "This is not a genome-representative panel. More genes are needed.",
    # The comparison that matters is IN-PANEL against NOT-IN-PANEL. Quoting the
    # genome-wide average alongside both made it read as a within-panel
    # contrast, which is not what it is.
    "one_line": "85.1% of captured genes have CpG-island promoters, against "
                "61.6% of the genes that were not captured. The panel does not "
                "span the promoter-class range, so more genes are needed.",
    "not_about": "This is about the promoter CLASS the panel samples, not about "
                 "individual genes being unusual. Every captured gene is fine; "
                 "the SET is not representative.",
    "measured": [
        # Verified against CURRENT_FINDINGS section 6c on 2026-08-17. The odds
        # ratio compares panel against NOT-IN-PANEL (61.6%), which is the right
        # contrast; against the genome-wide figure, which includes the panel, it
        # is 3.24.
        {"what": "CpG-island promoters", "panel": "85.1%", "genome": "63.8%",
         "not_in_panel": "61.6%",
         "note": "odds ratio 3.56 against non-panel genes (3.24 against "
                 "genome-wide, which includes the panel)"},
        {"what": "TATA-containing promoters", "panel": "1.3%", "genome": "2.8%",
         "not_in_panel": "2.9%",
         "note": "depleted 2.2-fold; scanned at -34 to -18 on hg38"},
    ],
    "cause": "TSS extraction then ATAC-accessibility selection in CD4. Accessible "
             "human promoters are overwhelmingly CpG-island and TATA-less, so the "
             "gating that chose the genes also chose the promoter class.",
    # The right framing is restriction of range, not circularity. Nothing feeds
    # back into anything; the panel simply does not span the axis. Calling it an
    # echo chamber invites the wrong rebuttal, because a reviewer can answer that
    # charge and still leave the real problem standing.
    "framing": "Restriction of range, not circularity. Nothing feeds back into "
               "anything: the panel does not span the promoter-class axis, so "
               "effects that depend on that contrast are attenuated before any "
               "test runs.",
    "scope": "Conclusions here apply to accessible, CpG-island-rich promoters in "
             "CD4+ T cells. They do not license statements about the human "
             "promoter repertoire as a whole.",
    # Which conclusions this threatens and which it does not. gene_desert
    # reaching AUC 0.806 with no feature selection proves the panel is not so
    # compressed that nothing separates.
    "affects": [
        {"claim": "the landscape is continuous (no clusters)", "status": "robust",
         "why": "restricted range would create spurious homogeneity, not destroy "
                "real clusters; and 18 raw / 19 shape components sit above the "
                "parallel-analysis noise ceiling, so the panel is far from "
                "architecturally uniform"},
        {"claim": "super-enhancers are displaced but not separated", "status": "robust",
         "why": "SE is about the enhancer landscape, not promoter class"},
        {"claim": "housekeeping genes have no architectural signature",
         "status": "QUALIFIED",
         "why": "in an 85%-CpG panel this compares CpG-island housekeeping genes "
                "against CpG-island everything-else; the contrast the hypothesis "
                "is about was selected out before the test ran"},
        # Narrowed 2026-08-17: it is the HOUSEKEEPING DISCRIMINANT that cannot be
        # asked, not core-promoter sequence classes in general. The Initiator set
        # has 80 panel genes and clears the floor.
        {"claim": "the TATA-vs-CpG housekeeping discriminant",
         "status": "CANNOT BE ASKED",
         "why": "23 panel genes are TATA-containing and only 4 are focused-TATA "
                "without a CpG island, both below the 25-gene floor; the two CpG "
                "sets cover 85% and 84% and exceed the 70% coverage ceiling"},
        {"claim": "core-promoter elements generally", "status": "PARTLY testable",
         "why": "coreprom_Inr has 80 panel genes (4%) and is usable. Inr is a "
                "core promoter element rather than the housekeeping "
                "discriminant, so it tests sequence-defined architecture without "
                "settling the housekeeping question"},
    ],
    # Verified 2026-08-17 against experiment_feature_search on the 85-feature
    # substrate: 0.812 with all features and no selection, 0.819 under greedy
    # search. The earlier 0.806 was the 73-feature run.
    "calibration": "The panel is not so compressed that nothing separates: "
                   "gene_desert reaches AUC 0.812 with no feature selection at "
                   "all. But that is a POSITIONAL variable. There is no "
                   "comparably clean calibration on a regulatory category that "
                   "is independent of chromatin assays, which is the sharpest "
                   "objection to the nulls and the second reason for the 20k "
                   "panel.",
    # Why this motivates data collection, which is the point of stating it.
    "consequence": "THIS IS THE LOAD-BEARING ARGUMENT FOR THE 20k PANEL, and it "
                   "is now a measurement rather than a convenience.",
    "why_more_data": [
        "Two of the four claim classes above are qualified or unanswerable "
        "purely because the panel does not span promoter classes. No amount of "
        "reanalysis fixes that; only more genes do.",
        "The only clean calibration is gene_desert at AUC 0.812, which is a "
        "POSITIONAL variable. A panel spanning promoter classes would allow "
        "calibration against a known REGULATORY distinction, which is the "
        "sharpest outstanding objection to every null result here.",
        "1,846 genes is also thin for representation learning (Aim 2), and "
        "three reference sets each cover over 80% of this panel, which would "
        "not happen genome-wide.",
    ],
}
