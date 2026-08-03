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
