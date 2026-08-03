"""Absolute path registry for mccprofiler-app.

Project policy: paths are hardcoded to ``/home/imm/grte4643/...`` deliberately.
Datasets span multiple mounts and relative paths break. Do not rewrite these to
relative paths.

This module is the ONLY place a filesystem path is defined. Nothing else in the
app hardcodes one. Pure stdlib so any python3 can run the audit, following the
precedent set by ``data_transfer/build_briefing.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# roots
# --------------------------------------------------------------------------

DPHIL = Path("/home/imm/grte4643/Documents/DPhil")
MCC = DPHIL / "Data_Exploration" / "MCC"

CD4 = MCC / "cd4_cleaned"
SCRIPTS = CD4 / "scripts_cleaned"
MCC_DATA = CD4 / "mcc"
ORCHID = MCC / "orchid"
LAB = DPHIL / "Lab" / "Protocol_20k"

AUDIT = SCRIPTS / "audit"
CONTINUOUS = AUDIT / "continuous_methods"
GW_AUDIT = AUDIT / "GW"
PROFILER = SCRIPTS / "mccprofiler"

# app-owned
APP_ROOT = Path(__file__).resolve().parents[2]
STORE = APP_ROOT / "store"


# --------------------------------------------------------------------------
# artefact registry
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Artefact:
    """One input the store builder reads.

    ``required`` artefacts abort the build if missing. Optional ones degrade a
    single view and are reported but tolerated.
    """

    key: str
    path: Path
    what: str
    required: bool = True
    phase: str = "P0"
    notes: str = ""


# The genome-wide panel is the default and the only one exposed initially.
# See PLAN.md trap #2 for why the 791 immune panel is held back.
PANEL_GW = "gw_cd4_1"

ARTEFACTS: list[Artefact] = [
    # ---- signal ---------------------------------------------------------
    Artefact(
        "clean_matrix_gw",
        SCRIPTS / "process/output_gw_cd4_1/scripts/outliers/results_combined_50bp_clean.pkl",
        "post-QC, post-outlier 8-channel signal matrix, 1846 x 40000",
        phase="P0",
        notes="4.7 GB. Read once by build_store, never by the server.",
    ),
    Artefact(
        "features_gw",
        PROFILER / "outputs_gw_cd4_1/features/features_zscored.pkl",
        "MCCProfiler feature matrix, 1846 x 91, z-scored",
        phase="P0",
    ),
    # ---- labels and coordinates -----------------------------------------
    Artefact(
        "archetype_labels",
        PROFILER / "outputs_gw_cd4_1/cluster/archetype_labels.tsv",
        "base k=4 archetypes for all 1,846 genes (HK 844 / ME 720 / sparse 261 / off 21)",
        phase="P0",
        notes=(
            "The k=4 file, NOT archetype_labels_k5.tsv which sits beside it and is "
            "blacklisted. Confirmed as the ARI reference by the 2026-07-29 handoff."
        ),
    ),
    Artefact(
        "me_subtypes_labels",
        GW_AUDIT / "output/me_subtypes_labels.tsv",
        "ME split for the 720 arch-ME genes only (constitutive 369 / effector 351)",
        phase="P0",
        notes=(
            "A SUBSET, not a full labelling — it covers only arch-ME. Overlay it on "
            "archetype_labels to get the 5 groups. Counts are the post-amount-"
            "correction 369/351; PROJECT_STATUS.md still quotes the stale 377/343."
        ),
    ),
    Artefact(
        "posteriors",
        GW_AUDIT / "output/posterior_membership.tsv",
        "per-gene archetype posteriors, confidence, entropy",
        phase="P5",
    ),
    Artefact(
        "gcpca_axes",
        CONTINUOUS / "gcpca_axes.npz",
        "gcPCA density-free coordinates, B1 and B2 backgrounds",
        phase="P5",
    ),
    # ---- named dimensions (the 08-03 reframe) ---------------------------
    Artefact(
        "dimension_names",
        CONTINUOUS / "dimension_names.tsv",
        "named axes: loadings, dominant class, per-anchor rho, variance %",
        phase="P3",
        notes="91-feature space. Do NOT join to dimension_reproducibility on PC index.",
    ),
    Artefact(
        "dimension_reproducibility",
        CONTINUOUS / "dimension_reproducibility.tsv",
        "per-axis cross-capture reproducibility, explained/reproduces flags",
        phase="P3",
        notes="63-shared-feature space. Namespaced separately. See PLAN.md trap #10.",
    ),
    Artefact(
        "structure_vs_noise",
        CONTINUOUS / "structure_vs_noise.tsv",
        "parallel analysis: real vs permuted-null dimension counts",
        phase="P3",
        notes="PC1 excess kurtosis field is UNTRIMMED (16.41). Trimmed is 1.01.",
    ),
    # ---- cohort view ----------------------------------------------------
    Artefact(
        "external_groups",
        CONTINUOUS / "external_group_stratification.tsv",
        "externally-defined gene groups vs size-matched random floor",
        phase="P4",
        notes="Display pct_retained and signal_over_random, never eta2.",
    ),
    Artefact(
        "anchor_table",
        SCRIPTS / "public/outputs/anchor_table.tsv",
        "gene-level biology anchors, 18,802 x 25",
        phase="P4",
        notes="Canonical. Has anchor_table.meta.json alongside recording per-column coverage.",
    ),
    Artefact(
        "gene_table",
        MCC / "collaboration/data/gene_table.tsv",
        "anchor_table superset, 18,802 x 42 — built for the Taipale work",
        phase="P4",
        notes=(
            "Preferred over anchor_table: adds gene_length (the uncontrolled confound "
            "flagged 07-31), gc_500bp/gc_2kb, in_gw_panel/in_immune_panel, "
            "viewpoint_chrom/pos, panel_n_peaks_all, and abs_vp_tss_offset which "
            "directly quantifies the viewpoint-vs-TSS trap."
        ),
    ),
    # ---- reproducibility demo -------------------------------------------
    Artefact(
        "cross_panel_reproducibility",
        CONTINUOUS / "cross_panel_reproducibility.tsv",
        "per-feature rho across the 116 twice-captured genes",
        phase="P8",
    ),
    # ---- confound panel -------------------------------------------------
    Artefact(
        "locus_confounds",
        CONTINUOUS / "locus_intrinsic_confounds.tsv",
        "probe-level confounds: GC, alignments, repeat length, density score",
        phase="P9",
    ),
    Artefact(
        "window_confounds",
        CONTINUOUS / "window_confounds.tsv",
        "window-level confounds incl. radial gradients (the key falsification)",
        phase="P9",
    ),
    Artefact(
        "window_sequence_properties",
        CONTINUOUS / "window_sequence_properties.tsv",
        "per-gene window GC, AT-run density, mappability over 12 distance bands",
        phase="P9",
        required=False,
    ),
    Artefact(
        "oligo_list",
        LAB / "oligo_design_pipeline/panel_ON_A/final_oligo_list.txt",
        "per-probe design properties: GC%, alignments, repeat_length, density_score",
        phase="P9",
        required=False,
        notes=(
            "Previously unused asset flagged 08-03. ON_A is the GW panel; note the "
            "sibling panel_ON_A_RUN1 and panel_ON_A_TEST dirs are NOT it."
        ),
    ),
    # ---- explanation layer ----------------------------------------------
    Artefact(
        "nested_baselines",
        CONTINUOUS / "nested_baselines.tsv",
        "91 features vs n_peaks vs magnitude basis, cross-validated",
        phase="P10",
        notes="The viability result. Leads the explanation layer.",
    ),
    # ---- other views ----------------------------------------------------
    Artefact(
        "ps_fits",
        CONTINUOUS / "ps_fits.tsv",
        "per-gene P(s) exponents: alpha both/near/far/asymmetry, fit R2",
        phase="P7",
    ),
    Artefact(
        "ps_binned_profiles",
        CONTINUOUS / "ps_binned_profiles.npy",
        "log-spaced binned P(s) curves per gene",
        phase="P7",
    ),
    Artefact(
        "ps_gene_ids",
        CONTINUOUS / "ps_gene_ids.json",
        "gene id order for the P(s) arrays",
        phase="P7",
    ),
    Artefact(
        "radar_profiles",
        CONTINUOUS / "archetype_radar_profiles.tsv",
        "archetype median + IQR bands per radar axis",
        phase="P7",
    ),
    Artefact(
        "insulation_targets",
        CONTINUOUS / "hichip_per_gene_targets.tsv",
        "per-gene CTCF HiChIP insulation, r=0.918 across replicates",
        phase="P9",
        required=False,
    ),
    # ---- BYOG (deferred) -------------------------------------------------
    Artefact(
        "frozen_bundle",
        ORCHID / "model/frozen_cd4_bundle.pkl",
        "frozen CD4 scaler, expected curve, feature order, centroids",
        phase="P13",
        required=False,
        notes="Bring-your-own-gene only. Proven by orchid at max |delta| = 0.0.",
    ),
]

# Legacy duplicates that must never be read in preference to the canonical copy.
# scripts_cleaned_legacy/ is dead per CLAUDE.md; anchor_table_v1.tsv is superseded.
STALE_DUPLICATES: dict[str, str] = {
    str(CD4 / "scripts_cleaned_legacy/benchmark_mcc_only/output/anchor_table.tsv"): (
        "legacy tree — use scripts_cleaned/public/outputs/anchor_table.tsv"
    ),
}

# --------------------------------------------------------------------------
# blacklist — build_store refuses to read these
# --------------------------------------------------------------------------

BLACKLIST: dict[str, str] = {
    "archetype_labels_k5.tsv": (
        "Encodes the dead HK/arch-poised split, which did not replicate on GW. "
        "The validated 5th group is the ME split in me_subtypes_labels.tsv."
    ),
}

# Viewpoint positions must come from the BED, never from parsing viewpoint_id.
VIEWPOINT_DIR = MCC_DATA / "viewpoints"


_BY_KEY: dict[str, Artefact] = {a.key: a for a in ARTEFACTS}


def by_key(key: str) -> Artefact:
    """Look up a registered artefact. Raises rather than returning None, so a
    typo fails at the call site instead of surfacing as a confusing None.path."""
    try:
        return _BY_KEY[key]
    except KeyError:
        raise KeyError(
            f"unknown artefact {key!r}; registered: {', '.join(sorted(_BY_KEY))}"
        ) from None


def search_paths() -> list[Path]:
    """Directories the audit globs when an artefact is missing at its recorded path."""
    return [SCRIPTS, CD4, MCC / "collaboration", LAB]
