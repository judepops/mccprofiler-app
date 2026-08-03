"""One line per feature, so the translator can select on features not just axes.

Why this exists. The model was given 91 bare feature names and no idea what any
of them measured, so it fell back on the five named PCs for everything. That is
the wrong default twice over: a PC name covers only 22 to 35 percent of its axis
(see scripts/diagnose_pc_names.py), and plenty of perfectly ordinary questions
name something a single feature measures exactly while no PC expresses it at
all. "Genes with one-tailed MCC" is the clean example. It is
`oe_asymmetry_mean_all` high, and there is no axis for it.

Descriptions come from `feature_geometry.describe()` where it has one, since
that is what the feature explainer already draws and the two must not diverge.
The 24 features it does not cover are written out below. Every line says what
HIGH means, because a feature filter takes a direction and a name alone does not
say which end is which.
"""

from __future__ import annotations

from .feature_geometry import describe

# The 24 with no geometry entry. Phrased as "high means ...", the same shape as
# the generated lines, so the model sees one consistent format.
EXTRA: dict[str, str] = {
    "bimodality_score": "How two-humped the contact distance distribution is. High "
        "means contacts fall into two separated distance regimes rather than one.",
    "contact_asymmetry": "Imbalance of signal between the two sides of the "
        "viewpoint. High means one-sided, a one-tailed profile.",
    "corr_oe_distance": "Correlation between O/E signal and genomic distance. High "
        "means enrichment grows with distance, unusual, most profiles decay.",
    "corr_raw_distance": "Correlation between raw signal and genomic distance. "
        "Normally strongly negative because contacts decay with distance.",
    "ctcf_signal_fraction": "Share of peak signal in CTCF-class peaks, O/E "
        "normalised. High means a CTCF-dominated contact profile.",
    "ctcf_signal_fraction_raw": "As above on raw signal, so it is not corrected for "
        "the distance decay.",
    "distance_to_nearest_peak_bp": "Distance from the viewpoint to the closest "
        "called peak. High means the nearest partner is far away.",
    "dominant_peak_distance_bp": "Distance from the viewpoint to the strongest peak. "
        "High means the main contact partner is distal.",
    "dominant_peak_isolation_bp": "Gap between the strongest peak and its nearest "
        "neighbouring peak. High means that peak stands alone.",
    "enhancer_signal_fraction_raw": "Share of raw peak signal in enhancer-class "
        "peaks. High means enhancer-dominated.",
    "mean_peak_gap_bp": "Average spacing between consecutive peaks. High means peaks "
        "are spread out; low means they are packed together.",
    "median_contact_distance": "Median distance of contact signal from the "
        "viewpoint. High means the profile sits far out, low means it hugs the bait.",
    "n_high_consensus_peaks_075": "Number of peaks called consistently across "
        "replicates at 0.75 consensus. High means many confidently reproducible "
        "partners.",
    "n_peaks_within_100kb": "Peak count inside 100 kb of the viewpoint. High means a "
        "locally busy neighbourhood.",
    "oe_distal_mean": "Mean O/E enrichment in the distal window. High means distal "
        "contacts exceed what distance alone predicts.",
    "oe_distance_kurtosis": "Peakedness of the O/E signal across distance. High means "
        "enrichment concentrates at a narrow band of distances.",
    "oe_max_cv": "Variability of peak O/E maxima. High means peaks differ a lot in "
        "strength; low means they are uniform.",
    "oe_proximal_max": "Strongest O/E enrichment near the viewpoint. High means an "
        "intense local contact.",
    "peak_dominance_index": "How much the top peak outweighs the rest. High means one "
        "peak dominates the profile.",
    "promoter_signal_fraction": "Share of peak signal in promoter-class peaks, O/E "
        "normalised. High means promoter-driven.",
    "promoter_signal_fraction_raw": "As above on raw signal, not distance-corrected.",
    "signal_entropy_distal": "Spread of signal across the distal region only. High "
        "means distal contacts are dispersed over many places rather than focused.",
    "signal_profile_kurtosis": "Peakedness of the whole profile. High means signal "
        "concentrates in a few sharp spikes.",
    "spacing_regularity": "How evenly spaced the peaks are. High means regular, "
        "periodic placement rather than clustered.",
}

# The single most important thing in this file. "Long-range enhancer contacts"
# is one condition about one set of peaks, not two conditions about a gene.
# Two marginal filters return genes that have long-range contacts SOMEWHERE and
# enhancers SOMEWHERE, which they need not be the same peaks, so the answer
# quietly stops being about long-range enhancers at all. 34 of the 91 features
# are element-class specific precisely so the joint question can be asked
# directly.
JOINT_HINTS = """\
JOINT CONDITIONS. Read this before combining filters.

When a question attaches a property to an ELEMENT CLASS (enhancer, promoter,
CTCF), that is ONE condition about those peaks, not two conditions about the
gene. Use the single feature ending in _enhancer / _promoter / _ctcf.

  "long-range enhancer interactions"
      RIGHT  max_distance_to_viewpoint_enhancer high        (one filter)
      WRONG  frac_far_distal high + n_peaks_enhancer high
             That returns genes with long-range contacts somewhere and
             enhancers somewhere. They need not be the same peaks, so the
             result stops being about long-range enhancers.

  "local CTCF contacts"     -> oe_local_enrichment_max_ctcf high, or
                               max_distance_to_viewpoint_ctcf low
  "distant promoter contacts" -> max_distance_to_viewpoint_promoter high
  "one strong enhancer"     -> raw_peak_max_max_enhancer high
  "broad enhancer contact"  -> oe_fwhm_bp_max_enhancer high
  "one-sided enhancer contacts" -> oe_asymmetry_mean_enhancer high
  "reproducible enhancer contacts" -> consensus_fraction_max_enhancer high
  "enhancer-dominated"      -> enhancer_signal_fraction_raw high

The rule generally: if the question says "<property> <element> contacts", look
for a feature whose name contains BOTH the property and the element, and prefer
that one feature over two. Only use two filters when the question genuinely has
two independent conditions, such as "long-range enhancer contacts in genes that
also have many CTCF sites".

Element-specific features exist for these properties:
  distance   max_distance_to_viewpoint_{enhancer,promoter,ctcf}
  count      n_peaks_{enhancer,promoter,ctcf}
  strength   raw_peak_max_max_*, oe_max_max_*, raw_log2_enrichment_max_*
  locality   oe_local_enrichment_max_*
  asymmetry  oe_asymmetry_{mean,max}_*      (one-tailed)
  tailedness oe_tailedness_{mean,max}_*
  width      oe_fwhm_bp_max_{enhancer,ctcf}
  agreement  consensus_fraction_max_{enhancer,ctcf}
"""

# Questions phrased in ordinary words that no axis expresses, mapped to the
# feature that answers them exactly. Listed separately because these are the
# cases where the model reliably reached for a PC and should not have.
PHRASE_HINTS = """\
Phrases that map to a FEATURE and to no axis at all:
  "one-tailed", "one-sided", "asymmetric", "all on one side"
      -> oe_asymmetry_mean_all high (or _ctcf / _enhancer / _promoter for one
         element class, contact_asymmetry for the raw version)
  "two-sided", "balanced", "symmetric"        -> oe_asymmetry_mean_all low
  "long tail", "heavy tailed", "trails off"   -> oe_tailedness_mean_all high
  "one dominant contact", "single partner"    -> peak_dominance_index high, or
                                                 frac_signal_in_top_peak high
  "spread out", "diffuse", "many partners"    -> signal_entropy high
  "sharp spike", "focal", "punctate"          -> signal_profile_kurtosis high
  "evenly spaced", "periodic peaks"           -> spacing_regularity high
  "two distance regimes", "bimodal"           -> bimodality_score high
  "isolated peak", "nothing nearby"           -> dominant_peak_isolation_bp high
  "busy locally", "crowded near the gene"     -> n_peaks_within_100kb high
  "how much signal", "contact amount"         -> total_mcc  (NOT PC1, which
                                                 correlates with it at +0.03)
"""


def _line(name: str) -> str:
    d = describe(name) or {}
    text = d.get("measures") or EXTRA.get(name)
    if not text:
        return f"  {name}"
    text = " ".join(str(text).split())
    if len(text) > 210:
        text = text[:207].rsplit(" ", 1)[0] + "..."
    return f"  {name}: {text}"


def build(features: list[str]) -> str:
    """The reference block injected into the prompt."""
    lines = "\n".join(_line(f) for f in sorted(features))
    return (
        "FEATURE REFERENCE. The query is built from these. Each is one measured "
        "quantity computed from the contact profile, so a feature filter asks for "
        "exactly the thing it names.\n\n"
        f"{JOINT_HINTS}\n"
        f"{PHRASE_HINTS}\n"
        f"All {len(features)} features, and what HIGH means:\n{lines}\n"
    )
