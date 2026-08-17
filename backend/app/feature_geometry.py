"""What each feature MEASURES, drawn on a real profile.

The feature table can say a gene is at the 99.9th percentile of
`bait_pileup_fraction` without ever showing what that sums. This module maps a
feature name to a drawable description, which region of the +/-1 Mb window it
reads, and how, so the app can shade that region on the gene's own trace and
alongside a contrasting gene.

Extends the spatial-zones schematic in
`audit/MCCProfiler/report/report.tex` from a static diagram of the cutoffs to
an overlay on actual signal.

**Honesty rule:** a feature that cannot be drawn faithfully gets `kind: None`
and a reason, rather than a plausible-looking overlay that misrepresents it.
The topology block is the clearest case, its adjacency is the complete graph
on active peaks (`topology.py:104`), so an arc diagram would imply measured
interaction structure that does not exist.
"""

from __future__ import annotations

# Distance bands, from mccprofiler.config.DISTANCE_BANDS. Mirrored, since every
# band is defined on |distance from viewpoint|.
BANDS = {
    "bait": (0, 2_000),
    "viewpoint_proximal": (0, 10_000),
    "local": (10_000, 50_000),
    "distal": (50_000, 250_000),
    "far_distal": (250_000, 1_000_000),
}


def _band(name: str, what: str, band: str) -> dict:
    lo, hi = BANDS[band]
    return {
        "kind": "band",
        "regions": [{"start_bp": lo, "end_bp": hi, "label": band, "mirrored": True}],
        "measures": what,
        "reads": f"the shaded region, as a fraction of the whole profile",
    }


# Only features whose geometry is unambiguous appear here. Everything else
# falls through to `describe()`'s honest default.
GEOMETRY: dict[str, dict] = {
    # -- G2 spatial allocation: the clearest case, a band fraction -----------
    "frac_viewpoint_proximal": _band(
        "frac_viewpoint_proximal",
        "Share of all contact signal falling within 10 kb of the viewpoint. High "
        "means bait-dominated, typical of silenced genes where almost everything "
        "sits at the viewpoint itself.",
        "viewpoint_proximal",
    ),
    "frac_local": _band(
        "frac_local",
        "Share of signal in the 10-50 kb ring: the immediate enhancer landscape "
        "around the promoter.",
        "local",
    ),
    "frac_distal": _band(
        "frac_distal",
        "Share of signal at 50-250 kb, mid-range contact, multi-enhancer "
        "architecture or a moderately extended domain.",
        "distal",
    ),
    "frac_far_distal": _band(
        "frac_far_distal",
        "Share of signal beyond 250 kb. Long-range contact, or noise, the two are "
        "not distinguished here.",
        "far_distal",
    ),
    "bait_pileup_fraction": _band(
        "bait_pileup_fraction",
        "Share of signal trapped in the +/-2 kb bait region. High means the capture "
        "dominates the profile and little signal is informative about regulation.",
        "bait",
    ),
    "distal_signal_density": {
        "kind": "band",
        "regions": [{"start_bp": 50_000, "end_bp": 1_000_000, "label": "distal+",
                     "mirrored": True}],
        "measures": "Mean signal per bin beyond 50 kb, log-compressed. Unlike the "
                    "band fractions this is a density, not a share, it does not "
                    "fall when the bait pile-up grows.",
        "reads": "average height across the shaded region",
    },
    "local_to_distal_ratio": {
        "kind": "two_band",
        "regions": [
            {"start_bp": 10_000, "end_bp": 50_000, "label": "local", "mirrored": True,
             "role": "numerator"},
            {"start_bp": 50_000, "end_bp": 250_000, "label": "distal", "mirrored": True,
             "role": "denominator"},
        ],
        "measures": "log1p(local) minus log1p(distal). Positive means signal is held "
                    "close; negative means it reaches out past 50 kb.",
        "reads": "the two shaded regions, compared",
    },

    # -- G1 intensity: thresholds on the height axis -------------------------
    "total_mcc": {
        "kind": "whole",
        "measures": "Sum of every bin. The crudest amount descriptor, and the one "
                    "most confounded by capture efficiency, it is log1p'd before "
                    "z-scoring and conditioned on downstream.",
        "reads": "the total area under the trace",
    },
    "max_mcc": {
        "kind": "threshold",
        "quantile": None,
        "measures": "The single tallest bin. Almost always at or beside the bait.",
        "reads": "the height of the tallest point",
    },
    "q90_mcc": {
        "kind": "threshold",
        "quantile": 0.90,
        "measures": "The 90th percentile of bin heights. High means much of the "
                    "window carries signal, not just the bait.",
        "reads": "the height at which 10% of bins are taller",
    },
    "q99_mcc": {
        "kind": "threshold",
        "quantile": 0.99,
        "measures": "The 99th percentile of bin heights, whether the brightest 1% "
                    "of the window stands above background.",
        "reads": "the height at which 1% of bins are taller",
    },

    # -- G4 shape ------------------------------------------------------------
    "signal_entropy": {
        "kind": "whole",
        "measures": "Shannon entropy of the signal treated as a distribution over "
                    "bins. High means contact is spread across the window; low "
                    "means it is concentrated in a few places.",
        "reads": "how evenly the trace is spread, not how tall it is",
    },
    "empty_band_fraction": {
        "kind": "empty",
        "measures": "Fraction of the window carrying no signal at all. The direct "
                    "complement of reach.",
        "reads": "the gaps, bins at zero",
    },
    "frac_signal_in_top_peak": {
        "kind": "top_peak",
        "measures": "Share of all signal in the single strongest peak. High means "
                    "one contact dominates the locus.",
        "reads": "the tallest peak, against everything else",
    },
    "gini_mcc": {
        "kind": "whole",
        "measures": "Gini coefficient of bin heights. 0 = every bin equal, 1 = all "
                    "signal in one bin. An inequality measure, not a location one.",
        "reads": "how unequally signal is distributed across bins",
    },

    # -- G3 distance moments -------------------------------------------------
    "mean_contact_distance": {
        "kind": "moment",
        "measures": "Signal-weighted mean distance from the viewpoint, with the "
                    "bait excluded so the pile-up does not dominate it.",
        "reads": "the centre of mass of the trace, bait removed",
    },
    "std_contact_distance": {
        "kind": "moment",
        "measures": "Spread of the signal-weighted distance distribution. Large "
                    "means contact at many scales at once.",
        "reads": "how far the trace spreads around its centre of mass",
    },
}

# Peak-class features share one geometry: highlight peaks of that class.
_CLASS_HINT = {
    "enhancer": "enhancer-class peaks",
    "ctcf": "CTCF-class peaks",
    "promoter": "promoter-class peaks",
    "all": "all called peaks",
}


def describe(name: str) -> dict:
    """Return a drawable description for a feature, or an honest refusal."""
    if name in GEOMETRY:
        return {"feature": name, "drawable": True, **GEOMETRY[name]}

    # Peak-class aggregations: drawable as "these peaks, highlighted".
    for cls, hint in _CLASS_HINT.items():
        if name.endswith(f"_{cls}") or name.startswith(f"n_peaks_{cls}"):
            return {
                "feature": name,
                "drawable": True,
                "kind": "peak_class",
                "element_class": None if cls == "all" else cls,
                "measures": f"Computed over {hint} only. The peak arm runs every "
                            f"aggregation once per element class, so this is the "
                            f"{cls} slice of it.",
                "reads": f"the highlighted {hint}",
            }

    # Topology: deliberately NOT drawn. See module docstring.
    if any(name.startswith(p) for p in
           ("mean_degree", "max_degree", "n_active_pairs", "frac_active_pairs",
            "n_isolates", "n_active_peaks")):
        return {
            "feature": name,
            "drawable": False,
            "kind": None,
            "measures": "Counts over the peak-peak graph.",
            "why_not_drawable":
                "The adjacency this is computed from connects every pair of active "
                "peaks unconditionally (topology.py:104), it is the complete graph, "
                "with no contact or distance criterion. So this reduces to an "
                "algebraic function of two counts, and drawing it as a network would "
                "imply interaction structure the data does not measure.",
        }

    return {
        "feature": name,
        "drawable": False,
        "kind": None,
        "measures": None,
        "why_not_drawable":
            "No unambiguous geometry defined for this feature yet. Rather than draw "
            "a guess, nothing is shown, see FEATURE_SPEC.md for the definition.",
    }
