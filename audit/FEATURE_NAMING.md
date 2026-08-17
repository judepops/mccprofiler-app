# Are the features named for what they measure?

Audit of all 73 features against their definitions in `mccprofiler/src/`, written
2026-08-16 after `frac_promoter_proximal` turned out not to be a promoter
feature. It is not the only one.

**Nothing here is a bug.** Every feature computes what its code says. The problem
is that several names invite a reading the code does not support, and at least
one of those misreadings has already reached written work.

---

## 1. The headline: `frac_promoter_proximal` is not about promoters

```python
DISTANCE_BANDS = {"promoter_proximal": (0, 10_000), ...}   # config.py
```

It is the fraction of signal within 10 kb of the **viewpoint**, measured as
`|bin - CENTER_BIN|`. Two separate problems:

1. **It is a distance ring, not an element class.** No promoter annotation is
   consulted. A gene with no promoter-classed peak at all still has a
   `frac_promoter_proximal`.
2. **The viewpoint is not the TSS.** The matrix is anchored on the midpoint of
   the captured region, a per-peak identifier. About 2% of genes differ by
   hundreds of kb (RERE by 292 kb). For those genes this measures a 10 kb ring
   around a point that is not the promoter.

**Rename to `frac_within_10kb_of_viewpoint`.** The honest name is unwieldy, which
is exactly why the misleading one survived.

## 2. Three distance words each mean two different things

This is the more dangerous class of problem, because the features sit side by
side in the same table and look comparable.

### "distal" is both a ring and a half-line

| feature | region | source |
|---|---|---|
| `frac_distal` | **50 to 250 kb** | `DISTANCE_BANDS["distal"]` |
| `local_to_distal_ratio` | **50 to 250 kb** | `_BAND_MASKS["distal"]` |
| `distal_signal_density` | **everything >= 50 kb** | `_DISTAL_MASK` |
| `signal_entropy_distal` | **everything >= 50 kb** | `_DISTAL_MASK` |
| `oe_distal_mean` | check before quoting | |

`distal_signal_density` includes the far-distal band; `frac_distal` excludes it.
They sound like the same region and are not. `local_to_distal_ratio` uses the
ring while `distal_signal_density` uses the half-line, in the same module.

### "proximal" is both 2 kb and 10 kb

| feature | region |
|---|---|
| `oe_proximal_max` | **+/- 2 kb**, the bait region (`BAIT_KB = 2`) |
| `bait_pileup_fraction` | **+/- 2 kb**, the bait region |
| `frac_promoter_proximal` | **0 to 10 kb** |

### "local" is both 10-50 kb and a 100 kb radius

| feature | region |
|---|---|
| `frac_local`, `local_to_distal_ratio` | **10 to 50 kb** (`DISTANCE_BANDS["local"]`) |
| `n_peaks_within_100kb` | **+/- 100 kb** (`LOCAL_RADIUS_KB = 100`) |

The feature name is honest here; the constant behind it is called `LOCAL_RADIUS`
while `local` elsewhere means 10-50 kb. `oe_local_enrichment_*` is a third use,
per-peak background enrichment, unrelated to either. **Check that one before
quoting it as anything to do with local contacts.**

### "band" is two different partitions

| | |
|---|---|
| `DISTANCE_BANDS` | 4 unequal rings: 10 / 50 / 250 / 1000 kb |
| `empty_band_fraction` | 20 equal 50 kb bands per side (`EMPTY_BAND_KB = 50`, `N_BANDS = 40`) |

`empty_band_fraction` is not the fraction of the four named bands that are empty.

## 3. The element-class features mean "any such element in the window"

`promoter_signal_fraction`, `n_peaks_promoter`, `raw_peak_max_max_promoter`,
`oe_local_enrichment_max_promoter`, `raw_log2_enrichment_max_promoter`,
`max_distance_to_viewpoint_promoter`.

Measured 2026-08-16: **only 5.8%** of promoter-classed peaks lie within 5 kb of
the gene's own TSS, median distance **115 kb**. So these are overwhelmingly
*other genes'* promoters. The names read as though they describe the gene's own
promoter and they do not.

That is not a defect for the promoter-assembly question, where other-promoter
contact is the quantity of interest (section 6d). It is a defect for anyone
reading `promoter_signal_fraction` as "how promoter-driven is this gene".

## 4. Names that are actively wrong about what they measure

- **`contact_asymmetry`.** Docstring says "positive means more signal
  downstream", but computed in *genome* coordinates, so downstream is upstream
  for half the genome and the gene-relative signal cancels. Strand-corrected
  mean is 13.7 SE from zero (section 6b). The name promises a directional
  quantity the implementation destroys.
- **`distance_to_viewpoint`** (in `annotated.tsv`, upstream of several
  features). **Unsigned.** Signed placement is `peak_midpoint - viewpoint_pos`.
- **The "topology" block.** After the 2026-08-16 removals it contains
  `spacing_regularity`, `corr_oe_distance`, `oe_max_cv` and `n_active_peaks`.
  None is topology. The genuinely topological features were the degenerate ones
  and they are gone. Renaming the block was already a carried-forward item from
  the 07-31 handoff.

## 5. Names that are fine

`total_mcc`, `max_mcc`, `q90_mcc`, `n_peaks_*`, `n_active_peaks`,
`n_high_consensus_peaks_075`, `consensus_fraction_*`, `raw_peak_max_*`,
`oe_max_*`, `gini_mcc`, `signal_entropy`, `frac_signal_in_top_peak`,
`peak_dominance_index`, `dominant_peak_distance_bp`,
`dominant_peak_isolation_bp`, `mean_peak_gap_bp`, `spacing_regularity`,
`bimodality_score`, `oe_distance_kurtosis`, `median_contact_distance`,
`std_contact_distance`, `corr_oe_distance`, `corr_raw_distance`.

Each says what it does. `bait_pileup_fraction` in particular is correctly named
after the bait rather than a biological region.

---

## What to do

**Do not bulk-rename.** Renaming a feature changes the correlation-prune order
and therefore which features survive, so a rename is a substrate change, not
cosmetics. The prune already demonstrably selects for the wrong thing.

Priority order:

1. **Fix the docstrings first**, especially `contact_asymmetry`, whose docstring
   is factually wrong rather than merely loose.
2. **Rename `frac_promoter_proximal` -> `frac_within_10kb_of_viewpoint`** at the
   next deliberate substrate change, and re-run the prune. This is the one that
   has actually misled.
3. **Pick one meaning for "distal"** and make `distal_signal_density` and
   `frac_distal` agree, or rename one of them.
4. **Rename the topology block** to what it now contains.
5. **Add a glossary to the briefing pack** giving the exact region for every
   distance word. Collaborators read the briefing, not `config.py`.

## The general lesson

Three of these problems come from the same source: **a distance from the
viewpoint was given a biological name.** "promoter_proximal", "local", "distal"
and "proximal" all sound like statements about regulatory architecture and are
all statements about base pairs from an arbitrary anchor. Given the viewpoint is
not the TSS, even the anchor is not what the name implies.

When adding features, name them after the measurement, not the interpretation.
`frac_within_10kb_of_viewpoint` cannot be misread; `frac_promoter_proximal`
already was.
