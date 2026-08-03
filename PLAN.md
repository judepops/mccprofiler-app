# mccprofiler-app — build plan

*Drafted 2026-08-03. Jude Popham, Davies + Yau labs, WIMM Oxford.*

An interactive app for exploring where a gene sits in the MCC regulatory
continuum: live contact profile, P(s) decay, feature radar against archetype
bands, and a confidence-graded archetype readout.

**Scope decisions taken 2026-08-03:**

| decision | choice |
|---|---|
| Gene coverage | Catalogue-first (1,846 GW genes precomputed), API shaped so bring-your-own-gene slots in later |
| Deployment | Localhost only. No data egress. No auth, no hosting story, unpublished MCC never leaves the workstation |
| Scale target | Design the data layer for the 20k panel now (HDF5 + index + multiscale), not for 1,846 |

---

## 1. The design rule that governs everything

**Position first. The label is a confidence-graded readout of position, never
the primary object.**

The project's own results forbid a category badge as the headline output:

- Gap statistic returns **k=1** on GW, with the gap declining monotonically
- HDBSCAN returns **one** cluster at min_cluster_size 25 and 50; BIRCH one at every k
- Dip test is **unimodal** on PC1–5 (min p 0.88–0.96)
- **44%** of active genes are probabilistic mixtures; the ME subtypes are only
  27–44% core
- AA median dominant weight is **0.56**, with only ~4% of genes "pure"

An app that answers "PDCD1 → arch-ME-effector" full stop would contradict the
thesis argument and reproduce exactly the dichotomisation the super-enhancer
critique targets (Altman & Royston 2006; Pott & Lieb 2015).

The honest version costs nothing extra and is a **better demo**, because all
three ingredients already exist on disk: gcPCA density-free coordinates,
per-gene posteriors with entropy, and AA mixture weights. A gene renders as a
*position* with a *mixture composition*; the label appears with its posterior
attached, and genes below the confidence threshold render as "mixture" rather
than being silently rounded to a corner.

This makes the app a live demonstration of the strongest measurement-theory
result in the project: the groups are **not discoverable** (gap k=1) but **are
reproducible once imposed** (Cohen's κ = 0.72 at k=3–4 across independent
captures of the same 116 genes).

Supporting citations for the explanation layer: Rousseeuw 1987 (a very low
average silhouette means seek an alternative to clustering as a model); Hennig
2015, *What are the true clusters?* (no unique true clusters; clustering becomes
scientific through transparent communication, not uniqueness).

---

## 1b. Revised 2026-08-03 — the axes are the product, not the archetypes

The 08-03 session named the real dimensions, and that supersedes the
archetype-centric framing above. **The primary coordinate system is the named
interpretable axes; the archetypes are one derived readout among several.**

`14_name_dimensions.py` → `dimension_names.tsv`. 19 components sit above the
95th-percentile noise ceiling and carry 78% of variance. The leading five are
interpretable:

| axis | reading | var % |
|---|---|---|
| PC1 | richness / reach vs emptiness — **this is the amount axis, label it as such** | 15.3 |
| PC2 | local vs long-range | 9.5 |
| PC3 | **promoter-driven vs enhancer-driven** | 8.2 |
| PC4 | enhancer vs CTCF composition | 6.3 |
| PC5 | concentrated vs dispersed | 5.9 |

**PC3 is the most quotable result in the project**: the classic
housekeeping-vs-developmental axis, recovered from contact architecture alone,
with no chromatin marks and no sequence. Placing a gene on PC3 live is a better
app than any label.

### The six candidate dimensions

Cross-referencing `dimension_reproducibility.tsv` (`explained=False` and
`reproduces=True`, excluding PC1 because it is amount) gives **PC5, PC6, PC9,
PC10, PC13, PC14** — summing to **20.3%** of variance, which matches the
handoff's "~20% of the real signal".

These are amount-independent, cross-capture reproducible, and not explained by
any biology anchor, by 3D insulation, by probe-level sequence or mappability, or
by window-level GC / AT-run density / mappability **including radial gradients**
(max |ρ| 0.069, r² 0.5% — and a radial gradient is the *only* artefact class that
could manufacture a shape feature, so this is the specific falsification and it
came back clean).

**Wording rule, from the handoff verbatim:** say "candidates for regulatory
variation not captured by existing assays". Never "unexplained" — that is
provisional by construction and Phase B is the test. This phrasing goes in the
UI string table, not just the docs.

### Display rule: numbers, never verdicts

The 08-03 decision removed binary survive/compromised classification because
PC5's verdict flipped twice as the threshold moved. **The app inherits this
rule**: show max |ρ| and r² against the largest association the test demonstrably
resolves. No pass/fail badges, no green ticks, no "confound-free" labels.

**PC5 carries a GC caveat that must render inline** wherever PC5 appears: most
GC-associated of the six in both independent tests, same sign, ~6% of variance at
probe level and ~2% at window level, versus <1% for the other five. Sign
agreement alone does not discriminate (4/6 agree, about chance) — what singles
PC5 out is ranking first in both by a 3× margin.

### Why the 91-feature substrate exists at all

`13_nested_baselines.py` answers the viability question and belongs in the
explanation layer, because it is the justification for the whole app:

| comparison | mean gain | targets won | p |
|---|---|---|---|
| 91 features vs `n_peaks` alone | +0.140 | **9/9** | 1.7e-13 |
| 91 vs the 11-feature magnitude basis | +0.034 | **9/9** | 7.5e-06 |

gnomAD LOEUF goes 0.031 → 0.210, carried almost entirely by shape. The nonlinear
arm was **uniformly worse**, which kills the "a linear model just couldn't reach
it" objection.

---

## 2. Architecture — the store is the contract

```
scripts_cleaned/  (science, cd4env)          mccprofiler/  (this repo)
     |                                              |
     |  build_store.py                              |
     |  reads 4.7 GB pickle + tables                |
     +---------------> store/ <--------------------+
                       profiles.h5                  server (mccapp env)
                       tables/*.parquet             never imports mccprofiler
                       manifest.json                never opens the 4.7 GB pickle
```

**The server never imports the science stack and never touches the raw pickle.**
It reads a purpose-built store. Three consequences worth having on purpose:

1. The API starts in milliseconds and has a tiny dependency surface.
2. `cd4env` cannot drift because the app needs a package version bumped.
3. Bring-your-own-gene later is a **worker process in cd4env that writes into
   the store** — not a rewrite of the server. This is the same frozen-bundle
   pattern orchid already proved (`orchid/code/_paths.py`, T1 reproduced CD4
   features bit-for-bit, max |Δ| = 0.0).

**Stack:** React + TypeScript + Vite + Tailwind. FastAPI + uvicorn. HDF5 (h5py)
for profiles, parquet for tables. Plot layer: canvas for the profile (40k points
downsampled), Observable Plot or D3 for radar/scatter, igv.js for the genomic
context panel.

**Cross-repo paths:** this repo sits at `cd4_cleaned/mccprofiler/`, the science
at `cd4_cleaned/scripts_cleaned/`. Use a `_paths.py` with absolute paths, per
project policy (hardcoded `/home/imm/grte4643/...` is deliberate — datasets span
mounts). Do not add `scripts_cleaned` to the app's import path.

---

## 3. Data layer

`build_store.py` runs in `cd4env`, reads the canonical artefacts, and emits a
versioned store. This is the bulk of the real engineering; the UI is the easy half.

### Profiles

Raw: 1,846 genes × 40,000 bins × 8 channels × float32 = **2.4 GB**. Not
servable per-request, and at 20k genes it is 26 GB.

Store as HDF5, **chunked per gene**, with a multiscale pyramid:

| level | bins | bp/bin | use |
|---|---|---|---|
| L0 | 40,000 | 50 | zoomed detail, read by range |
| L1 | 8,000 | 250 | mid zoom |
| L2 | 2,000 | 1,000 | default view |

Default L2 for all 8 channels is ~118 MB total across the panel — negligible.
Zoom requests read an L0 chunk range only for the visible window. This is what
makes 20k genes tractable rather than a rewrite.

Store **raw** and derive O/E and P(s)-residual server-side from the stored
expected curve, rather than storing three copies.

### Tables

All small enough to be parquet and loaded once at startup:

- features (1,846 × 91, raw + z-scored + panel percentile)
- archetype labels, posteriors, entropy, confidence
- gcPCA axes (B1 and B2 backgrounds), AA mixture weights
- P(s) fits: `alpha_both / near / far / asymmetry`, fit R², binned curves
- peaks with element class and signed distance
- viewpoint positions (from the viewpoint BED — **never** parsed from `viewpoint_id`)
- anchor table joined per gene (expression, tau, CpG, conservation, essentiality)
- **named dimensions** — `dimension_names.tsv` (loadings, dominant class, per-anchor ρ),
  `dimension_reproducibility.tsv` (namespaced separately, see trap #10),
  `structure_vs_noise.tsv`
- **confounds** — `locus_intrinsic_confounds.tsv`, `window_confounds.tsv`,
  `window_sequence_properties.tsv`, plus per-probe design properties from
  `panel_ON_A/final_oligo_list.txt`
- **baselines** — `nested_baselines.tsv` for the explanation layer
- cross-panel reproducibility, external-group stratification, HiChIP per-gene
  insulation targets

Everything in `audit/continuous_methods/` except `hichip_insulation_bins.tsv.gz`
(37 MB) is under 600 KB, so the whole tabular layer is trivial to load at startup.

### manifest.json — provenance is not optional here

The store records source path, mtime, git commit of `scripts_cleaned`, panel id,
gene count, and build date for **every** input. There are known-stale artefacts
in the tree (see §7) and a silent pickup would put wrong numbers in front of a
supervisor. The build script **refuses to run** if it resolves a blacklisted
input.

### Phase 0 includes a path audit

I have not verified every input path below still exists at the name recorded in
the handoffs. `00_audit_paths.py` resolves and stats each one, and reports
missing/ambiguous before anything is built:

| artefact | recorded location |
|---|---|
| clean matrix (GW) | `process/output_gw_cd4_1/scripts/outliers/results_combined_50bp_clean.pkl` |
| features (GW) | `mccprofiler/outputs_gw_cd4_1/features/features_zscored.pkl` |
| ME-split labels | `audit/GW/output/me_subtypes_labels.tsv` |
| posteriors | `audit/GW/output/posterior_membership.tsv` |
| gcPCA axes | `audit/continuous_methods/gcpca_axes.npz` |
| P(s) fits | `audit/continuous_methods/ps_fits.tsv`, `ps_binned_profiles.npy` |
| peaks | `annotated.tsv` (GW, 40,987 peaks) |
| viewpoints | `cd4_cleaned/mcc/viewpoints/<source>_viewpoints_final.tsv` |
| radar reference | `audit/Archetype_Tests/07_radar.py` |
| frozen bundle (BYOG) | `orchid/model/frozen_cd4_bundle.pkl` |
| anchor table | path needs confirming — recorded under `Cowork/` but 07-29 reports 1,813/1,846 GW genes resolved locally |

---

## 4. API surface

Shaped now so BYOG is additive, not a refactor. A submitted gene returns a
synthetic id that every `/genes/{id}` endpoint then serves unchanged.

```
GET  /api/panels                              panels + default + provenance
GET  /api/genes?panel=&q=&limit=              search / autocomplete
GET  /api/genes/{id}                          metadata, viewpoint, label + posterior,
                                              coordinates, mixture weights, entropy
GET  /api/genes/{id}/profile
       ?channel=mcc&mode=raw|oe|psresid&level=&start=&end=
GET  /api/genes/{id}/peaks                    element class, signed distance
GET  /api/genes/{id}/features                 91 features, z-scores, percentiles
GET  /api/genes/{id}/ps                       binned curve + alpha fits + regime break
GET  /api/genes/{id}/radar?axes=discriminative|curated
GET  /api/panel/embedding?method=gcpca|aa|pca all-gene 2D coords for the map
GET  /api/panel/archetypes                    centroids, IQR bands, n, prose
GET  /api/bigwig/{source}/{file}              range-request passthrough for igv.js

POST /api/submit                              (phase 2) -> job id -> synthetic gene id
```

---

## 5. Views

1. **Viewpoint profile** — the core, custom-built. ±1 Mb aligned at bin 20,000,
   toggle raw / O/E / P(s)-residual, channel selector, the four feature distance
   bands shaded at their real cutoffs (0–10 kb, 10–50 kb, 50–250 kb, 250 kb–1 Mb),
   peaks marked and coloured by element class. Brush to zoom pulls L0.
2. **Dimension profile — the primary readout (revised 08-03).** The gene's
   position on each named axis as a percentile bar against the panel
   distribution, PC1 explicitly labelled *amount*, PC3 given prominence as
   promoter-driven vs enhancer-driven. Each axis expands to show its top ±
   loading features and its confound numbers.
3. **Continuum map** — all genes in PC or gcPCA space, axis pair selectable
   (default PC2 × PC3, which is the interpretable plane, not PC1 × PC2 which is
   amount-dominated). Selected gene highlighted; archetype regions as soft
   density contours **not** hard boundaries. Colour by posterior confidence so
   the 44% mixture body is visually obvious.
4. **Archetype readout** — demoted to a secondary panel. Posterior bars +
   entropy + AA mixture weights. Below-threshold genes read "mixture", not a corner.
5. **Cohort view (added 08-03 from FEEDBACK.md) — the strongest single view.**
   Pick a built-in reference set or paste a gene list; see that set's distribution
   on each named axis against the panel. This is a live demo of the best
   non-circular result in the project: externally-defined groups separate
   architecture at **5.9–41.1×** over a size-matched random floor and retain
   ≥94% of it within both density and insulation strata. Unlike the archetype
   version it needs no within/pooled-ratio argument, because the groups were not
   defined from the features.

   Three design constraints, all load-bearing:
   - **Report coverage prominently.** The panel is 1,846 of ~20,000 genes, so a
     pasted list of 40 collaborator hits may match 4. Show "7 of your 40 genes
     are in the panel" and degrade honestly rather than plotting 4 points as if
     they were 40. Built-in sets are the primary path until BYOG lands.
   - **Apply the ≥25-gene filter.** Of 23 reference sets, `Roadmap_CTCF_bound`
     has 1 GW gene and `ChromHMM_polycomb` has 7. Enrichments over those are
     noise and must not render with the same confidence as Eisenberg_HK's 633.
   - **Display `pct_retained` and `signal_over_random`, not η².** Pooled η² runs
     0.005–0.023; the absolute value is not the claim, the ratio is. Same lesson
     as the archetype circularity caveat.

   Keep `gene_desert_bottomQ_density` visible as the **positive control** — it is
   defined from gene density and is the one group that collapses under density
   stratification (11% retained) while surviving insulation (100%). It is what
   makes the other six rows meaningful.
6. **Ranked lists** — the inverse of the cohort view: given an axis, which genes
   sit at the extremes. Nearly free once the store exists, and it is how a user
   finds a gene worth looking at when they do not already have one in mind.
7. **Export** — gene × 91 features + axis percentiles + coordinates + archetype
   proximity, for one gene, a list, or the whole panel. The CellProfiler analogy
   only holds if the output is portable; CellProfiler's deliverable *is* a
   feature table. Minimum bar for "tool, not demo".
8. **Reproducibility demo** — the 116 genes measured twice in independent
   captures, both profiles overlaid and both positions on the map. Makes κ=0.72
   and the 0.511 technical noise floor visceral instead of a caption. Cheap:
   `cross_panel_reproducibility.tsv` already holds all 63 features at n=116
   (median ρ 0.752, `frac_local` 0.984, asymmetry features <0.3), so it is a
   table read, not an analysis. **Label it**: the immune side is the stale 791
   baseline missing PDCD1 + 15 genes, so the overlap may shift after the re-run.
   Acceptable for a technical-noise demo, not for a biological claim.
9. **Radar** — gene trace over archetype IQR bands. Port from `07_radar.py`
   (`archetype_radar_profiles.tsv` already on disk). **Keep the bands.** They
   overlap almost completely and that honesty is the point.
10. **P(s)** — log-log with fitted α, the ~100 kb regime break (α_near 0.750,
   α_far 1.291), gene's α against the panel distribution.
11. **Case-study views (added 08-03)** — a better-formatted replacement for the
   `audit/MCCProfiler/output_leiden_k7/J1_case_studies` figures, whose content is
   right but which are unreadable: 2541x4909 px matplotlib dumps with 8 pt type,
   a legend over the data, and a "y-axis capped at 36" note wedged into a corner.
   Four parts, all interactive rather than baked:
   - **Peak annotations on the profile** — coloured by element class, radius on a
     sqrt scale so *area* tracks `peak_max` (a linear radius exaggerates tall
     peaks roughly quadratically). Click selects the nearest peak.
   - **Peak detail** — the tallest-peak zoom panel: FWHM, sharpness, enrichment,
     consensus fraction, class.
   - **Topology arcs** — the peak-peak graph that `mean_degree` and
     `frac_active_pairs` summarise, drawn as arcs. *Not yet built.*
   - **Feature table** — 91 features grouped by block, with **percentile leading
     and z secondary**. The original showed rank as `#45/791`, which needs mental
     arithmetic; "95th percentile" is directly readable. One stated colour scale
     rather than the original's per-cell reference scheme.
12. **Gene comparison** — two genes side by side: profiles, axis positions,
   distance in feature space. The natural follow-up to a lookup.
13. **Confound panel (new, 08-03)** — per dimension, the probe-level and
   window-level association numbers from `locus_intrinsic_confounds.tsv` and
   `window_confounds.tsv`, including the radial-gradient result. Numbers and
   effect sizes only, per §1b. Also exposes the per-probe design properties
   (GC%, alignment count, repeat length, density score) that exist for
   **1846/1846** GW genes in `panel_ON_A/final_oligo_list.txt` — zero fetch, and
   previously unused.
14. **Genomic context** — igv.js, all 8 channels, GENCODE. Deferred, **not cut**
    (see decision D2 below).
15. **Explanation layer** — woven through, not a separate About page. Leads with
    the nested-baseline justification (§1b), then the structure-vs-noise evidence.
    Reuse the `html_dphil/dphil_html.py` house style so it matches the existing
    explainers.
16. **SE comparison** (optional, high rhetorical value) — a gene's SE status
    against its continuum position. **Display-only. Never touches the feature path.**

### The Lab page (added 08-03)

A second page, deliberately outside the product, for exploratory and
methodological views: diagnostics, checks, and anything whose caveats are too
heavy for a page a collaborator might skim. Nothing on it should be quoted
without its caveat attached.

First occupant is the **reproducibility view** — the 116 genes captured in both
panels, per-feature agreement plus a paired scatter against the identity line.
It is the only direct measurement of technical noise available (same gene, two
independent captures, same pipeline) and it pre-empts "how do I know this is
real" better than prose can. It reproduces the handoff exactly: median rho
0.752, 40/63 above 0.7, one below 0.3, with `frac_local` (0.984) at one end and
the `oe_asymmetry_*` family (0.27-0.37) at the other.

It lives in the Lab rather than the app because the immune side is the stale 791
baseline, so it is defensible as a technical demonstration and not as biology.

### Two decisions recorded against FEEDBACK.md

**D1 — the therapeutic framing ships with its effect size or not at all.**
FEEDBACK.md proposes attaching "different CRISPR strategies and redundancy
expectations" to PC3. PC3 carries 8.2% of variance and correlates with tau at
ρ = −0.221; the best R² against any anchor anywhere in the project is 0.077.
Hanging a therapeutic implication on that reproduces the overclaim pattern the
08-03 session retracted twice under challenge, and "the sentence that makes a
clinician care" is a rhetorical goal, not a scientific one. **Resolution:** keep
it, in the explanation layer, on PC3, framed as *motivation* with the number
rendered adjacent — never as an implication of the result.

**D2 — igv.js is deferred, not cut.** FEEDBACK.md argues UCSC/IGV already do
this well. That is an argument against *building a browser*, which igv.js already
satisfies — it is an embed, and the only infrastructure it needs is one
range-serving endpoint. It also has a function the feedback misses: the
viewpoint-anchored ±1 Mb profile is an unfamiliar view carrying no genomic
coordinates, and the first question any viewer asks is "where actually is this?".
igv.js is the orientation anchor that makes the primary view trustworthy.
**Resolution:** demote below cohort and export, keep in scope.

---

## 6. Why not UCSC

Considered and rejected for the primary view, for three independent reasons:

- **It cannot do the primary view at all.** UCSC and IGV are coordinate-anchored;
  MCCProfiler's object is *viewpoint-anchored*, every gene aligned at bin 20,000.
  That alignment is the entire reason genes are comparable. No genome browser
  renders it.
- **All three UCSC routes require the bigWigs on a public URL.** The REST API
  (`api.genome.ucsc.edu`) serves only UCSC-hosted tracks; track hubs and
  `hgRenderTracks` both need UCSC's servers to fetch the files. That is a
  data-release decision on unpublished pre-Transfer data, not a technical choice
  — and it contradicts the localhost-only decision above.
- **`hgRenderTracks` returns a PNG.** Static, not "dynamic live".

**igv.js instead** for the genomic-context panel: reads bigWig over HTTP range
requests, so FastAPI serves the local files with `Accept-Ranges` and nothing
leaves the machine. UCSC survives only as an outbound deep link, which needs no
hosting and no API.

---

## 7. Data traps

1. **Never read `archetype_labels_k5.tsv`.** It encodes the dead HK/arch-poised
   split, which did not replicate on GW. The validated 5th group is the ME split
   (`me_subtypes_labels.tsv`). Blacklisted in the store builder.
2. **Default to GW (1,846).** The 791 immune baseline is missing 16 genes
   including PDCD1 from the chromosome-edge bug and the re-run is still
   outstanding across four handoffs. If the immune panel is exposed at all, label
   it ascertainment-biased (AMI 5.7× vs GW; SigClust ~22 SD at matched n and
   features) and stale.
3. **Viewpoint ≠ TSS.** Load positions via `io.py::load_viewpoints`. Never parse
   the numeric suffix of `viewpoint_id` — it is a per-peak midpoint. For ~2% of
   genes the offset is hundreds of kb (RERE: 292 kb), which renders as subtly
   wrong rather than obviously broken.
4. **`distance_to_viewpoint` is unsigned.** `oe_matrix` *is* strand-flipped while
   `peaks_df.distance_to_viewpoint` is signed-genomic — use `resolve_peak_bin()`.
   Getting this wrong puts peaks on the wrong side of the profile.
5. **ATAC gates, it is never a clustering feature.** The explanation layer must
   say this correctly or it undercuts the independence argument. Clustering reads
   only `combined_matrix`.
6. **Asymmetry features are the least reproducible** (ρ < 0.3 across independent
   captures of the same genes, vs median 0.752). Flag or omit in the radar.
7. **The technical noise floor is real**: gene-to-itself distance 5.36 vs
   gene-to-other 10.49, ratio 0.511. About half a typical between-gene distance
   is technical. Worth surfacing in the UI rather than hiding.
8. **Use trimmed statistics, not raw (08-03).** The `meeting_2026-07-31.html`
   deck overclaims with "5.2× more extreme genes", skew −2.50, excess kurtosis
   16.4 — all outlier-inflated. Trimmed values are **3.7×, −0.84, 1.01**. Note
   `structure_vs_noise.tsv` stores the *untrimmed* PC1 excess kurtosis (16.41),
   so reading that field naively reproduces the overclaim.
9. **SigClust: quote z = −11.86 (trimmed), not −5.38 or −5.26.** Trimming makes
   the result *stronger* — the outliers were inflating the null's variance, so
   the departure from Gaussian is a property of the bulk distribution, not a
   tail artefact.
10. **`dimension_names.tsv` and `dimension_reproducibility.tsv` are keyed on PC
    number but computed in different feature spaces** — PC1 variance reads 15.27
    in one and 14.17 in the other (91-feature vs 63-shared). **Do not join them
    on PC index without recording which space each came from.** The store builder
    must namespace them.
11. **Staging is dead (08-03).** The 2026-07-27 two-stage active/silenced design
    was overturned — active-only loses in all four configurations, and on the
    immune panel it destroys signal (−1.63 SD). Do not build a staged view.

---

## 8. Build phases

| phase | deliverable |
|---|---|
| **P0** | Path audit; `build_store.py`; manifest + provenance guard; blacklist enforcement |
| **P1** | FastAPI skeleton, gene + panel endpoints, React/Vite/Tailwind shell, gene search |
| **P2** | Viewpoint profile plot — the core view, raw/OE toggle, bands, peaks, brush-zoom |
| **P3** | Dimension profile — the primary readout |
| **P4** | **Cohort view** with coverage reporting and the ≥25-gene filter |
| **P5** | Continuum map + archetype readout with posteriors and mixture state |
| **P6** | Export + ranked lists |
| **P7** | Radar with IQR bands; feature table; P(s) panel |
| **P8** | Reproducibility demo (116 twice-captured genes) |
| **P9** | Confound panel |
| **P10** | Explanation layer in `html_dphil` house style |
| **P11** | igv.js genomic context + bigWig range passthrough |
| **P12** | Gene comparison; SE panel (optional) |
| **P13** | BYOG worker in cd4env writing to the store (deferred, API already shaped) |

Commit at phase boundaries. PR at the end of P2 (first genuinely usable slice)
unless asked sooner.

---

## 9. Environments

| env | role |
|---|---|
| `cd4env` | `build_store.py` and the future BYOG worker. Has the science stack + mccprofiler |
| `mccapp` (new) | server only — fastapi, uvicorn, h5py, pyarrow. Deliberately thin |
| node/npm | frontend. `node_modules` gitignored |

Keeping these separate is what stops the app from becoming a reason to touch
`cd4env`.

---

## 10. Open questions

- **Repo rename pending.** The project is named **`mccprofiler-app`** (decided
  2026-08-03) to disambiguate it from the feature package at
  `scripts_cleaned/mccprofiler/` that it imports from. The GitHub repo and local
  directory are still `mccprofiler` — rename both before the first PR:
  `gh repo rename mccprofiler-app` (GitHub redirects the old URL), then move the
  local dir and `git remote set-url`.
- **Whether to expose the immune panel at all**, given trap #2.
- **Whether the archetype prose is written or generated.** The four archetypes
  need one honest paragraph each; these are supervisor-facing and should be
  Jude's words, not mine.
