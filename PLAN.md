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
2. **Continuum map** — all genes in gcPCA space, selected gene highlighted,
   archetype regions as soft density contours **not** hard boundaries. Colour by
   posterior confidence so the 44% mixture body is visually obvious.
3. **Archetype readout** — posterior bars + entropy + AA mixture weights.
   Below-threshold genes read "mixture", not a corner.
4. **Radar** — gene trace over archetype IQR bands. Port from `07_radar.py`.
   **Keep the bands.** They overlap almost completely and that honesty is the point.
5. **P(s)** — log-log with fitted α, the ~100 kb regime break (α_near 0.750,
   α_far 1.291), gene's α against the panel distribution.
6. **Feature table** — 91 features with panel percentile and block grouping.
7. **Genomic context** — igv.js, all 8 channels, GENCODE. Secondary view.
8. **Explanation layer** — woven through, not a separate About page. Reuse the
   `html_dphil/dphil_html.py` house style so it matches the existing explainers.
9. **SE comparison** (optional, high rhetorical value) — a gene's SE status
   against its continuum position. **Display-only. Never touches the feature path.**

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

---

## 8. Build phases

| phase | deliverable |
|---|---|
| **P0** | Path audit; `build_store.py`; manifest + provenance guard; blacklist enforcement |
| **P1** | FastAPI skeleton, gene + panel endpoints, React/Vite/Tailwind shell, gene search |
| **P2** | Viewpoint profile plot — the core view, raw/OE toggle, bands, peaks, brush-zoom |
| **P3** | Continuum map + archetype readout with posteriors and mixture state |
| **P4** | Radar with IQR bands; feature table |
| **P5** | P(s) panel |
| **P6** | igv.js genomic context + bigWig range passthrough |
| **P7** | Explanation layer in `html_dphil` house style |
| **P8** | SE comparison panel (optional) |
| **P9** | BYOG worker in cd4env writing to the store (deferred, API already shaped) |

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
