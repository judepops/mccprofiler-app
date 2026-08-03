# mccprofiler-app

Interactive explorer for where a gene sits in the MCC regulatory continuum —
viewpoint-aligned contact profile, named architecture axes, P(s) decay, cohort
comparison, and a confidence-graded archetype readout.

Human CD4+ T cells, Micro Capture-C. Jude Popham, Davies + Yau labs, WIMM Oxford.

See [PLAN.md](PLAN.md) for the design, the data traps, and the phase order.

## The one rule

**Position first. The label is a confidence-graded readout of position, never
the primary object.** The gap statistic returns k=1, HDBSCAN returns one
cluster, the dip test is unimodal, and 44% of active genes are mixtures — so a
bare category badge would contradict the science it is meant to display. Genes
render as coordinates with a mixture composition; labels carry their posterior.

Corollary: **every figure is data-backed.** The API serves numbers and the
frontend draws them. No pre-rendered image is ever served — a static PNG cannot
be brushed, zoomed, filtered or recoloured, and those interactions are the point.

## Layout

```
backend/
  app/
    paths.py          absolute path registry — the only place a path is defined
    store_schema.py   the store contract (channels, pyramid, tables)
    store.py          read-side access; the server touches nothing else
    main.py           FastAPI
  scripts/
    00_audit_paths.py resolve and check every input, before building
    build_store.py    read the 4.4 GB pickle once, emit the store
store/                built locally, gitignored
```

The store is the contract. `build_store.py` runs in `cd4env` and may import
`mccprofiler`; **the server never does**, never opens the 4.4 GB pickle, and does
not know where the science tree lives. Bring-your-own-gene later becomes a
`cd4env` worker writing into the store, not a server rewrite.

## Setup

```bash
# 1. check every input resolves (pure stdlib, any python3)
python3 backend/scripts/00_audit_paths.py

# 2. build the store — reads the clean pickle, ~70 s, writes ~1.6 GB
conda activate cd4env
python backend/scripts/build_store.py

# 3. serve
conda activate mccapp
cd backend && uvicorn app.main:app --reload --port 8000
```

`mccapp` is deliberately thin: fastapi, uvicorn, h5py, pandas, pyarrow, numpy.
Keeping it separate from `cd4env` is what stops the app becoming a reason to
touch the pipeline environment.

## Store

| | |
|---|---|
| Panel | `gw_cd4_1`, 1,846 genes after QC and outlier removal |
| Channels | 8 — `mcc` is the clustering substrate; `atac` **gates membership and is never a clustering feature**; the rest are validation |
| Profiles | `profiles.h5`, 3-level pyramid at 50 bp / 250 bp / 1 kb, one gene per chunk, float32, gzip |
| Tables | parquet, largest is 18,802 × 42 |
| Provenance | `manifest.json` stamps every input with size, mtime and the `scripts_cleaned` commit |

Genes are anchored on the **experimental viewpoint midpoint, not the canonical
TSS** — for ~2% of genes these differ by hundreds of kb. Viewpoint positions are
loaded from the BED, never parsed from `viewpoint_id`.

## Endpoints

```
GET /api/health                      store version, panel, provenance
GET /api/provenance                  every input + the blacklist + channel roles
GET /api/genes?q=&group=&limit=      search
GET /api/genes/{gene}                label, posterior, mixture, P(s) fit
GET /api/genes/{gene}/profile        ?channel=&mode=raw|oe&level=&start_bp=&end_bp=
GET /api/dimensions                  named axes + reproducibility (separate spaces)
GET /api/cohorts                     external gene groups vs random floor
```
