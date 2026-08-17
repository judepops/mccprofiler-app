"""FastAPI app, serves the store, nothing else.

Localhost only by decision (PLAN.md scope table): unpublished MCC never leaves
the workstation. No auth, because there is no network surface to protect.

Every endpoint returns data, never a rendered image. The frontend draws.

    conda activate mccapp
    uvicorn app.main:app --reload --port 8000    # from backend/
"""

from __future__ import annotations

from functools import lru_cache
from typing import Literal

import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.middleware.cors import CORSMiddleware

from . import env as _env
from . import store_schema as S

# Before anything reads credentials.
_LOADED_ENV = _env.load()
from . import translate as _translate
from .feature_geometry import describe as describe_feature
from .feature_reference import build as build_feature_reference
from .query import QueryError, query_schema, run as run_query
from .store import StoreMissing, get_store

# Feature distance bands, mirrored from mccprofiler.config.DISTANCE_BANDS so the
# plot shades the same cutoffs the features are actually defined on.
DISTANCE_BANDS = {
    "viewpoint_proximal": (0, 10_000),
    "local": (10_000, 50_000),
    "distal": (50_000, 250_000),
    "far_distal": (250_000, 1_000_000),
}

app = FastAPI(
    title="mccprofiler-app",
    description="Gene position in the MCC regulatory continuum.",
    version="0.1.0",
)

# Vite dev server. Localhost only, see module docstring.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    # POST is required: /api/ask, /api/query and /api/cohorts/compare all take a
    # body. This read GET only, which passed every curl test (curl sends no
    # preflight) while failing in the browser, where the OPTIONS preflight 400s
    # and fetch reports the generic "Failed to fetch".
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)


def _distinct_axes(x: str, y: str) -> None:
    """Reject x == y with a 400 instead of dying with a 500.

    Selecting the same column twice makes `emb[["gene_id", x, y]]` produce a
    DUPLICATE column, after which `df[x]` returns a DataFrame rather than a
    Series and the float conversion raises. The frontend now swaps rather than
    duplicating, but a hand-typed URL should get a usable message and a plot
    against itself is a diagonal line carrying no information anyway.
    """
    if x == y:
        raise HTTPException(
            400,
            f"x and y are both {x!r}; a plane needs two different axes "
            f"(plotting an axis against itself gives a diagonal line).",
        )


def _embeddings_all(s_):
    """Raw PCs and the amount-corrected sPCs in one frame, joined on gene_id.

    Added 2026-08-16. Several endpoints read `embeddings` directly and so were
    blind to the shape axes the app now displays by default: the enrichment grid
    answered "unknown axes 'spc1' or 'spc2'" for the page's own default plane.
    Anything that resolves an axis key should go through here.
    """
    emb = s_.table("embeddings")
    if s_.has("shape_embeddings"):
        sh = s_.table("shape_embeddings")
        emb = emb.merge(sh, on="gene_id", how="left")
    return emb


def store():
    try:
        return get_store()
    except StoreMissing as e:
        raise HTTPException(503, str(e)) from None


def _clean(obj):
    """NaN is not valid JSON. Convert to None so the client sees a gap, not a 500."""
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if np.isnan(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


# ---------------------------------------------------------------------------
# meta
# ---------------------------------------------------------------------------


@app.get("/api/health")
def health():
    s = store()
    m = s.manifest
    return {
        "ok": True,
        "store_version": m["store_version"],
        "built": m["built"],
        "panel": m["panel"],
        "n_genes": len(s.genes),
        "channels": s.channels,
        "scripts_cleaned_commit": m.get("scripts_cleaned_commit"),
    }


@app.get("/api/panel-bias")
def panel_bias():
    """What the panel is and is not. Rendered at the top of every page.

    Served rather than hardcoded in the frontend so the numbers have one home,
    and so the per-claim status table travels with them: the bias threatens some
    conclusions and not others, and a bare "the panel is biased" banner invites
    a reader to discount everything or nothing.
    """
    return S.PANEL_BIAS


@app.get("/api/provenance")
def provenance():
    """Every input with size, mtime and the commit it was built from.

    Exposed rather than buried: several artefacts in this tree have stale twins
    (see the blacklist), so a viewer should be able to check what they are
    looking at.
    """
    m = store().manifest
    return {
        "built": m["built"],
        "scripts_cleaned_commit": m.get("scripts_cleaned_commit"),
        "app_commit": m.get("app_commit"),
        "inputs": m["inputs"],
        "blacklist": m["blacklist"],
        "channel_roles": S.CHANNEL_ROLE,
    }


# ---------------------------------------------------------------------------
# genes
# ---------------------------------------------------------------------------


@app.get("/api/genes")
def list_genes(
    q: str | None = Query(None, description="symbol prefix or substring"),
    group: str | None = None,
    limit: int = Query(50, le=500),
):
    s = store()
    g = s.genes
    if q:
        u = q.upper()
        starts = g[g["gene_symbol"].str.upper().str.startswith(u)]
        contains = g[g["gene_symbol"].str.upper().str.contains(u, regex=False)]
        g = pd.concat([starts, contains]).drop_duplicates(subset="gene_id")
    if group:
        g = g[g["group"] == group]

    cols = [c for c in ("gene_id", "gene_symbol", "group", "max_posterior",
                        "confidence_class", "viewpoint_chrom", "viewpoint_pos")
            if c in g.columns]
    return {"n": int(len(g)), "genes": _clean(g[cols].head(limit).to_dict("records"))}


@app.get("/api/genes/{gene}")
def get_gene(gene: str):
    s = store()
    gid = s.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")

    row = s.genes[s.genes["gene_id"] == gid].iloc[0].to_dict()
    symbol = row["gene_symbol"]
    mix = s.posterior_mix(symbol)

    # The label is a readout of position, never the primary object (PLAN.md §1).
    # `is_mixture` is surfaced at the top level so a client cannot render a bare
    # badge without also having the caveat to hand.
    canonical = row.get("group")
    display = S.ARCHETYPE_DISPLAY.get(canonical, {})

    out = {
        "gene_id": gid,
        "gene_symbol": symbol,
        "viewpoint": {"chrom": row.get("viewpoint_chrom"), "pos": row.get("viewpoint_pos")},
        "archetype": {
            "group": canonical,
            "display": display.get("display"),
            "architecture": display.get("architecture"),
            "caveat": S.CONTINUUM_CAVEAT,
            # Travels with the label, like the caveat, so a client cannot render
            # a confident badge without the staleness to hand.
            "provisional": True,
            "provisional_note": S.ARCHETYPE_PROVISIONAL,
            "max_posterior": row.get("max_posterior"),
            "confidence_class": row.get("confidence_class"),
            "is_core": bool(row.get("is_core", False)),
            "is_mixture": not bool(row.get("is_core", False)),
            "entropy": row.get("assignment_entropy"),
            "mixture": mix,
            "is_qc_class": row.get("group") == "arch-off",
        },
    }

    if s.has("ps_fits"):
        ps = s.table("ps_fits")
        key = "gene_id" if "gene_id" in ps.columns else ps.columns[0]
        hit = ps[ps[key].astype(str) == gid]
        if not hit.empty:
            out["ps"] = _clean(hit.iloc[0].to_dict())
            # A bare exponent is unreadable. alpha_far = 0.204 looks like a
            # finding until you know the panel median is 1.291, at which point it
            # is a near-flat outlier and probably a noise floor rather than
            # unusually long-range contact. Same median-plus-percentile pattern
            # as the contact table.
            out["ps_reference"] = _ps_reference()

    return _clean(out)


@lru_cache(maxsize=1)
def _ps_panel() -> dict:
    """Panel distribution of every P(s) exponent, computed once."""
    s_ = store()
    if not s_.has("ps_fits"):
        return {}
    ps = s_.table("ps_fits")
    out = {}
    for c in ps.columns:
        if c in ("gene_id", "symbol_key"):
            continue
        v = pd.to_numeric(ps[c], errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size:
            out[c] = v
    return out


def _ps_reference() -> dict:
    """Median, p10, p90 and this-gene percentile for each exponent."""
    panel = _ps_panel()
    return {
        c: {
            "median": round(float(np.median(v)), 3),
            "p10": round(float(np.quantile(v, 0.10)), 3),
            "p90": round(float(np.quantile(v, 0.90)), 3),
        }
        for c, v in panel.items()
    }


@app.get("/api/ps-percentile/{gene}")
def ps_percentile(gene: str):
    """Where one gene sits in the panel for each P(s) exponent."""
    s_ = store()
    gid = s_.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")
    ps = s_.table("ps_fits")
    key = "gene_id" if "gene_id" in ps.columns else ps.columns[0]
    hit = ps[ps[key].astype(str) == gid]
    if hit.empty:
        raise HTTPException(404, "no P(s) fit for this gene")
    row = hit.iloc[0]
    panel = _ps_panel()
    return _clean({
        c: round(float((v < float(row[c])).mean() * 100.0), 1)
        for c, v in panel.items()
        if pd.notna(row.get(c))
    })


@app.get("/api/genes/{gene}/profile")
def get_profile(
    gene: str,
    channel: str = S.PRIMARY_CHANNEL,
    mode: Literal["raw", "oe"] = "raw",
    level: int = S.DEFAULT_LEVEL,
    start_bp: int | None = None,
    end_bp: int | None = None,
):
    """Viewpoint-aligned contact profile.

    Distances are signed bp from the experimental viewpoint at the window
    centre. This is NOT a genomic-coordinate view, the alignment is what makes
    genes comparable, and no genome browser can reproduce it.
    """
    s = store()
    gid = s.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")
    try:
        p = s.profile(gid, channel=channel, level=level, mode=mode,
                      start_bp=start_bp, end_bp=end_bp)
    except (KeyError, ValueError) as e:
        raise HTTPException(400, str(e)) from None

    p["bands"] = [
        {"name": name, "start_bp": lo, "end_bp": hi}
        for name, (lo, hi) in DISTANCE_BANDS.items()
    ]
    p["channel_role"] = S.CHANNEL_ROLE.get(channel)
    return _clean(p)


# ---------------------------------------------------------------------------
# dimensions, the primary coordinate system (PLAN.md §1b)
# ---------------------------------------------------------------------------


@app.get("/api/dimensions")
def dimensions():
    s = store()
    dims = s.table("dimensions")
    rep = s.table("dimension_reproducibility") if s.has("dimension_reproducibility") else None

    out = _clean(dims.to_dict("records"))

    # Deliberately NOT joined on PC index: the two tables are computed in
    # different feature spaces (91 vs 63-shared) and their variance shares
    # disagree. Returned side by side, each labelled with its space.
    return {
        "feature_space": "full91",
        "dimensions": out,
        "reproducibility": {
            "feature_space": "shared63",
            "note": "Computed in the 63-shared-feature space. Not aligned to the "
                    "above by index, variance shares differ between spaces.",
            "rows": _clean(rep.to_dict("records")) if rep is not None else [],
        },
        "display_rule": "Report max |rho| and r2. No survive/compromised verdicts, "
                        "a threshold-dependent verdict is a claim that has to be "
                        "defended, and PC5's flipped twice on threshold choice.",
    }


@app.get("/api/genes/{gene}/peaks")
def get_peaks(gene: str):
    """Pre-called peaks for one gene, positioned by SIGNED offset.

    The source table's `distance_to_viewpoint` is unsigned; offset_bp is
    computed in the store as peak_midpoint - viewpoint_pos and verified against
    the unsigned column, so peaks land on the correct side of the profile.
    """
    s = store()
    gid = s.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")
    if not s.has("peaks"):
        return {"peaks": [], "note": "peaks table not built"}

    pk = s.table("peaks")
    hit = pk[pk["symbol_key"] == gid.upper()].copy()
    hit = hit.sort_values("peak_max", ascending=False)
    return {
        "gene_id": gid,
        "n": int(len(hit)),
        "by_class": {k: int(v) for k, v in hit["re"].value_counts().items()},
        "colors": S.ELEMENT_COLOR,
        "peaks": _clean(hit.to_dict("records")),
    }


@app.get("/api/genes/{gene}/features")
def get_features(gene: str):
    """All 91 features for one gene, grouped by block, with panel percentile.

    Percentile is the point: a raw z means little on its own, but "97th
    percentile of the panel on far-distal fraction" is directly readable.
    """
    s = store()
    gid = s.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")
    if not s.has("features"):
        raise HTTPException(503, "features table not built")

    f = s.table("features")
    hit = f[f["symbol_key"] == gid.upper()]
    if hit.empty:
        raise HTTPException(404, f"no features for {gene!r}")

    blocks: dict[str, dict] = {}
    for _, r in hit.iterrows():
        code, desc = S.feature_block(str(r["feature"]))
        b = blocks.setdefault(code, {"block": code, "description": desc, "features": []})
        b["features"].append({
            "name": r["feature"],
            "z": float(r["z"]),
            "percentile": float(r["percentile"]),
        })

    for b in blocks.values():
        b["features"].sort(key=lambda x: -abs(x["z"]))

    order = [c for c, _, _ in S.FEATURE_BLOCKS] + ["other"]
    return {
        "gene_id": gid,
        "n_features": int(len(hit)),
        "blocks": [blocks[c] for c in order if c in blocks],
    }


@app.get("/api/dimensions/scree")
def scree(space: Literal["corrected", "raw"] = "corrected"):
    """Variance per component against a permuted-noise ceiling.

    The ceiling is PCA on a per-feature-permuted copy: the variance a component
    of this size explains when there is provably nothing to find. Components
    above it are the "real dimensions".

    `space` must match whatever the caller is showing loadings for. It defaults
    to the corrected axes, like /api/dimensions/{pc}/loadings: the panel that
    consumes both was previously drawing a RAW scree beside CORRECTED loadings,
    so "PC1, 17.23%" sat next to sPC1's features. Mixing the two is worse than
    either alone because nothing on screen says they disagree.
    """
    s_ = store()
    tbl = "shape_scree" if space == "corrected" else "pc_scree"
    if not s_.has(tbl):
        tbl = "pc_scree"
    if not s_.has(tbl):
        raise HTTPException(503, "no scree table built")
    df = s_.table(tbl)
    corrected = tbl == "shape_scree"
    return _clean({
        "space": "amount-corrected" if corrected else "raw (provenance only)",
        "prefix": "sPC" if corrected else "PC",
        "n_above_noise": int(df["above_noise"].sum()),
        "method": "parallel analysis against a per-feature-permuted null",
        "note": "Components above the ceiling carry structure that survives "
                "destroying all feature-feature covariance while keeping every "
                "marginal distribution intact."
                + (" Magnitude is projected out before this PCA, so no component "
                   "here carries overall signal." if corrected else
                   " These are the raw components; PC1 and PC3 both carry "
                   "magnitude on the current substrate."),
        "rows": df.to_dict("records"),
    })


@lru_cache(maxsize=4)
def _feature_reference(features: tuple[str, ...]) -> str:
    """Cached: the feature list is fixed for the life of the store.

    Takes a tuple, not a list. lru_cache hashes its arguments, so a list here
    raises `unhashable type` on every call.
    """
    return build_feature_reference(list(features))


def _annotate_cohorts(s_, genes: list[dict]) -> list[dict]:
    """Attach external set membership to an already-selected gene list.

    Annotation, never selection. Super-enhancer membership is included because
    seeing it here is exactly the post-hoc comparison it is licensed for, and
    excluding it would hide the one label the project most wants to argue
    against.
    """
    if not genes or not s_.has("cohort_membership"):
        return genes
    cm = s_.table("cohort_membership")
    by_gene: dict[str, list[str]] = {}
    for sym, grp in zip(cm["symbol_key"], cm["group"]):
        by_gene.setdefault(sym, []).append(grp)
    out = []
    for g in genes:
        key = (g.get("gene_symbol") or "").upper()
        out.append({**g, "in_sets": sorted(by_gene.get(key, []))})
    return out


def _overlap(d: float) -> float:
    """Overlapping coefficient of two equal-variance normals separated by d.

    The number that matters for reading these panels. Every external set clears
    significance against 1,846 genes, so a p-value says almost nothing; how much
    the two distributions actually overlap says whether you could ever see the
    difference by looking.
    """
    from math import erf, sqrt
    return 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(d) / (2.0 * sqrt(2.0)))))


@app.get("/api/enrichment/grid")
def enrichment_grid(
    x: str = "pc2",
    y: str = "pc3",
    bins: int = Query(14, ge=6, le=30),
    min_n: int = Query(8, ge=3),
    n_perm: int = Query(200, ge=0, le=1000),
):
    """Where each external reference set sits on the map, as small multiples.

    Post-hoc by construction and therefore legitimate: the embedding is built
    from MCC features alone, and these sets are painted on afterwards. They
    never entered the coordinates. That is the same licence the super-enhancer
    comparison has, and the opposite of selecting genes by a label and then
    describing their architecture.

    Three things make the difference between a figure and a misleading figure,
    all of them here.

    `min_n` floors the count per cell. A 2D embedding is dense in the middle and
    sparse at the edges, so without a floor the strongest apparent enrichments
    are cells holding two genes.

    `n_perm` builds a LABEL-permutation null, which is the control this figure
    actually needs. The UMAP null already in the app asks whether the embedding
    has structure; it says nothing about whether an enrichment is real. Shuffling
    set membership across genes, keeping the set size fixed, gives the enrichment
    magnitude that arises by chance on this exact grid. `threshold` is the 95th
    percentile of the per-permutation maximum |log2|, so it is corrected for
    scanning many cells rather than being a per-cell p-value.

    Enrichment is log2(observed rate / panel base rate), with a half-count added
    so an empty cell is a finite number rather than negative infinity.
    """
    s_ = store()
    if not s_.has("embeddings") or not s_.has("cohort_membership"):
        raise HTTPException(503, "embeddings or cohort_membership not built")

    _distinct_axes(x, y)
    emb = _embeddings_all(s_)
    if x not in emb.columns or y not in emb.columns:
        raise HTTPException(404, f"unknown axes {x!r} or {y!r}")

    genes = s_.genes[["gene_id", "symbol_key"]].merge(
        emb[["gene_id", x, y]], on="gene_id", how="left"
    ).dropna(subset=[x, y])

    xs = genes[x].to_numpy(float)
    ys = genes[y].to_numpy(float)

    # EQUAL-WIDTH edges over the robust range, deliberately not percentile edges.
    # Percentile edges give every cell the same gene count, which is better
    # behaved statistically but renders the plot in rank space: the shape of the
    # cloud is destroyed by construction and the picture no longer looks like
    # the map it claims to describe. Since the grid is drawn on top of the real
    # scatter, the cells must correspond to real coordinates. The 1st to 99th
    # percentile range keeps a handful of extreme genes from stretching the grid
    # to the point where everything lands in one cell.
    def edges(v: np.ndarray) -> np.ndarray:
        lo, hi = np.percentile(v, [1, 99])
        if hi <= lo:
            lo, hi = v.min(), max(v.max(), v.min() + 1e-9)
        return np.linspace(lo, hi, bins + 1)

    xe = edges(xs)
    ye = edges(ys)
    ix = np.clip(np.digitize(xs, xe[1:-1]), 0, len(xe) - 2)
    iy = np.clip(np.digitize(ys, ye[1:-1]), 0, len(ye) - 2)
    nx, ny = len(xe) - 1, len(ye) - 1
    cell = ix * ny + iy
    total = np.bincount(cell, minlength=nx * ny).astype(float)

    cm = s_.table("cohort_membership")
    sizes = cm.groupby("group").size()
    usable = sorted(sizes[sizes >= S.MIN_GROUP_N].index)

    strat = {}
    if s_.has("external_groups"):
        eg = s_.table("external_groups")
        for g, sub in eg.groupby("group"):
            r = sub.sort_values("pct_retained").iloc[0]
            strat[g] = {"stratifier": r["stratifier"],
                        "pct_retained": float(r["pct_retained"]),
                        "signal_over_random": float(r["signal_over_random"])}

    sym = genes["symbol_key"].to_numpy()
    rng = np.random.default_rng(0)
    n_genes = len(genes)
    panels = []

    # All dimensions above the noise ceiling, per-dimension z-scored so each
    # weighs equally in the centroid distance. Used for the displacement test.
    #
    # The space follows the axes being displayed. This was previously ALWAYS the
    # raw PCs, so a grid drawn on sPC1 x sPC2 reported `strongest_dim: pc3` and
    # a displacement measured in a space the viewer was not looking at. Worse,
    # the raw space carries magnitude (PC1 r 0.505, PC3 r 0.519), so part of any
    # displacement measured there is "this set has more signal" - which is the
    # confound every Aim 3 result corrects for. CURRENT_FINDINGS section 4 is
    # computed in the corrected space, so the app now matches it.
    use_shape = x.startswith("spc") or y.startswith("spc")
    scree_tbl = "shape_scree" if (use_shape and s_.has("shape_scree")) else "pc_scree"
    pre = "spc" if scree_tbl == "shape_scree" else "pc"
    D, dim_cols = None, []
    if s_.has(scree_tbl):
        sc = s_.table(scree_tbl)
        dim_cols = [f"{pre}{int(r.pc)}" for r in sc.itertuples()
                    if r.above_noise and f"{pre}{int(r.pc)}" in emb.columns]
        if dim_cols:
            full = s_.genes[["gene_id"]].merge(
                emb[["gene_id"] + dim_cols], on="gene_id", how="left"
            ).set_index("gene_id").loc[genes["gene_id"]].to_numpy(float)
            D = (full - full.mean(0)) / full.std(0)

    for g in usable:
        members = set(cm.loc[cm["group"] == g, "symbol_key"])
        hit = np.fromiter((s in members for s in sym), bool, n_genes)
        k = hit.sum()
        if k < S.MIN_GROUP_N:
            continue
        base = k / n_genes
        obs = np.bincount(cell[hit], minlength=nx * ny).astype(float)

        with np.errstate(divide="ignore", invalid="ignore"):
            rate = (obs + 0.5) / (total + 1.0)
            log2 = np.log2(rate / base)
        log2[total < min_n] = np.nan

        # Label-permutation null, corrected for scanning the whole grid.
        threshold = None
        if n_perm:
            maxima = np.empty(n_perm)
            for b in range(n_perm):
                pick = rng.choice(n_genes, size=k, replace=False)
                o = np.bincount(cell[pick], minlength=nx * ny).astype(float)
                r = (o + 0.5) / (total + 1.0)
                l = np.abs(np.log2(r / base))
                l[total < min_n] = np.nan
                maxima[b] = np.nanmax(l) if np.isfinite(l).any() else 0.0
            threshold = float(np.quantile(maxima, 0.95))

        # Cell-wise testing only sees PATCHES. A set that varies smoothly across
        # the map spreads its signal over many cells, each below threshold, and
        # would be reported as "nothing" when it is in fact strongly positioned.
        # So also measure the simplest gradient: how far the members sit along
        # each axis relative to everyone else, as a standardised difference.
        shift = {}
        for key, vals in ((x, xs), (y, ys)):
            a, b = vals[hit], vals[~hit]
            sd = vals.std()
            shift[key] = round(float((a.mean() - b.mean()) / sd), 3) if sd else 0.0

        # Displacement in the FULL space, not just the two axes on screen. A set
        # can be flat on this plane and strongly positioned on PC7, and the
        # patch and gradient numbers would both miss it. This is the honest
        # single answer to "is this set positioned at all".
        #
        # Reported next to an overlap percentage on purpose. With 1,846 genes
        # every set clears significance, so the p-value is not the interesting
        # part; how much the two distributions actually overlap is.
        disp = None
        if D is not None:
            obs = np.linalg.norm(D[hit].mean(0))
            null = np.empty(min(n_perm, 200) or 1)
            for b in range(len(null)):
                null[b] = np.linalg.norm(D[rng.choice(n_genes, size=k, replace=False)].mean(0))
            sd = null.std()
            best_j = int(np.abs(D[hit].mean(0)).argmax())
            a_, b_ = D[hit, best_j], D[~hit, best_j]
            sp = np.sqrt(((k - 1) * a_.var(ddof=1)
                          + (n_genes - k - 1) * b_.var(ddof=1)) / (n_genes - 2))
            cd = float((a_.mean() - b_.mean()) / sp) if sp else 0.0
            disp = {
                "z": round(float((obs - null.mean()) / sd), 1) if sd else None,
                "p": round(float((np.sum(null >= obs) + 1) / (len(null) + 1)), 4),
                "strongest_dim": dim_cols[best_j],
                "cohens_d": round(cd, 2),
                "overlap_pct": round(float(_overlap(cd) * 100), 1),
            }

        # Distributions along the displayed axes, members against the rest.
        # Added 2026-08-16 because "displaced but not separated" is the central
        # claim and two heavily overlapping curves show it in a way no effect
        # size does: a reader who sees d = -0.31 has to be told that means 88%
        # overlap, whereas a reader who sees the curves has already understood.
        # Shared bin edges over the pooled range so the two are comparable, and
        # densities rather than counts because the groups differ ~10-fold in n.
        dists = {}
        for ax in (x, y):
            if ax not in emb.columns:
                continue
            v = emb.set_index("gene_id").loc[genes["gene_id"], ax].to_numpy(float)
            ok_v = np.isfinite(v)
            if ok_v.sum() < 20:
                continue
            edges = np.histogram_bin_edges(v[ok_v], bins=28)
            hin, _ = np.histogram(v[hit & ok_v], bins=edges, density=True)
            hout, _ = np.histogram(v[~hit & ok_v], bins=edges, density=True)
            a_v, b_v = v[hit & ok_v], v[~hit & ok_v]
            sp_v = np.sqrt((a_v.var() + b_v.var()) / 2)
            d_v = float((a_v.mean() - b_v.mean()) / sp_v) if sp_v else 0.0
            dists[ax] = {
                "edges": [round(float(e), 4) for e in edges],
                "members": [round(float(h), 5) for h in hin],
                "rest": [round(float(h), 5) for h in hout],
                "member_mean": round(float(a_v.mean()), 3),
                "rest_mean": round(float(b_v.mean()), 3),
                "cohens_d": round(d_v, 2),
                "overlap_pct": round(float(_overlap(d_v) * 100), 1),
            }

        finite = log2[np.isfinite(log2)]
        panels.append({
            "group": g,
            "n": int(k),
            "axis_shift": shift,
            # Per-axis member vs rest densities; see the comment above.
            "distributions": dists,
            # Displacement in the full above-noise space, with the overlap that
            # displacement actually corresponds to. Both, always: every set
            # clears significance here, so the effect size is the honest number.
            "displacement": disp,
            # Indices into `points`, so each tile can draw the real scatter with
            # this set highlighted instead of an abstract heatmap. The picture
            # is the thing people read; the statistics go on top of it.
            "members": [int(i) for i in np.flatnonzero(hit)],
            "cells": [None if not np.isfinite(v) else round(float(v), 3)
                      for v in log2],
            "max_abs": float(np.abs(finite).max()) if finite.size else 0.0,
            # How many cells clear the permutation threshold. Zero is a real and
            # useful answer: it says this set is spread across the map.
            "n_above_null": (int((np.abs(finite) > threshold).sum())
                             if threshold is not None else None),
            "null_threshold": threshold,
            "is_super_enhancer": g == "dbSUPER_CD4_SE_TSS_pm50kb",
            "is_positive_control": g == "gene_desert_bottomQ_density",
            "stratification": strat.get(g),
        })

    # Sorted by the larger of the two signals, patch or gradient, so a set that
    # is strongly positioned but smoothly so is not buried at the bottom.
    def strength(p):
        sh = max(abs(v) for v in p["axis_shift"].values()) if p["axis_shift"] else 0
        return -(max((p["n_above_null"] or 0) / 10.0, sh))

    panels.sort(key=strength)

    return _clean({
        # ALL_AXIS_*, not PC_*: the grid defaults to the corrected plane, and
        # PC_LABELS has no entry for "spc1" so it fell through to the bare key
        # with poles=None, leaving the tiles unlabelled and unoriented.
        "x_axis": {"key": x, "label": S.ALL_AXIS_LABELS.get(x, x),
                   "poles": S.ALL_AXIS_POLES.get(x)},
        "y_axis": {"key": y, "label": S.ALL_AXIS_LABELS.get(y, y),
                   "poles": S.ALL_AXIS_POLES.get(y)},
        "nx": nx, "ny": ny, "min_n": min_n, "n_perm": n_perm,
        "n_genes": int(n_genes),
        "cell_totals": [int(v) for v in total],
        # The shared point cloud, in real axis units, drawn faintly on every
        # tile so each small multiple IS the map rather than an abstraction of
        # it. Rounded to 3dp: this is 1,846 pairs and full float precision
        # would triple the payload for no visible difference.
        "points": [[round(float(a), 3), round(float(b), 3)] for a, b in zip(xs, ys)],
        "extent": {"x0": float(xe[0]), "x1": float(xe[-1]),
                   "y0": float(ye[0]), "y1": float(ye[-1])},
        "is_umap": x.startswith("umap") or y.startswith("umap"),
        "is_null": x.startswith("umap_null") or y.startswith("umap_null"),
        "panels": panels,
        "caveats": {
            "post_hoc": "The embedding is built from MCC features alone. These sets "
                        "were painted on afterwards and never entered the "
                        "coordinates, which is what makes the overlap worth reading.",
            "null": "Colour is log2 of the observed rate over the panel base rate. "
                    "The outlined cells are the ones exceeding a label-permutation "
                    "null, membership shuffled across genes at fixed set size, taken "
                    "at the 95th percentile of the per-shuffle MAXIMUM so it is "
                    "corrected for scanning every cell. Sets with no outlined cells "
                    "are spread across the map, which is a result, not a failure.",
            "robustness": "Trust the gradient over the patch count. Patch counts "
                          "depend on where the cell boundaries fall, and moving from "
                          "percentile to equal-width edges changed several of them "
                          "by a cell or two, while every gradient stayed identical "
                          "because it does not use the grid at all.",
            "gradient": "Cell-wise testing only finds patches. A set that varies "
                        "smoothly across the map spreads its signal thinly and can "
                        "clear no single cell while still being strongly positioned, "
                        "so each panel also reports the standardised shift of its "
                        "members along each axis. Read both: patches and gradients "
                        "are different claims and this figure can show either.",
            "atac": "ATAC gates which genes and peaks exist at all, so sets derived "
                    "from accessibility are not independent of the coordinates.",
            "density": "gene_desert_bottomQ_density is the positive control for the "
                       "gene-density confound: it retains only 11% of its effect "
                       "under density stratification while insulation-based effects "
                       "survive. Read every panel with its stratification figure.",
            "umap": S.UMAP_CAVEAT,
        },
    })


@app.get("/api/embedding/loadings")
def plane_loadings(x: str = S.DEFAULT_PLANE[0], y: str = S.DEFAULT_PLANE[1],
                   top: int = 8):
    """Loading vectors for the two displayed axes, for a biplot overlay.

    Ranked by length in the plane, sqrt(lx^2 + ly^2), not by either component
    alone. A feature can load hard on x and not at all on y, and ranking by x
    would draw it as a long arrow that says nothing about the vertical spread
    the viewer is looking at.

    Returns `available: False` rather than an error for UMAP axes. UMAP has no
    loadings at all: it is a nonlinear embedding with no linear map back to
    features, so an arrow over it would be an invention.
    """
    s_ = store()
    # Shape axes must read shape_loadings. Reading pc_loadings for an "spc" axis
    # would draw arrows from a different coordinate system onto the plot, which
    # is worse than drawing none.
    shape = x.startswith("spc") or y.startswith("spc")
    tbl = "shape_loadings" if shape else "pc_loadings"
    if not s_.has(tbl):
        raise HTTPException(503, f"{tbl} not built")

    def pc_index(a: str) -> int | None:
        pre = "spc" if shape else "pc"
        return int(a[len(pre):]) if a.startswith(pre) and a[len(pre):].isdigit() else None

    ix, iy = pc_index(x), pc_index(y)
    if ix is None or iy is None:
        return {
            "available": False,
            "reason": "Loadings exist only for the principal components. UMAP is "
                      "nonlinear and has no linear map back to the features, so "
                      "there is no arrow to draw.",
        }

    df = s_.table("pc_loadings")
    lx = df[df["pc"] == ix].set_index("feature")["loading"]
    ly = df[df["pc"] == iy].set_index("feature")["loading"]
    if lx.empty or ly.empty:
        raise HTTPException(404, f"no loadings for {x} or {y}")

    both = pd.DataFrame({"x": lx, "y": ly}).dropna()
    both["length"] = np.hypot(both["x"], both["y"])
    sel = both.sort_values("length", ascending=False).head(top)

    return _clean({
        "available": True,
        "x_axis": x,
        "y_axis": y,
        "n_features": int(len(both)),
        "vectors": [
            {"feature": f, "x": float(r.x), "y": float(r.y), "length": float(r.length)}
            for f, r in sel.iterrows()
        ],
        "note": "Arrow direction is the direction in which that feature increases. "
                "Length is how strongly the feature loads on this plane. PCA sign "
                "is arbitrary, so read the contrast between opposite arrows, not "
                "the absolute orientation.",
    })


@app.get("/api/dimensions/{pc}/loadings")
def loadings(pc: int, top: int = 15, space: Literal["corrected", "raw"] = "corrected"):
    """Feature loadings for one component, the evidence for its name.

    A named axis is an interpretation of its loadings. Serving the name without
    them would ask the reader to take the interpretation on trust.

    `space` defaults to the amount-corrected axes, matching the display. Pass
    "raw" for the provenance PCs; the two are different components and their
    numbering does NOT correspond, so never compare sPC3 with PC3.
    """
    s_ = store()
    tbl = "shape_loadings" if space == "corrected" else "pc_loadings"
    if not s_.has(tbl):
        tbl = "pc_loadings"
    if not s_.has(tbl):
        raise HTTPException(503, "no loadings table built")
    df = s_.table(tbl)
    hit = df[df["pc"] == pc]
    if hit.empty:
        raise HTTPException(404, f"no loadings for PC{pc}")

    ranked = hit.reindex(hit["loading"].abs().sort_values(ascending=False).index)
    sel = ranked.head(top).sort_values("loading", ascending=False)

    # Scree must come from the same space as the loadings, or the variance share
    # printed beside a component would belong to a different component.
    scree_tbl = "shape_scree" if tbl == "shape_loadings" else "pc_scree"
    scree_row = None
    if s_.has(scree_tbl):
        sc = s_.table(scree_tbl)
        r = sc[sc["pc"] == pc]
        if not r.empty:
            scree_row = r.iloc[0].to_dict()

    # The label must come from the same space as the loadings. Serving shape
    # loadings under a raw PC label is exactly the two-coordinate-system
    # confusion this endpoint was changed to remove.
    shape_space = tbl == "shape_loadings"
    key = f"{'spc' if shape_space else 'pc'}{pc}"
    return _clean({
        "pc": pc,
        "space": "amount-corrected" if shape_space else "raw (provenance only)",
        "label": S.ALL_AXIS_LABELS.get(key, key.upper()),
        "poles": S.ALL_AXIS_POLES.get(key),
        "variance_pct": scree_row.get("variance_pct") if scree_row else None,
        "above_noise": bool(scree_row.get("above_noise")) if scree_row else None,
        "n_features": int(len(hit)),
        "loadings": sel.to_dict("records"),
    })


@app.get("/api/embedding")
def embedding(
    x: str = S.DEFAULT_EMBEDDING_AXES[0],
    y: str = S.DEFAULT_EMBEDDING_AXES[1],
    highlight: str | None = None,
):
    """All genes as 2D coordinates, for the continuum map.

    Any pair of stored axes can be crossed. PC1 is the amount axis, so PC1 x PC2
    largely re-sorts genes by signal depth; PC2 x PC3 is the interpretable plane
    and is the default.
    """
    s = store()
    if not s.has("embeddings"):
        raise HTTPException(503, "embeddings table not built")

    # Both coordinate systems are served from one endpoint, keyed by the axis
    # prefix: "spc*" is the amount-corrected space (the default and the one the
    # analysis uses), "pc*"/"umap*" are the raw space kept for the toggle. They
    # are joined on gene_id so a caller can cross a shape axis with a UMAP axis
    # if they insist, though that is rarely meaningful.
    emb = _embeddings_all(s)

    for ax in (x, y):
        if ax not in emb.columns:
            raise HTTPException(
                400,
                f"unknown axis {ax!r}; available: "
                f"{', '.join(c for c in emb.columns if c not in ('gene_id', 'symbol_key'))}",
            )

    # `region` is the display taxonomy; `group` is the superseded 2026-07-21
    # labelling, kept so provenance joins still work.
    gcols = ["gene_id", "gene_symbol", "group", "max_posterior", "is_core"]
    gcols += [c for c in ("region", "top_weight") if c in s.genes.columns]
    g = s.genes[gcols]
    df = emb[["gene_id", x, y]].merge(g, on="gene_id", how="left")

    hl = s.resolve(highlight) if highlight else None
    _pts = df.rename(columns={x: "x", y: "y"})

    return {
        "x_axis": {"key": x, "label": S.ALL_AXIS_LABELS.get(x, x),
                   "poles": S.ALL_AXIS_POLES.get(x)},
        "y_axis": {"key": y, "label": S.ALL_AXIS_LABELS.get(y, y),
                   "poles": S.ALL_AXIS_POLES.get(y)},
        "space": ("amount-corrected" if x.startswith("spc") and y.startswith("spc")
                  else "raw (provenance only)" if x.startswith("pc") and y.startswith("pc")
                  else "mixed"),
        "space_note": S.SHAPE_SPACE_NOTE,
        "axes_available": [
            {"key": c, "label": S.ALL_AXIS_LABELS.get(c, c),
             "space": "corrected" if c.startswith("spc") else "raw"}
            for c in emb.columns if c not in ("gene_id", "symbol_key")
        ],
        "highlight": hl,
        "is_umap": x.startswith("umap") or y.startswith("umap"),
        "is_null": x.startswith("umap_null") or y.startswith("umap_null"),
        "umap_caveat": S.UMAP_CAVEAT,
        "continuum_caveat": S.CONTINUUM_CAVEAT,
        "n": int(len(df)),
        # Filter against the RENAMED frame, not `df`. Checking `df.columns` here
        # silently dropped "x" and "y" (they are still spc1/spc2 at that point),
        # so every point came back without coordinates and the map rendered
        # empty. A membership test against the wrong frame fails quietly, which
        # is why this looked like a frontend bug.
        "points": _clean(
            _pts[[c for c in ("gene_id", "gene_symbol", "x", "y", "group",
                              "region", "top_weight", "max_posterior", "is_core")
                  if c in _pts.columns]].to_dict("records")
        ),
    }


# ---------------------------------------------------------------------------
# lab, exploratory views, not part of the main app
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def _repro_per_feature() -> pd.DataFrame:
    """Per-feature cross-capture rho, computed from `reproducibility_pairs`.

    REPOINTED 2026-08-17. This used to read the `reproducibility` table, which is
    a straight copy of
    `scripts_cleaned/audit/continuous_methods/cross_panel_reproducibility.tsv`,
    dated 2026-07-31 and never regenerated. By today that file is 17 days stale
    and **24 of its 63 features no longer exist in the substrate**, including
    every `oe_asymmetry_*` and `oe_tailedness_*`, which were deleted on
    2026-08-16 precisely BECAUSE they failed a reproducibility threshold.

    So the published median of 0.752 was computed over a set containing 14
    features removed for being unreproducible, and it understates the current
    substrate. Recomputed here on the live pairs:

        stale table, all 63 features        median rho 0.752
        the 24 retired features              median rho 0.654
        stale table, 39 surviving features   median rho 0.833
        CURRENT pairs, 41 features           median rho 0.819
                                             33/41 above 0.7, ZERO below 0.3

    The "one feature below 0.3" quoted in the app CLAUDE.md was in the retired
    set. No unreproducible feature remains.

    `reproducibility_pairs` is rebuilt from the immune pickle on every store
    build, so it cannot go stale the way a copied tsv can, and all 41 of its
    features are live. That is why it is now the single source.

    Spearman by hand rather than via scipy: the thin server env has no scipy, and
    rho on ranks is a two-line computation.
    """
    s_ = store()
    pairs = s_.table("reproducibility_pairs")
    rows = []
    for feat, g in pairs.groupby("feature"):
        a = g["gw"].rank().to_numpy(float)
        b = g["immune"].rank().to_numpy(float)
        rho = float(np.corrcoef(a, b)[0, 1]) if len(g) > 2 else float("nan")
        pr = float(np.corrcoef(g["gw"].to_numpy(float),
                               g["immune"].to_numpy(float))[0, 1]) if len(g) > 2 else float("nan")
        rows.append({"feature": feat, "spearman_rho": rho,
                     "pearson_r": pr, "n": int(len(g))})
    return (pd.DataFrame(rows)
            .sort_values("spearman_rho", ascending=False)
            .reset_index(drop=True))


@app.get("/api/lab/reproducibility")
def lab_reproducibility(feature: str | None = None):
    """The 116 genes captured in both panels, feature by feature.

    The only direct measurement of technical reproducibility available: same
    gene, two independent captures, same pipeline. Makes the 0.511 noise floor
    and kappa = 0.72 concrete rather than quoted.

    LAB ONLY. The immune side is the stale 791 baseline (missing PDCD1 + 15),
    so this is a technical demonstration, never a biological claim.
    """
    s_ = store()
    if not s_.has("reproducibility_pairs"):
        raise HTTPException(503, "reproducibility pairs not built")

    per_feature = _repro_per_feature()
    out = {
        "caveat": "Immune side is the stale 791 baseline, missing PDCD1 and 15 other "
                  "genes from the chromosome-edge bug. Valid as a technical "
                  "reproducibility demonstration; not a biological claim. The overlap "
                  "may shift after the re-run.",
        "n_genes": int(per_feature["n"].max()) if "n" in per_feature else None,
        "n_features": int(len(per_feature)),
        "median_rho": float(per_feature["spearman_rho"].median()),
        "n_above_0_7": int((per_feature["spearman_rho"] >= 0.7).sum()),
        "n_below_0_3": int((per_feature["spearman_rho"] < 0.3).sum()),
        "note": "Computed live from the paired values, not read from the "
                "2026-07-31 tsv, which is 17 days stale and holds 24 features "
                "that no longer exist. The asymmetry and tailedness features "
                "that used to sit at the bottom of this list were removed on "
                "2026-08-16 for failing reproducibility, which is why the median "
                "is now 0.819 rather than the published 0.752 and nothing is "
                "left below rho 0.3.",
        "source": "reproducibility_pairs (rebuilt every store build)",
        "per_feature": _clean(per_feature.to_dict("records")),
    }

    if feature and s_.has("reproducibility_pairs"):
        pairs = s_.table("reproducibility_pairs")
        hit = pairs[pairs["feature"] == feature]
        if hit.empty:
            raise HTTPException(404, f"no paired values for {feature!r}")
        out["pairs"] = {
            "feature": feature,
            "points": _clean(hit[["symbol_key", "gw", "immune"]].to_dict("records")),
        }
    return _clean(out)


# NOTE: panel.csv is declared BEFORE {gene}.csv. FastAPI matches routes in
# declaration order, so the parameterised route would otherwise capture
# "panel" as a gene name and 404.
@app.get("/api/export/panel.csv")
def export_panel():
    """The whole panel: gene, label, posterior, and every stored coordinate."""
    s_ = store()
    cols = [c for c in ("gene_id", "gene_symbol", "viewpoint_chrom", "viewpoint_pos",
                        "group", "max_posterior", "confidence_class", "is_core")
            if c in s_.genes.columns]
    df = s_.genes[cols].copy()
    if s_.has("embeddings"):
        emb = _embeddings_all(s_).drop(columns=["symbol_key"], errors="ignore")
        df = df.merge(emb, on="gene_id", how="left")
    return Response(
        content=df.to_csv(index=False),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="mccprofiler_panel.csv"'},
    )


@app.get("/api/export/{gene}.csv")
def export_gene(gene: str):
    """One gene, everything known about it, as CSV.

    The CellProfiler analogy only holds if the output is portable, its
    deliverable IS a feature table. This is the line between a demo and a tool.
    """
    s_ = store()
    gid = s_.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")

    rows: list[dict] = []
    row = s_.genes[s_.genes["gene_id"] == gid].iloc[0]
    for k in ("gene_id", "gene_symbol", "viewpoint_chrom", "viewpoint_pos",
              "group", "max_posterior", "confidence_class"):
        if k in row:
            rows.append({"section": "gene", "name": k, "value": row[k]})

    if s_.has("features"):
        f = s_.table("features")
        for _, r in f[f["symbol_key"] == gid.upper()].iterrows():
            rows.append({"section": "feature", "name": r["feature"],
                         "value": r["z"], "percentile": r["percentile"]})

    if s_.has("embeddings"):
        e = s_.table("embeddings")
        hit = e[e["symbol_key"] == gid.upper()]
        if not hit.empty:
            for c in hit.columns:
                if c in ("gene_id", "symbol_key"):
                    continue
                rows.append({"section": "coordinate", "name": c,
                             "value": float(hit.iloc[0][c])})

    df = pd.DataFrame(rows)
    csv = df.to_csv(index=False)
    return Response(
        content=csv,
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{gid}_mccprofiler.csv"'},
    )


# ---------------------------------------------------------------------------
# gene finder, the LLM translates, pandas retrieves
# ---------------------------------------------------------------------------


def _vocabulary(s_) -> dict:
    """What a query may refer to. Also the model's entire view of the data."""
    # Corrected axes FIRST and raw PCs after, because this list is also the
    # model's entire view of the coordinate system. Advertising raw PCs first
    # would invite queries against axes that carry magnitude (PC1 r 0.505,
    # PC3 r 0.519), i.e. filters that partly select on how much signal a gene
    # has. The ask endpoint currently passes no axes at all, by decision, but
    # the deterministic /api/query path and the UI both read this.
    emb = s_.table("embeddings") if s_.has("embeddings") else None
    raw = [c for c in (emb.columns if emb is not None else [])
           if c not in ("gene_id", "symbol_key")]
    sh = s_.table("shape_embeddings") if s_.has("shape_embeddings") else None
    shape = [c for c in (sh.columns if sh is not None else [])
             if c not in ("gene_id", "symbol_key")]
    axes = shape + raw
    cohorts: list[str] = []
    if s_.has("cohort_membership"):
        counts = s_.table("cohort_membership")["group"].value_counts()
        cohorts = sorted(counts[counts >= S.MIN_GROUP_N].index.tolist())
    groups = sorted(x for x in s_.genes["group"].dropna().unique())
    features = sorted(s_.table("features")["feature"].unique()) if s_.has("features") else []
    return {"axes": axes, "cohorts": cohorts, "groups": groups, "features": features}


@app.get("/api/vocabulary")
def vocabulary():
    """The query vocabulary, so the UI can build a query without the model."""
    s_ = store()
    v = _vocabulary(s_)
    return {
        **v,
        "axis_labels": {a: S.ALL_AXIS_LABELS.get(a, a) for a in v["axes"]},
        "axis_poles": {a: S.ALL_AXIS_POLES[a] for a in v["axes"]
                       if a in S.ALL_AXIS_POLES},
        "axis_space": {a: ("corrected" if a.startswith("spc") else "raw")
                       for a in v["axes"]},
        "space_note": S.SHAPE_SPACE_NOTE,
        "archetype_labels": {g: S.ARCHETYPE_DISPLAY.get(g, {}).get("display", g)
                             for g in v["groups"]},
    }


@app.post("/api/query")
def execute_query(query: dict):
    """Run a query DSL object. Deterministic; needs no API key.

    This is the retrieval path. /api/ask only adds a translation step in front
    of it, so anything the model can ask for, the UI can ask for directly.
    """
    try:
        return _clean(run_query(query, store()))
    except QueryError as e:
        raise HTTPException(400, str(e)) from None


@app.get("/api/translate/status")
def translate_status():
    """Whether the plain-English box will work, and which model backs it."""
    return _translate.status()


@app.post("/api/ask")
def ask(body: dict):
    """Translate a question into a query, then run it.

    The model never sees the data, only the question and a schema naming the
    available axes and features. The parsed query comes back with the results so
    the translation can be checked and corrected: a misread question shows up as
    a wrong query, not a wrong gene list.

    Selection is on contact architecture ONLY. External reference sets are not
    in the schema the model receives, so it cannot select on them even if the
    question invites it. That is the super-enhancer rule generalised: borrowed
    labels are post-hoc validation, and a tool whose claim is that architecture
    is more informative than those labels cannot use them to choose its genes.
    Membership is attached to the results afterwards, as annotation.
    """
    question = (body or {}).get("question", "").strip()
    if not question:
        raise HTTPException(400, "no question given")

    s_ = store()
    v = _vocabulary(s_)

    # No axis poles, because there are no axis filters. See below.
    pole_lines = ""
    # Empty axis AND cohort lists, so neither filter type exists in the schema
    # rather than being merely discouraged. A model cannot emit what it has no
    # vocabulary for, which is a stronger guarantee than an instruction.
    #
    # Axes are gone because a principal component is a mixture: the concept each
    # one is named after carries only 22 to 35 percent of it
    # (scripts/diagnose_pc_names.py), so "PC2 low" silently selects on several
    # things at once where `frac_far_distal` high selects on one. The 91
    # features say what they mean, including 34 that are element-class specific,
    # which is what lets "long-range enhancer contacts" be asked as the single
    # joint condition it actually is.
    schema = _translate.flat_schema([], [], v["groups"], v["features"])

    try:
        query = _translate.translate(
            question, schema, pole_lines, _feature_reference(tuple(v["features"]))
        )
    except _translate.TranslationUnavailable as e:
        raise HTTPException(503, str(e)) from None
    except _translate.TranslationFailed as e:
        raise HTTPException(502, f"translation failed: {e}") from None

    # A wholly unanswerable question translates to nothing executable, and no
    # filters means every gene passes, so the honest "cannot answer" would
    # render as all 1,846 genes matched. Refuse instead of returning the panel.
    #
    # Must check `rank` as well as `filters`. A ranking with a declared proxy
    # (no filters, `unsupported` set to say a TAD was approximated by contact
    # reach) is a perfectly good answer, and testing filters alone threw it away.
    if (not query.get("filters") and not query.get("rank")
            and (query.get("unsupported") or "").strip()):
        return _clean({
            "question": question,
            "query": query,
            "interpretation": query.get("interpretation"),
            "unsupported": query["unsupported"].strip(),
            "provider": _translate.active_provider(),
            "n_matched": 0,
            "steps": [],
            "genes": [],
            "note": "Nothing in this question maps to something the store holds, "
                    "so no filter was built. See what is missing above.",
        })

    try:
        result = run_query(query, s_)
    except QueryError as e:
        # A query that cannot execute is still worth showing: it is the evidence
        # of what the model misread.
        return _clean({"question": question, "query": query, "error": str(e)})

    # Post-hoc only. These sets played no part in choosing the genes above, which
    # is what makes them worth reading: the overlap is an observation about an
    # architecture-selected list, not a property it was selected for.
    result = dict(result)
    result["genes"] = _annotate_cohorts(s_, result.get("genes", []))

    return _clean({
        "question": question,
        "query": query,
        "interpretation": query.get("interpretation"),
        "reasoning": query.get("reasoning") or None,
        "annotation_note": "External set membership is shown for the returned genes "
                           "but played no part in selecting them. Selection used "
                           "contact architecture only.",
        # Set when part of the question asks for something the store does not
        # hold, TADs and insulation being the usual cases. Shown rather than
        # swallowed, because the failure mode worth avoiding is a confident
        # answer to a question the data cannot address.
        "unsupported": (query.get("unsupported") or "").strip() or None,
        "provider": _translate.active_provider(),
        **result,
        "disclaimer": "The model translated your question into the query shown. It "
                      "never saw the data. Check the query, and if it misread you, "
                      "edit it and re-run.",
    })


@app.get("/api/feature/{name}")
def feature_explainer(name: str, level: int = 2):
    """What a feature measures, drawn on two real profiles.

    Returns the geometry to overlay plus the panel's highest- and lowest-scoring
    gene on this feature, with their profiles, so the feature is shown as a
    contrast between two real traces rather than defined in prose.
    """
    s_ = store()
    if not s_.has("features"):
        raise HTTPException(503, "features table not built")

    f = s_.table("features")
    sub = f[f["feature"] == name]
    if sub.empty:
        raise HTTPException(404, f"unknown feature {name!r}")

    geom = describe_feature(name)

    # Exemplars: the extremes make the geometry legible in a way the median
    # never does.
    hi = sub.loc[sub["z"].idxmax()]
    lo = sub.loc[sub["z"].idxmin()]

    examples = []
    for role, row in (("high", hi), ("low", lo)):
        gid = s_.resolve(str(row["symbol_key"]))
        if gid is None:
            continue
        try:
            prof = s_.profile(gid, level=level, mode="raw")
        except (KeyError, ValueError):
            continue
        g = s_.genes[s_.genes["gene_id"] == gid].iloc[0]
        examples.append({
            "role": role,
            "gene_id": gid,
            "gene_symbol": g["gene_symbol"],
            "group": g.get("group"),
            "z": float(row["z"]),
            "percentile": float(row["percentile"]),
            "profile": prof,
        })

    return _clean({
        **geom,
        "distribution": {
            "min": float(sub["z"].min()),
            "max": float(sub["z"].max()),
            "median": float(sub["z"].median()),
        },
        "examples": examples,
        "caveat": "These are the panel extremes, chosen to make the geometry "
                  "visible. Most genes sit between them, the feature is a "
                  "continuous quantity, not a two-way split.",
    })


@app.get("/api/ranked")
def ranked(
    axis: str = S.DEFAULT_PLANE[0],
    direction: Literal["top", "bottom", "both"] = "both",
    limit: int = Query(25, le=200),
    include_qc: bool = False,
):
    """Given an axis, which genes sit at its extremes?

    The inverse of the cohort view, and how you find a gene worth looking at
    when you do not already have one in mind. Percentile is reported alongside
    the raw score because a score is meaningless without the distribution.

    `include_qc` is FALSE by default, and this matters. The 21 `arch-off` genes
    are a QC class, near-empty signal with every descriptor at floor, and they
    are not biology. Measured 2026-08-16: on sPC1 they span -27.3 to -3.3 while
    every other gene in the panel sits between -8.7 and +9.1, so they fill the
    entire negative tail and the view returns twelve QC failures instead of the
    genes a reader is looking for. Set include_qc=true to see them; the response
    always reports how many were withheld so the filtering is never silent.

    This is a DISPLAY filter only. It does not address the deeper problem that
    those same 21 genes carry 25.2% of sPC1's variance and measurably rotate the
    components (|cos| 0.67 to 0.90 when they are dropped before the PCA), which
    is a question about the coordinate system itself, not about this endpoint.
    """
    s_ = store()
    if not s_.has("embeddings"):
        raise HTTPException(503, "embeddings table not built")
    emb = _embeddings_all(s_)
    if axis not in emb.columns:
        raise HTTPException(
            400,
            f"unknown axis {axis!r}; available: "
            f"{', '.join(c for c in emb.columns if c not in ('gene_id', 'symbol_key'))}",
        )

    df = emb[["gene_id", axis]].merge(
        s_.genes[["gene_id", "gene_symbol", "group", "max_posterior"]],
        on="gene_id", how="left",
    )
    n_qc = int((df["group"] == "arch-off").sum())
    if not include_qc:
        df = df[df["group"] != "arch-off"]
    # Percentile AFTER the filter, so it describes the distribution shown.
    df["percentile"] = df[axis].rank(pct=True) * 100
    df = df.rename(columns={axis: "value"})

    def take(asc: bool) -> list[dict]:
        return _clean(df.sort_values("value", ascending=asc).head(limit).to_dict("records"))

    poles = S.ALL_AXIS_POLES.get(axis)
    out = {
        "axis": axis,
        "label": S.ALL_AXIS_LABELS.get(axis, axis),
        "poles": poles,
        "n": int(len(df)),
        "qc_excluded": 0 if include_qc else n_qc,
        "qc_note": (None if include_qc else
                    f"{n_qc} arch-off (QC) genes withheld. They have near-empty "
                    f"signal with every descriptor at floor, and they occupy the "
                    f"whole negative tail of the leading axes, so including them "
                    f"returns QC failures rather than genes. Pass include_qc=true "
                    f"to show them."),
        "caveat": "Extremes of a continuum, not a category. A gene at the 99th "
                  "percentile differs from one at the 95th by degree, and roughly "
                  "half a typical between-gene distance is technical noise.",
    }
    if direction in ("top", "both"):
        out["top"] = take(asc=False)
        out["top_pole"] = poles["pos"] if poles else None
    if direction in ("bottom", "both"):
        out["bottom"] = take(asc=True)
        out["bottom_pole"] = poles["neg"] if poles else None
    return out


@app.get("/api/explain")
def explain():
    """The argument the app is making, with its numbers.

    Leads with the nested-baseline result because that is the justification for
    the 91-feature substrate existing at all, without it, `n_peaks` would do.
    """
    s_ = store()
    out: dict = {
        "why_these_features": {
            "claim": "The 91-feature substrate beats counting peaks on every "
                     "target tested.",
            "detail": "Nested, cross-validated, fixed cell state, identical folds, "
                      "regularisation tuned inside each training fold. The 9/9 "
                      "consistency matters more than the size: there is no target "
                      "where the extra 80 features are dead weight. The nonlinear "
                      "arm was uniformly WORSE, which kills the objection that the "
                      "features carry more and a linear model simply could not "
                      "reach it.",
        },
        "why_not_clusters": {
            "claim": "The landscape is a continuum, and that is a positive result "
                     "rather than a failure to cluster.",
            "detail": "SigClust rejects a single Gaussian (z = -11.86 trimmed) with "
                      "the permuted control clean at p = 1.000, while the dip test "
                      "finds no multimodality. Two tests with opposite assumptions. "
                      "The gap statistic returns k=1 with the gap declining "
                      "monotonically, HDBSCAN returns one cluster, and 44% of active "
                      "genes are mixtures.",
        },
        "why_name_regions_at_all": {
            "claim": "The groups are not discoverable, but they are reproducible "
                     "once imposed.",
            "detail": "Using the 116 genes captured in BOTH panels, cluster once on "
                      "the pooled matrix, then assign each gene's two independent "
                      "captures separately, Cohen's kappa peaks at 0.72 for k=3-4. "
                      "That separates two claims usually conflated: the groups "
                      "cannot be found from the data's density, but once defined "
                      "they are reproducible measurements. Only the first failed.",
            "citations": ["Rousseeuw 1987", "Hennig 2015, What are the true clusters?",
                          "Altman & Royston 2006", "Pott & Lieb 2015"],
        },
        "what_it_is_not": {
            "claim": "Effect sizes are small in absolute terms.",
            "detail": "The best R-squared against any biology anchor anywhere is "
                      "0.077, and partial distance correlation shows roughly 80% of "
                      "the association is mediated by genomic position. Continuous "
                      "coordinates capture marginally more biology than discrete "
                      "labels; neither captures much beyond where a gene sits.",
        },
        "noise_floor": {
            "claim": "About half a typical between-gene distance is technical.",
            "detail": "Gene-to-itself distance 5.36 vs gene-to-other 10.49 across "
                      "repeat captures, ratio 0.511. This bounds how well anything "
                      "downstream can perform.",
        },
        "atac_is_not_a_feature": {
            "claim": "ATAC gates membership; it is never a clustering input.",
            "detail": "Clustering reads only the MCC channel. ATAC still determines "
                      "results indirectly by gating which genes and which peaks "
                      "enter, so ATAC-based validation is NOT independent.",
        },
    }

    if s_.has("nested_baselines"):
        nb = s_.table("nested_baselines")
        out["why_these_features"]["table"] = _clean(nb.to_dict("records"))
    if s_.has("dimensions"):
        out["dimensions"] = _clean(s_.table("dimensions").to_dict("records"))
    if s_.has("external_groups"):
        out["external_validation"] = _clean(
            s_.table("external_groups").to_dict("records"))

    out["provenance"] = {
        "panel": s_.manifest["panel"],
        "n_genes": len(s_.genes),
        "built": s_.manifest["built"],
        "scripts_cleaned_commit": s_.manifest.get("scripts_cleaned_commit"),
    }
    return _clean(out)


@app.get("/api/archetypes")
def archetypes():
    """Named regions of the continuum, described by architecture.

    Display names are derived from the top discriminating features and make no
    functional claim. The pipeline's own labels borrow biological categories
    they do not track, see ARCHETYPE_DISPLAY for the measurements.
    """
    s = store()
    counts = s.genes["group"].value_counts().to_dict()
    rows = []
    for canonical, meta in S.ARCHETYPE_DISPLAY.items():
        n = int(counts.get(canonical, 0))
        sub = s.genes[s.genes["group"] == canonical]
        rows.append({
            "canonical": canonical,
            "display": meta["display"],
            "architecture": meta["architecture"],
            "top_features": meta["top_features"],
            "is_qc": meta.get("is_qc", False),
            "n": n,
            "n_core": int(sub["is_core"].sum()) if "is_core" in sub else None,
            "n_mixture": int((~sub["is_core"]).sum()) if "is_core" in sub else None,
        })
    return {
        "caveat": S.CONTINUUM_CAVEAT,
        "naming_note": "Names describe contact architecture, not gene function. "
                       "The pipeline's canonical ids are retained for provenance.",
        "archetypes": _clean(rows),
    }


@app.get("/api/cohorts/compare")
def cohort_compare(
    group: str | None = None,
    symbols: str | None = None,
    axes: str = "pc1,pc2,pc3,pc4,pc5",
):
    """Where does a gene set sit on each axis, against the rest of the panel?

    Either a built-in `group`, or a pasted `symbols` list. For a pasted list the
    response leads with coverage, because the panel is 1,846 of ~20,000 genes,
    a collaborator's 40 hits may match 4, and plotting 4 points as though they
    were 40 would be the wrong answer delivered confidently.
    """
    s = store()
    if not s.has("embeddings"):
        raise HTTPException(503, "embeddings table not built")

    emb = s.table("embeddings")
    axis_keys = [a.strip() for a in axes.split(",") if a.strip()]
    for a in axis_keys:
        if a not in emb.columns:
            raise HTTPException(400, f"unknown axis {a!r}")

    coverage = None
    if symbols:
        wanted = [w.strip().upper() for w in symbols.replace("\n", ",").split(",") if w.strip()]
        member = set(emb.loc[emb["symbol_key"].isin(wanted), "symbol_key"])
        coverage = {
            "requested": len(set(wanted)),
            "in_panel": len(member),
            "matched": sorted(member)[:200],
            "missing": sorted(set(wanted) - member)[:50],
            "note": "The panel is 1,846 of ~20,000 genes. Genes not captured cannot "
                    "be placed; they are not absent from the biology.",
        }
        label = "pasted list"
    elif group:
        if not s.has("cohort_membership"):
            raise HTTPException(503, "cohort_membership not built")
        cm = s.table("cohort_membership")
        member = set(cm.loc[cm["group"] == group, "symbol_key"])
        if not member:
            raise HTTPException(404, f"no panel genes in group {group!r}")
        label = group
    else:
        raise HTTPException(400, "pass either group= or symbols=")

    if len(member) < S.MIN_GROUP_N:
        note = (f"Only {len(member)} panel genes, below the {S.MIN_GROUP_N}-gene "
                f"floor. Treat any apparent difference as noise.")
    else:
        note = None

    in_set = emb["symbol_key"].isin(member)
    rows = []
    for a in axis_keys:
        v_in = emb.loc[in_set, a].dropna()
        v_out = emb.loc[~in_set, a].dropna()
        if len(v_in) < 2 or len(v_out) < 2:
            continue
        # Pooled-SD standardised difference. Reported as an effect size, not a
        # p-value: with n in the hundreds almost anything reaches significance.
        sd = np.sqrt(((len(v_in) - 1) * v_in.var() + (len(v_out) - 1) * v_out.var())
                     / max(len(v_in) + len(v_out) - 2, 1))
        d = (v_in.mean() - v_out.mean()) / sd if sd > 0 else 0.0
        rows.append({
            "axis": a,
            "label": S.ALL_AXIS_LABELS.get(a, a),
            "n_in": int(len(v_in)),
            "cohen_d": float(d),
            "in_quartiles": [float(v_in.quantile(q)) for q in (0.25, 0.5, 0.75)],
            "out_quartiles": [float(v_out.quantile(q)) for q in (0.25, 0.5, 0.75)],
        })

    rows.sort(key=lambda r: -abs(r["cohen_d"]))
    return _clean({
        "label": label,
        "n_in_panel": int(len(member)),
        "coverage": coverage,
        "small_set_warning": note,
        "axes": rows,
    })


@app.get("/api/cohorts")
def cohorts():
    """Externally-defined gene groups, the strongest non-circular evidence.

    These groups were not defined from the features, so unlike the archetypes
    they need no within/pooled-ratio argument to be interpretable.
    """
    s = store()
    if not s.has("cohort_membership"):
        raise HTTPException(503, "cohort_membership not built")

    counts = s.table("cohort_membership")["group"].value_counts()

    # Stratification statistics exist for only some sets, the seven the audit
    # ran. Attach them where present rather than restricting the list to them.
    strat: dict[str, dict] = {}
    if s.has("external_groups"):
        eg = s.table("external_groups")
        for name, sub in eg.groupby("group"):
            strat[str(name)] = {
                str(r["stratifier"]): {
                    "pct_retained": float(r["pct_retained"]),
                    "signal_over_random": float(r["signal_over_random"]),
                }
                for _, r in sub.iterrows()
            }

    rows = []
    for name, n in counts.items():
        rows.append({
            "group": str(name),
            "n_in_panel": int(n),
            "usable": bool(n >= S.MIN_GROUP_N),
            "is_positive_control": name == "gene_desert_bottomQ_density",
            "is_super_enhancer": "SE" in str(name) or "dbSUPER" in str(name),
            "stratification": strat.get(str(name)),
        })
    rows.sort(key=lambda r: -r["n_in_panel"])

    return _clean({
        "min_group_n": S.MIN_GROUP_N,
        "n_offered": sum(r["usable"] for r in rows),
        "n_filtered_out": sum(not r["usable"] for r in rows),
        "positive_control": "gene_desert_bottomQ_density",
        "notes": {
            "why_external": "These groups were not defined from the MCC features, so "
                            "unlike the archetypes they need no within/pooled-ratio "
                            "argument to be interpretable.",
            "positive_control": "gene_desert_bottomQ_density is defined FROM gene "
                                "density and is the one group that collapses under "
                                "density stratification (11% retained) while surviving "
                                "insulation (100%). That the method detects a genuinely "
                                "positional group as positional is what makes the other "
                                "rows meaningful.",
            "display": "Report pct_retained and signal_over_random, not eta-squared. "
                       "Pooled eta-squared runs 0.005-0.023; the ratio is the claim, "
                       "not the absolute value.",
            "super_enhancer": "Super-enhancer sets are shown for post-hoc comparison "
                              "ONLY. They are never inputs to the features or the "
                              "clustering, that would be circular.",
        },
        "rows": rows,
    })

@app.get("/api/taxonomy")
def taxonomy():
    """Reach x composition regions, their radar profiles, and the mixture stats.

    The regions are a RESOLUTION CHOICE on a continuum, not discovered clusters:
    Leiden returns a single community below resolution 0.4, HDBSCAN returns zero
    across sixteen conditions, and on the 116 twice-captured genes Leiden finds
    no communities at all. What justifies naming them is that both levels
    reproduce across independent captures (reach ARI 0.741, composition 0.616)
    and are seed-stable (0.977).

    The structure is deliberately ASYMMETRIC. A CTCF-dominated group exists only
    among far-reaching genes; forcing a third mid-range group yields four genes
    and collapses stability to 0.796. That matches CTCF biology, since CTCF loops
    are long-range structural contacts.
    """
    s_ = store()
    if not s_.has("taxonomy"):
        raise HTTPException(503, "taxonomy table not built")
    tax = s_.table("taxonomy")
    wcols = [c for c in tax.columns if c.startswith("w_")]
    regions = [c[2:] for c in wcols]

    counts = tax["region"].value_counts()
    rows = []
    for r in regions:
        sub = tax[tax["region"] == r]
        rows.append({
            "region": r,
            "reach": r.split("-")[0],
            "composition": r.split("-", 1)[1],
            "n": int(counts.get(r, 0)),
            "mean_top_weight": round(float(sub["top_weight"].mean()), 3) if len(sub) else None,
        })

    radar = []
    if s_.has("taxonomy_radar"):
        rd = s_.table("taxonomy_radar")
        for axis, g in rd.groupby("axis"):
            radar.append({
                "axis": axis,
                "feature": g["feature"].iloc[0],
                "values": {r["region"]: float(r["median"]) for _, r in g.iterrows()},
                "q25": {r["region"]: float(r["q25"]) for _, r in g.iterrows()},
                "q75": {r["region"]: float(r["q75"]) for _, r in g.iterrows()},
            })

    return _clean({
        "regions": rows,
        "radar": radar,
        "radar_note": "Panel percentile of the amount-corrected feature. 50 is "
                      "the panel median, so no axis can mean 'more signal'.",
        # Computed on the FITTED genes only. The 9 under-evidenced genes carry
        # top_weight = NaN since the 2026-08-17 re-fit, and `NaN < 0.5` is False,
        # so including them would have counted every unlabelled gene as EXCEEDING
        # the threshold and reported 99.3% instead of 99.78%. Two decimals for the
        # same reason build_taxonomy.py now uses them: `.0%` rounded 99.78% up to
        # "100%", which overstates the load-bearing claim of this whole taxonomy.
        "mixture": {
            "n_fitted": int(tax["top_weight"].notna().sum()),
            "median_top_weight": round(float(tax["top_weight"].median()), 3),
            "pct_below_half": round(
                float((tax.loc[tax["top_weight"].notna(), "top_weight"] < 0.5)
                      .mean() * 100), 2),
            "n_at_or_above_half": int((tax["top_weight"] >= 0.5).sum()),
            "n_under_evidenced": int(tax["region"].isna().sum()),
            "uniform": round(1 / len(regions), 3),
        },
        "caveat": S.CONTINUUM_CAVEAT,
        # How these regions were arrived at, served so the derivation is visible
        # in the app rather than buried in a script docstring. Every number here
        # was measured, and two of the steps exist because an earlier version
        # was wrong.
        "derivation": [
            {"step": 1, "name": "Amount-corrected substrate",
             "detail": "Magnitude basis projected out of all 85 features, then "
                       "the 8 basis features dropped, leaving 77. Every "
                       "component correlates 0.000 with overall signal by "
                       "construction.",
             "result": "77 features x 1,846 genes"},
            {"step": 2, "name": "Split 1: reach",
             "detail": "KMeans k=2. This is sPC1, the largest axis: mean "
                       "contact distance, peak gap and far-distal fraction.",
             "result": "mid-range 1,074 / far-reaching 772",
             "evidence": "seed stability 0.998 (mid) and 0.965 (far); "
                         "cross-capture ARI 0.741"},
            {"step": 3, "name": "Split 2: composition, WITHIN each reach half",
             "detail": "k chosen per half, not forced symmetric. Composition is "
                       "orthogonal to reach (1-4% of its variance on sPC1), "
                       "which is why flat clustering could never find it: reach "
                       "dominates the distance metric and swamps it.",
             "result": "mid-range k=2, far-reaching k=3",
             "evidence": "mid k=3 yields a 4-gene group and drops stability "
                         "0.998 -> 0.796; far k=3 stays balanced at 238 and "
                         "0.950. Cross-capture ARI 0.616"},
            {"step": 4, "name": "Naming by architecture",
             "detail": "Each community assigned the element class it is richest "
                       "in, by one-to-one (Hungarian) assignment. No external "
                       "gene set is used: the old pipeline named its HK cluster "
                       "by best Eisenberg Fisher p, which made any later "
                       "'enriched for housekeeping' claim circular.",
             "result": "ctcf / enhancer / promoter per half"},
            {"step": 5, "name": "Memberships, not labels",
             "detail": "Softmax over distance to each region centroid. Every "
                       "gene in the panel has a top weight below 0.5, so a hard "
                       "label would assert a belonging the data does not "
                       "support.",
             "result": "median top weight 0.28 against 0.20 for uniform"},
        ],
        "k_evidence": [
            {"half": "mid-range", "n": 1074, "k": 2, "smallest": 519,
             "stability": 0.998, "chosen": True},
            {"half": "mid-range", "n": 1074, "k": 3, "smallest": 4,
             "stability": 0.796, "chosen": False},
            {"half": "mid-range", "n": 1074, "k": 4, "smallest": 3,
             "stability": 0.853, "chosen": False},
            {"half": "far-reaching", "n": 772, "k": 2, "smallest": 333,
             "stability": 0.965, "chosen": False},
            {"half": "far-reaching", "n": 772, "k": 3, "smallest": 238,
             "stability": 0.950, "chosen": True},
            {"half": "far-reaching", "n": 772, "k": 4, "smallest": 92,
             "stability": 0.916, "chosen": False},
        ],
    })


# Reach bands for the per-gene evidence grid. The taxonomy's reach split is
# 50-250 kb against >250 kb (frac_distal vs frac_far_distal); the <50 kb row is
# shown because it is where most signal usually is and its absence would make
# the shares unreadable.
def _annotate_mixture(mix, grid: dict) -> list[dict]:
    """Per-gene membership over the 3 x 3 band-by-class grid.

    COMPLETE BY CONSTRUCTION, 2026-08-17. Nine cells (prox / mid / far x
    promoter / enhancer / ctcf) cover 100% of every gene's peak signal, so there
    is no residual and no "outside the scheme" component to explain.

    Why not the clustered regions. Those are the panel-level grouping and their
    reach level is RELATIVE (contained vs extended); a gene's dominant absolute
    band matches its clustered reach label only 20% of the time. Using them to
    describe one gene produced statements that contradicted themselves, e.g.
    BCCIP at "24% far-enhancer" beside "no peaks of this kind". The clustered
    regions are still reported for that gene, separately, as the group it was
    assigned to -- a different and legitimate question.

    Composition, by contrast, agrees: the clustered class matches the dominant
    class for 63% of genes (promoter 72%, enhancer 57%, ctcf 51%) against 33% by
    chance. That dimension of the clustering is sound and is unchanged.
    """
    rows = {r["band"]: r for r in (grid.get("rows") or [])}
    if not rows or grid.get("no_peaks_called"):
        return [{
            "region": "no peaks called",
            "weight": 1.0,
            "weight_within_described": None,
            "own_signal_pct": None,
            "supported": True,
            "note": "no peaks were called for this gene, so there is nothing to "
                    "describe; its MCC signal is normal and the annotation is "
                    "what is missing",
        }]
    out = []
    for band in ("prox", "mid", "far"):
        row = rows.get(band) or {}
        for comp in ("promoter", "enhancer", "ctcf"):
            v = row.get(comp)
            share = 0.0 if v is None else float(v)
            out.append({
                "region": f"{band}-{comp}",
                "band_label": _BAND_LABEL[band],
                "weight": round(share / 100.0, 4),
                "weight_within_described": round(share / 100.0, 4),
                "own_signal_pct": round(share, 1),
                "supported": bool(share >= 1.0),
            })
    out.sort(key=lambda r: -r["weight"])
    return out


def _panel_reference(rows: list[dict]) -> dict:
    """Panel median, this gene's percentile, and the zero-share, per cell."""
    ref = _panel_cell_medians()
    pct = {}
    for r in rows:
        for c in ("promoter", "enhancer", "ctcf"):
            key = f"{r['band']}-{c}"
            v = r.get(c)
            d = ref["_dists"].get(key)
            if v is None or d is None or not len(d):
                continue
            # Strictly-below, so a gene sitting in the zero pile does not get
            # credit for the whole tie.
            pct[key] = round(float((d < float(v)).mean() * 100), 0)
    return {
        "cells": ref["cells"],
        "bands": ref["bands"],
        "zero_pct": ref["zero_pct"],
        "percentile": pct,
        "note": "Median is the clearer reference for the prox and mid cells. In "
                "the far row 41-63% of genes sit at zero, so a percentile there "
                "is inflated by ties and should be read with the zero-share.",
    }


@lru_cache(maxsize=1)
def _panel_cell_medians() -> dict:
    """Panel median share for each (band x class) cell, plus band totals.

    Added 2026-08-17. A cell reading 40% is uninterpretable on its own: for
    prox-enhancer the panel median is 17.4%, so 40% is the 85th percentile, but
    for prox-promoter the median is 17.7% and for far-enhancer it is 0.0%. Every
    number in the gene table now carries its reference.

    Cached: it is a whole-panel statistic and does not change between requests.
    """
    s_ = get_store()
    pk = s_.table("peaks").merge(s_.genes[["symbol_key"]], on="symbol_key", how="inner")
    absd = pk["offset_bp"].abs()
    band = pd.cut(absd, [0, 50_000, 250_000, 10 ** 9],
                  labels=["prox", "mid", "far"])
    pk = pk.assign(band=band)
    C = pk.pivot_table(index="symbol_key", columns=["band", "re"],
                       values="peak_max", aggfunc="sum", observed=True).fillna(0)
    C = C.div(C.sum(axis=1), axis=0) * 100
    cells, bands, dists, zero = {}, {}, {}, {}
    for b in ("prox", "mid", "far"):
        tot = None
        for c in ("promoter", "enhancer", "ctcf"):
            col = C[(b, c)] if (b, c) in C.columns else None
            key = f"{b}-{c}"
            cells[key] = round(float(col.median()), 1) if col is not None else 0.0
            # Kept so a percentile can be computed per gene, and so the share of
            # genes at exactly zero can be reported. That share matters: 63% of
            # genes have NO far-range enhancer contact, so any nonzero value
            # ranks above 63% and a percentile there reads as more notable than
            # it is. The median is the better reference for the six prox/mid
            # cells (7-24% zeros); the percentile earns its place in the far row.
            dists[key] = (np.sort(col.to_numpy(float))
                          if col is not None else np.zeros(1))
            zero[key] = (round(float((col == 0).mean() * 100), 0)
                         if col is not None else 100.0)
            tot = col if tot is None else tot + col
        bands[b] = round(float(tot.median()), 1) if tot is not None else 0.0
    return {"cells": cells, "bands": bands, "_dists": dists, "zero_pct": zero}


def _beyond_composition(grid: dict) -> dict:
    """Element-class composition of contacts BEYOND 50 kb.

    The clustering keys on reach x composition, so this is the slice of the
    evidence table that explains an assignment. Every gene's largest single cell
    is usually proximal, which made the assigned group look like it contradicted
    the table: a gene reading "mostly promoter, within 50 kb" could be assigned
    extended-ctcf because 26% of its signal is CTCF beyond 50 kb.
    """
    rows = {r["band"]: r for r in (grid.get("rows") or [])}
    out, total = {}, 0.0
    for cls in ("promoter", "enhancer", "ctcf"):
        v = sum(float((rows.get(b) or {}).get(cls) or 0.0) for b in ("mid", "far"))
        out[cls] = v
        total += v
    if total <= 0:
        return {"total_pct": 0.0, "shares": {}, "dominant": None}
    return {
        "total_pct": round(total, 1),
        "shares": {k: round(100.0 * v / total, 1) for k, v in out.items()},
        "dominant": max(out, key=out.get),
    }


_REACH_BANDS = [(0, 50_000, "prox"),
                (50_000, 250_000, "mid"),
                (250_000, 10 ** 9, "far")]
_BAND_LABEL = {"prox": "<50 kb", "mid": "50-250 kb", "far": ">250 kb"}


def _reach_composition_grid(s_, gid: str) -> dict:
    """Share of a gene's peak signal by reach band x element class."""
    if not s_.has("peaks"):
        return {}
    pk = s_.table("peaks")
    sym = s_.genes.set_index("gene_id")["symbol_key"].get(gid)
    g = pk[pk["symbol_key"] == sym]
    if g.empty:
        # NO PEAKS CALLED. Distinct from "peaks exist but outside the named
        # bands", and it must not be treated as missing-and-therefore-ignorable:
        # these genes have normal MCC signal (PADI4 total 260, max bin 34) and
        # normal CTCF ChIP, so this is an annotation gap, not an empty locus.
        return {
            "n_peaks": 0,
            "rows": [{"band": b[2], "n_peaks": 0,
                      "promoter": None, "enhancer": None, "ctcf": None}
                     for b in _REACH_BANDS],
            "total_signal": 0.0,
            "no_peaks_called": True,
            "described_pct": 0.0,
            "proximal_pct": None,
            "weakly_described": True,
            "coverage_note": "No peaks were called for this gene, so there is no "
                             "element-class or reach evidence at all. Its MCC "
                             "signal is normal; the annotation is missing. The "
                             "taxonomy cannot describe it.",
        }
    absd = g["offset_bp"].abs()
    total = float(g["peak_max"].sum())
    rows = []
    for lo, hi, name in _REACH_BANDS:
        sel = g[(absd >= lo) & (absd < hi)]
        cells = {}
        for cls in ("promoter", "enhancer", "ctcf"):
            v = float(sel[sel["re"] == cls]["peak_max"].sum())
            cells[cls] = round(100.0 * v / total, 1) if total else 0.0
        rows.append({"band": name, "n_peaks": int(len(sel)), **cells})
    # HOW MUCH OF THIS GENE THE TAXONOMY ACTUALLY DESCRIBES.
    #
    # The reach levels are the 50-250 kb and >250 kb bands, so any signal inside
    # 50 kb is outside both and is not described by the scheme at all. That is
    # not a corner case: 55.5% of a typical gene's peak signal is proximal, and
    # 497 of 1,843 genes have over 70% there. LCK has 87%, which is why three of
    # its five region weights had no supporting peaks -- its architecture is
    # almost entirely in a zone the taxonomy does not name.
    #
    # Reporting this turns "LCK is 20% far-enhancer" from a false statement into
    # a bounded one: the taxonomy describes 13% of LCK, and within that 13% it
    # leans as stated.
    # All three bands are named now, so coverage is 100% by construction and
    # there is no "outside the scheme". `proximal_pct` is kept because it is
    # genuinely informative: the panel median is ~55%.
    prox = rows[0]["promoter"] + rows[0]["enhancer"] + rows[0]["ctcf"]
    described = 100.0
    return {
        "n_peaks": int(len(g)),
        "total_signal": round(total, 1),
        "rows": rows,
        "classes": ["promoter", "enhancer", "ctcf"],
        "proximal_pct": round(prox, 1),
        "described_pct": described,
        # Panel reference for every cell and band, so a number can be read as
        # high or typical rather than just large.
        "panel_median": _panel_reference(rows),
        "coverage_note": (
            f"The taxonomy's reach levels cover 50-250 kb and >250 kb, so it "
            f"describes {described}% of this gene's peak signal; {round(prox,1)}% "
            f"lies inside 50 kb, which the scheme does not name. Panel median "
            f"for proximal signal is ~55%."
        ),
        "weakly_described": bool(described < 25),
    }


@app.get("/api/taxonomy/map")
def taxonomy_map(x: str = S.DEFAULT_PLANE[0], y: str = S.DEFAULT_PLANE[1]):
    """Every gene's position on the chosen plane, coloured by nearest region.

    Served as parallel arrays rather than objects: 1,846 points x 3 numbers is
    a fifth the payload of a list of dicts, and this is redrawn on every axis
    change.

    Works over UMAP as well as the sPC planes, deliberately. The regions were
    defined in the 77-dimensional corrected space, NOT on any 2D embedding, so
    seeing them on two different projections is a check rather than a
    tautology: if they only look coherent on the plane they were drawn in, they
    are an artefact of that plane.
    """
    s_ = store()
    if not s_.has("taxonomy"):
        raise HTTPException(503, "taxonomy table not built")
    _distinct_axes(x, y)
    emb = _embeddings_all(s_)
    for ax in (x, y):
        if ax not in emb.columns:
            raise HTTPException(400, f"unknown axis {ax!r}")
    tax = s_.table("taxonomy")[["gene_id", "region", "top_weight"]]
    # The 9 under-evidenced genes carry region = NaN since the 2026-08-17 re-fit.
    # They are dropped from the MAP because a point with no region has no colour
    # and no membership to size it by, and plotting them uncoloured would invite
    # the reading that they are a sixth region. They remain visible on the gene
    # page, which states plainly that the taxonomy does not describe them.
    #
    # Dropped explicitly rather than left to `sorted()`, which raised
    # TypeError: '<' not supported between 'float' and 'str' on the NaN.
    tax = tax[tax["region"].notna()]
    df = emb[["gene_id", x, y]].merge(tax, on="gene_id", how="inner")
    g = s_.genes[["gene_id", "gene_symbol"]]
    df = df.merge(g, on="gene_id", how="left")
    regions = sorted(df["region"].dropna().unique())
    idx = {r: i for i, r in enumerate(regions)}
    return _clean({
        "x_axis": {"key": x, "label": S.ALL_AXIS_LABELS.get(x, x),
                   "poles": S.ALL_AXIS_POLES.get(x)},
        "y_axis": {"key": y, "label": S.ALL_AXIS_LABELS.get(y, y),
                   "poles": S.ALL_AXIS_POLES.get(y)},
        "regions": regions,
        "n": int(len(df)),
        "xs": [round(float(v), 4) for v in df[x]],
        "ys": [round(float(v), 4) for v in df[y]],
        "r": [int(idx[v]) for v in df["region"]],
        "w": [round(float(v), 3) for v in df["top_weight"]],
        "symbols": df["gene_symbol"].astype(str).tolist(),
        "note": "Regions were defined in the 77-dimensional corrected space, "
                "not on this plane. Opacity is the gene's top membership, so "
                "faint points are genes that sit between regions, which is most "
                "of the panel.",
    })


@lru_cache(maxsize=1)
def _tvd_quantiles() -> dict:
    """Panel distribution of total variation distance from the uniform floor.

    TVD is the natural single-number summary of "does this gene's profile say
    anything": 0 means the weights are exactly uniform and carry no information,
    1 means a hard label. It is also exactly the total above-floor excess divided
    by 100, since the weights sum to 100 and the deviations therefore sum to
    zero, so the positive and negative deviations are equal in magnitude. That
    identity is what lets one bar carry both which regions a gene leans toward
    and how much lean there is in total.

    Served rather than hardcoded so the summary bar's scale follows the data.
    """
    s_ = store()
    if not s_.has("taxonomy"):
        return {}
    tax = s_.table("taxonomy")
    wcols = [c for c in tax.columns if c.startswith("w_")]
    if not wcols:
        return {}
    Wm = tax[wcols].to_numpy(float)
    k = Wm.shape[1]
    tvd = 0.5 * np.abs(Wm - 1.0 / k).sum(1)
    tvd = tvd[np.isfinite(tvd)]
    return {
        "median": round(float(np.median(tvd)), 3),
        "p90": round(float(np.quantile(tvd, 0.90)), 3),
        "max": round(float(tvd.max()), 3),
    }


@lru_cache(maxsize=1)
def _nearest_distance_quantiles() -> dict:
    """Panel distribution of the nearest-centroid distance.

    Needed because the softmax has a SECOND failure mode, opposite to the uniform
    floor and found the same way (2026-08-17, from SACS). The weight uses
    exp(-d^2/tau), so a gene far from every centroid gets a confident-looking
    profile from a small RELATIVE gap: SACS sits 31.15 to 34.39 from the five
    centroids, a 10% spread, but in d^2 that is a gap of 212, which the softmax
    turns into 62.9% for the nearest region. It is not near that region, it is
    far from all of them.

    This is not an isolated case. Across the panel, corr(nearest distance, top
    weight) = +0.436, and 29 of the 54 genes reaching any weight of 40% sit in
    the top distance decile, so more than half the confident assignments come
    from a tenth of the panel selected for being unlike everything.
    """
    s_ = store()
    if not s_.has("taxonomy"):
        return {}
    tax = s_.table("taxonomy")
    dcols = [c for c in tax.columns if c.startswith("d_")]
    if not dcols:
        return {}
    near = tax[dcols].to_numpy(float).min(1)
    near = near[np.isfinite(near)]
    return {
        "median": round(float(np.median(near)), 2),
        "p90": round(float(np.quantile(near, 0.90)), 2),
        "values": near,
    }


def _nearest_pct(near: float | None) -> float | None:
    """This gene's nearest-centroid distance as a panel percentile."""
    q = _nearest_distance_quantiles()
    vals = q.get("values")
    if near is None or vals is None or not len(vals):
        return None
    return round(float((vals < near).mean() * 100), 1)


@app.get("/api/genes/{gene}/taxonomy")
def gene_taxonomy(gene: str):
    """One gene's membership across the regions, as a blend.

    Never a single label. Every gene in the panel has a top weight below 0.5,
    so naming one region for a gene would overstate what the data supports.
    """
    s_ = store()
    gid = s_.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")
    if not s_.has("taxonomy"):
        raise HTTPException(503, "taxonomy table not built")
    tax = s_.table("taxonomy")
    row = tax[tax["gene_id"].astype(str) == gid]
    if row.empty:
        raise HTTPException(404, f"no taxonomy row for {gene!r}")
    row = row.iloc[0]
    wcols = [c for c in tax.columns if c.startswith("w_")]
    mix = sorted(((c[2:], float(row[c])) for c in wcols), key=lambda x: -x[1])

    # THE DERIVATION, so a weight can be checked rather than trusted.
    # weight_i = exp(-d_i^2 / tau) / sum_j exp(-d_j^2 / tau)
    # Showing 0.38 without showing the distance that produced it is the same
    # failure as showing an archetype name without its loadings.
    grid_for_assignment = _reach_composition_grid(s_, gid)
    dcols = {c[2:]: float(row[c]) for c in tax.columns if c.startswith("d_")}
    tau = float(row["tau"]) if "tau" in tax.columns else None
    near = min(dcols.values()) if dcols else None

    # THE UNIFORM FLOOR. A softmax over K regions cannot go below 1/K for a gene
    # that resembles nothing in particular, so a raw percentage is unreadable:
    # with K=5, "20%" is not one fifth of the gene's architecture, it is exactly
    # no evidence either way. Measured on the panel 2026-08-17, the softmax has
    # almost no room to move: nearest-centroid distance median 7.29 against
    # furthest 9.67, a spread of only 31% of the nearest distance, because in 77
    # dimensions all five centroids sit at similar distance from any point. The
    # consequence is that 63% of genes have their ENTIRE profile within 10 pp of
    # uniform and only 3% reach any weight of 40%.
    #
    # This surfaced from a real misreading: a gene 79% within 50 kb and 0% beyond
    # 250 kb showed 20% extended-enhancer and 19% extended-ctcf, which looks like
    # a contradiction and is in fact the floor (and, at 19%, marginally below it).
    # So every weight is reported as a DEVIATION from the floor as well as a
    # share, and a profile that never leaves the floor is flagged as carrying no
    # information rather than being drawn as five confident-looking bars.
    uniform_pct = 100.0 / len(mix) if mix else None

    steps = []
    if dcols and tau:
        for r, w in mix:
            d = dcols.get(r)
            if d is None:
                continue
            steps.append({
                "region": r,
                "distance": round(d, 3),
                # How much further than the closest region, which is what the
                # softmax actually responds to.
                "excess_over_nearest": round(d - near, 3),
                "exp_term": round(float(pow(2.718281828, -(d * d) / tau)), 5),
                "weight": round(w, 4),
                # Signed distance from the uniform floor, in percentage points.
                # Positive means the gene is closer to this region than chance;
                # negative is evidence AGAINST, which a raw share cannot express.
                "deviation_pp": (round(w * 100 - uniform_pct, 1)
                                 if uniform_pct is not None else None),
            })

    return _clean({
        "gene_id": gid,
        "nearest_region": row["region"],
        "top_weight": float(row["top_weight"]),
        "mixedness": float(row["mixedness"]),
        # EXCLUDED FROM THE FIT, 2026-08-17. Nine genes have fewer than
        # config.ACTIVE_NPEAKS_MIN = 3 called peaks, so there is no architecture
        # to compare against a centroid. They are marked rather than hidden: a
        # gene that simply renders nothing looks like a loading failure, and one
        # of the nine (SACS, 1 peak) still has a contact table to draw, which
        # without this flag would sit there implying a description exists.
        "under_evidenced": bool(row.get("under_evidenced", False)),
        "under_evidenced_note": (
            "Fewer than 3 called peaks, so this gene was excluded from the "
            "taxonomy fit and carries no region. Its MCC signal is normal; the "
            "peak annotation is what is missing. Left in, these genes took the "
            "three highest memberships in the panel, because the softmax squares "
            "distance and a gene with no evidence sits far from every centroid."
            if bool(row.get("under_evidenced", False)) else None),
        # Each weight is cross-checked against the gene's OWN peaks. A softmax
        # over five regions distributes weight everywhere, so a gene can carry a
        # substantial weight for a region it has no evidence of: LCK has ZERO
        # peaks beyond 250 kb and a 20% far-enhancer weight; GATA3 has ZERO
        # promoter-classed signal and a 17% far-promoter weight. The weight is
        # not wrong -- it correctly reports proximity to a centroid in 77
        # dimensions -- but read as biology it is an overstatement, so it is
        # marked rather than left to be misread.
        "mixture": _annotate_mixture(mix, _reach_composition_grid(s_, gid)),
        # How decisive the assignment is. Measured panel-wide: median margin
        # over the runner-up is 8.1%, and 59% of genes have a runner-up within
        # 10%. So a single group label, shown alone, overstates for most genes.
        "assignment": {
            "group": row["region"],
            "runner_up": (steps[1]["region"] if len(steps) > 1 else None),
            "margin_pct": (round((steps[1]["distance"] - steps[0]["distance"])
                                 / steps[0]["distance"] * 100, 1)
                           if len(steps) > 1 and steps[0]["distance"] else None),
            "close_call": (bool(len(steps) > 1 and steps[0]["distance"]
                                and (steps[1]["distance"] - steps[0]["distance"])
                                / steps[0]["distance"] * 100 < 10)),
            "panel_median_margin_pct": 8.1,
            # What the clustering keys on is reach x composition, so the
            # composition of contacts BEYOND the proximal band is the part of the
            # table that explains the assignment. Without this the assignment
            # looks like it contradicts the biggest cell, which is usually
            # proximal for every gene.
            "beyond_50kb": _beyond_composition(grid_for_assignment),
        },
        # Whether this gene's profile says anything at all. Panel references are
        # measured, not asserted: see the uniform-floor comment above.
        "flatness": {
            "uniform_pct": round(uniform_pct, 1) if uniform_pct else None,
            "n_regions": len(mix),
            # Total variation distance from uniform. 0 = the profile carries no
            # information, 1 = a hard label.
            "tvd": (round(0.5 * sum(abs(w - 1.0 / len(mix)) for _, w in mix), 3)
                    if mix else None),
            "panel_median_tvd": _tvd_quantiles().get("median", 0.119),
            "panel_p90_tvd": _tvd_quantiles().get("p90"),
            "panel_max_tvd": _tvd_quantiles().get("max"),
            "max_deviation_pp": (round(max(abs(w * 100 - uniform_pct)
                                           for _, w in mix), 1)
                                 if mix and uniform_pct else None),
            # A profile that never leaves the floor is not a five-way blend, it
            # is an absence of evidence, and it is the common case.
            "near_uniform": (bool(mix and uniform_pct
                                  and max(abs(w * 100 - uniform_pct)
                                          for _, w in mix) <= 10)),
            "panel_pct_near_uniform": 63,
            "panel_pct_above_40": 3,
            # THE OPPOSITE FAILURE. See _nearest_distance_quantiles: a gene far
            # from every centroid gets a confident profile out of a small
            # relative gap, because the softmax squares the distance. A high
            # top weight here means "unlike everything" rather than "clearly
            # this one", and it must not be read as a decisive assignment.
            "nearest_distance": round(near, 2) if near is not None else None,
            "nearest_distance_pct": _nearest_pct(near),
            "panel_median_nearest_distance": _nearest_distance_quantiles().get("median"),
            "distance_inflated": (
                bool(near is not None
                     and _nearest_distance_quantiles().get("p90") is not None
                     and near >= _nearest_distance_quantiles()["p90"])),
            "distance_note": "corr(nearest distance, top weight) is +0.436 "
                             "across the panel, and 29 of the 54 genes reaching "
                             "any weight of 40% sit in the top distance decile. "
                             "So over half the confident-looking assignments "
                             "come from a tenth of the panel selected for being "
                             "unlike everything, not for resembling one group.",
            "note": "A softmax over 5 regions cannot fall below 20% for a gene "
                    "that resembles nothing in particular, so 20% is the floor "
                    "rather than a fifth of the architecture. The distances "
                    "leave little room to move: nearest-centroid distance is "
                    "7.29 against a furthest of 9.67 at the panel median, "
                    "because in 77 dimensions every centroid is far from every "
                    "point. 63% of genes have their whole profile within 10 pp "
                    "of the floor.",
        },
        "derivation": {
            "formula": "weight_i = exp(-d_i^2 / tau) / sum_j exp(-d_j^2 / tau)",
            "tau": tau,
            "tau_note": "Temperature, set to the median nearest-centroid squared "
                        "distance across the panel so the scale is a property of "
                        "the data. Set it from the within-region spread instead "
                        "and the softmax flattens to near-uniform, which is what "
                        "the first version did.",
            "nearest_distance": round(near, 3) if near is not None else None,
            "steps": steps,
            "reading": "Distances are in the 77-dimensional amount-corrected "
                       "space. In that many dimensions everything is far from "
                       "everything, so what matters is the GAP between the "
                       "nearest region and the rest, not the absolute distance.",
        },
        # THE EVIDENCE GRID. The taxonomy's two axes are reach and composition,
        # and both are directly visible in this gene's own peaks: reach is where
        # they sit in the window, composition is what class they are. So the
        # membership is a summary of two things you can point at, and this grid
        # is that summary before any modelling.
        #
        # It also says things the blend cannot. GATA3 carries a 17% far-promoter
        # weight while having ZERO promoter-classed signal: a small weight means
        # "not especially far from that centroid in 77 dimensions", NOT evidence
        # of that architecture. Showing both lets a reader see the disagreement
        # instead of inheriting it.
        "evidence": _reach_composition_grid(s_, gid),
        "evidence_note": "Share of this gene's peak signal by reach band and "
                         "element class. This is the raw material the membership "
                         "summarises. A region can carry weight while its "
                         "element class shows 0% here, which means the gene is "
                         "merely not far from that centroid, not that it has "
                         "that architecture.",
        "note": "A blend, not a label. Top weight is below 0.5 for every gene "
                "in the panel, so the nearest region is a coordinate, not a "
                "category.",
    })


# ---------------------------------------------------------------------------
# compare suggestions
# ---------------------------------------------------------------------------

def _shape_matrix(s_):
    """Every gene's coordinates in the amount-corrected space, plus an index.

    All retained sPCs, not the two on screen. A pair of genes can sit on top of
    each other in the sPC1 x sPC2 plane and be far apart in the space, so
    picking neighbours from the displayed plane would suggest comparisons that
    look adjacent and are not.
    """
    emb = _embeddings_all(s_)
    cols = [c for c in emb.columns if c.startswith("spc")]
    cols.sort(key=lambda c: int(c[3:]))
    sub = emb[["gene_id"] + cols].dropna()
    return sub["gene_id"].astype(str).to_numpy(), sub[cols].to_numpy(float), cols


@lru_cache(maxsize=1)
def _median_nn_distance() -> float:
    """Median nearest-neighbour distance across the panel, for scale.

    Computed once. A suggested neighbour at 4.1 means nothing until you know
    whether the typical gene's nearest neighbour is at 2 or at 8.
    """
    _, X, _ = _shape_matrix(store())
    # 1,846 x 1,846 in float32 is ~13 MB, so the full matrix is cheaper than
    # any clever alternative and this runs once per process.
    g2 = (X ** 2).sum(1)
    d2 = g2[:, None] + g2[None, :] - 2.0 * (X @ X.T)
    np.fill_diagonal(d2, np.inf)
    return round(float(np.median(np.sqrt(np.maximum(d2.min(1), 0.0)))), 2)


@app.get("/api/genes/{gene}/compare-suggestions")
def compare_suggestions(gene: str, limit_per_group: int = 3):
    """Which genes are worth comparing this one against, and why.

    A free-text box assumes the reader already has a second gene in mind. Usually
    they do not, and the comparisons that actually teach something are not the
    ones that come to mind: they are the ones that hold one axis fixed and vary
    the other. The taxonomy has exactly two axes, so this offers

      nearest        the closest architecture in the space, the "is this gene
                     unusual at all" question
      isolate_reach  same composition, most different reach
      isolate_comp   same reach, most different composition
      exemplar       the gene sitting closest to each region centroid, so the
                     reader can see what the label is supposed to look like,
                     including the runner-up region for a close call
      contrast       the furthest gene in the panel, the upper bound on how
                     different two genes here get

    Distances are reported so the reader can weigh them against the noise floor:
    across the 116 twice-captured genes a gene sits 5.36 from ITSELF and 10.49
    from another gene, on the 91-feature substrate. A suggested pair closer than
    that separation is closer than the measurement error, which is worth knowing
    before reading a difference as biology.
    """
    s_ = store()
    gid = s_.resolve(gene)
    if gid is None:
        raise HTTPException(404, f"gene {gene!r} not in the panel")

    ids, X, cols = _shape_matrix(s_)
    pos = {g: i for i, g in enumerate(ids)}
    if gid not in pos:
        raise HTTPException(503, "no shape coordinates for this gene")
    i0 = pos[gid]
    d = np.sqrt(((X - X[i0]) ** 2).sum(1))

    genes = s_.table("genes")
    genes = genes.assign(gene_id=genes["gene_id"].astype(str))
    meta = genes.set_index("gene_id")

    # NEVER SUGGEST A GENE WITH TOO LITTLE EVIDENCE TO COMPARE AGAINST.
    #
    # Caught immediately on the first run: "Maximum contrast" returned TTPAL,
    # which has ZERO called peaks. It is the furthest gene in the panel because
    # it has no data, not because its architecture differs, and offering it as
    # the exemplar of difference is the worst possible suggestion.
    #
    # Nine genes are affected: 3 with no peaks at all (ASB9, PADI4, TTPAL) and 6
    # with one or two. They are also the panel's most confident-looking
    # assignments (median top weight 48.7% at zero peaks and 36.4% at one or two,
    # against 28.3% panel-wide) precisely because a gene with no evidence lands
    # far from every centroid and the softmax squares distance. SACS has ONE
    # peak, and its entire architecture description rests on it.
    #
    # The threshold is config.ACTIVE_NPEAKS_MIN = 3, the pipeline's own existing
    # floor for calling a gene active, rather than a number invented here.
    MIN_PEAKS = 3
    n_peaks = (s_.table("peaks").groupby("symbol_key").size()
               if s_.has("peaks") else None)
    thin: set[str] = set()
    if n_peaks is not None:
        counts = meta["symbol_key"].map(n_peaks).fillna(0)
        thin = set(counts[counts < MIN_PEAKS].index.astype(str))
    me = meta.loc[gid] if gid in meta.index else None
    my_reach = str(me["reach"]) if me is not None else None
    my_comp = str(me["composition"]) if me is not None else None

    def pack(idx: int, why: str) -> dict:
        g = str(ids[idx])
        m = meta.loc[g] if g in meta.index else None
        return {
            "gene_symbol": (str(m["gene_symbol"]) if m is not None else g),
            "gene_id": g,
            "region": (str(m["region"]) if m is not None else None),
            "distance": round(float(d[idx]), 2),
            "why": why,
        }

    # Drop the gene itself by identity rather than by position. Slicing [1:] was
    # wrong for a thin PRIMARY gene (SACS, 1 peak): the filter removes it from
    # the order, so position 0 is already another gene and [1:] silently skipped
    # the true nearest neighbour.
    order = np.array([j for j in np.argsort(d)
                      if str(ids[j]) != gid and str(ids[j]) not in thin])
    groups: list[dict] = []

    nearest = [pack(int(j), "closest architecture in the panel")
               for j in order[:limit_per_group]]
    if nearest:
        groups.append({
            "key": "nearest",
            "label": "Most similar",
            "note": "Nearest in the amount-corrected space across all retained "
                    "components. If even the closest gene is far, this gene is "
                    "architecturally unusual.",
            "genes": nearest,
        })

    # One axis at a time. This is the comparison that isolates a cause, and it is
    # the one a reader will not think to construct.
    if my_reach and my_comp:
        same_comp = np.array([
            j for j in order
            if j != i0 and str(ids[j]) in meta.index
            and str(meta.loc[str(ids[j])]["composition"]) == my_comp
            and str(meta.loc[str(ids[j])]["reach"]) != my_reach
        ])
        if same_comp.size:
            groups.append({
                "key": "isolate_reach",
                "label": f"Same composition ({my_comp}), other reach",
                "note": "Holds element-class composition fixed and flips reach, "
                        "so a difference between the two is about distance, not "
                        "about what the contacts are.",
                "genes": [pack(int(j), f"{my_comp}, opposite reach")
                          for j in same_comp[:limit_per_group]],
            })

        same_reach = np.array([
            j for j in order
            if j != i0 and str(ids[j]) in meta.index
            and str(meta.loc[str(ids[j])]["reach"]) == my_reach
            and str(meta.loc[str(ids[j])]["composition"]) != my_comp
        ])
        if same_reach.size:
            groups.append({
                "key": "isolate_comp",
                "label": f"Same reach ({my_reach}), other composition",
                "note": "Holds reach fixed and varies composition, so a "
                        "difference is about what the contacts are rather than "
                        "how far they go.",
                "genes": [pack(int(j), f"{my_reach}, different composition")
                          for j in same_reach[:limit_per_group]],
            })

    # Region exemplars, drawn from the taxonomy's own distance columns rather
    # than recomputed, so the app cannot disagree with the assignment it shows.
    if s_.has("taxonomy"):
        tax = s_.table("taxonomy")
        tax = tax.assign(gene_id=tax["gene_id"].astype(str))
        trow = tax[tax["gene_id"] == gid]
        # Exemplars are drawn from this table too, so it needs the same evidence
        # filter: an "exemplar" with two peaks teaches nothing.
        pool_tax = tax[~tax["gene_id"].isin(thin)]
        dcols = [c for c in tax.columns if c.startswith("d_")]
        if not trow.empty and dcols:
            mine = trow.iloc[0]
            ranked = sorted(dcols, key=lambda c: float(mine[c]))
            picks = []
            for rank, c in enumerate(ranked[:2]):
                region = c[2:]
                cand = pool_tax[pool_tax["gene_id"] != gid].nsmallest(1, c)
                if cand.empty:
                    continue
                best = cand.iloc[0]
                bid = str(best["gene_id"])
                m = meta.loc[bid] if bid in meta.index else None
                picks.append({
                    "gene_symbol": (str(m["gene_symbol"]) if m is not None else bid),
                    "gene_id": bid,
                    "region": region,
                    "distance": (round(float(d[pos[bid]]), 2)
                                 if bid in pos else None),
                    "why": ("clearest example of your assigned group"
                            if rank == 0 else
                            "clearest example of the runner-up group"),
                })
            if picks:
                groups.append({
                    "key": "exemplar",
                    "label": "Group exemplars",
                    "note": "The gene sitting closest to each region centroid. "
                            "Comparing against these shows what the label is "
                            "meant to look like, which matters most when the "
                            "assignment is a close call.",
                    "genes": picks,
                })

    if order.size > 1:
        groups.append({
            "key": "contrast",
            "label": "Maximum contrast",
            "note": "The furthest gene in the panel. An upper bound on how "
                    "different two genes get here.",
            "genes": [pack(int(order[-1]), "furthest gene in the panel")],
        })

    return _clean({
        "gene_id": gid,
        "n_components": len(cols),
        "distance_note": "Euclidean distance in the amount-corrected space over "
                         f"{len(cols)} components. For scale, repeat captures of "
                         "the SAME gene sit 5.36 apart against 10.49 for two "
                         "different genes, so a pair closer than a few units "
                         "differs by less than the measurement error.",
        # The reference that makes "closest gene: 4.1" readable. Without it a
        # distance is a number with no scale attached.
        "panel_median_nearest": _median_nn_distance(),
        "groups": groups,
    })
