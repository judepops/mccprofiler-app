"""FastAPI app — serves the store, nothing else.

Localhost only by decision (PLAN.md scope table): unpublished MCC never leaves
the workstation. No auth, because there is no network surface to protect.

Every endpoint returns data, never a rendered image. The frontend draws.

    conda activate mccapp
    uvicorn app.main:app --reload --port 8000    # from backend/
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.middleware.cors import CORSMiddleware

from . import store_schema as S
from .store import StoreMissing, get_store

# Feature distance bands, mirrored from mccprofiler.config.DISTANCE_BANDS so the
# plot shades the same cutoffs the features are actually defined on.
DISTANCE_BANDS = {
    "promoter_proximal": (0, 10_000),
    "local": (10_000, 50_000),
    "distal": (50_000, 250_000),
    "far_distal": (250_000, 1_000_000),
}

app = FastAPI(
    title="mccprofiler-app",
    description="Gene position in the MCC regulatory continuum.",
    version="0.1.0",
)

# Vite dev server. Localhost only — see module docstring.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


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

    return _clean(out)


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
    centre. This is NOT a genomic-coordinate view — the alignment is what makes
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
# dimensions — the primary coordinate system (PLAN.md §1b)
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
                    "above by index — variance shares differ between spaces.",
            "rows": _clean(rep.to_dict("records")) if rep is not None else [],
        },
        "display_rule": "Report max |rho| and r2. No survive/compromised verdicts — "
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
def scree():
    """Variance per component against a permuted-noise ceiling.

    The ceiling is PCA on a per-feature-permuted copy: the variance a component
    of this size explains when there is provably nothing to find. Components
    above it are the "real dimensions". This store reproduces the 19 that
    structure_vs_noise.tsv reports.
    """
    s_ = store()
    if not s_.has("pc_scree"):
        raise HTTPException(503, "pc_scree not built")
    df = s_.table("pc_scree")
    return _clean({
        "n_above_noise": int(df["above_noise"].sum()),
        "method": "parallel analysis against a per-feature-permuted null",
        "note": "Components above the ceiling carry structure that survives "
                "destroying all feature-feature covariance while keeping every "
                "marginal distribution intact.",
        "rows": df.to_dict("records"),
    })


@app.get("/api/dimensions/{pc}/loadings")
def loadings(pc: int, top: int = 15):
    """Feature loadings for one component — the evidence for its name.

    A named axis is an interpretation of its loadings. Serving the name without
    them would ask the reader to take the interpretation on trust.
    """
    s_ = store()
    if not s_.has("pc_loadings"):
        raise HTTPException(503, "pc_loadings not built")
    df = s_.table("pc_loadings")
    hit = df[df["pc"] == pc]
    if hit.empty:
        raise HTTPException(404, f"no loadings for PC{pc}")

    ranked = hit.reindex(hit["loading"].abs().sort_values(ascending=False).index)
    sel = ranked.head(top).sort_values("loading", ascending=False)

    scree_row = None
    if s_.has("pc_scree"):
        sc = s_.table("pc_scree")
        r = sc[sc["pc"] == pc]
        if not r.empty:
            scree_row = r.iloc[0].to_dict()

    return _clean({
        "pc": pc,
        "label": S.PC_LABELS.get(f"pc{pc}", f"PC{pc}"),
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
    emb = s.table("embeddings")

    for ax in (x, y):
        if ax not in emb.columns:
            raise HTTPException(
                400,
                f"unknown axis {ax!r}; available: "
                f"{', '.join(c for c in emb.columns if c not in ('gene_id', 'symbol_key'))}",
            )

    g = s.genes[["gene_id", "gene_symbol", "group", "max_posterior", "is_core"]]
    df = emb[["gene_id", x, y]].merge(g, on="gene_id", how="left")

    hl = s.resolve(highlight) if highlight else None

    return {
        "x_axis": {"key": x, "label": S.PC_LABELS.get(x, x),
                   "poles": S.PC_POLES.get(x)},
        "y_axis": {"key": y, "label": S.PC_LABELS.get(y, y),
                   "poles": S.PC_POLES.get(y)},
        "axes_available": [
            {"key": c, "label": S.PC_LABELS.get(c, c)}
            for c in emb.columns if c not in ("gene_id", "symbol_key")
        ],
        "highlight": hl,
        "is_umap": x.startswith("umap") or y.startswith("umap"),
        "is_null": x.startswith("umap_null") or y.startswith("umap_null"),
        "umap_caveat": S.UMAP_CAVEAT,
        "continuum_caveat": S.CONTINUUM_CAVEAT,
        "n": int(len(df)),
        "points": _clean(
            df.rename(columns={x: "x", y: "y"})[
                ["gene_id", "gene_symbol", "x", "y", "group", "max_posterior", "is_core"]
            ].to_dict("records")
        ),
    }


# ---------------------------------------------------------------------------
# lab — exploratory views, not part of the main app
# ---------------------------------------------------------------------------


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
    if not s_.has("reproducibility"):
        raise HTTPException(503, "reproducibility table not built")

    per_feature = s_.table("reproducibility").sort_values("spearman_rho", ascending=False)
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
        "note": "Least reproducible are the asymmetry features — candidates for "
                "removal, and the reason the profile view treats them cautiously.",
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
        emb = s_.table("embeddings").drop(columns=["symbol_key"], errors="ignore")
        df = df.merge(emb, on="gene_id", how="left")
    return Response(
        content=df.to_csv(index=False),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="mccprofiler_panel.csv"'},
    )


@app.get("/api/export/{gene}.csv")
def export_gene(gene: str):
    """One gene, everything known about it, as CSV.

    The CellProfiler analogy only holds if the output is portable — its
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


@app.get("/api/explain")
def explain():
    """The argument the app is making, with its numbers.

    Leads with the nested-baseline result because that is the justification for
    the 91-feature substrate existing at all — without it, `n_peaks` would do.
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
            "detail": "Using the 116 genes captured in BOTH panels — cluster once on "
                      "the pooled matrix, then assign each gene's two independent "
                      "captures separately — Cohen's kappa peaks at 0.72 for k=3-4. "
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
    they do not track — see ARCHETYPE_DISPLAY for the measurements.
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
    response leads with coverage, because the panel is 1,846 of ~20,000 genes —
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
        note = (f"Only {len(member)} panel genes — below the {S.MIN_GROUP_N}-gene "
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
            "label": S.PC_LABELS.get(a, a),
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
    """Externally-defined gene groups — the strongest non-circular evidence.

    These groups were not defined from the features, so unlike the archetypes
    they need no within/pooled-ratio argument to be interpretable.
    """
    s = store()
    if not s.has("cohort_membership"):
        raise HTTPException(503, "cohort_membership not built")

    counts = s.table("cohort_membership")["group"].value_counts()

    # Stratification statistics exist for only some sets — the seven the audit
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
                              "clustering — that would be circular.",
        },
        "rows": rows,
    })
