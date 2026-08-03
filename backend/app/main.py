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
from fastapi import FastAPI, HTTPException, Query
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
    out = {
        "gene_id": gid,
        "gene_symbol": symbol,
        "viewpoint": {"chrom": row.get("viewpoint_chrom"), "pos": row.get("viewpoint_pos")},
        "archetype": {
            "group": row.get("group"),
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


@app.get("/api/cohorts")
def cohorts():
    """Externally-defined gene groups — the strongest non-circular evidence.

    These groups were not defined from the features, so unlike the archetypes
    they need no within/pooled-ratio argument to be interpretable.
    """
    s = store()
    df = s.table("external_groups")
    usable = df[df["passes_min_n"]] if "passes_min_n" in df.columns else df
    return {
        "min_group_n": S.MIN_GROUP_N,
        "n_offered": int(usable["group"].nunique()),
        "n_filtered_out": int(df["group"].nunique() - usable["group"].nunique()),
        "positive_control": "gene_desert_bottomQ_density",
        "note": "gene_desert_bottomQ_density is the positive control: defined FROM "
                "gene density, it is the one group that collapses under density "
                "stratification while surviving insulation. It is what makes the "
                "other rows meaningful.",
        "rows": _clean(usable.to_dict("records")),
    }
