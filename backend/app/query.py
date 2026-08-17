"""The query DSL, and the deterministic filter that executes it.

This module is the retrieval half of the natural-language gene finder, and it
knows nothing about language models. A query is a plain object; executing it is
plain pandas. The same query always returns the same genes.

That separation is the whole design. The model's only job is to turn a sentence
into one of these objects; if it misreads the sentence, the result is a
*visibly wrong query* the user can see and correct, not a plausible gene list
with no provenance. Everything here works with no API key, the UI can build
and run queries directly.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from . import store_schema as S


class QueryError(ValueError):
    """A query that cannot be executed, phrased for the user."""


# ---------------------------------------------------------------------------
# schema
# ---------------------------------------------------------------------------

# Handed to the model as the structured-output schema, and used here to validate
# whatever comes back. Deliberately small: every filter maps to one unambiguous
# operation over the store, so there is nothing for a translation to be vague
# about.
def query_schema(axes: list[str], cohorts: list[str], groups: list[str],
                 features: list[str]) -> dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["filters"],
        "properties": {
            "filters": {
                "type": "array",
                "description": "Conditions combined with AND.",
                "items": {
                    "anyOf": [
                        {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["type", "axis", "direction"],
                            "properties": {
                                "type": {"const": "axis"},
                                "axis": {"enum": axes},
                                "direction": {"enum": ["high", "low"]},
                                "percentile": {
                                    "type": "number",
                                    "description": "Cutoff, default 25. 'high' keeps "
                                                   "genes above the (100-p)th "
                                                   "percentile; 'low' keeps genes "
                                                   "below the pth.",
                                },
                            },
                        },
                        {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["type", "cohort"],
                            "properties": {
                                "type": {"const": "cohort"},
                                "cohort": {"enum": cohorts},
                            },
                        },
                        {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["type", "group"],
                            "properties": {
                                "type": {"const": "archetype"},
                                "group": {"enum": groups},
                            },
                        },
                        {
                            "type": "object",
                            "additionalProperties": False,
                            "required": ["type", "feature", "direction"],
                            "properties": {
                                "type": {"const": "feature"},
                                "feature": {"enum": features},
                                "direction": {"enum": ["high", "low"]},
                                "percentile": {"type": "number"},
                            },
                        },
                    ]
                },
            },
            "sort": {
                "type": "object",
                "additionalProperties": False,
                "required": ["axis", "direction"],
                "properties": {
                    "axis": {"enum": axes},
                    "direction": {"enum": ["high", "low"]},
                },
            },
            "interpretation": {
                "type": "string",
                "description": "One plain sentence restating the query, for the user "
                               "to check. Name the direction explicitly, e.g. 'genes "
                               "in the bottom 25% of PC2 (long-range)'.",
            },
        },
    }


# ---------------------------------------------------------------------------
# execution
# ---------------------------------------------------------------------------


def describe(f: dict, poles: dict[str, dict[str, str]]) -> str:
    """Render one filter as the sentence the user checks the translation against."""
    t = f.get("type")
    if t == "axis":
        ax, d = f["axis"], f["direction"]
        p = f.get("percentile", 25)
        pole = (poles.get(ax) or {}).get("pos" if d == "high" else "neg")
        end = f", {pole}" if pole else ""
        side = f"top {p:g}%" if d == "high" else f"bottom {p:g}%"
        return f"{ax.upper()} in the {side}{end}"
    if t == "cohort":
        return f"member of {f['cohort']}"
    if t == "archetype":
        disp = S.ARCHETYPE_DISPLAY.get(f["group"], {}).get("display", f["group"])
        return f"in the {disp} region"
    if t == "feature":
        d = f["direction"]
        p = f.get("percentile", 25)
        side = f"top {p:g}%" if d == "high" else f"bottom {p:g}%"
        return f"{f['feature']} in the {side}"
    raise QueryError(f"unknown filter type {t!r}")


def run(query: dict, store) -> dict[str, Any]:
    """Execute a query against the store. Pure pandas, no model involved."""
    genes = store.genes
    emb = store.table("embeddings") if store.has("embeddings") else None
    feats = store.table("features") if store.has("features") else None
    cm = store.table("cohort_membership") if store.has("cohort_membership") else None

    keep = pd.Series(True, index=genes.index)
    steps: list[dict] = []

    for f in query.get("filters", []):
        t = f.get("type")
        before = int(keep.sum())

        if t == "axis":
            if emb is None:
                raise QueryError("axis filters need the embeddings table")
            ax = f["axis"]
            if ax not in emb.columns:
                raise QueryError(f"unknown axis {ax!r}")
            vals = genes["gene_id"].map(emb.set_index("gene_id")[ax])
            pct = vals.rank(pct=True) * 100
            p = float(f.get("percentile", 25))
            keep &= (pct >= 100 - p) if f["direction"] == "high" else (pct <= p)

        elif t == "cohort":
            if cm is None:
                raise QueryError("cohort filters need the cohort_membership table")
            members = set(cm.loc[cm["group"] == f["cohort"], "symbol_key"])
            if not members:
                raise QueryError(f"no panel genes in cohort {f['cohort']!r}")
            keep &= genes["symbol_key"].isin(members)

        elif t == "archetype":
            keep &= genes["group"] == f["group"]

        elif t == "feature":
            if feats is None:
                raise QueryError("feature filters need the features table")
            sub = feats[feats["feature"] == f["feature"]]
            if sub.empty:
                raise QueryError(f"unknown feature {f['feature']!r}")
            pct = genes["symbol_key"].map(sub.set_index("symbol_key")["percentile"])
            p = float(f.get("percentile", 25))
            keep &= (pct >= 100 - p) if f["direction"] == "high" else (pct <= p)

        else:
            raise QueryError(f"unknown filter type {t!r}")

        steps.append({
            "filter": f,
            "reads_as": describe(f, S.PC_POLES),
            "before": before,
            "after": int(keep.sum()),
        })

    hits = genes[keep].copy()

    # ---- ranking ----------------------------------------------------------
    # "Genes with X" is a question about degree, not membership, and a
    # percentile cut answers it badly: one cut at the top 25% can only ever
    # return a quarter of the panel however specific the question was. Ranking
    # on a composite of the features that express X returns the MOST X, so
    # selectivity comes from the question rather than from the threshold.
    #
    # Features are already z-scored, so a signed mean is a well-posed composite.
    # Weights are honoured but default to 1: an unweighted mean is the honest
    # default when nothing justifies saying one component matters more.
    rank = query.get("rank") or []
    contributions: list[dict] = []
    if rank and feats is not None:
        wide = feats.pivot(index="gene_id", columns="feature", values="z")
        comp = pd.Series(0.0, index=hits["gene_id"])
        used = 0.0
        for r in rank:
            name = r.get("feature")
            if name not in wide.columns:
                raise QueryError(f"unknown feature {name!r}")
            sign = -1.0 if r.get("direction") == "low" else 1.0
            w = float(r.get("weight", 1.0))
            vals = hits["gene_id"].map(wide[name]).fillna(0.0).values
            comp += sign * w * vals
            used += abs(w)
            contributions.append({
                "feature": name,
                "direction": r.get("direction", "high"),
                "weight": w,
            })
        if used:
            comp /= used
        hits["score"] = comp.values

        # Per-feature percentiles for every ranked gene, and a flag for any
        # criterion the gene FAILS.
        #
        # A weighted sum lets a gene compensate: asked for "short-range,
        # promoter-dense, no CTCF", ZNF226 ranked third on 98th-percentile
        # promoter peaks and 100th-percentile promoter signal while sitting at
        # the 6th percentile for short-range. That is a defensible hit and a
        # hidden trade-off, and hiding it is the problem. The scorecard makes
        # each gene's weakest criterion visible instead of averaging it away.
        pct = feats.pivot(index="gene_id", columns="feature", values="percentile")
        cards = []
        for gid in hits["gene_id"]:
            row = {}
            for c in contributions:
                v = pct[c["feature"]].get(gid)
                if v is None or pd.isna(v):
                    continue
                # "low" asked for => a LOW percentile is the good outcome
                good = (100 - v) if c["direction"] == "low" else v
                row[c["feature"]] = {
                    "percentile": round(float(v), 1),
                    "meets": round(float(good), 1),
                }
            weak = sorted(row.items(), key=lambda kv: kv[1]["meets"])
            cards.append({
                "criteria": row,
                "weakest": weak[0][0] if weak else None,
                "weakest_meets": weak[0][1]["meets"] if weak else None,
                "fails_any": bool(weak and weak[0][1]["meets"] < 25),
            })
        hits["scorecard"] = cards
        hits = hits.sort_values("score", ascending=False)
        steps.append({
            "reads_as": "ranked by " + ", ".join(
                f"{c['feature']} {c['direction']}" for c in contributions
            ),
            "before": int(len(hits)),
            "after": int(len(hits)),
        })

    sort = query.get("sort")
    if not rank and sort and emb is not None and sort.get("axis") in emb.columns:
        hits["_sort"] = hits["gene_id"].map(emb.set_index("gene_id")[sort["axis"]])
        hits = hits.sort_values("_sort", ascending=sort["direction"] == "low")
        hits = hits.drop(columns="_sort")

    cols = [c for c in ("gene_id", "gene_symbol", "group", "max_posterior",
                        "viewpoint_chrom", "viewpoint_pos", "score", "scorecard")
            if c in hits.columns]

    limit = int(query.get("limit", 25 if rank else 50))
    return {
        "n_matched": int(len(hits)),
        "ranked": bool(rank),
        "rank_by": contributions or None,
        "genes": hits[cols].head(limit).to_dict("records"),
        "steps": steps,
        # Small results are the expected outcome of a specific question, not an
        # error, but they are also where over-reading is easiest, so say so.
        # Only for filters: a ranking always returns its top N, so a short list
        # there means nothing.
        "note": (
            "Fewer than 5 genes match. That is a specific question, not necessarily "
            "a meaningful group, with 1,846 genes, narrow conjunctions land on "
            "handfuls by chance."
            if not rank and len(hits) < 5 else None
        ),
    }
