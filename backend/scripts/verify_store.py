#!/usr/bin/env python3
"""Re-check the store's invariants after a build.

Every silent-data bug in this project so far — the wrong anchor table, the
all-NaN posterior merge, the case-mismatched `Cxorfnn` joins, the label subset
that covered 39% of the panel — produced a store that loaded fine and served
wrong numbers. None raised. Each was caught only because someone noticed a
total was off by eleven, or a column was NaN.

This turns that noticing into an assertion. It reads the built store and checks
the things that were wrong before, plus the published numbers the store must
reproduce.

    conda activate mccapp   # or cd4env
    python backend/scripts/verify_store.py

Exit code 0 = all checks pass, 1 = at least one failed. Safe to run in a loop
after every build.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from app import paths as P  # noqa: E402
from app import store_schema as S  # noqa: E402
from app.store import Store, StoreMissing  # noqa: E402

# Numbers the store must reproduce, each traceable to a handoff or an audit
# output. A mismatch means either the store is wrong or the science moved --
# both are worth stopping for.
EXPECTED = {
    "n_genes": 1846,
    "n_features": 91,
    "label_counts": {
        "arch-HK": 844,
        "arch-ME-constitutive": 369,
        "arch-ME-effector": 351,
        "arch-sparse": 261,
        "arch-off": 21,
    },
    "n_core": 1040,                  # 56.3%, vs FINDINGS.md's "56% core"
    "repro_genes": 116,              # genes captured in both panels
    "repro_features": 63,            # shared features
    "repro_median_rho": 0.752,       # 2026-07-31 handoff
    "pca_variance_top5": [15.27, 9.47, 8.21, 6.33, 5.88],  # dimension_names.tsv
    "n_real_dimensions": 19,         # structure_vs_noise.tsv
    "usable_cohorts": 21,            # of 23 reference sets, >= 25 panel genes
}

PASS, FAIL, WARN = "ok  ", "FAIL", "warn"
_results: list[tuple[str, str, str]] = []


def check(name: str, ok: bool, detail: str = "", warn_only: bool = False) -> bool:
    _results.append((PASS if ok else (WARN if warn_only else FAIL), name, detail))
    return ok


def approx(a: float, b: float, tol: float = 0.01) -> bool:
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------


def verify(s: Store) -> None:
    # The store's own join guards raise before verify can inspect anything —
    # which is correct (fail at the earliest point), but a traceback is a poor
    # report. Catch it and record it as the failed check it is.
    try:
        genes = s.genes
    except Exception as e:  # noqa: BLE001
        check("store loads", False, f"{type(e).__name__}: {e}")
        return
    check("store loads", True)

    # -- shape ---------------------------------------------------------------
    check("gene count", len(genes) == EXPECTED["n_genes"],
          f"{len(genes)} (expected {EXPECTED['n_genes']})")

    # -- joins: the failure mode that produced NaN and dropped 11 genes -------
    for col in ("group", "max_posterior", "viewpoint_chrom", "viewpoint_pos"):
        if col not in genes.columns:
            check(f"column {col}", False, "missing")
            continue
        n_missing = int(genes[col].isna().sum())
        check(f"{col} complete", n_missing == 0, f"{n_missing} null")

    # The Cxorfnn bug: totals were 11 short and nothing raised.
    if "group" in genes.columns:
        total = int(genes["group"].value_counts().sum())
        check("labels cover every gene", total == len(genes),
              f"{total} labelled of {len(genes)}")

        counts = genes["group"].value_counts().to_dict()
        for grp, want in EXPECTED["label_counts"].items():
            got = int(counts.get(grp, 0))
            check(f"  {grp}", got == want, f"{got} (expected {want})")

    # -- posteriors: the all-NaN merge ---------------------------------------
    if "is_core" in genes.columns:
        n_core = int(genes["is_core"].sum())
        check("core gene count", n_core == EXPECTED["n_core"],
              f"{n_core} = {100 * n_core / len(genes):.1f}% "
              f"(expected {EXPECTED['n_core']}, ~56%)")

    # -- features ------------------------------------------------------------
    if s.has("features"):
        f = s.table("features")
        n_feat = f["feature"].nunique()
        check("feature count", n_feat == EXPECTED["n_features"],
              f"{n_feat} (expected {EXPECTED['n_features']})")
        per_gene = f.groupby("symbol_key").size()
        check("every gene has every feature",
              bool((per_gene == n_feat).all()),
              f"min {per_gene.min()}, max {per_gene.max()}")
        # Percentiles are per-feature ranks, so each feature must span 0-100.
        span = f.groupby("feature")["percentile"].agg(["min", "max"])
        check("percentiles are per-feature",
              bool((span["min"] < 5).all() and (span["max"] > 95).all()),
              "each feature spans the panel")

    # -- peaks: the unsigned-distance trap -----------------------------------
    if s.has("peaks"):
        pk = s.table("peaks")
        half = (S.N_BINS_L0 * S.BIN_SIZE_L0) // 2
        out = int((pk["offset_bp"].abs() > half).sum())
        check("all peaks inside the window", out == 0, f"{out} outside +/-1 Mb")
        check("peaks are signed both ways",
              bool((pk["offset_bp"] < 0).any() and (pk["offset_bp"] > 0).any()),
              f"{int((pk['offset_bp'] < 0).sum())} left, "
              f"{int((pk['offset_bp'] > 0).sum())} right")
        classes = set(pk["re"].unique())
        check("element classes", classes == {"enhancer", "ctcf", "promoter"},
              ", ".join(sorted(classes)))

    # -- embeddings and PCA --------------------------------------------------
    if s.has("embeddings"):
        emb = s.table("embeddings")
        check("embedding covers every gene", len(emb) == len(genes),
              f"{len(emb)} rows")
        nan_cols = [c for c in emb.columns if emb[c].isna().any()]
        check("no NaN coordinates", not nan_cols, ", ".join(nan_cols[:5]))

    if s.has("pc_scree"):
        sc = s.table("pc_scree")
        got = [round(float(v), 2) for v in sc["variance_pct"].head(5)]
        check("PCA reproduces dimension_names.tsv",
              all(approx(a, b, 0.02) for a, b in zip(got, EXPECTED["pca_variance_top5"])),
              f"{got} vs {EXPECTED['pca_variance_top5']}")
        n_real = int(sc["above_noise"].sum())
        check("real dimensions above noise",
              n_real == EXPECTED["n_real_dimensions"],
              f"{n_real} (expected {EXPECTED['n_real_dimensions']})")

    if s.has("pc_loadings") and s.has("features"):
        ld = s.table("pc_loadings")
        n_per_pc = ld.groupby("pc").size()
        check("loadings cover every feature",
              bool((n_per_pc == EXPECTED["n_features"]).all()),
              f"{n_per_pc.min()}-{n_per_pc.max()} per component")

    # -- reproducibility -----------------------------------------------------
    if s.has("reproducibility"):
        rp = s.table("reproducibility")
        check("shared features", len(rp) == EXPECTED["repro_features"],
              f"{len(rp)} (expected {EXPECTED['repro_features']})")
        if "n" in rp.columns:
            check("twice-captured genes", int(rp["n"].max()) == EXPECTED["repro_genes"],
                  f"{int(rp['n'].max())} (expected {EXPECTED['repro_genes']})")
        med = float(rp["spearman_rho"].median())
        check("median reproducibility rho",
              approx(med, EXPECTED["repro_median_rho"], 0.005),
              f"{med:.3f} (expected {EXPECTED['repro_median_rho']})")

    # -- cohorts -------------------------------------------------------------
    if s.has("cohort_membership"):
        cm = s.table("cohort_membership")
        counts = cm["group"].value_counts()
        usable = int((counts >= S.MIN_GROUP_N).sum())
        check("usable reference sets", usable == EXPECTED["usable_cohorts"],
              f"{usable} of {counts.size} clear the {S.MIN_GROUP_N}-gene floor")
        unknown = set(cm["symbol_key"]) - set(genes["symbol_key"])
        check("cohort members are panel genes", not unknown,
              f"{len(unknown)} unknown" if unknown else "")

    # -- profiles ------------------------------------------------------------
    try:
        chans = s.channels
        check("channel count", len(chans) == len(S.CHANNELS),
              f"{len(chans)}: {', '.join(chans)}")
        check("mcc channel present", S.PRIMARY_CHANNEL in chans)

        gid = str(genes.iloc[0]["gene_id"])
        for lvl in sorted(S.LEVELS):
            p = s.profile(gid, level=lvl)
            check(f"profile L{lvl}", p["n"] == S.level_bins(lvl),
                  f"{p['n']} bins @ {p['bp_per_bin']} bp")
        # O/E must differ from raw, or the expected curve was not applied.
        raw = np.array(s.profile(gid, mode="raw")["values"])
        oe = np.array(s.profile(gid, mode="oe")["values"])
        check("O/E differs from raw", not np.allclose(raw, oe),
              "expected curve applied")
    except Exception as e:  # noqa: BLE001
        check("profiles readable", False, str(e))

    # -- provenance ----------------------------------------------------------
    m = s.manifest
    check("manifest records the panel", m.get("panel") == P.PANEL_GW, str(m.get("panel")))
    check("manifest records a commit", bool(m.get("scripts_cleaned_commit")),
          str(m.get("scripts_cleaned_commit")))
    blacklisted = [k for k, v in m.get("inputs", {}).items()
                   if v.get("present") and Path(v["path"]).name in P.BLACKLIST]
    check("no blacklisted input was read", not blacklisted, ", ".join(blacklisted))


def main() -> int:
    try:
        s = Store(P.STORE)
    except StoreMissing as e:
        print(f"\n{e}\n")
        return 1

    print("\nmccprofiler-app — store verification\n")
    verify(s)

    width = max(len(n) for _, n, _ in _results) + 2
    for status, name, detail in _results:
        print(f"  {status}  {name:<{width}} {detail}")

    failed = [n for st, n, _ in _results if st == FAIL]
    print()
    if failed:
        print(f"  FAIL — {len(failed)} of {len(_results)} checks failed: "
              f"{', '.join(failed[:6])}")
    else:
        print(f"  PASS — {len(_results)} checks")
    print()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
