#!/usr/bin/env python3
"""Re-check the store's invariants after a build.

Every silent-data bug in this project so far, the wrong anchor table, the
all-NaN posterior merge, the case-mismatched `Cxorfnn` joins, the label subset
that covered 39% of the panel, produced a store that loaded fine and served
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
    "n_features": 85,      # 91 -> 73 (S2.1/S2.1b) -> 85 (NaN fix), 2026-08-16
    "label_counts": {
        "arch-HK": 844,
        "arch-ME-constitutive": 369,
        "arch-ME-effector": 351,
        "arch-sparse": 261,
        "arch-off": 21,
    },
    "n_core": 1040,                  # 56.3%, vs FINDINGS.md's "56% core"
    "repro_genes": 116,              # genes captured in both panels
    # RE-BASELINED 2026-08-17, and the reason matters because CLAUDE.md warns
    # against re-baselining silently.
    #
    # These two checks used to read the `reproducibility` table, which is a copy
    # of audit/continuous_methods/cross_panel_reproducibility.tsv, dated
    # 2026-07-31 and never regenerated. 24 of its 63 features no longer exist in
    # the substrate, including every oe_asymmetry_* and oe_tailedness_*, which
    # were deleted on 2026-08-16 BECAUSE they failed a reproducibility threshold.
    #
    # So both checks were passing for a reason unrelated to the thing being
    # checked: the file had 63 rows and a median of 0.752 because nobody had
    # touched it. Same failure shape as the CORS bug and the vacuous
    # `tsc --noEmit`, and a third instance of "if a check has never once failed,
    # verify the verifier".
    #
    # Now computed from `reproducibility_pairs`, which is rebuilt from the immune
    # pickle on every store build and whose 41 features are all live. The median
    # RISES, because the features removed for being unreproducible were dragging
    # it down:
    #     stale, all 63                median rho 0.752   <- the old expectation
    #     the 24 retired features      median rho 0.654
    #     stale, 39 survivors          median rho 0.833
    #     current pairs, 41 features   median rho 0.819   33/41 above 0.7
    "repro_features": 41,            # live shared features, from the pairs table
    "repro_median_rho": 0.819,       # recomputed 2026-08-17, was 0.752 stale
    # Re-baselined 2026-08-16 after removing 4 degenerate topology features and
    # 16 unreproducible per-peak moments, and protecting n_active_peaks.
    # NOTE: these no longer match audit/continuous_methods/dimension_names.tsv,
    # which was computed on the 91-feature substrate and is now STALE. The check
    # below compares against these constants, not against that file. Regenerate
    # dimension_names.tsv upstream before quoting it anywhere.
    "pca_variance_top5": [15.42, 9.66, 8.81, 5.56, 5.36],
    "n_real_dimensions": 18,         # 19 -> 14 -> 18; structure_vs_noise.tsv stale
    "usable_cohorts": 21,            # of 23 reference sets, >= 25 panel genes
    # Amount-corrected space, added 2026-08-16 and now the display default.
    # Magnitude is projected out before PCA, so these components differ from the
    # raw PCs entirely: sPC3 is NOT PC3 and the numbering does not correspond.
    "shape_variance_top5": [10.00, 8.35, 7.34, 5.59, 5.25],
    # Renamed 2026-08-17: mid/far -> contained/extended. The reach axis is a
    # RELATIVE position, and the absolute-band names matched a gene's dominant
    # band only 20% of the time.
    # RE-FIT 2026-08-17 excluding 9 under-evidenced genes (<3 called peaks:
    # ASB9, CD209, CEP15, DEPDC1, NCALD, PADI4, SACS, TLR4, TTPAL). They carried
    # the three highest top weights in the panel while having zero to two peaks,
    # because the membership softmax squares distance and a gene with no evidence
    # sits far from every centroid. Excluded from the FIT and given region = NA.
    #
    # Counts moved by more than the nine removed: reach went 1074/772 to
    # 1119/718, so 36 other genes crossed the reach boundary. That is consistent
    # with the measured substrate sensitivity (section 6g) and is a reminder the
    # boundary is soft, not that something went wrong.
    #
    # Previous, for provenance: contained-promoter 555, contained-enhancer 519,
    # extended-enhancer 281, extended-ctcf 253, extended-promoter 238.
    "taxonomy_regions": {
        "contained-promoter": 586, "contained-enhancer": 533,
        "extended-promoter": 252, "extended-enhancer": 236,
        "extended-ctcf": 230,
    },
    "n_under_evidenced": 9,          # region = NA by design, not a join failure
    "shape_basis_n": 8,              # MAG_OVERALL lists 11; 2 degenerate, and
                                     # raw_peak_max_max_all now loses the prune
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
    # The store's own join guards raise before verify can inspect anything,
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

    # -- amount-corrected space, the DISPLAY DEFAULT since 2026-08-16 --------
    # The load-bearing check is the last one. The whole reason this space is the
    # display is that no component carries amount; if that stops being true the
    # app is back to plotting "this gene has more signal" as if it were shape,
    # which is the failure the switch was made to prevent.
    if s.has("shape_scree"):
        sc = s.table("shape_scree")
        got = [round(float(v), 2) for v in sc["variance_pct"].head(5)]
        check("shape PCA variance",
              all(approx(a, b, 0.02) for a, b in zip(got, EXPECTED["shape_variance_top5"])),
              f"{got} vs {EXPECTED['shape_variance_top5']}")

    if s.has("shape_basis"):
        nb = len(s.table("shape_basis"))
        check("magnitude basis resolves to 9 of MAG_OVERALL's 11",
              nb == EXPECTED["shape_basis_n"],
              f"{nb} (mean_degree and mean_degree_raw no longer exist)")

    if s.has("shape_embeddings") and s.has("features"):
        se = s.table("shape_embeddings")
        check("shape embeddings cover the panel", len(se) == EXPECTED["n_genes"],
              f"{len(se)} rows")
        W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
        basis = list(s.table("shape_basis")["feature"]) if s.has("shape_basis") else []
        basis = [b for b in basis if b in W.columns]
        if basis:
            M = W[basis].to_numpy(float)
            amount = (M - M.mean(0)).mean(1)
            sub = se.set_index("gene_id").reindex(W.index)
            worst = max(abs(float(np.corrcoef(sub[f"spc{i}"], amount)[0, 1]))
                        for i in range(1, 11) if f"spc{i}" in sub.columns)
            check("no shape component carries amount", worst < 0.05,
                  f"max |r| over the first 10 = {worst:.4f}")

    # -- taxonomy ------------------------------------------------------------
    # The load-bearing check is the last one. If any gene ever exceeds a top
    # weight of 0.5 for most of the panel, the regions have stopped being areas
    # of a continuum and the whole framing needs revisiting.
    if s.has("taxonomy"):
        tx = s.table("taxonomy")
        check("taxonomy covers every gene", len(tx) == EXPECTED["n_genes"],
              f"{len(tx)} rows")
        counts = tx["region"].value_counts().to_dict()
        check("taxonomy region counts",
              counts == EXPECTED["taxonomy_regions"],
              f"{counts}")
        # The under-evidenced rows carry NaN weights by design, so they are
        # excluded from the weight checks rather than allowed to fail them.
        ue = tx["under_evidenced"].fillna(False).astype(bool) \
            if "under_evidenced" in tx.columns else pd.Series(False, index=tx.index)
        check("under-evidenced genes are marked, not missing",
              int(ue.sum()) == EXPECTED["n_under_evidenced"],
              f"{int(ue.sum())} marked (expected {EXPECTED['n_under_evidenced']})")
        check("under-evidenced genes carry no region",
              int(tx.loc[ue, "region"].notna().sum()) == 0,
              f"{int(tx.loc[ue, 'region'].notna().sum())} wrongly labelled")

        fit = tx[~ue]
        wcols = [c for c in tx.columns if c.startswith("w_")]
        sums = fit[wcols].sum(axis=1)
        check("memberships sum to 1", bool(((sums - 1).abs() < 1e-3).all()),
              f"min {sums.min():.4f}, max {sums.max():.4f}")
        pct = float((fit["top_weight"] < 0.5).mean() * 100)
        check("regions are areas of a continuum, not clusters", pct > 95,
              f"{pct:.2f}% of genes have top weight < 0.5 "
              f"({int((fit['top_weight'] >= 0.5).sum())} exceptions)")

    if "region" in genes.columns:
        # Exactly the 9 under-evidenced genes may lack a region. Zero would mean
        # the exclusion silently failed; more would mean a join dropped rows, the
        # failure mode that lost 11 Cxorfnn genes once before.
        n_na = int(genes["region"].isna().sum())
        check("only under-evidenced genes lack a region",
              n_na == EXPECTED["n_under_evidenced"],
              f"{n_na} without a region (expected exactly "
              f"{EXPECTED['n_under_evidenced']})")

    # -- reproducibility -----------------------------------------------------
    # Computed from `reproducibility_pairs`, NOT the stale `reproducibility`
    # copy. See the EXPECTED block above for why.
    if s.has("reproducibility_pairs"):
        pairs = s.table("reproducibility_pairs")
        per_feat = pairs.groupby("feature").apply(
            lambda g: pd.Series({
                "spearman_rho": float(np.corrcoef(g["gw"].rank(),
                                                  g["immune"].rank())[0, 1]),
                "n": int(len(g)),
            }),
            include_groups=False,
        )
        check("shared features (live)",
              len(per_feat) == EXPECTED["repro_features"],
              f"{len(per_feat)} (expected {EXPECTED['repro_features']})")
        check("twice-captured genes",
              int(per_feat["n"].max()) == EXPECTED["repro_genes"],
              f"{int(per_feat['n'].max())} (expected {EXPECTED['repro_genes']})")
        med = float(per_feat["spearman_rho"].median())
        check("median reproducibility rho",
              approx(med, EXPECTED["repro_median_rho"], 0.005),
              f"{med:.3f} (expected {EXPECTED['repro_median_rho']})")
        # The 2026-08-16 removal of the asymmetry/tailedness families was
        # justified on a rho >= 0.70 threshold, so nothing unreproducible should
        # survive. This asserts the removal actually worked rather than trusting
        # that it did.
        n_bad = int((per_feat["spearman_rho"] < 0.3).sum())
        check("no unreproducible feature survives", n_bad == 0,
              f"{n_bad} features below rho 0.3")

        # The stale copy is still in the store for provenance. Assert it is NOT
        # being mistaken for current: if someone regenerates the tsv, this fires
        # and the numbers above should be revisited together.
        if s.has("reproducibility"):
            stale = s.table("reproducibility")
            live = set(s.table("features")["feature"].unique())
            retired = sorted(set(stale["feature"]) - live)
            check("stale reproducibility copy is unchanged (provenance only)",
                  len(retired) == 24,
                  f"{len(retired)} of {len(stale)} features retired; "
                  f"regenerate the tsv or drop it from paths.py")

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

    print("\nmccprofiler-app, store verification\n")
    verify(s)

    width = max(len(n) for _, n, _ in _results) + 2
    for status, name, detail in _results:
        print(f"  {status}  {name:<{width}} {detail}")

    failed = [n for st, n, _ in _results if st == FAIL]
    print()
    if failed:
        print(f"  FAIL, {len(failed)} of {len(_results)} checks failed: "
              f"{', '.join(failed[:6])}")
    else:
        print(f"  PASS, {len(_results)} checks")
    print()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
