#!/usr/bin/env python3
"""P1, build the app store from the canonical pipeline artefacts.

Runs in `cd4env` (needs numpy, pandas, h5py, pyarrow). Reads the 4.4 GB clean
pickle once, writes a purpose-built store, and records provenance for every
input. The server then reads only the store.

    conda activate cd4env
    python backend/scripts/build_store.py                # profiles + tables
    python backend/scripts/build_store.py --tables-only  # skip the pickle
    python backend/scripts/build_store.py --levels 1,2   # skip native-resolution

Refuses to read blacklisted artefacts (see app/paths.py BLACKLIST) rather than
silently producing a store that disagrees with the current science.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from app import paths as P  # noqa: E402
from app import store_schema as S  # noqa: E402

# The builder is the ONE place allowed to import the science package, so the
# viewpoint convention cannot drift between pipeline and app. The server never
# does this. Mirrors orchid's thin-wrapper-import pattern.
sys.path.insert(0, str(P.PROFILER / "src"))


def log(msg: str) -> None:
    print(f"  {msg}", flush=True)


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------


def git_commit(repo: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def stamp(path: Path) -> dict:
    if path is None or not path.exists():
        return {"path": str(path), "present": False}
    st = path.stat()
    return {
        "path": str(path),
        "present": True,
        "size": st.st_size,
        "mtime": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds"),
    }


def assert_not_blacklisted(path: Path) -> None:
    name = Path(path).name
    if name in P.BLACKLIST:
        raise SystemExit(
            f"\nREFUSING to read blacklisted artefact: {path}\n  {P.BLACKLIST[name]}\n"
        )


# ---------------------------------------------------------------------------
# profiles
# ---------------------------------------------------------------------------


def downsample(row: np.ndarray, factor: int) -> np.ndarray:
    """Mean-pool a 1D profile by an integer factor."""
    if factor == 1:
        return row
    return row.reshape(-1, factor).mean(axis=1)


def build_profiles(levels: list[int], batch: int) -> dict:
    import h5py
    import pickle

    pkl = P.by_key("clean_matrix_gw").path
    assert_not_blacklisted(pkl)

    log(f"loading {pkl.name} ({pkl.stat().st_size / 1024**3:.1f} GB), this takes a few minutes")
    with open(pkl, "rb") as f:
        data = pickle.load(f)

    gene_ids = [str(g) for g in data["gene_ids"]]
    source_labels = [str(s) for s in data["source_labels"]]
    n_genes = len(gene_ids)
    log(f"{n_genes} genes, {len(S.CHANNELS)} channels")

    present = [k for k in S.CHANNELS if k in data]
    missing = [k for k in S.CHANNELS if k not in data]
    if missing:
        log(f"NOTE: channels absent from pickle, skipped: {', '.join(missing)}")

    P.STORE.mkdir(parents=True, exist_ok=True)
    h5_path = P.STORE / "profiles.h5"

    with h5py.File(h5_path, "w") as h5:
        h5.create_dataset(
            S.H5_GENE_IDS,
            data=np.array(gene_ids, dtype=h5py.string_dtype("utf-8")),
        )

        for key in present:
            chan = S.CHANNELS[key]
            mat = data[key]
            if mat.shape[1] != S.N_BINS_L0:
                raise SystemExit(f"{key}: expected {S.N_BINS_L0} bins, got {mat.shape[1]}")

            # Panel-mean expected curve, for server-side O/E. Computed at native
            # resolution on float64 then stored float32, matching expected.py.
            expected = np.asarray(mat, dtype=np.float64).mean(axis=0)
            h5.create_dataset(
                f"{S.H5_EXPECTED}/{chan}",
                data=expected.astype(np.float32),
                compression=S.COMPRESSION,
                compression_opts=S.COMPRESSION_OPTS,
            )

            for lvl in levels:
                factor = S.LEVELS[lvl]
                nbins = S.level_bins(lvl)
                dset = h5.create_dataset(
                    f"{S.H5_PROFILES}/{chan}/L{lvl}",
                    shape=(n_genes, nbins),
                    dtype=np.float32,
                    chunks=(1, nbins),          # one gene per chunk
                    compression=S.COMPRESSION,
                    compression_opts=S.COMPRESSION_OPTS,
                    shuffle=True,
                )
                for start in range(0, n_genes, batch):
                    stop = min(start + batch, n_genes)
                    block = np.asarray(mat[start:stop], dtype=np.float32)
                    if factor > 1:
                        block = block.reshape(stop - start, nbins, factor).mean(axis=2)
                    dset[start:stop] = block
                log(f"{chan} L{lvl}: {n_genes} x {nbins} @ {S.level_bp(lvl)} bp")

            del mat
            data[key] = None  # release as we go; the pickle holds ~4.7 GB

    size = h5_path.stat().st_size / 1024**2
    log(f"profiles.h5 written, {size:.0f} MB")

    return {
        "gene_ids": gene_ids,
        "source_labels": source_labels,
        "n_genes": n_genes,
        "channels": [S.CHANNELS[k] for k in present],
        "levels": levels,
        "h5_mb": round(size, 1),
    }


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------


def read_table(path: Path, **kw) -> pd.DataFrame:
    assert_not_blacklisted(path)
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, **kw)
    return add_symbol_key(df)


def add_symbol_key(df: pd.DataFrame) -> pd.DataFrame:
    """Add an uppercase join key alongside gene_symbol.

    The pipeline's cluster outputs uppercase gene symbols while the signal
    pickle preserves original casing, so the 11 `Cxorfnn` genes, the only
    symbols with interior lowercase, silently failed to join. Every
    symbol-keyed merge in the app uses `symbol_key`, never `gene_symbol`.
    """
    if "gene_symbol" in df.columns and "symbol_key" not in df.columns:
        df = df.copy()
        df["symbol_key"] = df["gene_symbol"].astype(str).str.upper()
    return df


def build_shape_embeddings(fm: pd.DataFrame, fnames: list[str], seed: int = 0):
    """Amount-corrected coordinates. THIS IS THE DISPLAY DEFAULT.

    Why this exists, added 2026-08-16. The app was plotting raw PCs while every
    Aim 3 result was computed on the amount-corrected substrate, so the report
    carried two coordinate systems and readers had to know which was which. That
    is the same confusion behind the earlier sPC/PC mix-ups.

    It also fixes a naming problem rather than papering over it. On the
    73-feature substrate the raw PC1 correlates 0.505 with the magnitude basis
    and PC3 correlates 0.519, so "amount" is smeared across two components and
    neither can be honestly named. The app's default PC2xPC3 plane had been
    chosen *because* PC1 was the amount axis, and that stopped being true.

    Here magnitude is projected out first, so every component correlates 0.000
    with it by construction and the axes are named for what they measure:

        sPC1  reach            max_distance_to_viewpoint, mean_peak_gap
        sPC2  enhancer vs promoter
        sPC3  CTCF vs promoter
        sPC4  spread vs bait-piled
        sPC5  mid-range vs both extremes

    The magnitude basis is MAG_OVERALL, imported rather than redefined so there
    is one definition. It lists 11 features and resolves to 9 here, because
    mean_degree and mean_degree_raw were dropped as degenerate on 2026-08-16;
    the resolved basis is written to the `shape_basis` table for provenance.

    Raw PCs remain available in `embeddings` / `pc_loadings` / `pc_scree` for
    provenance and for the toggle, and `verify_store.py` still checks them.
    """
    import sys as _sys
    _sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                        "cd4_cleaned/scripts_cleaned/audit/GW/scripts")
    from _shape import MAG_OVERALL
    from sklearn.decomposition import PCA

    mag = [m for m in MAG_OVERALL if m in fnames]
    keep = [f for f in fnames if f not in mag]
    log(f"shape space: {len(keep)} features after removing a {len(mag)}-feature "
        f"magnitude basis (MAG_OVERALL lists {len(MAG_OVERALL)})")

    A = np.column_stack([np.ones(len(fm)), fm[mag].to_numpy(dtype=np.float64)])
    X = fm[keep].to_numpy(dtype=np.float64)
    beta, *_ = np.linalg.lstsq(A, X, rcond=None)
    X = X - A @ beta
    X = (X - X.mean(0)) / X.std(0)

    pca = PCA(n_components=min(30, X.shape[1]), random_state=seed).fit(X)
    scores = pca.transform(X)
    var = pca.explained_variance_ratio_ * 100
    log("shape variance %: " + ", ".join(f"sPC{i+1} {v:.2f}" for i, v in enumerate(var[:5])))

    # Assert the correction actually worked. If any component still carries
    # amount, the display would repeat the problem it was built to remove.
    amount = fm[mag].to_numpy(dtype=np.float64)
    amount = (amount - amount.mean(0)).mean(1)
    worst = max(abs(float(np.corrcoef(scores[:, i], amount)[0, 1]))
                for i in range(min(10, scores.shape[1])))
    log(f"max |r| between any of the first 10 shape components and magnitude: {worst:.4f}")
    if worst > 0.05:
        raise RuntimeError(f"shape space still carries magnitude (|r| = {worst:.3f})")

    rng_null = np.random.default_rng(seed)
    Xn = X.copy()
    for j in range(Xn.shape[1]):
        rng_null.shuffle(Xn[:, j])
    null_var = PCA(n_components=pca.n_components_, random_state=seed) \
        .fit(Xn).explained_variance_ratio_ * 100
    log(f"shape parallel analysis: {int((var > null_var).sum())} components above noise")

    emb = pd.DataFrame({"gene_id": fm["gene_id"].astype(str).tolist()})
    for i in range(scores.shape[1]):
        emb[f"spc{i + 1}"] = scores[:, i]
    load = pd.DataFrame([{"pc": i + 1, "feature": f, "loading": float(v)}
                         for i in range(pca.components_.shape[0])
                         for f, v in zip(keep, pca.components_[i])])
    scree = pd.DataFrame({"pc": np.arange(1, len(var) + 1), "variance_pct": var,
                          "cumulative_pct": np.cumsum(var), "noise_pct": null_var,
                          "above_noise": var > null_var})
    basis = pd.DataFrame({"feature": mag, "role": "magnitude_basis"})
    return emb, load, scree, basis


def build_embeddings(fm: pd.DataFrame, fnames: list[str], seed: int = 0
                     ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-gene 2D coordinates: PCA, gcPCA, UMAP, and a UMAP null.

    PCA is the default in the UI because it is linear and preserves distances,
    so a continuum renders as a continuum.

    UMAP is included because people expect it, but it optimises a local-neighbour
    objective and is well documented to produce apparent clusters from continuous
    data. Rather than caveat that in prose, the same UMAP is fitted to a
    per-feature-permuted copy of the matrix, every marginal preserved exactly,
    all joint structure destroyed, so there is provably nothing to find. If the
    null looks similarly clustered, the artefact is self-evident. This mirrors
    the null used by audit/Archetype_Tests/12_structure_evidence.py.
    """
    from sklearn.decomposition import PCA

    X = fm[fnames].to_numpy(dtype=np.float64)
    gids = fm["gene_id"].astype(str).tolist()

    pca = PCA(n_components=min(30, X.shape[1]), random_state=seed).fit(X)
    scores = pca.transform(X)
    var = pca.explained_variance_ratio_ * 100
    log("PCA variance %: " + ", ".join(f"PC{i+1} {v:.2f}" for i, v in enumerate(var[:5])))

    # Loadings and a scree, so the axis names can be checked rather than trusted.
    # A named axis is an interpretation of its loadings; showing the name without
    # the loadings asks the reader to take it on faith.
    load_rows = []
    for i in range(pca.components_.shape[0]):
        for j, fname in enumerate(fnames):
            load_rows.append({"pc": i + 1, "feature": fname,
                              "loading": float(pca.components_[i, j])})
    loadings = pd.DataFrame(load_rows)

    # Parallel-analysis noise ceiling: PCA on a per-feature-permuted copy gives
    # the variance a component of this size explains when there is nothing to
    # find. Components above it are the "real dimensions" the 08-03 session
    # counted (19 on this substrate).
    rng_null = np.random.default_rng(seed)
    Xnull = X.copy()
    for j in range(Xnull.shape[1]):
        rng_null.shuffle(Xnull[:, j])
    null_var = PCA(n_components=pca.n_components_, random_state=seed) \
        .fit(Xnull).explained_variance_ratio_ * 100
    n_real = int((var > null_var).sum())
    log(f"parallel analysis: {n_real} components above the permuted-noise ceiling")

    scree = pd.DataFrame({
        "pc": np.arange(1, len(var) + 1),
        "variance_pct": var,
        "cumulative_pct": np.cumsum(var),
        "noise_pct": null_var,
        "above_noise": var > null_var,
    })

    out = pd.DataFrame({"gene_id": gids})
    for i in range(scores.shape[1]):
        out[f"pc{i + 1}"] = scores[:, i]

    # gcPCA density-free coordinates, if the audit produced them.
    try:
        g = np.load(P.by_key("gcpca_axes").path, allow_pickle=True)
        gc_ids = [str(s) for s in g["gene_ids"]]
        order = {gid: i for i, gid in enumerate(gc_ids)}
        idx = [order.get(gid) for gid in gids]
        if all(i is not None for i in idx):
            for key, tag in (("B1_density_explained", "b1"),
                             ("B2_marginal_within_density", "b2")):
                arr = g[key][idx]
                out[f"gcpca_{tag}_1"] = arr[:, 0]
                out[f"gcpca_{tag}_2"] = arr[:, 1]
            log("gcPCA density-free axes attached (B1, B2)")
        else:
            log("WARNING: gcPCA gene order does not cover the panel; axes skipped")
    except Exception as e:  # noqa: BLE001 - optional input, never fatal
        log(f"gcPCA axes unavailable ({e}); skipped")

    import umap

    def fit_umap(mat: np.ndarray, tag: str) -> None:
        emb = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                        random_state=seed).fit_transform(mat)
        out[f"{tag}_1"] = emb[:, 0]
        out[f"{tag}_2"] = emb[:, 1]

    fit_umap(X, "umap")
    log("UMAP fitted on the real matrix")

    rng = np.random.default_rng(seed)
    Xn = X.copy()
    for j in range(Xn.shape[1]):
        rng.shuffle(Xn[:, j])          # permute each feature independently
    fit_umap(Xn, "umap_null")
    log("UMAP fitted on the per-feature-permuted null (marginals kept, joint destroyed)")

    out["symbol_key"] = out["gene_id"].str.upper()
    return out, loadings, scree


def build_repro_pairs(fm_gw: pd.DataFrame, fnames_gw: list[str]) -> pd.DataFrame:
    """The 116 genes captured in BOTH panels, feature by feature, in each.

    This is the only direct measurement of technical reproducibility available:
    the same gene, two independent captures, the same pipeline. It is what makes
    kappa = 0.72 and the 0.511 noise floor concrete rather than quoted.

    LAB PAGE ONLY. The immune side is the stale 791 baseline missing PDCD1 and
    15 others, so this is defensible as a technical demonstration and NOT as a
    biological claim. The overlap may shift after the 791 re-run.
    """
    import pickle as _pickle

    im_path = P.PROFILER / "outputs/features/features_zscored.pkl"
    if not im_path.exists():
        log("immune features absent; reproducibility pairs skipped")
        return pd.DataFrame(columns=["symbol_key", "feature", "gw", "immune"])

    with open(im_path, "rb") as f:
        im = _pickle.load(f)
    fm_im, fnames_im = im["feature_matrix"], list(im["feature_names"])

    shared = [f for f in fnames_gw if f in set(fnames_im)]
    g = fm_gw.copy()
    g["symbol_key"] = g["gene_id"].astype(str).str.upper()
    i = fm_im.copy()
    i["symbol_key"] = i["gene_id"].astype(str).str.upper()

    both = sorted(set(g["symbol_key"]) & set(i["symbol_key"]))
    log(f"reproducibility: {len(both)} genes in both panels, {len(shared)} shared features")

    gi = g.set_index("symbol_key").loc[both, shared]
    ii = i.set_index("symbol_key").loc[both, shared]

    rows = []
    for feat in shared:
        for sym in both:
            rows.append({"symbol_key": sym, "feature": feat,
                         "gw": float(gi.at[sym, feat]), "immune": float(ii.at[sym, feat])})
    return pd.DataFrame(rows)


def build_cohorts(fm: pd.DataFrame) -> pd.DataFrame:
    """Per-gene membership of the externally-defined reference sets.

    Reuses `benchmarking.loaders.build_reference_sets` verbatim rather than
    re-deriving the sets, so the app's cohorts are the same objects the audit
    scripts test. Membership is restricted to panel genes here, the sets are
    genome-wide, and a group with 5,391 members genome-wide may have far fewer
    in a 1,846-gene panel, which is the number that governs whether it is worth
    displaying at all.
    """
    sys.path.insert(0, str(P.SCRIPTS / "benchmarking/src"))
    from benchmarking.loaders import (build_reference_sets, load_anchor,
                                      load_bio_labels, load_dice, load_lambert)

    anchor = load_anchor()
    refs = build_reference_sets(anchor, load_lambert()[0], load_dice(), load_bio_labels())

    syms = [str(g).upper() for g in fm["gene_id"]]
    rows = []
    for name, members in sorted(refs.items()):
        hit = [s for s in syms if s in members]
        rows.extend({"group": name, "symbol_key": s} for s in hit)

    df = pd.DataFrame(rows)
    counts = df["group"].value_counts()
    usable = int((counts >= S.MIN_GROUP_N).sum())
    log(f"cohorts: {len(refs)} reference sets, {usable} with >= {S.MIN_GROUP_N} "
        f"panel genes ({len(refs) - usable} too small to offer)")
    return df


def build_genes(gene_ids: list[str], source_labels: list[str]) -> pd.DataFrame:
    """Gene index with viewpoint positions loaded from the BED, never parsed
    from viewpoint_id (PLAN.md trap #3)."""
    from mccprofiler.io import load_viewpoints

    vps = load_viewpoints(gene_ids, source_labels, bed_dir=P.VIEWPOINT_DIR)
    log(f"viewpoints resolved for {len(vps)}/{len(gene_ids)} genes from the BEDs")

    rows = []
    for gid, src in zip(gene_ids, source_labels):
        suffix = f"_{src}"
        symbol = gid[: -len(suffix)] if gid.endswith(suffix) else gid
        chrom, pos = vps.get(gid, (None, None))
        rows.append(
            {
                "gene_id": gid,
                "gene_symbol": symbol,
                "source": src,
                "viewpoint_chrom": chrom,
                "viewpoint_pos": pos,
            }
        )
    return add_symbol_key(pd.DataFrame(rows))


def build_tables(gene_ids: list[str] | None, source_labels: list[str] | None) -> dict:
    out_dir = P.STORE / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, dict] = {}

    def write(name: str, df: pd.DataFrame, source: Path | None = None) -> None:
        path = out_dir / f"{name}.parquet"
        df.to_parquet(path, index=False)
        written[name] = {
            "rows": len(df),
            "cols": len(df.columns),
            "source": stamp(source) if source else None,
        }
        log(f"{name}: {len(df)} x {len(df.columns)}")

    # gene index -------------------------------------------------------------
    # The taxonomy region is joined on here so every view can colour by it.
    # `group` is KEPT unchanged for provenance and for joining back to the
    # 2026-07-21 labels the audit outputs reference; it is no longer the display
    # taxonomy, having been fit on the 91-feature substrate and named partly by
    # Eisenberg overlap.
    if gene_ids is not None:
        g = build_genes(gene_ids, source_labels)
        tax_p = P.by_key("taxonomy_labels").path
        if tax_p.exists():
            tx = read_table(tax_p)[["gene_id", "region", "reach", "composition",
                                    "top_weight", "mixedness"]]
            g["gene_id"] = g["gene_id"].astype(str)
            tx["gene_id"] = tx["gene_id"].astype(str)
            before = len(g)
            g = g.merge(tx, on="gene_id", how="left")
            if len(g) != before:
                raise SystemExit("taxonomy join changed the gene count")
            log(f"genes: region attached, {int(g['region'].isna().sum())} without one")
        write("genes", g)

    # named dimensions, the 08-03 reframe, kept in separate namespaces because
    # the two tables are keyed on PC index but computed in different feature
    # spaces (PLAN.md trap #10).
    dim = P.by_key("dimension_names").path
    write("dimensions", read_table(dim), dim)

    rep = P.by_key("dimension_reproducibility").path
    df_rep = read_table(rep)
    df_rep = df_rep.rename(columns={"PC": "component_index"})
    df_rep["feature_space"] = "shared63"
    write("dimension_reproducibility", df_rep, rep)

    # peaks ------------------------------------------------------------------
    # Signed placement is computed here, once, from peak_midpoint minus the
    # viewpoint. `distance_to_viewpoint` in the source table is UNSIGNED and
    # using it would put half the peaks on the wrong side of the profile.
    if gene_ids is not None:
        pk_path = P.by_key("annotated_peaks").path
        pk = read_table(pk_path, low_memory=False)
        pk.columns = [c.lower().replace(" ", "_") for c in pk.columns]
        pk["symbol_key"] = pk["viewpoint"].astype(str).str.upper()

        vp = build_genes(gene_ids, source_labels)[["symbol_key", "viewpoint_pos"]]
        pk = pk.merge(vp, on="symbol_key", how="inner")
        pk["offset_bp"] = pk["peak_midpoint"].astype("int64") - pk["viewpoint_pos"].astype("int64")

        # Sanity: the unsigned source column must equal |our signed offset|.
        # A mismatch means the viewpoint join is wrong, which would be invisible
        # on a plot but wrong everywhere.
        delta = (pk["offset_bp"].abs() - pk["distance_to_viewpoint"].abs()).abs()
        bad = int((delta > 1).sum())
        if bad:
            log(f"WARNING: {bad}/{len(pk)} peaks disagree with distance_to_viewpoint "
                f"(max {delta.max():.0f} bp), check the viewpoint join")
        else:
            log(f"peak placement verified against unsigned distance for all {len(pk)}")

        # Only peaks inside the +/-1 Mb window can be drawn.
        in_window = pk["offset_bp"].abs() <= (S.N_BINS_L0 * S.BIN_SIZE_L0) // 2
        log(f"peaks: {len(pk)} mapped, {int((~in_window).sum())} outside the window")
        keep = ["symbol_key", "viewpoint", "viewpoint_id", "chromosome", "start", "end",
                "peak_midpoint", "offset_bp", "peak_max", "peak_size", "sharpness",
                "log2_enrichment", "consensus_fraction", "re"]
        write("peaks", pk.loc[in_window, [c for c in keep if c in pk.columns]], pk_path)

    # features ---------------------------------------------------------------
    # 91 z-scored features plus panel percentile, for the feature table view.
    import pickle as _pickle
    fpath = P.by_key("features_gw").path
    with open(fpath, "rb") as f:
        F = _pickle.load(f)
    fm = F["feature_matrix"].copy()
    fnames = list(F["feature_names"])
    idcol = [c for c in fm.columns if c not in fnames]
    fm = fm.rename(columns={idcol[0]: "gene_id"}) if idcol else fm
    long = fm.melt(id_vars="gene_id", value_vars=fnames,
                   var_name="feature", value_name="z")
    # Panel percentile per feature, what "where does this gene sit" means.
    long["percentile"] = long.groupby("feature")["z"].rank(pct=True) * 100
    long["symbol_key"] = long["gene_id"].astype(str).str.upper()
    write("features", long, fpath)
    log(f"features: {len(fnames)} features x {fm.shape[0]} genes (long form)")

    # embeddings -------------------------------------------------------------
    emb, loadings, scree = build_embeddings(fm, fnames)
    write("embeddings", emb)
    write("pc_loadings", loadings)
    write("pc_scree", scree)

    # reach x composition taxonomy -------------------------------------------
    # Replaces the flat archetypes for display. Membership columns (w_*) are the
    # point: every gene in the panel has a top weight below 0.5, so a hard label
    # would assert a belonging the data does not support.
    for key, name in (("taxonomy_labels", "taxonomy"),
                      ("taxonomy_radar", "taxonomy_radar")):
        pth = P.by_key(key).path
        if pth.exists():
            write(name, read_table(pth), pth)
        else:
            log(f"WARNING: {key} missing, taxonomy views will be unavailable")

    # amount-corrected ("shape") embeddings, the DISPLAY DEFAULT ---------------
    semb, sload, sscree, sbasis = build_shape_embeddings(fm, fnames)
    write("shape_embeddings", semb)
    write("shape_loadings", sload)
    write("shape_scree", sscree)
    write("shape_basis", sbasis)

    # reproducibility pairs (lab page) ---------------------------------------
    write("reproducibility_pairs", build_repro_pairs(fm, fnames))

    # cohort membership ------------------------------------------------------
    write("cohort_membership", build_cohorts(fm))

    # cohort view ------------------------------------------------------------
    grp = P.by_key("external_groups").path
    df_grp = read_table(grp)
    df_grp["passes_min_n"] = df_grp["n_in_group"] >= S.MIN_GROUP_N
    write("external_groups", df_grp, grp)

    # anchors, gene_table is the 42-col superset, preferred over anchor_table
    anc = P.by_key("gene_table").path
    write("anchors", read_table(anc, low_memory=False), anc)

    # reproducibility demo ---------------------------------------------------
    cpr = P.by_key("cross_panel_reproducibility").path
    write("reproducibility", read_table(cpr), cpr)

    # confounds, numbers only, no verdicts (PLAN.md §1b display rule) --------
    for key, name in (("locus_confounds", "confounds_locus"),
                      ("window_confounds", "confounds_window")):
        art = P.by_key(key)
        if art.path.exists():
            write(name, read_table(art.path), art.path)

    # explanation layer ------------------------------------------------------
    nb = P.by_key("nested_baselines").path
    write("nested_baselines", read_table(nb), nb)

    # P(s) -------------------------------------------------------------------
    ps = P.by_key("ps_fits").path
    write("ps_fits", read_table(ps), ps)

    # radar ------------------------------------------------------------------
    rad = P.by_key("radar_profiles").path
    write("radar", read_table(rad), rad)

    # labels -----------------------------------------------------------------
    # Two files compose the final labelling. archetype_labels.tsv covers all
    # 1,846 genes at k=4; me_subtypes_labels.tsv covers ONLY the 720 arch-ME
    # genes. Using the ME file alone would label 39% of the panel and silently
    # drop the rest, it is a subset, not a labelling.
    base_p = P.by_key("archetype_labels").path
    me_p = P.by_key("me_subtypes_labels").path
    base = read_table(base_p)
    me = read_table(me_p)

    labels = base.merge(
        me.drop(columns=["gene_symbol"]), on="symbol_key", how="left",
        validate="one_to_one",
    )
    labels["group"] = labels["me_subtype"].fillna("").where(
        labels["me_subtype"].notna(), labels["archetype"]
    )
    labels["group"] = labels.apply(
        lambda r: f"arch-ME-{r['me_subtype'].split('-')[-1].lower()}"
        if pd.notna(r["me_subtype"]) else r["archetype"],
        axis=1,
    )
    # arch-off is a QC class (near-empty signal), not a biological archetype.
    labels["is_biological"] = labels["group"] != "arch-off"

    n_me = int(labels["me_subtype"].notna().sum())
    if n_me != (base["archetype"] == "arch-ME").sum():
        raise SystemExit(
            f"ME split covers {n_me} genes but archetype_labels has "
            f"{(base['archetype'] == 'arch-ME').sum()} arch-ME, merge is wrong"
        )
    log(f"labels composed: {len(labels)} genes, {n_me} carry an ME subtype")
    log("  " + ", ".join(f"{k} {v}" for k, v in labels["group"].value_counts().items()))
    write("labels", labels, base_p)

    post = P.by_key("posteriors").path
    df_post = read_table(post)
    required = {"gene_symbol", "max_posterior", "assignment_entropy", "confidence_class"}
    missing = required - set(df_post.columns)
    if missing:
        raise SystemExit(
            f"posteriors table is missing {sorted(missing)}; got {list(df_post.columns)}"
        )
    # The pipeline already bins confidence as core(>0.8) / leaning(0.5-0.8) etc.
    # Derive the boolean from the same threshold rather than inventing a second
    # convention that could drift from confidence_class.
    df_post["is_core"] = df_post["max_posterior"] >= S.CORE_POSTERIOR_MIN
    n_core = int(df_post["is_core"].sum())
    log(f"posteriors: {n_core}/{len(df_post)} core (>{S.CORE_POSTERIOR_MIN}), "
        f"{len(df_post) - n_core} mixtures")
    write("posteriors", df_post, post)

    return written


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables-only", action="store_true",
                    help="skip the 4.4 GB pickle and rebuild tables only")
    ap.add_argument("--levels", default="0,1,2",
                    help="pyramid levels to write (default 0,1,2)")
    ap.add_argument("--batch", type=int, default=128,
                    help="genes per write block")
    args = ap.parse_args()

    levels = sorted({int(x) for x in args.levels.split(",") if x.strip()})
    for lvl in levels:
        if lvl not in S.LEVELS:
            raise SystemExit(f"unknown level {lvl}; known: {sorted(S.LEVELS)}")

    started = datetime.now(timezone.utc)
    print("\nmccprofiler-app, build store\n")

    profiles_meta: dict = {}
    gene_ids = source_labels = None

    if args.tables_only:
        log("--tables-only: skipping profiles")
        genes_pq = P.STORE / "tables" / "genes.parquet"
        if genes_pq.exists():
            g = pd.read_parquet(genes_pq)
            gene_ids, source_labels = list(g["gene_id"]), list(g["source"])
    else:
        profiles_meta = build_profiles(levels, args.batch)
        gene_ids = profiles_meta["gene_ids"]
        source_labels = profiles_meta["source_labels"]

    print()
    tables = build_tables(gene_ids, source_labels)

    manifest = {
        "store_version": S.STORE_VERSION,
        "built": started.isoformat(timespec="seconds"),
        "duration_s": round((datetime.now(timezone.utc) - started).total_seconds(), 1),
        "panel": P.PANEL_GW,
        "scripts_cleaned_commit": git_commit(P.SCRIPTS),
        "app_commit": git_commit(P.APP_ROOT),
        "profiles": {k: v for k, v in profiles_meta.items()
                     if k not in {"gene_ids", "source_labels"}},
        "tables": tables,
        "inputs": {a.key: stamp(a.path) for a in P.ARTEFACTS},
        "blacklist": P.BLACKLIST,
    }
    (P.STORE / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\n  store at {P.STORE}")
    print(f"  manifest.json written, {manifest['duration_s']}s\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
