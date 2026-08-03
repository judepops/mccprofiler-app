"""Read-side access to the store.

The server touches nothing else. It does not import mccprofiler, does not open
the 4.4 GB pickle, and does not know where the science tree lives — everything
comes from `store/`, built by backend/scripts/build_store.py.

Tables are small (largest is 18,802 x 42) so they load once at startup. Profiles
live in HDF5 and are read per gene, one chunk per read.
"""

from __future__ import annotations

import json
import threading
from functools import cached_property
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from . import store_schema as S


class StoreMissing(RuntimeError):
    """Raised when the store has not been built yet."""


class Store:
    def __init__(self, root: Path):
        self.root = Path(root)
        self._h5_lock = threading.Lock()
        self._h5: h5py.File | None = None

        if not (self.root / "manifest.json").exists():
            raise StoreMissing(
                f"no store at {self.root}\n"
                "build it first:  conda activate cd4env && "
                "python backend/scripts/build_store.py"
            )

    # -- provenance ---------------------------------------------------------

    @cached_property
    def manifest(self) -> dict:
        return json.loads((self.root / "manifest.json").read_text())

    # -- tables -------------------------------------------------------------

    @cached_property
    def _tables(self) -> dict[str, pd.DataFrame]:
        out: dict[str, pd.DataFrame] = {}
        tdir = self.root / "tables"
        for path in sorted(tdir.glob("*.parquet")):
            out[path.stem] = pd.read_parquet(path)
        return out

    def table(self, name: str) -> pd.DataFrame:
        try:
            return self._tables[name]
        except KeyError:
            raise KeyError(
                f"no table {name!r}; present: {', '.join(sorted(self._tables))}"
            ) from None

    def has(self, name: str) -> bool:
        return name in self._tables

    # -- gene index ---------------------------------------------------------

    @staticmethod
    def _merge_checked(left: pd.DataFrame, right: pd.DataFrame, on: str,
                       what: str, min_coverage: float = 0.5) -> pd.DataFrame:
        """Left-join, and refuse a join that matched almost nothing.

        A mis-keyed merge produces an all-NaN column rather than an error, which
        is how `confidence` reached the API as NaN for every gene. Coverage is
        asserted so the failure is loud.
        """
        merged = left.merge(right, on=on, how="left", validate="one_to_one")
        added = [c for c in right.columns if c != on]
        if added:
            coverage = merged[added[0]].notna().mean()
            if coverage < min_coverage:
                raise ValueError(
                    f"{what}: join on {on!r} matched only {coverage:.1%} of rows — "
                    f"check the key. right columns: {list(right.columns)}"
                )
        return merged

    @cached_property
    def genes(self) -> pd.DataFrame:
        """One row per gene, joined with labels and posteriors."""
        g = self.table("genes").copy()

        if self.has("labels"):
            lab = self.table("labels")
            cols = [c for c in ("gene_symbol", "archetype", "group", "is_biological")
                    if c in lab.columns]
            g = self._merge_checked(g, lab[cols], "gene_symbol", "labels")

        if self.has("posteriors"):
            post = self.table("posteriors")
            cols = [c for c in ("gene_symbol", "max_posterior", "assignment_entropy",
                                "confidence_class", "is_core", "soft_argmax")
                    if c in post.columns]
            g = self._merge_checked(g, post[cols], "gene_symbol", "posteriors")
            # A left-join can introduce NaN, which turns a bool column to object
            # and makes `~col` raise. Genes with no posterior are not core.
            if "is_core" in g.columns:
                g["is_core"] = g["is_core"].fillna(False).astype(bool)

        return g

    def posterior_mix(self, gene_symbol: str) -> dict[str, float]:
        """Per-archetype posteriors for one gene — the mixture, not the label."""
        post = self.table("posteriors")
        row = post.loc[post["gene_symbol"] == gene_symbol]
        if row.empty:
            return {}
        pcols = [c for c in post.columns if c.startswith("p_")]
        return {c[2:]: float(row.iloc[0][c]) for c in pcols}

    @cached_property
    def _gene_row(self) -> dict[str, int]:
        """gene_id -> row index in the HDF5 profile datasets."""
        with self._open() as h5:
            ids = [s.decode() if isinstance(s, bytes) else str(s)
                   for s in h5[S.H5_GENE_IDS][:]]
        return {gid: i for i, gid in enumerate(ids)}

    def resolve(self, query: str) -> str | None:
        """Accept a gene_id or a bare symbol; return the canonical gene_id."""
        if query in self._gene_row:
            return query
        hits = self.genes.loc[self.genes["gene_symbol"].str.upper() == query.upper(),
                              "gene_id"]
        return str(hits.iloc[0]) if len(hits) else None

    # -- profiles -----------------------------------------------------------

    def _open(self) -> h5py.File:
        return h5py.File(self.root / "profiles.h5", "r")

    @cached_property
    def channels(self) -> list[str]:
        with self._open() as h5:
            return sorted(h5[S.H5_PROFILES].keys())

    def profile(
        self,
        gene_id: str,
        channel: str = S.PRIMARY_CHANNEL,
        level: int = S.DEFAULT_LEVEL,
        mode: str = "raw",
        start_bp: int | None = None,
        end_bp: int | None = None,
    ) -> dict:
        """Return one gene's profile.

        Distances are signed bp relative to the viewpoint, which sits at the
        centre of the window. `mode='oe'` divides by the panel-mean expected
        curve with the same epsilon floor the pipeline uses, so the app's O/E
        agrees with mccprofiler's.
        """
        row = self._gene_row.get(gene_id)
        if row is None:
            raise KeyError(f"unknown gene_id {gene_id!r}")
        if level not in S.LEVELS:
            raise ValueError(f"unknown level {level}; known: {sorted(S.LEVELS)}")

        nbins = S.level_bins(level)
        bp = S.level_bp(level)
        half = nbins // 2

        # window -> bin slice, clamped
        lo, hi = 0, nbins
        if start_bp is not None:
            lo = max(0, int(np.floor(start_bp / bp)) + half)
        if end_bp is not None:
            hi = min(nbins, int(np.ceil(end_bp / bp)) + half)
        if hi <= lo:
            raise ValueError("empty window")

        with self._h5_lock, self._open() as h5:
            grp = h5[f"{S.H5_PROFILES}/{channel}"]
            if f"L{level}" not in grp:
                raise KeyError(f"channel {channel!r} has no level {level}")
            values = grp[f"L{level}"][row, lo:hi].astype(np.float64)
            expected_l0 = h5[f"{S.H5_EXPECTED}/{channel}"][:]

        if mode == "oe":
            factor = S.LEVELS[level]
            exp = (expected_l0.reshape(-1, factor).mean(axis=1)
                   if factor > 1 else expected_l0)
            values = values / np.maximum(exp[lo:hi], S.OE_EPS)
        elif mode != "raw":
            raise ValueError(f"unknown mode {mode!r}; use 'raw' or 'oe'")

        return {
            "gene_id": gene_id,
            "channel": channel,
            "mode": mode,
            "level": level,
            "bp_per_bin": bp,
            "start_bp": (lo - half) * bp,
            "end_bp": (hi - half) * bp,
            "n": int(hi - lo),
            "values": values.tolist(),
        }


_store: Store | None = None


def get_store(root: Path | None = None) -> Store:
    global _store
    if _store is None:
        from .paths import STORE
        _store = Store(root or STORE)
    return _store
