"""Which panel genes have no peaks called but clearly should have?

Found 2026-08-17 while checking why LCK-like genes get misleading taxonomy
memberships. Three genes (PADI4, TTPAL, ASB9) have ZERO rows in annotated.tsv
despite:

  * being present in the viewpoint list for gw_cd4_1,
  * carrying total MCC signal at or above the panel median,
  * and showing local O/E enrichment comparable to, or higher than, genes that
    DID get peaks called. TTPAL reaches O/E max 59.4 with 73 bins above O/E 5,
    against 20.8-46.5 and 32-52 for matched genes that were called.

They are also absent from `unfiltered.tsv` (132,447 rows), which rules out every
downstream threshold: consensus, coverage, IDR. Whatever drops them happens
between the viewpoint list and the peak output, and it happens silently.

`config.py` documents exactly this failure mode for the viewpoint join: "a
mismatch silently drops all peaks for that source rather than raising". The
viewpoint file also holds 2,000 rows against 1,846 panel genes, so 154 entries
are already unaccounted for.

This script does the systematic version: every panel gene, ranked by how much
enrichment it has relative to how many peaks it got, so the audit is not limited
to the three that happened to surface.

WHY IT MATTERS BEYOND THREE GENES. A gene with no peaks has every element-class
and per-class feature at a default, which is how absence became an extreme value
in the first place (see the NaN fix, CURRENT_FINDINGS section 0). It also makes
the taxonomy undescribable for that gene: the app now reports "no peaks called"
at 100% rather than inventing a membership.

Run in cd4env. Writes audit/scripts/output/missing_peaks.tsv.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import warnings
import numpy as np
import pandas as pd
from app.store import get_store

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)

PEAKS = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
         "scripts_cleaned/process/output_resgw_10000/utils/peaks/annotated.tsv")
VIEWPOINTS = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
              "mcc/viewpoints/gw_cd4_1_viewpoints_final.tsv")

s = get_store()
genes = s.genes[["gene_id", "gene_symbol", "symbol_key"]].copy()

pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()
n_peaks = pk.groupby("symbol_key").size()

vp = pd.read_csv(VIEWPOINTS, sep="\t", header=None)
vp.columns = ["chrom", "start", "end", "viewpoint"] + \
             [f"c{i}" for i in range(4, vp.shape[1])]
vp["symbol_key"] = vp["viewpoint"].astype(str).str.upper()
print(f"{len(genes)} panel genes, {len(vp)} viewpoints, "
      f"{pk['symbol_key'].nunique()} genes with >=1 peak")

# Enrichment per gene, computed from the stored profiles so this is independent
# of the peak caller entirely.
rows = []
for gid, sym, key in genes.itertuples(index=False):
    try:
        raw = np.asarray(s.profile(gid, channel="mcc", level=2, mode="raw")["values"], float)
        oe = np.asarray(s.profile(gid, channel="mcc", level=2, mode="oe")["values"], float)
    except Exception:                                              # noqa: BLE001
        continue
    raw = raw[np.isfinite(raw)]
    oe = oe[np.isfinite(oe)]
    if raw.size == 0:
        continue
    nz = raw[raw > 0]
    rows.append({
        "gene_symbol": sym,
        "symbol_key": key,
        "n_peaks": int(n_peaks.get(key, 0)),
        "in_viewpoint_list": bool((vp["symbol_key"] == key).any()),
        "total_signal": round(float(raw.sum()), 1),
        "max_bin": round(float(raw.max()), 2),
        "oe_max": round(float(oe.max()) if oe.size else 0.0, 1),
        "bins_oe_gt3": int((oe > 3).sum()),
        "bins_oe_gt5": int((oe > 5).sum()),
        "frac_nonzero": round(float((raw > 0).mean()), 3),
    })
d = pd.DataFrame(rows)

# The reference distribution: what enrichment looks like among genes that WERE
# called. A gene with no peaks whose enrichment sits inside this range is a
# calling failure, not a quiet locus.
called = d[d.n_peaks > 0]
ref = {c: (called[c].quantile(0.10), called[c].median())
       for c in ("total_signal", "oe_max", "bins_oe_gt5")}
print("\nreference (genes WITH peaks): 10th percentile / median")
for c, (q10, med) in ref.items():
    print(f"  {c:<16}{q10:>10.1f}{med:>10.1f}")

d["above_p10_signal"] = d.total_signal >= ref["total_signal"][0]
d["above_p10_oe"] = d.oe_max >= ref["oe_max"][0]
d["above_p10_bins"] = d.bins_oe_gt5 >= ref["bins_oe_gt5"][0]
d["n_criteria_met"] = (d.above_p10_signal.astype(int)
                       + d.above_p10_oe.astype(int)
                       + d.above_p10_bins.astype(int))

# SHOULD have peaks: none called, yet enrichment inside the called range.
suspect = d[(d.n_peaks == 0) & (d.n_criteria_met >= 2)].copy()
suspect = suspect.sort_values("bins_oe_gt5", ascending=False)

print("\n" + "=" * 88)
print("GENES WITH ZERO PEAKS THAT SHOULD HAVE THEM")
print("=" * 88)
if len(suspect):
    print(f"{'gene':<10}{'total':>8}{'oe max':>8}{'O/E>5':>7}{'O/E>3':>7}"
          f"{'nonzero':>9}{'in VP list':>12}{'criteria':>10}")
    for r in suspect.itertuples():
        print(f"{r.gene_symbol:<10}{r.total_signal:>8.0f}{r.oe_max:>8.1f}"
              f"{r.bins_oe_gt5:>7}{r.bins_oe_gt3:>7}{r.frac_nonzero:>9.2f}"
              f"{str(r.in_viewpoint_list):>12}{r.n_criteria_met:>8}/3")
else:
    print("  none")

# Near-misses: very few peaks relative to strong enrichment. Same mechanism,
# partially expressed, and these are the ones nobody would ever notice.
low = d[(d.n_peaks > 0) & (d.n_peaks <= 3)
        & (d.bins_oe_gt5 >= ref["bins_oe_gt5"][1])].copy()
low = low.sort_values("bins_oe_gt5", ascending=False)
print("\n" + "=" * 88)
print("UNDER-CALLED: 1-3 peaks despite above-MEDIAN enrichment")
print("=" * 88)
if len(low):
    print(f"{'gene':<10}{'peaks':>7}{'total':>8}{'oe max':>8}{'O/E>5':>7}")
    for r in low.head(25).itertuples():
        print(f"{r.gene_symbol:<10}{r.n_peaks:>7}{r.total_signal:>8.0f}"
              f"{r.oe_max:>8.1f}{r.bins_oe_gt5:>7}")
    print(f"  ... {len(low)} genes total")
else:
    print("  none")

# Viewpoints that never became panel genes at all: the 154-row discrepancy.
orphan = vp[~vp["symbol_key"].isin(set(genes.symbol_key))]
print("\n" + "=" * 88)
print("VIEWPOINTS THAT ARE NOT PANEL GENES")
print("=" * 88)
print(f"  {len(orphan)} of {len(vp)} viewpoints have no gene in the store.")
print(f"  Of those, {int(orphan['symbol_key'].isin(set(pk.symbol_key)).sum())} DO "
      f"have peaks called, so they were captured and processed but dropped from "
      f"the panel later (QC or outlier removal).")

d["should_have_peaks"] = d.index.isin(suspect.index)
d["under_called"] = d.index.isin(low.index)
d.to_csv(OUT / "missing_peaks.tsv", sep="\t", index=False)
orphan[["chrom", "start", "end", "viewpoint"]].to_csv(
    OUT / "viewpoints_not_in_panel.tsv", sep="\t", index=False)

print("\n" + "=" * 88)
print(f"  {len(suspect)} genes should have peaks and have none")
print(f"  {len(low)} genes are under-called (1-3 peaks, above-median enrichment)")
print(f"  {len(orphan)} viewpoints never became panel genes")
print(f"\nwrote {OUT / 'missing_peaks.tsv'}")
print(f"wrote {OUT / 'viewpoints_not_in_panel.tsv'}")
