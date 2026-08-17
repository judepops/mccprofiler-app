"""Is a membership percentage a claim, or an artefact of the softmax?

Written 2026-08-17 from a reader's question rather than a test. A gene with 79%
of its peak signal within 50 kb and ZERO beyond 250 kb was displayed carrying
20% extended-enhancer and 19% extended-ctcf, which reads as a contradiction.

It is not a contradiction and not a bug in the weights. With K = 5 regions, a
softmax cannot fall below 20% for a gene that resembles nothing in particular,
so 20% is the FLOOR: no evidence either way, and below it, evidence against. The
display was presenting the floor as though it were a share.

This script measures how much room the softmax actually has, and then finds the
opposite failure, which is worse because it looks like a strong result:

    weight_i = exp(-d_i^2 / tau) / sum_j exp(-d_j^2 / tau)

The distance is SQUARED. So a gene far from every centroid gets a confident
profile out of a small RELATIVE gap. That means the highest memberships in the
panel belong to the genes least like anything, i.e. exactly the genes a reader
would pick out as clean examples.

Nothing here affects the clustering. Assignment is argmin over distances and
never touches the softmax. This is entirely about what may be reported.

Run with the mccapp interpreter from the app root:
    /home/imm/grte4643/miniconda3/envs/mccapp/bin/python \
        audit/scripts/diagnose_softmax_floor.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAX = ROOT / "store" / "tables" / "taxonomy.parquet"
GENES = ROOT / "store" / "tables" / "genes.parquet"

t = pd.read_parquet(TAX)
dcols = sorted(c for c in t.columns if c.startswith("d_"))
wcols = sorted(c for c in t.columns if c.startswith("w_"))
D = t[dcols].to_numpy(float)
W = t[wcols].to_numpy(float) * 100.0
K = W.shape[1]
uniform = 100.0 / K
sym = (pd.read_parquet(GENES)
       .assign(gene_id=lambda x: x.gene_id.astype(str))
       .set_index("gene_id")["gene_symbol"])
symbols = t["gene_id"].astype(str).map(sym).fillna(t["gene_id"].astype(str))

near, far = D.min(1), D.max(1)
top, low = W.max(1), W.min(1)
# Total variation distance from uniform: 0 = the profile says nothing at all,
# 1 = a hard label. The single most useful summary of "does this mean anything".
tvd = 0.5 * np.abs(W / 100.0 - 1.0 / K).sum(1)

print("=" * 78)
print(f"HOW MUCH ROOM DOES THE SOFTMAX HAVE?   K = {K}, floor = {uniform:.1f}%")
print("=" * 78)
print(f"  nearest-centroid distance   median {np.median(near):6.2f}"
      f"   range {near.min():.2f}-{near.max():.2f}")
print(f"  furthest-centroid distance  median {np.median(far):6.2f}")
print(f"  spread as a fraction of the nearest distance:"
      f" median {np.median((far - near) / near) * 100:.0f}%")
print("\n  In 77 dimensions every centroid is far from every point, so the")
print("  distances cannot separate. This is concentration of measure, not a")
print("  badly chosen temperature: retuning tau trades one failure mode below")
print("  for the other, it cannot remove both.")

print("\n" + "=" * 78)
print("FAILURE 1: THE FLOOR. Most profiles carry no information.")
print("=" * 78)
print(f"  top weight       median {np.median(top):5.1f}%   (uniform {uniform:.1f}%)")
print(f"  lowest weight    median {np.median(low):5.1f}%")
print(f"  top minus lowest median {np.median(top - low):5.1f} pp")
print(f"  TVD from uniform median {np.median(tvd):.3f}"
      f"   p90 {np.quantile(tvd, 0.9):.3f}   max {tvd.max():.3f}")
flat = np.abs(W - uniform).max(1) <= 10
print(f"\n  ENTIRE profile within 10 pp of the floor: {int(flat.sum())} of "
      f"{len(W)}  ({flat.mean() * 100:.0f}%)")
print(f"  any weight >= 40%:                       {int((top >= 40).sum())} "
      f"({(top >= 40).mean() * 100:.0f}%)")

print("\n" + "=" * 78)
print("FAILURE 2: DISTANCE-INFLATED CONFIDENCE. The dangerous one.")
print("=" * 78)
print("  The weight squares the distance, so a gene far from EVERY centroid")
print("  gets a confident profile from a small relative gap.")
print(f"\n  corr(nearest distance, top weight) = {np.corrcoef(near, top)[0, 1]:+.3f}")
hi = near >= np.quantile(near, 0.90)
print(f"  top weight, top distance decile: median {np.median(top[hi]):.1f}%"
      f"   vs {np.median(top[~hi]):.1f}% for the rest")
print(f"  of the {int((top >= 40).sum())} genes reaching any weight >= 40%, "
      f"{int((top[hi] >= 40).sum())} are in the top distance decile")
print("\n  So more than half the confident-looking assignments come from a")
print("  tenth of the panel selected for being unlike everything. Those are")
print("  precisely the genes a reader would quote as clean examples.")

print("\n  worked example, the most confident-looking genes in the panel:")
print(f"  {'gene':<12}{'top w':>8}{'d near':>9}{'d far':>8}"
      f"{'rel gap':>9}{'d^2 gap':>10}{'pct':>7}")
for j in np.argsort(-top)[:8]:
    rel = (far[j] - near[j]) / near[j] * 100
    d2 = far[j] ** 2 - near[j] ** 2
    pct = (near < near[j]).mean() * 100
    print(f"  {symbols.iloc[j]:<12}{top[j]:7.1f}%{near[j]:9.2f}{far[j]:8.2f}"
          f"{rel:8.0f}%{d2:10.0f}{pct:6.0f}%")
print("\n  Read the last two columns together: a 10% gap in distance becomes a")
print("  gap of hundreds once squared, and the softmax reports that as")
print("  certainty. The `pct` column is the distance percentile: near 100 means")
print("  the gene is far from all five regions, not close to one.")

print("\n" + "=" * 78)
print("HOW THE APP NOW REPORTS THIS")
print("=" * 78)
p90 = np.quantile(near, 0.90)
print(f"  floor              {uniform:.1f}%  (weights shown as signed deviations)")
print(f"  near_uniform       max |deviation| <= 10 pp        -> {int(flat.sum())} genes")
print(f"  distance_inflated  nearest distance >= p90 = {p90:.2f}"
      f"  -> {int(hi.sum())} genes")
for g in ("SACS", "GATA3", "BCCIP", "LCK"):
    m = symbols == g
    if not m.any():
        continue
    j = int(np.flatnonzero(m.to_numpy())[0])
    tag = ("distance-inflated" if near[j] >= p90
           else "flat" if flat[j] else "informative")
    print(f"    {g:<8}d {near[j]:6.2f} (p{(near < near[j]).mean() * 100:4.1f})"
          f"  max dev {np.abs(W[j] - uniform).max():5.1f} pp"
          f"  tvd {tvd[j]:.3f}   -> {tag}")

print("\n  The clustering is untouched: assignment is argmin over distances and")
print("  never uses the softmax. Only what may be REPORTED has changed.")
