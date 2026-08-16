"""A health check on every feature in the substrate.

Written 2026-08-16 after four separate feature problems were found one at a time
by accident. Each was a different failure mode, and each was caught only because
something downstream looked odd:

  DEGENERATE      mean_degree and friends were exact algebraic functions of two
                  peak counts, found because mean_degree was PC1's second
                  strongest correlate and that seemed too good.
  UNREPRODUCIBLE  the per-peak moment families sat at rho 0.27 to 0.68, found
                  while asking which components were trustworthy.
  MIS-ORIENTED    contact_asymmetry is computed in genome coordinates, so a
                  gene-relative signal cancels across strands. Found while
                  sourcing strand for an unrelated test.
  CONFOUNDED      nothing in the substrate had ever been checked against gene
                  length, which turned out to explain a headline result.

This runs all four checks over every feature at once, so the next one is found
deliberately rather than by luck.

  1. PREDICTABILITY   R2 of each feature from all the others. A feature at
                      R2 > 0.95 is a deterministic function of the rest and adds
                      nothing. This is the check that would have caught
                      mean_degree without anyone noticing it in a loadings plot.
  2. REPRODUCIBILITY  cross-panel rho on the twice-captured genes.
  3. CONFOUNDING      correlation with gene length, expression, CpG, gene
                      density, and chr19 membership.
  4. DISTRIBUTION     skew, kurtosis, zero-inflation and saturation. A feature
                      that is 90% zeros or pinned at a bound is not measuring
                      what its name says across most of the panel.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_predict
from app.store import get_store

warnings.filterwarnings("ignore")
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
TSS_BED = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/"
           "01_tss_regions.bed")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
feats = list(W.columns)
X = W.to_numpy(float)
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
print(f"{len(feats)} features, {len(W)} genes\n")

# --- 2. reproducibility ----------------------------------------------------
rep = s.table("reproducibility_pairs")
rho = {}
for f, sub in rep.groupby("feature"):
    if len(sub) > 10 and sub["gw"].std() > 0 and sub["immune"].std() > 0:
        rho[f] = float(np.corrcoef(sub["gw"].rank(), sub["immune"].rank())[0, 1])

# --- 3. confounders --------------------------------------------------------
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")
bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed = bed.drop_duplicates("symbol_key").set_index("symbol_key")

conf = pd.DataFrame(index=range(len(sym)))
for name, col in [("length", "log10_gene_length"), ("expr", "gtex_whole_blood_tpm"),
                  ("cpg", "cpg_oe_ratio"), ("density", "gene_density_500kb")]:
    v = pd.Series(sym).map(gt[col]).to_numpy(float)
    if name == "expr":
        v = np.log1p(v)
    conf[name] = np.where(np.isfinite(v), v, np.nanmedian(v))
conf["chr19"] = (pd.Series(sym).map(bed["chrom"]) == "chr19").astype(float).to_numpy()
C = conf.to_numpy(float)
C = (C - C.mean(0)) / np.where(C.std(0) > 0, C.std(0), 1)

# --- 1. predictability from the other features -----------------------------
print("computing leave-one-feature-out predictability (this is the slow part)...")
kf = KFold(5, shuffle=True, random_state=0)
pred_r2, conf_r2 = [], []
for j in range(len(feats)):
    y = X[:, j]
    others = np.delete(X, j, axis=1)
    p = cross_val_predict(LinearRegression(), others, y, cv=kf)
    pred_r2.append(max(0.0, 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    p2 = cross_val_predict(LinearRegression(), C, y, cv=kf)
    conf_r2.append(max(0.0, 1 - ((y - p2) ** 2).sum() / ((y - y.mean()) ** 2).sum()))

# --- 4. distribution -------------------------------------------------------
raw_path = Path("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
                "scripts_cleaned/mccprofiler/outputs_gw_cd4_1/features/features_raw.pkl")
raw = pd.read_pickle(raw_path)
R = raw["feature_matrix"] if isinstance(raw, dict) else raw

d = pd.DataFrame({
    "feature": feats,
    "predictable_R2": np.round(pred_r2, 3),
    "rho": [round(rho[f], 2) if f in rho else np.nan for f in feats],
    "conf_R2": np.round(conf_r2, 3),
})
for j, f in enumerate(feats):
    v = R[f].to_numpy(float) if f in R.columns else np.full(len(R), np.nan)
    v = v[np.isfinite(v)]
    d.loc[j, "zero_frac"] = round(float((v == 0).mean()), 2) if len(v) else np.nan
    d.loc[j, "skewness"] = round(float(pd.Series(v).skew()), 1) if len(v) else np.nan

print("\n" + "=" * 92)
print("1. PREDICTABLE FROM THE OTHER FEATURES  (R2 > 0.90 adds little; > 0.95 is a derivative)")
print("=" * 92)
hi = d.sort_values("predictable_R2", ascending=False)
for r in hi.head(10).itertuples():
    flag = "  <-- DERIVATIVE" if r.predictable_R2 > .95 else ("  <-- redundant" if r.predictable_R2 > .90 else "")
    print(f"   {r.predictable_R2:.3f}  {r.feature:<38}{flag}")
print(f"   ... {int((d.predictable_R2 > .95).sum())} above 0.95, "
      f"{int((d.predictable_R2 > .90).sum())} above 0.90, median {d.predictable_R2.median():.2f}")

print("\n" + "=" * 92)
print("2. REPRODUCIBILITY  (rho on the twice-captured genes)")
print("=" * 92)
meas = d.dropna(subset=["rho"]).sort_values("rho")
for r in meas.head(8).itertuples():
    flag = "  <-- below 0.70" if r.rho < .70 else ""
    print(f"   {r.rho:.2f}  {r.feature:<38}{flag}")
print(f"   ... {int((meas.rho < .70).sum())} of {len(meas)} measured features below 0.70, "
      f"median {meas.rho.median():.2f}   ({len(d) - len(meas)} not measured twice)")

print("\n" + "=" * 92)
print("3. CONFOUNDED  (R2 from length, expression, CpG, density, chr19)")
print("=" * 92)
cf = d.sort_values("conf_R2", ascending=False)
for r in cf.head(8).itertuples():
    flag = "  <-- >25% is gene properties" if r.conf_R2 > .25 else ""
    print(f"   {r.conf_R2:.3f}  {r.feature:<38}{flag}")
print(f"   ... {int((d.conf_R2 > .25).sum())} above 0.25, {int((d.conf_R2 > .10).sum())} above 0.10, "
      f"median {d.conf_R2.median():.3f}")

print("\n" + "=" * 92)
print("4. DISTRIBUTION  (zero-inflated or extreme skew)")
print("=" * 92)
dd = d.dropna(subset=["zero_frac"])
bad = dd[(dd.zero_frac > .5) | (dd.skewness.abs() > 5)].sort_values("zero_frac", ascending=False)
for r in bad.head(10).itertuples():
    why = []
    if r.zero_frac > .5:
        why.append(f"{r.zero_frac:.0%} zeros")
    if abs(r.skewness) > 5:
        why.append(f"skew {r.skewness:+.0f}")
    print(f"   {r.feature:<38}{', '.join(why)}")
print(f"   ... {int((dd.zero_frac > .5).sum())} features over 50% zeros, "
      f"{int((dd.skewness.abs() > 5).sum())} with |skew| > 5")

out = Path(__file__).resolve().parent / "output" / "feature_health.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(out, sep="\t", index=False)

print("\n" + "=" * 92)
print("FEATURES FAILING MORE THAN ONE CHECK")
print("=" * 92)
d["fails"] = ((d.predictable_R2 > .90).astype(int)
              + (d.rho.fillna(1) < .70).astype(int)
              + (d.conf_R2 > .25).astype(int)
              + (d.zero_frac.fillna(0) > .5).astype(int))
mult = d[d.fails >= 2].sort_values("fails", ascending=False)
if len(mult):
    for r in mult.itertuples():
        print(f"   {r.feature:<38} predictable {r.predictable_R2:.2f}  "
              f"rho {r.rho if pd.notna(r.rho) else float('nan'):.2f}  "
              f"conf {r.conf_R2:.2f}  zeros {r.zero_frac:.0%}")
else:
    print("   none")
print(f"\nwrote {out}")
