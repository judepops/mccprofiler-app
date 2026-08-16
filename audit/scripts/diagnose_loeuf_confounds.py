"""Does the constraint result survive gene length, expression and CpG density?

The forward plan makes constraint the most promising positive result: LOEUF is
predicted at r = 0.210 by all 91 features and **0.191 by the amount-corrected
shape substrate alone**, so 91% of the signal is shape rather than magnitude.
Nothing else in the panel behaves like that, and Phase B is built on it.

Three confounders have never been controlled anywhere in this project:

  GENE LENGTH   LOEUF is a depletion statistic, so it scales with coding
                sequence length by construction: longer genes accumulate more
                expected loss-of-function variants and give more power to
                detect depletion. Length ALSO changes contact features
                mechanically, because a longer gene occupies more of the
                plus-or-minus 1 Mb window. This was flagged as a gap in the
                feature review weeks ago and never addressed.
  EXPRESSION    correlates with constraint and with contact amount.
  CpG DENSITY   same.

So the question is whether the 0.191 is architecture or is length wearing an
architectural costume. Nested, cross-validated, matching the methodology of
audit/continuous_methods/nested_baselines.tsv:

  1. confounders only          length, expression, CpG
  2. shape only                the amount-corrected substrate
  3. confounders + shape       does shape ADD anything?
  4. shape, confounders partialled out of BOTH sides

Model 3 versus 1 is the test that matters. If shape adds nothing over the
confounders, Phase B is not worth building and an hour here has saved a phase.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.model_selection import KFold
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
SEED = 0
N_SPLITS = 5
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
cols = ["gnomad_loeuf", "log10_gene_length", "gtex_whole_blood_tpm", "cpg_oe_ratio"]
m = pd.DataFrame({"symbol_key": sym}).merge(
    gt[["symbol_key"] + cols].drop_duplicates("symbol_key"), on="symbol_key", how="left")

ok = m[cols].notna().all(axis=1).to_numpy()
print(f"{ok.sum()} of {len(m)} genes have LOEUF and all three confounders\n")

y = m.loc[ok, "gnomad_loeuf"].to_numpy(float)
C = m.loc[ok, ["log10_gene_length", "gtex_whole_blood_tpm", "cpg_oe_ratio"]].to_numpy(float)
C[:, 1] = np.log1p(C[:, 1])                 # TPM is heavy-tailed
C = (C - C.mean(0)) / C.std(0)

mag = [x for x in MAG_OVERALL if x in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
SHAPE = PCA(n_components=18, random_state=0).fit_transform(X - A @ beta)[ok]
SHAPE = (SHAPE - SHAPE.mean(0)) / SHAPE.std(0)

print("How much do the confounders relate to LOEUF and to each other?")
names = ["log10 gene length", "log blood TPM", "CpG o/e"]
for i, nm in enumerate(names):
    print(f"   r(LOEUF, {nm:<18}) = {np.corrcoef(y, C[:, i])[0, 1]:+.3f}")
print(f"   r(length, CpG)              = {np.corrcoef(C[:, 0], C[:, 2])[0, 1]:+.3f}")


def cv_r(Xd, yd, seed=SEED):
    """Cross-validated Pearson r between prediction and truth."""
    kf = KFold(N_SPLITS, shuffle=True, random_state=seed)
    pred = np.zeros_like(yd, dtype=float)
    for tr, te in kf.split(Xd):
        mdl = RidgeCV(alphas=np.logspace(-2, 4, 25)).fit(Xd[tr], yd[tr])
        pred[te] = mdl.predict(Xd[te])
    return float(np.corrcoef(pred, yd)[0, 1])


def resid(target, on):
    return target - LinearRegression().fit(on, target).predict(on)


rows = [
    ("1. confounders only (length, expression, CpG)", cv_r(C, y)),
    ("2. shape only (amount-corrected, 18 comps)", cv_r(SHAPE, y)),
    ("3. confounders + shape", cv_r(np.column_stack([C, SHAPE]), y)),
    ("4. shape, confounders partialled from both",
     cv_r(np.column_stack([resid(SHAPE[:, j], C) for j in range(SHAPE.shape[1])]),
          resid(y, C))),
]
print("\n" + "=" * 74)
print("NESTED TEST: does shape predict constraint beyond length/expression/CpG?")
print("=" * 74)
for lab, r in rows:
    print(f"  {lab:<48} r = {r:+.3f}")

add = rows[2][1] - rows[0][1]
print(f"\n  shape adds {add:+.3f} over the confounders alone")
print(f"  after partialling, shape retains r = {rows[3][1]:+.3f} "
      f"against {rows[1][1]:+.3f} unadjusted "
      f"({100 * rows[3][1] / rows[1][1]:.0f}% retained)")

# Which single confounder does the damage?
print("\nleave-one-confounder-in, to see which one matters:")
for i, nm in enumerate(names):
    r = cv_r(np.column_stack([resid(SHAPE[:, j], C[:, [i]]) for j in range(SHAPE.shape[1])]),
             resid(y, C[:, [i]]))
    print(f"   partialling {nm:<18} shape r = {r:+.3f}")

print("\nModel 3 versus 1 is the test. If shape adds nothing over the confounders,")
print("the constraint result is length and expression, and Phase B is not worth")
print("building on it.")
