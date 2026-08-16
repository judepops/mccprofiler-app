"""How much of the 91-feature space do simple gene properties explain?

The predictive justification for Aim 1 is retracted: contact shape adds nothing
to gene length, expression and CpG density when predicting existing annotations
(`diagnose_loeuf_confounds.py`). That answers "do the features predict better",
and the answer is no.

This asks the different and more appropriate question: **do the features CONTAIN
something those properties do not?** A tool can fail to beat simple covariates at
predicting known labels while still measuring a quantity nobody has measured
before. Those are separate claims and only the second is what a new assay
readout is for.

Reported as variance explained, per feature and over the whole matrix:

  R2 of each feature on [log10 gene length, log blood TPM, CpG o/e]
  the same on the retained principal components
  and, for scale, the same against the 11-feature magnitude basis, which is a
  property of the assay rather than of the gene

Low numbers mean the features are not a re-description of gene properties. High
numbers would mean the substrate is largely a proxy for length and expression
and the whole project rests on a re-parameterisation.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_predict
from _shape import MAG_OVERALL
from app.store import get_store

warnings.filterwarnings("ignore")
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
CONF = ["log10_gene_length", "gtex_whole_blood_tpm", "cpg_oe_ratio"]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
m = (pd.DataFrame({"symbol_key": sym})
     .merge(gt[["symbol_key"] + CONF].drop_duplicates("symbol_key"),
            on="symbol_key", how="left"))
ok = m[CONF].notna().all(axis=1).to_numpy()
C = m.loc[ok, CONF].to_numpy(float)
C[:, 1] = np.log1p(C[:, 1])
C = (C - C.mean(0)) / C.std(0)
print(f"{ok.sum()} of {len(m)} genes have all three properties\n")

feats = list(W.columns)
X = W.to_numpy(float)[ok]
MAG = W[[x for x in MAG_OVERALL if x in W.columns]].to_numpy(float)[ok]


def var_explained(Y, P, cv=True):
    """Share of Y's variance explained by predictors P, cross-validated."""
    out = []
    kf = KFold(5, shuffle=True, random_state=0)
    for j in range(Y.shape[1]):
        y = Y[:, j]
        if y.std() == 0:
            out.append(0.0)
            continue
        pred = (cross_val_predict(LinearRegression(), P, y, cv=kf) if cv
                else LinearRegression().fit(P, y).predict(P))
        out.append(max(0.0, 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    return np.array(out)


r2 = var_explained(X, C)
print("=" * 78)
print("VARIANCE OF THE 91 FEATURES EXPLAINED BY LENGTH, EXPRESSION AND CpG")
print("=" * 78)
print(f"  mean R2 across features      {r2.mean():.3f}   ({r2.mean() * 100:.1f}%)")
print(f"  median R2                    {np.median(r2):.3f}")
print(f"  total variance explained     "
      f"{(r2 * X.var(0)).sum() / X.var(0).sum():.3f}   "
      f"({100 * (r2 * X.var(0)).sum() / X.var(0).sum():.1f}%)")
print(f"  features above R2 0.10       {(r2 > .10).sum()} of {len(feats)}")
print(f"  features above R2 0.25       {(r2 > .25).sum()} of {len(feats)}")

order = np.argsort(r2)[::-1]
print("\n  most explained by gene properties:")
for i in order[:6]:
    print(f"     {r2[i]:.3f}  {feats[i]}")
print("  least explained:")
for i in order[-4:]:
    print(f"     {r2[i]:.3f}  {feats[i]}")

# For scale: how much does the assay's own magnitude basis explain?
keep = [c for c in feats if c not in MAG_OVERALL and c not in DEGEN]
Xs = W[keep].to_numpy(float)[ok]
r2_mag = var_explained(Xs, MAG)
print(f"\n  for scale, the 11-feature magnitude basis explains "
      f"{(r2_mag * Xs.var(0)).sum() / Xs.var(0).sum():.1%} of the non-magnitude features")

# And the retained components
A = np.column_stack([np.ones(ok.sum()), MAG])
beta, *_ = np.linalg.lstsq(A, Xs, rcond=None)
D = PCA(n_components=18, random_state=0).fit_transform(Xs - A @ beta)
r2_pc = var_explained(D, C)
print("\n" + "=" * 78)
print("THE SAME, ON THE RETAINED COMPONENTS (magnitude-corrected substrate)")
print("=" * 78)
for j in range(8):
    bar = "#" * int(r2_pc[j] * 100)
    print(f"  shape-PC{j + 1:<3} R2 {r2_pc[j]:.3f}  {bar}")
print(f"  mean over 18 components      {r2_pc.mean():.3f}   ({r2_pc.mean() * 100:.1f}%)")
print(f"  components above R2 0.10     {(r2_pc > .10).sum()} of 18")

print("\n" + "=" * 78)
print(f"HEADLINE: gene length, expression and CpG density explain "
      f"{100 * (r2 * X.var(0)).sum() / X.var(0).sum():.0f}% of the feature space")
print("=" * 78)
