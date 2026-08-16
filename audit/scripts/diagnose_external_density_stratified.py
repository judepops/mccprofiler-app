"""Does the Aim 3 positive control survive gene-density stratification?

Addendum 3 established that chromatin-state categories separate on the
amount-corrected substrate (bivalent d = -0.58, silenced d = -0.73) while
super-enhancers do not (p = 0.18). The super-enhancer null is only interpretable
if those positives are real, because they are what demonstrates the test can
detect separation at all.

But `gene_desert_bottomQ_density` remains the LARGEST displacement even after
amount correction (z = 20.8), so genomic context is uncorrected and is the
strongest thing in the space. Silenced genes are enriched in gene-poor regions.
So "silenced genes separate on shape" could still be "silenced genes live in
gene deserts", which is the amount trap one level over.

Two independent checks:

  1. RESIDUALISE. Add gene density to the confounder basis alongside the
     11-feature MAG_OVERALL magnitude basis, and re-run displacement.
  2. STRATIFY. Compute the effect within gene-density deciles and pool. This
     asks the same question without assuming the density relationship is linear,
     and gives a retention percentage directly comparable to the
     101.8% archetype retention already reported.

If both retain, the positive control holds and the super-enhancer null stands.
If either collapses, Aim 3 is uninterpretable and must be withdrawn.

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
from _shape import MAG_OVERALL

from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
N_PERM = 2000
N_DECILE = 10
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
FOCUS = ["ChromHMM_bivalent", "Roadmap_silenced", "Lambert_TF",
         "GWAS_immune_hot", "Eisenberg_HK", "dbSUPER_CD4_SE_TSS_pm50kb",
         "gene_desert_bottomQ_density"]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
genes = s.genes.set_index("gene_id").loc[W.index]
sym = genes["symbol_key"].to_numpy()

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
dens = (pd.Series(sym, name="symbol_key").to_frame()
        .merge(gt[["symbol_key", "gene_density_500kb"]], on="symbol_key", how="left")
        ["gene_density_500kb"].to_numpy(float))
have = np.isfinite(dens)
print(f"gene density available for {have.sum()} of {len(dens)} genes")
dens_filled = np.where(have, dens, np.nanmedian(dens))


def substrate(extra_confounders=None):
    """MAG_OVERALL-corrected shape, optionally also residualised on extras."""
    mag = [m for m in MAG_OVERALL if m in W.columns]
    keep = [c for c in W.columns if c not in mag and c not in DEGEN]
    cols = [np.ones(len(W)), W[mag].to_numpy(float)]
    if extra_confounders is not None:
        cols.append(extra_confounders)
    A = np.column_stack(cols)
    X = W[keep].to_numpy(float)
    beta, *_ = np.linalg.lstsq(A, X, rcond=None)
    R = X - A @ beta
    p = PCA(n_components=18, random_state=0).fit(R)
    D = p.transform(R)
    return (D - D.mean(0)) / D.std(0)


def displace(D, hit, n_perm=N_PERM):
    k = int(hit.sum())
    n = len(D)
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(n_perm)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    return dict(n=k, z=(obs - null.mean()) / null.std(),
                p=(np.sum(null >= obs) + 1) / (n_perm + 1),
                d=float((a_.mean() - b_.mean()) / sp) if sp else 0.0, dim=j)


def stratified_d(D, hit, strata, j):
    """Effect on component j, computed within strata and pooled by weight.

    Makes no assumption that the density relationship is linear, which the
    residualisation does.
    """
    num = den = 0.0
    for sname in np.unique(strata):
        m = strata == sname
        a_, b_ = D[m & hit, j], D[m & ~hit, j]
        if len(a_) < 3 or len(b_) < 3:
            continue
        sp = np.sqrt(((len(a_) - 1) * a_.var(ddof=1) + (len(b_) - 1) * b_.var(ddof=1))
                     / max(len(a_) + len(b_) - 2, 1))
        if sp <= 0:
            continue
        w = len(a_)
        num += w * (a_.mean() - b_.mean()) / sp
        den += w
    return num / den if den else np.nan


D_amt = substrate()
D_both = substrate(dens_filled.reshape(-1, 1))
decile = pd.qcut(dens_filled, N_DECILE, labels=False, duplicates="drop")

cm = s.table("cohort_membership")
rows = []
for g in FOCUS:
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, len(sym))
    if hit.sum() < S.MIN_GROUP_N:
        continue
    a = displace(D_amt, hit)
    b = displace(D_both, hit)
    sd = stratified_d(D_amt, hit, decile, a["dim"])
    rows.append(dict(group=g, n=a["n"],
                     d_amt=a["d"], p_amt=a["p"],
                     d_both=b["d"], p_both=b["p"],
                     d_strat=sd,
                     retained=100 * abs(sd) / abs(a["d"]) if a["d"] else np.nan))

df = pd.DataFrame(rows)
print("\n" + "=" * 92)
print("GENE-DENSITY CHECK ON THE AIM 3 POSITIVE CONTROLS")
print("=" * 92)
print(f"{'set':<34}{'n':>5}{'d (amt-corr)':>14}{'d (+density)':>14}"
      f"{'d (stratified)':>16}{'retained':>10}")
for r in df.itertuples():
    print(f"{r.group:<34}{r.n:>5}{r.d_amt:>13.2f}{r.d_both:>14.2f}"
          f"{r.d_strat:>16.2f}{r.retained:>9.0f}%")

print("\np-values, amount-corrected vs amount+density-corrected:")
for r in df.itertuples():
    flag = "" if r.p_both < 0.05 else "   <- no longer significant"
    print(f"  {r.group:<34} {r.p_amt:.4f}  ->  {r.p_both:.4f}{flag}")

print("\nReading: if the positives (bivalent, silenced) retain most of their")
print("effect under BOTH residualisation and stratification, the Aim 3 positive")
print("control holds and the super-enhancer null is interpretable. If they")
print("collapse, the separation was genomic context and Aim 3 must be withdrawn.")
