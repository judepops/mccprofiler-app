"""WHERE in the space does each external set find its displacement?

The effect size says how far a set is displaced. This asks something different
and, for the super-enhancer null, more legible: which component does it lean on?

A set with real architectural structure should find its maximal displacement in
a well-measured, high-variance component. A set with no architectural signature
has no preferred direction, so its largest lean lands wherever noise happens to
be largest, which is a low-variance, poorly-reproducing component.

That converts "d = -0.16" into "no preferred direction anywhere in the
well-measured part of the space", which is a qualitative statement, easier to
draw, and harder to argue with than an effect size close to zero.

Reports, for every reference set on the magnitude-corrected substrate:
  - the component carrying its largest displacement, and that component's rank
  - the variance that component explains
  - the share of that component's loading on features that reproduce (rho > 0.7)
  - the effect size and permutation p, for reference

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
TRUST = 0.7
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()

rep = s.table("reproducibility_pairs")
rho = {}
for f, sub in rep.groupby("feature"):
    if len(sub) > 10 and sub["gw"].std() > 0 and sub["immune"].std() > 0:
        rho[f] = float(np.corrcoef(sub["gw"].rank(), sub["immune"].rank())[0, 1])

mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
pca = PCA(n_components=18, random_state=0).fit(X - A @ beta)
D = pca.transform(X - A @ beta)
D = (D - D.mean(0)) / D.std(0)
n = len(D)

# Trust per component, computed once.
trust = []
for j in range(D.shape[1]):
    load = pd.Series(pca.components_[j], index=keep)
    m = load[[f for f in keep if f in rho]]
    w = (m ** 2) / (m ** 2).sum()
    trust.append(float((w * (pd.Series({f: rho[f] for f in m.index}) > TRUST)).sum()))

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
sets = sorted(sizes[sizes >= S.MIN_GROUP_N].index)

rows = []
for g in sets:
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in members for v in sym), bool, n)
    k = int(hit.sum())
    if k < S.MIN_GROUP_N:
        continue
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    rows.append(dict(group=g, n=k, axis=j + 1,
                     var_pct=pca.explained_variance_ratio_[j] * 100,
                     trust=trust[j] * 100,
                     d=float((a_.mean() - b_.mean()) / sp) if sp else 0.0,
                     p=(np.sum(null >= obs) + 1) / (N_PERM + 1)))

df = pd.DataFrame(rows).sort_values("axis")
print("=" * 88)
print("WHERE EACH SET FINDS ITS DISPLACEMENT, magnitude-corrected substrate")
print("=" * 88)
print(f"{'set':<34}{'n':>6}{'axis':>7}{'var%':>8}{'trust':>8}{'d':>8}{'p':>9}")
for r in df.itertuples():
    mark = "  <-- SE" if "dbSUPER" in r.group else ""
    print(f"{r.group:<34}{r.n:>6}{'sPC' + str(r.axis):>7}{r.var_pct:>7.1f}%"
          f"{r.trust:>7.0f}%{r.d:>8.2f}{r.p:>9.4f}{mark}")

sig = df[df["p"] < 0.05]
print(f"\nOf the {len(sig)} sets that clear p < 0.05:")
print(f"  {int((sig['axis'] <= 4).sum())} find their displacement in the top 4 components")
print(f"  median axis rank {sig['axis'].median():.0f}, "
      f"median trust {sig['trust'].median():.0f}%")
ns = df[df["p"] >= 0.05]
if len(ns):
    print(f"Of the {len(ns)} that do not:")
    print(f"  median axis rank {ns['axis'].median():.0f}, "
          f"median trust {ns['trust'].median():.0f}%")

se = df[df["group"].str.contains("dbSUPER")]
if len(se):
    r = se.iloc[0]
    worse = int((df["axis"] > r["axis"]).sum())
    print(f"\nSuper-enhancers: axis sPC{int(r['axis'])} "
          f"({r['var_pct']:.1f}% of variance, {r['trust']:.0f}% trusted). "
          f"Only {worse} of {len(df)} sets lean on a lower-ranked component.")

print("\nRank correlation between significance and axis rank:")
print(f"  spearman(-log10 p, axis rank) = "
      f"{np.corrcoef((-np.log10(df['p'])).rank(), df['axis'].rank())[0, 1]:+.2f}")
print("  Negative means: the more displaced a set is, the higher-variance the")
print("  component it displaces along. Sets with no architecture lean on noise.")
