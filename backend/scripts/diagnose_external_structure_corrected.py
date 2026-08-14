"""Does ANY external set separate once amount is properly removed?

`diagnose_external_structure.py` ran on the raw principal components. B1 then
established that PC1 is the amount axis (r = 0.623 with the MAG_OVERALL basis),
and PC3 is amount-contaminated at -0.45. So every displacement measured there
sits on an amount-loaded axis, and the ordering of the sets tracks how
amount-loaded their strongest axis is. Roadmap_silenced, the largest
displacement, loads on PC1, which reduces "silenced genes are displaced" to
"silenced genes have less signal", true by construction and silent about
architecture.

That makes Aim 3 untested rather than established. This script is the test.

Same procedure as before, but on the amount-corrected substrate: every feature
regressed on the full 11-feature MAG_OVERALL basis and replaced by its residual,
with the basis features themselves dropped, exactly as
`audit/GW/scripts/_shape.py::corrected_shape` does. Also drops the four
degenerate topology features confirmed in A2, three of which are in the
magnitude basis anyway.

Two outcomes, both decisive:

  Some sets still separate  -> those are the real positive control, and the
                              super-enhancer comparison against them is valid.
  Nothing separates         -> the displacement analysis only ever said "amount
                              varies between gene sets", and Aim 3 as written
                              must be withdrawn.

Also tests a mechanism-based housekeeping definition (ribosomal protein genes)
against the expression-based Eisenberg list, since the objection "housekeeping
genes form promoter assemblies so they should separate" is really an objection
about whether Eisenberg is a good proxy for that mechanism.

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
from _shape import MAG_OVERALL

from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
N_PERM = 2000
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
genes = s.genes.set_index("gene_id").loc[W.index]
sym = genes["symbol_key"].to_numpy()

mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
Xr = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, Xr, rcond=None)
Xs = Xr - A @ beta                      # amount-corrected shape substrate

print(f"magnitude basis {len(mag)} features, degenerate dropped {len(DEGEN)}, "
      f"shape substrate {Xs.shape[1]} features")

# Retain the same number of components as the raw analysis for comparability,
# capped by a parallel-analysis ceiling computed on this substrate.
rng0 = np.random.default_rng(0)
Xn = Xs.copy()
for j in range(Xn.shape[1]):
    rng0.shuffle(Xn[:, j])
p_real = PCA(n_components=min(30, Xs.shape[1]), random_state=0).fit(Xs)
p_null = PCA(n_components=min(30, Xs.shape[1]), random_state=0).fit(Xn)
n_real = int((p_real.explained_variance_ratio_ > p_null.explained_variance_ratio_).sum())
print(f"components above the noise ceiling on this substrate: {n_real}")

D = p_real.transform(Xs)[:, :n_real]
D = (D - D.mean(0)) / D.std(0)
n = len(D)

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
sets = sorted(sizes[sizes >= S.MIN_GROUP_N].index)

# A mechanism-based housekeeping set, to test whether the Eisenberg label rather
# than the biology is what fails to separate.
rp = np.array([bool(x) and (x.startswith("RPL") or x.startswith("RPS")) for x in sym])
extra = {"ribosomal_protein_genes (mechanism)": rp}
print(f"ribosomal protein genes in panel: {rp.sum()}")


def displacement(hit):
    k = int(hit.sum())
    if k < S.MIN_GROUP_N:
        return None
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[hit].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(D[hit].mean(0)).argmax())
    a_, b_ = D[hit, j], D[~hit, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    d = float((a_.mean() - b_.mean()) / sp) if sp else 0.0
    from math import erf, sqrt
    ov = 2 * (1 - 0.5 * (1 + erf(abs(d) / (2 * sqrt(2)))))
    return dict(n=k, z=(obs - null.mean()) / null.std(), p=(np.sum(null >= obs) + 1) / (N_PERM + 1),
                dim=f"sPC{j + 1}", d=d, overlap=ov * 100)


rows = []
for g in sets:
    members = set(cm.loc[cm["group"] == g, "symbol_key"])
    r = displacement(np.fromiter((v in members for v in sym), bool, n))
    if r:
        rows.append(dict(group=g, **r))
for lab, hit in extra.items():
    r = displacement(hit)
    if r:
        rows.append(dict(group=lab, **r))

df = pd.DataFrame(rows).sort_values("z", ascending=False)
print("\n" + "=" * 82)
print("DISPLACEMENT ON THE AMOUNT-CORRECTED SHAPE SUBSTRATE")
print("=" * 82)
print(f"{'set':<38}{'n':>6}{'z':>8}{'p':>9}{'dim':>7}{'d':>8}{'overlap':>9}")
for r in df.itertuples():
    star = "***" if r.p < 0.001 else ("**" if r.p < 0.01 else ("*" if r.p < 0.05 else ""))
    print(f"{r.group:<38}{r.n:>6}{r.z:>8.1f}{r.p:>9.4f}{r.dim:>7}"
          f"{r.d:>8.2f}{r.overlap:>8.0f}% {star}")

sig = df[df["p"] < 0.05]
print(f"\n{len(sig)} of {len(df)} sets still displaced after amount correction "
      f"(was 21 of 21 on the raw components).")
big = df[df["d"].abs() >= 0.5]
print(f"{len(big)} reach |d| >= 0.5. Largest |d| = {df['d'].abs().max():.2f}, "
      f"which leaves {df.loc[df['d'].abs().idxmax(), 'overlap']:.0f}% overlap.")

for name in ["Roadmap_silenced", "Eisenberg_HK", "dbSUPER_CD4_SE_TSS_pm50kb",
             "ribosomal_protein_genes (mechanism)"]:
    hit = df[df["group"] == name]
    if len(hit):
        r = hit.iloc[0]
        print(f"  {name:<38} z {r.z:6.1f}  d {r.d:+.2f}  {r.overlap:.0f}% overlap")
