"""The own-TSS vs other-locus split, applied to EVERY element class and EVERY set.

`diagnose_promoter_assembly.py` tested one class (promoter) against one set
(Eisenberg_HK) and found something. Testing one combination and reporting it is
the selective-search error that `experiment_feature_search.py` exists to catch,
so this runs the full grid and corrects for it.

Three element classes x three measures (own / other / total) x every usable set.
Every p-value is Benjamini-Hochberg corrected across the whole grid, and the
pre-specified test (Eisenberg_HK x promoter x other) is reported separately
because it was specified before the grid was run and does not owe multiplicity.

WHY THE SPLIT IS DEFINED ON THE GENE'S OWN TSS, NOT THE VIEWPOINT. The MCC
matrix is anchored on the midpoint of the captured region, which is a per-peak
identifier and NOT the canonical TSS; about 2% of genes differ by hundreds of kb
(RERE by 292 kb). Splitting on the viewpoint would therefore mislabel exactly the
genes where the distinction matters most. TSS comes from `01_tss_regions.bed`.

Run in cd4env. About ten minutes.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from _shape import MAG_OVERALL
from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)

PEAKS = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
         "scripts_cleaned/process/output_resgw_10000/utils/peaks/annotated.tsv")
TSS_BED = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/"
           "01_tss_regions.bed")
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
OWN_KB = 5
N_PERM = 2000
CLASSES = ["promoter", "enhancer", "ctcf"]

pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()
bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key").set_index("symbol_key")
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sk_all = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()

pk["own_tss"] = pk["symbol_key"].map(bed["tss"])
pk = pk[np.isfinite(pk["own_tss"])].copy()
pk["d_own"] = (pk["peak_midpoint"] - pk["own_tss"]).abs()
pk["is_own"] = pk["d_own"] <= OWN_KB * 1000

tot = pk.groupby("symbol_key")["peak_max"].sum()
cols = {}
for cls in CLASSES:
    sub = pk[pk["re"] == cls]
    o = sub[sub["is_own"]].groupby("symbol_key")["peak_max"].sum()
    t = sub[~sub["is_own"]].groupby("symbol_key")["peak_max"].sum()
    cols[f"{cls}_own"] = o
    cols[f"{cls}_other"] = t
F = pd.DataFrame({"tot": tot})
for k, v in cols.items():
    F[k] = v
F = F.fillna(0.0)
for cls in CLASSES:
    F[f"{cls}_own"] /= F["tot"].replace(0, np.nan)
    F[f"{cls}_other"] /= F["tot"].replace(0, np.nan)
    F[f"{cls}_total"] = F[f"{cls}_own"] + F[f"{cls}_other"]
MEAS = [f"{c}_{w}" for c in CLASSES for w in ("own", "other", "total")]

V = pd.DataFrame({m: pd.Series(sk_all).map(F[m]).to_numpy(float) for m in MEAS})
ok = V.notna().all(axis=1).to_numpy()
V, sk = V[ok].reset_index(drop=True), sk_all[ok]
Wm = W.reset_index(drop=True)[ok].reset_index(drop=True)

print(f"{len(sk)} genes, {len(MEAS)} measures")
for cls in CLASSES:
    sub = pk[pk["re"] == cls]
    print(f"  {cls:<10} {len(sub):>6} peaks, {sub['is_own'].mean():>6.1%} at the "
          f"gene's own TSS (within {OWN_KB} kb)")

# Full control set, as in section 6d: magnitude, density, length, expression, CpG.
A = [np.ones(len(sk))] + [Wm[c].to_numpy(float) for c in MAG_OVERALL if c in Wm.columns]
for col, lg in [("gene_density_500kb", 0), ("log10_gene_length", 0),
                ("gtex_whole_blood_tpm", 1), ("cpg_oe_ratio", 0)]:
    v = pd.Series(sk).map(gt[col]).to_numpy(float)
    if lg:
        v = np.log1p(v)
    A.append(np.where(np.isfinite(v), v, np.nanmedian(v)))
A = np.nan_to_num(np.column_stack(A), nan=0.0)

Y = V.to_numpy(float)
beta, *_ = np.linalg.lstsq(A, Y, rcond=None)
R = Y - A @ beta
R = (R - R.mean(0)) / R.std(0)

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
SETS = {}
for g in sorted(sizes.index):
    mem = set(cm.loc[cm["group"] == g, "symbol_key"])
    m = np.fromiter((v in mem for v in sk), bool, len(sk))
    if S.MIN_GROUP_N <= m.sum() <= S.MAX_GROUP_COVERAGE * len(sk):
        SETS[g] = m
print(f"  {len(SETS)} usable sets -> {len(SETS) * len(MEAS)} tests\n")

rows = []
idx = np.arange(len(sk))
for g, m in SETS.items():
    for j, meas in enumerate(MEAS):
        r = R[:, j]
        d = (r[m].mean() - r[~m].mean()) / np.sqrt((r[m].var() + r[~m].var()) / 2)
        nul = np.empty(N_PERM)
        for i in range(N_PERM):
            p = rng.permutation(idx)[:m.sum()]
            q = np.setdiff1d(idx, p)
            nul[i] = (r[p].mean() - r[q].mean()) / np.sqrt(
                (r[p].var() + r[q].var()) / 2)
        rows.append(dict(group=g, measure=meas, n=int(m.sum()), d=d,
                         p=max((np.abs(nul) >= abs(d)).mean(), 1 / N_PERM)))

df = pd.DataFrame(rows).sort_values("p")
mm = len(df)
df["rank"] = np.arange(1, mm + 1)
# Benjamini-Hochberg. The cumulative minimum runs from the LARGEST p downward,
# not from the smallest upward: q_i = min over j >= i of (p_j * m / j). A forward
# cummin on p-ascending data drags every later q down to the smallest early
# value and reports essentially everything as significant.
df["q"] = (df["p"] * mm / df["rank"])[::-1].cummin()[::-1].clip(upper=1.0)

out = Path(__file__).resolve().parent / "output" / "own_vs_other_all.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)

print("=" * 88)
print(f"SURVIVING BENJAMINI-HOCHBERG FDR 0.05 ACROSS ALL {mm} TESTS")
print("=" * 88)
sig = df[df["q"] < 0.05]
if len(sig):
    print(f"{'set':<30}{'measure':<18}{'n':>5}{'d':>8}{'p':>9}{'q':>9}")
    for r in sig.itertuples():
        print(f"{r.group:<30}{r.measure:<18}{r.n:>5}{r.d:>+8.3f}"
              f"{r.p:>9.4f}{r.q:>9.4f}")
else:
    print("  nothing survives correction")

print("\n" + "=" * 88)
print("THE PRE-SPECIFIED TEST, which does not owe multiplicity")
print("=" * 88)
ps = df[(df.group == "Eisenberg_HK") & (df.measure == "promoter_other")]
if len(ps):
    r = ps.iloc[0]
    print(f"  Eisenberg_HK x promoter_other   d = {r['d']:+.3f}  "
          f"p = {r['p']:.4f}  (BH q = {r['q']:.4f})")
    print(f"  Its q is shown only for context. It was specified before this grid")
    print(f"  was run, so the p stands on its own; the q answers the different")
    print(f"  question of whether it would survive as a discovery FROM the grid.")

print("\n" + "=" * 88)
print("HOW MUCH OF THE GRID IS NOISE")
print("=" * 88)
print(f"  tests run                : {mm}")
print(f"  nominally p < 0.05       : {(df.p < 0.05).sum()}  "
      f"(expected by chance {0.05 * mm:.1f})")
print(f"  surviving BH q < 0.05    : {(df.q < 0.05).sum()}")
print("\n  own-TSS measures should be mostly null. If they are not, the split is")
print("  picking up viewpoint placement rather than biology:")
for cls in CLASSES:
    o = df[df.measure == f"{cls}_own"]
    print(f"     {cls}_own    {(o.p < 0.05).sum():>2} of {len(o)} sets nominal, "
          f"{(o.q < 0.05).sum()} after correction")
print(f"\nwrote {out}")
