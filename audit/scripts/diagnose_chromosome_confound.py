"""Is a gene set's displacement explained by which chromosome its genes sit on?

Found 2026-08-16 while testing the Lambert DBD subfamilies, and it is a
confound class that had not been controlled at all.

Splitting Lambert TFs by DNA-binding domain showed the entire TF effect coming
from C2H2 zinc fingers (|d| 0.75) with non-C2H2 TFs at nothing (p = 0.52). But
C2H2-ZF genes are 38% chr19 in this panel against 6% elsewhere, because of the
KRAB-ZNF arrays, and chr19 membership BY ITSELF displaces at |d| = 0.63,
p = 0.0005. Other chromosomes give 0.23 to 0.26.

So chr19 has distinctive contact architecture and any set enriched there
inherits it. Two ways to ask whether a set survives, and they answer different
questions, so both are reported:

  EXCLUDE   drop chr19 genes and re-test. Asks whether the set separates
            among genes elsewhere in the genome.
  CONTROL   add chromosome as a one-hot covariate alongside the magnitude
            basis. Asks whether the set differs from its own chromosomal
            neighbours.

A set that fails EXCLUDE but passes CONTROL is displaced relative to its
neighbours while its overall displacement is substantially locus-driven. That
is the Lambert_TF case and it is why Lambert_TF is no longer the primary
positive control.

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
DEGEN = ["mean_degree", "mean_degree_raw", "frac_active_pairs", "n_isolates_raw"]
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
chrom = (pd.DataFrame({"symbol_key": sym})
         .merge(gt[["symbol_key", "tss_chrom"]], on="symbol_key", how="left")
         ["tss_chrom"].to_numpy())

mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
X = W[keep].to_numpy(float)
onehot = pd.get_dummies(pd.Series(chrom)).to_numpy(float)[:, 1:]


def substrate(with_chrom: bool):
    cols = [np.ones(len(W)), W[mag].to_numpy(float)]
    if with_chrom:
        cols.append(onehot)
    A = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(A, X, rcond=None)
    D = PCA(n_components=18, random_state=0).fit_transform(X - A @ beta)
    return (D - D.mean(0)) / D.std(0)


def disp(D, mask, sub=None):
    d_ = D if sub is None else D[sub]
    m_ = mask if sub is None else mask[sub]
    n = len(d_)
    k = int(m_.sum())
    if k < S.MIN_GROUP_N:
        return None
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(d_[m_].mean(0)))
    null = np.array([np.linalg.norm(d_[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(d_[m_].mean(0)).argmax())
    a_, b_ = d_[m_, j], d_[~m_, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    return k, abs((a_.mean() - b_.mean()) / sp), (np.sum(null >= obs) + 1) / (N_PERM + 1)


D_plain, D_chrom = substrate(False), substrate(True)
not19 = chrom != "chr19"
cm = s.table("cohort_membership")

print("=" * 96)
print("CHROMOSOME AS A CONFOUND")
print("=" * 96)
print("chromosome membership alone, on the magnitude-corrected substrate:")
for c in sorted(set(chrom[pd.notna(chrom)]), key=lambda x: -(chrom == x).sum())[:6]:
    r = disp(D_plain, chrom == c)
    if r:
        flag = "   <-- outlier" if r[1] > 0.4 else ""
        print(f"   {c:<8} n={r[0]:>4}  |d| {r[1]:.2f}  p={r[2]:.4f}{flag}")

print(f"\n{'set':<30}{'all genes':>22}{'excluding chr19':>24}"
      f"{'chrom controlled':>22}{'chr19':>8}")
sizes = cm.groupby("group").size()
for g in sorted(sizes[sizes >= S.MIN_GROUP_N].index):
    mem = set(cm.loc[cm["group"] == g, "symbol_key"])
    hit = np.fromiter((v in mem for v in sym), bool, len(sym))
    if hit.sum() > S.MAX_GROUP_COVERAGE * len(sym):
        continue
    a, b, c_ = disp(D_plain, hit), disp(D_plain, hit, not19), disp(D_chrom, hit)
    f = lambda r: f"|d|{r[1]:.2f} p{r[2]:.4f}" if r else "below floor"
    share = 100 * (chrom[hit] == "chr19").mean()
    warn = "  <-- locus-driven" if (a and b and a[2] < .05 <= b[2]) else ""
    print(f"{g:<30}{f(a):>22}{f(b):>24}{f(c_):>22}{share:>7.0f}%{warn}")

print("\nEXCLUDE asks whether the set separates elsewhere in the genome.")
print("CONTROL asks whether it differs from its own chromosomal neighbours.")
print("A set failing EXCLUDE but passing CONTROL is substantially locus-driven.")
