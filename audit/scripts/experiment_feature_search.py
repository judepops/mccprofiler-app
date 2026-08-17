"""How well can ANY feature combination separate the external gene sets?

The question is worth asking and the naive form is not answerable. With 73
features and 19 sets, a search over subsets will always return something that
separates well, because it is selecting on the outcome. "Subset X separates
super-enhancers at AUC 0.8" is not a finding unless you know what subset X would
have achieved against labels that mean nothing.

So every search here is run twice, identically:

    REAL      the true set membership
    NULL      the same membership, shuffled, same set size, same procedure

The interesting quantity is the GAP. A real ceiling of 0.80 against a null
ceiling of 0.78 means the search found almost nothing and the 0.80 is selection.
A real ceiling of 0.80 against a null of 0.62 means there is something there.

Three searches, from cheap and unbiased to expensive and heavily selected:

  1. ALL FEATURES        no selection at all, the honest baseline
  2. RANDOM SUBSETS      many random subsets at each size, best and median.
                         Shows how much the ceiling rises with search breadth
                         rather than with signal.
  3. GREEDY FORWARD      forward selection maximising cross-validated AUC, the
                         most aggressive search and therefore the one whose null
                         matters most.

Separation is measured as cross-validated AUC of a linear discriminant, which
answers "can a linear rule tell this set from the rest" directly, rather than
the centroid-displacement statistic used elsewhere. AUC 0.5 is chance; the
overlap figures quoted in CURRENT_FINDINGS correspond to AUC in the 0.55 to 0.65
range.

Run in cd4env. Takes roughly ten minutes.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import warnings
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from _shape import MAG_OVERALL
from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")
SEED = 0
N_RANDOM = 120           # random subsets per size
SIZES = [2, 4, 8, 16, 32]
GREEDY_K = 12            # forward-selection steps
CV = 4

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sym = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
feats = list(W.columns)

# Amount-corrected substrate, so the search cannot simply rediscover magnitude.
mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in feats if c not in mag]
A = np.column_stack([np.ones(len(W)), W[mag].to_numpy(float)])
X = W[keep].to_numpy(float)
beta, *_ = np.linalg.lstsq(A, X, rcond=None)
X = X - A @ beta
X = (X - X.mean(0)) / X.std(0)
n, p = X.shape
print(f"{n} genes, {p} amount-corrected features ({len(mag)} magnitude features removed)\n")

cm = s.table("cohort_membership")
sizes = cm.groupby("group").size()
SETS = {}
for g in sorted(sizes.index):
    if not (S.MIN_GROUP_N <= sizes[g] <= S.MAX_GROUP_COVERAGE * n):
        continue
    mem = set(cm.loc[cm["group"] == g, "symbol_key"])
    SETS[g] = np.fromiter((v in mem for v in sym), bool, n)
print(f"{len(SETS)} sets within the size floor and coverage ceiling\n")


def auc(cols, y, seed=SEED):
    """Cross-validated AUC of a linear discriminant on these columns."""
    if len(cols) == 0:
        return 0.5
    cvv = StratifiedKFold(CV, shuffle=True, random_state=seed)
    try:
        pr = cross_val_predict(LDA(), X[:, cols], y, cv=cvv, method="predict_proba")[:, 1]
        return float(roc_auc_score(y, pr))
    except Exception:                                            # noqa: BLE001
        return 0.5


def shuffled(y, rng):
    z = np.zeros(n, bool)
    z[rng.choice(n, size=int(y.sum()), replace=False)] = True
    return z


def greedy(y, k=GREEDY_K, seed=SEED):
    """Forward selection maximising cross-validated AUC."""
    chosen, best = [], 0.5
    remaining = list(range(p))
    trail = []
    for _ in range(k):
        scores = [(auc(chosen + [j], y, seed), j) for j in remaining]
        sc, j = max(scores)
        if sc <= best + 1e-4:
            break
        best, chosen = sc, chosen + [j]
        remaining.remove(j)
        trail.append(best)
    return chosen, best, trail


rng = np.random.default_rng(SEED)
rows = []
print(f"{'set':<30}{'all':>7}{'rand best':>11}{'greedy':>9}   "
      f"{'NULL all':>9}{'NULL rand':>11}{'NULL greedy':>12}{'gap':>7}")
print("-" * 98)

for g, y in SETS.items():
    y_null = shuffled(y, rng)

    a_all, n_all = auc(list(range(p)), y), auc(list(range(p)), y_null)

    best_rand = best_rand_null = 0.5
    for size in SIZES:
        for _ in range(N_RANDOM // len(SIZES)):
            cols = rng.choice(p, size=min(size, p), replace=False).tolist()
            best_rand = max(best_rand, auc(cols, y))
            best_rand_null = max(best_rand_null, auc(cols, y_null))

    _, a_greedy, _ = greedy(y)
    _, n_greedy, _ = greedy(y_null)

    gap = a_greedy - n_greedy
    rows.append(dict(group=g, n=int(y.sum()), all=a_all, rand=best_rand, greedy=a_greedy,
                     null_all=n_all, null_rand=best_rand_null, null_greedy=n_greedy, gap=gap))
    flag = "  <--" if gap > .10 else ""
    print(f"{g:<30}{a_all:>7.3f}{best_rand:>11.3f}{a_greedy:>9.3f}   "
          f"{n_all:>9.3f}{best_rand_null:>11.3f}{n_greedy:>12.3f}{gap:>+7.3f}{flag}")

d = pd.DataFrame(rows)
out = Path(__file__).resolve().parent / "output" / "feature_search.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
d.to_csv(out, sep="\t", index=False)

print("\n" + "=" * 98)
print("WHAT THE SEARCH ITSELF MANUFACTURES")
print("=" * 98)
print(f"  using all features, no selection : real {d['all'].median():.3f}  "
      f"null {d['null_all'].median():.3f}   gap {d['all'].median() - d['null_all'].median():+.3f}")
print(f"  best of {N_RANDOM} random subsets  : real {d['rand'].median():.3f}  "
      f"null {d['null_rand'].median():.3f}   gap {d['rand'].median() - d['null_rand'].median():+.3f}")
print(f"  greedy forward selection         : real {d['greedy'].median():.3f}  "
      f"null {d['null_greedy'].median():.3f}   gap {d['greedy'].median() - d['null_greedy'].median():+.3f}")
print(f"\n  search inflation on NULL labels: {d['null_all'].median():.3f} with no selection "
      f"rises to {d['null_greedy'].median():.3f} under greedy search.")
print(f"  That rise is entirely artefact. Any real result must clear it.")

print("\n" + "=" * 98)
print("SETS WHERE THE SEARCH BEATS ITS OWN NULL BY MORE THAN 0.10 AUC")
print("=" * 98)
win = d[d.gap > .10].sort_values("gap", ascending=False)
if len(win):
    for r in win.itertuples():
        print(f"   {r.group:<30} greedy {r.greedy:.3f} vs null {r.null_greedy:.3f}  "
              f"gap {r.gap:+.3f}  (n={r.n})")
else:
    print("   none")
se = d[d.group.str.contains("dbSUPER")]
if len(se):
    r = se.iloc[0]
    print(f"\n   super-enhancers: greedy {r['greedy']:.3f} vs null {r['null_greedy']:.3f}, "
          f"gap {r['gap']:+.3f}")
print(f"\nwrote {out}")
