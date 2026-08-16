"""The checks that were run interactively and never saved as scripts.

Written 2026-08-16 so that changing the feature space does not lose them. Each
of these changed a conclusion when it was first run, and each would otherwise
have to be reconstructed from a commit message.

  1. PANEL COVERAGE. Which reference sets exceed MAX_GROUP_COVERAGE and so
     cannot give an interpretable set-versus-rest contrast.
  2. TAXONOMY ROBUSTNESS. Leave-one-out, reassignment of the arguable sets, and
     a random-relabelling null for the DIRECT versus OUTPUT contrast.
  3. LAMBERT DBD SPLIT. Whether the transcription-factor signal is uniform
     across DNA-binding-domain families or comes from one of them.
  4. chr19 VERSUS DENSITY. Whether the chromosome effect is gene density
     wearing a chromosome label.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/chromprofiler/benchmarking/src")

import warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
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
cm = s.table("cohort_membership")
n_panel = len(W)
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
meta = pd.DataFrame({"symbol_key": sym}).merge(
    gt[["symbol_key", "tss_chrom", "gene_density_500kb"]].drop_duplicates("symbol_key"),
    on="symbol_key", how="left")
chrom = meta["tss_chrom"].to_numpy()
dens = meta["gene_density_500kb"].to_numpy(float)
dens = np.where(np.isfinite(dens), dens, np.nanmedian(dens))

mag = [m for m in MAG_OVERALL if m in W.columns]
keep = [c for c in W.columns if c not in mag and c not in DEGEN]
X = W[keep].to_numpy(float)


def substrate(extra=None):
    cols = [np.ones(n_panel), W[mag].to_numpy(float)]
    if extra is not None:
        cols.append(extra)
    A = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(A, X, rcond=None)
    D = PCA(n_components=18, random_state=0).fit_transform(X - A @ beta)
    return (D - D.mean(0)) / D.std(0)


def disp(D, mask):
    n = len(D)
    k = int(mask.sum())
    if k < S.MIN_GROUP_N:
        return None
    rng = np.random.default_rng(0)
    obs = float(np.linalg.norm(D[mask].mean(0)))
    null = np.array([np.linalg.norm(D[rng.choice(n, size=k, replace=False)].mean(0))
                     for _ in range(N_PERM)])
    j = int(np.abs(D[mask].mean(0)).argmax())
    a_, b_ = D[mask, j], D[~mask, j]
    sp = np.sqrt(((k - 1) * a_.var(ddof=1) + (n - k - 1) * b_.var(ddof=1)) / (n - 2))
    return k, abs((a_.mean() - b_.mean()) / sp), (np.sum(null >= obs) + 1) / (N_PERM + 1)


# --- 1. panel coverage -----------------------------------------------------
print("=" * 78)
print(f"1. PANEL COVERAGE  (ceiling {S.MAX_GROUP_COVERAGE:.0%}, floor {S.MIN_GROUP_N})")
print("=" * 78)
counts = cm.groupby("group").size().sort_values(ascending=False)
over = []
for g, v in counts.items():
    if v < S.MIN_GROUP_N:
        continue
    cov = v / n_panel
    if cov > S.MAX_GROUP_COVERAGE:
        over.append(g)
        print(f"   {g:<34}{v:>6}  {cov:>5.0%} of panel, rest = {n_panel - v}  EXCLUDED")
print(f"   {len(over)} of {(counts >= S.MIN_GROUP_N).sum()} sets exceed the ceiling")

# --- 2. taxonomy robustness ------------------------------------------------
D = substrate()
DEFN_DIRECT = {"Lambert_TF", "CpG_island_promoter", "phastCons_2kb_topQ",
               "GWAS_immune_hot", "GWAS_total_topQ", "GWAS_immune_any",
               "gnomAD_pLI_topQ", "Roadmap_silenced", "ChromHMM_bivalent",
               "ChromHMM_active_TSS"}
EXCLUDE = {"dbSUPER_CD4_SE_TSS_pm50kb", "gene_desert_bottomQ_density"}
eff = {}
for g in counts[counts >= S.MIN_GROUP_N].index:
    if g in EXCLUDE or counts[g] > S.MAX_GROUP_COVERAGE * n_panel:
        continue
    mem = set(cm.loc[cm["group"] == g, "symbol_key"])
    r = disp(D, np.fromiter((v in mem for v in sym), bool, n_panel))
    if r:
        eff[g] = r[1]
A = [v for g, v in eff.items() if g in DEFN_DIRECT]
B = [v for g, v in eff.items() if g not in DEFN_DIRECT]
base = mannwhitneyu(A, B, alternative="greater")[1]
print("\n" + "=" * 78)
print("2. TAXONOMY ROBUSTNESS, DIRECT vs OUTPUT")
print("=" * 78)
print(f"   baseline  DIRECT n={len(A)} med {np.median(A):.2f} | "
      f"OUTPUT n={len(B)} med {np.median(B):.2f} | p = {base:.4f}")
worst = 0.0
for g in eff:
    a2 = [v for k, v in eff.items() if k in DEFN_DIRECT and k != g]
    b2 = [v for k, v in eff.items() if k not in DEFN_DIRECT and k != g]
    worst = max(worst, mannwhitneyu(a2, b2, alternative="greater")[1])
print(f"   leave-one-out, worst case p = {worst:.4f}")
ARG = ["CpG_island_promoter", "phastCons_2kb_topQ", "gnomAD_pLI_topQ", "bio_dev_TF_k3"]
flip = {g: (g not in DEFN_DIRECT) for g in ARG if g in eff}
a3 = [v for g, v in eff.items() if (g in DEFN_DIRECT) != flip.get(g, False)]
b3 = [v for g, v in eff.items() if (g in DEFN_DIRECT) == flip.get(g, False)]
if a3 and b3:
    print(f"   all arguable sets reassigned, p = "
          f"{mannwhitneyu(a3, b3, alternative='greater')[1]:.4f}")
rng = np.random.default_rng(0)
ks = list(eff)
hits = sum(mannwhitneyu([eff[k] for k in rng.permutation(ks)[:len(A)]],
                        [eff[k] for k in rng.permutation(ks)[len(A):]],
                        alternative="greater")[1] <= base for _ in range(2000))
print(f"   random relabellings reaching the observed p: {hits}/2000 = {hits / 2000:.3f}")

# --- 3. Lambert DBD split --------------------------------------------------
print("\n" + "=" * 78)
print("3. LAMBERT DBD SPLIT: is the TF signal uniform across families?")
print("=" * 78)
try:
    from benchmarking.loaders import load_lambert
    _, dbd = load_lambert()
    c2h2 = {g for g, d in dbd.items() if isinstance(d, str) and d == "C2H2 ZF"}
    allTF = set(cm.loc[cm["group"] == "Lambert_TF", "symbol_key"])
    for lab, mem in [("Lambert_TF (all)", allTF), ("  C2H2 ZF", allTF & c2h2),
                     ("  non-C2H2", allTF - c2h2)]:
        r = disp(D, np.fromiter((v in mem for v in sym), bool, n_panel))
        print(f"   {lab:<22}" + (f"n={r[0]:>4}  |d| {r[1]:.2f}  p={r[2]:.4f}"
                                 if r else "below the floor"))
except Exception as e:                                        # noqa: BLE001
    print(f"   unavailable: {e}")

# --- 4. chr19 versus density ----------------------------------------------
print("\n" + "=" * 78)
print("4. chr19: a locus effect, or gene density wearing a chromosome label?")
print("=" * 78)
is19 = chrom == "chr19"
print(f"   density: chr19 median {np.median(dens[is19]):.1f} vs "
      f"{np.median(dens[~is19]):.1f} elsewhere, "
      f"d = {(dens[is19].mean() - dens[~is19].mean()) / dens.std():+.2f}")
for lab, D_ in [("magnitude only", D), ("magnitude + density", substrate(dens.reshape(-1, 1)))]:
    r = disp(D_, is19)
    print(f"   chr19 {lab:<24} |d| {r[1]:.2f}  p={r[2]:.4f}")
j = int(np.abs(D[is19].mean(0)).argmax())
dec = pd.qcut(dens, 10, labels=False, duplicates="drop")
num = den = 0.0
for q in np.unique(dec):
    m_ = dec == q
    a_, b_ = D[m_ & is19, j], D[m_ & ~is19, j]
    if len(a_) < 3 or len(b_) < 3:
        continue
    sp = np.sqrt(((len(a_) - 1) * a_.var(ddof=1) + (len(b_) - 1) * b_.var(ddof=1))
                 / max(len(a_) + len(b_) - 2, 1))
    if sp > 0:
        num += len(a_) * (a_.mean() - b_.mean()) / sp
        den += len(a_)
base19 = disp(D, is19)[1]
print(f"   chr19 {'density-decile stratified':<24} |d| {abs(num / den):.2f}  "
      f"(retains {100 * abs(num / den) / base19:.0f}%)")
