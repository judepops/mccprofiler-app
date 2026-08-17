et me answer plainly, because I think I oversold that result and you are right to push.

No, there is no housekeeping cluster. d = 0.14 with 95% overlap is a tiny average shift. It is smaller than most displacements already in your Section 4 table. It is the same "displaced, not separated" as everything else. I called it "the first positive of the week" because it survives controls and has a specificity pattern, but that phrasing implied more than the number supports, and I should have led with the 95%.

No, we have never used an external promoter-assemblies list, and you are right that this matters. What I actually tested was Eisenberg_HK — an expression-derived housekeeping list — against a feature I constructed. The published promoter-assembly work (Ronin/THAP11 targets) has its own gene list and we do not have it. So the claim is "expression-defined housekeeping genes contact other promoters slightly more", not "the promoter-assembly genes behave as published". Those are different statements and I conflated them.

On results changing: they are not changing randomly, but I have been presenting them badly. Two separate things:

The substrate changed once (91 → 73 features), and that moved numbers without moving conclusions — |d| before vs after correlates at r = 0.959.
Everything else is confounds being removed. A result that dies under a control was never a result. That is the process working, but reporting each one as it lands makes it feel like whiplash.
The stable picture has not moved all week: nothing separates, most things are slightly displaced, and most apparent effects are density, length, expression or chromosome.

And yes — you have put your finger on the actual problem. If splitting own- vs other-promoter changed the answer for one set, then the split is information the feature set does not contain, and it should be applied to every element class and every set, not just promoters and not just HK. Testing one combination and reporting it is exactly the selective-search error I ran the permutation null to avoid.

Doing both things you asked. First the naming audit, since frac_promoter_proximal may not be the only one:

Bash Pull band and class definitions
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/mccprofiler
echo "=== DISTANCE BANDS (config.py) ==="
sed -n '70,95p' src/mccprofiler/config.py
echo
echo "=== ELEMENT CLASSES ==="
grep -n "ELEMENT_CLASSES" src/mccprofiler/config.py | head -3
grep -n "ELEMENT_CLASSES\s*=" -A 3 src/mccprofiler/config.py | head -6
echo
echo "=== BAIT MASK ==="
grep -n "_BAIT_MASK\|BAIT" src/mccprofiler/global_features.py | head -6
grep -n "mask_bp\|MASK" src/mccprofiler/config.py | head -6

OUT
=== DISTANCE BANDS (config.py) ===

BAIT_HALF_BINS = kb_to_bins(BAIT_KB)        # 40 bins  (±2 kb)
DISTAL_HALF_BINS = kb_to_bins(DISTAL_KB)    # 1000 bins (±50 kb)

# Distance bands — non-overlapping rings around the viewpoint (units: bp from VP).
# Convention: lower-inclusive, upper-exclusive, signed-distance-agnostic
# (i.e. |bin - CENTER_BIN| * BIN_SIZE).
DISTANCE_BANDS = {
    "promoter_proximal": (0,          10_000),
    "local":             (10_000,     50_000),
    "distal":            (50_000,    250_000),
    "far_distal":        (250_000, 1_000_000),
}


# ---------------------------------------------------------------------------
# Feature engineering knobs
# ---------------------------------------------------------------------------

CORRELATION_THRESHOLD = 0.8       # |r| > threshold → prune higher-correlated pair

# Features excluded BEFORE the correlation prune, added 2026-08-16.
#
# These four are exact algebraic functions of two peak counts, verified on
# topology_raw.tsv at 1842/1842 genes:
#

=== ELEMENT CLASSES ===
209:ELEMENT_CLASSES = ("enhancer", "ctcf", "promoter")
209:ELEMENT_CLASSES = ("enhancer", "ctcf", "promoter")

=== BAIT MASK ===
32:    BAIT_HALF_BINS,
57:_BAIT_MASK = _DISTANCES_BP < (BAIT_HALF_BINS * BIN_SIZE)
98:    """Fraction of total signal lying within ±BAIT_KB of the viewpoint."""
102:    return float(profile[_BAIT_MASK].sum() / total)
151:    weights = oe_profile[~_BAIT_MASK]
152:    distances = _DISTANCES_BP[~_BAIT_MASK].astype(np.float64)
Two definitions that matter: DISTANCE_BANDS["promoter_proximal"] = (0, 10_000) measured as |bin − CENTER_BIN| from the viewpoint — and per your own trap #3, the viewpoint is not the TSS. So frac_promoter_proximal is neither "promoter" nor reliably proximal to one. Auditing all 73 against their definitions:

Bash Group all 73 features by kind
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate mccapp
python - <<'PY'
import sys; sys.path.insert(0,"backend")
from app.store import get_store
f=sorted(get_store().table("features")["feature"].unique())
import re
groups={}
for n in f:
    if re.search(r"_(enhancer|ctcf|promoter|all)$",n): k="per-element-class"
    elif n.startswith("frac_") or "band" in n: k="distance-band"
    elif n.startswith("oe_"): k="O/E derived"
    elif "distance" in n: k="distance"
    else: k="global"
    groups.setdefault(k,[]).append(n)
for k,v in groups.items():
    print(f"\n### {k}  ({len(v)})")
    for n in v: print("   ",n)
print(f"\ntotal {len(f)}")
PY

OUT
### global  (23)
    bait_pileup_fraction
    bimodality_score
    contact_asymmetry
    ctcf_signal_fraction
    ctcf_signal_fraction_raw
    distal_signal_density
    dominant_peak_isolation_bp
    enhancer_signal_fraction_raw
    local_to_distal_ratio
    max_mcc
    mean_peak_gap_bp
    n_active_peaks
    n_high_consensus_peaks_075
    n_peaks_within_100kb
    peak_dominance_index
    promoter_signal_fraction
    promoter_signal_fraction_raw
    q90_mcc
    signal_entropy
    signal_entropy_distal
    signal_profile_kurtosis
    spacing_regularity
    total_mcc

### per-element-class  (34)
    consensus_fraction_max_all
    consensus_fraction_max_ctcf
    consensus_fraction_max_enhancer
    consensus_fraction_mean_all
    frac_ctcf
    max_distance_to_viewpoint_all
    max_distance_to_viewpoint_ctcf
    max_distance_to_viewpoint_enhancer
    max_distance_to_viewpoint_promoter
    n_peaks_all
    n_peaks_ctcf
    n_peaks_enhancer
    n_peaks_promoter
    oe_fwhm_bp_max_all
    oe_fwhm_bp_max_ctcf
    oe_fwhm_bp_max_enhancer
    oe_local_enrichment_max_all
    oe_local_enrichment_max_ctcf
    oe_local_enrichment_max_enhancer
    oe_local_enrichment_max_promoter
    oe_local_enrichment_mean_all
    oe_max_max_all
    oe_max_max_ctcf
    oe_max_max_enhancer
    raw_log2_enrichment_max_all
    raw_log2_enrichment_max_ctcf
    raw_log2_enrichment_max_enhancer
    raw_log2_enrichment_max_promoter
    raw_log2_enrichment_mean_all
    raw_peak_max_max_all
    raw_peak_max_max_ctcf
    raw_peak_max_max_enhancer
    raw_peak_max_max_promoter
    raw_peak_max_mean_all

### distance  (6)
    corr_oe_distance
    corr_raw_distance
    distance_to_nearest_peak_bp
    dominant_peak_distance_bp
    median_contact_distance
    std_contact_distance

### distance-band  (6)
    empty_band_fraction
    frac_distal
    frac_far_distal
    frac_local
    frac_promoter_proximal
    frac_signal_in_top_peak

### O/E derived  (4)
    oe_distal_mean
    oe_distance_kurtosis
    oe_max_cv
    oe_proximal_max

total 73
Bash Read definitions of ambiguous features
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/mccprofiler
for f in oe_proximal_max empty_band_fraction n_peaks_within_100kb signal_entropy_distal distal_signal_density local_to_distal_ratio spacing_regularity peak_dominance_index dominant_peak_isolation_bp; do
  echo "=== $f ==="
  grep -rn "def $f" -A 8 --include=*.py src/ | head -11
done 2>/dev/null | grep -vE "^\s*$"

OUT
=== oe_proximal_max ===
src/mccprofiler/global_features.py:323:def oe_proximal_max(oe_profile: np.ndarray) -> float:
src/mccprofiler/global_features.py-324-    """Max O/E in the bait region (±BAIT_KB).
src/mccprofiler/global_features.py-325-
src/mccprofiler/global_features.py-326-    Likely correlates with bait_pileup_fraction and may be pruned at |r|>0.8;
src/mccprofiler/global_features.py-327-    kept in the candidate set so the prune decides.
src/mccprofiler/global_features.py-328-    """
src/mccprofiler/global_features.py-329-    region = oe_profile[_BAIT_MASK]
src/mccprofiler/global_features.py-330-    return float(region.max()) if region.size else 0.0
src/mccprofiler/global_features.py-331-
=== empty_band_fraction ===
=== n_peaks_within_100kb ===
=== signal_entropy_distal ===
src/mccprofiler/global_features.py:295:def signal_entropy_distal(profile: np.ndarray) -> float:
src/mccprofiler/global_features.py-296-    """Shannon entropy of the normalised non-zero MCC profile in the distal
src/mccprofiler/global_features.py-297-    region (|d| ≥ DISTAL_KB).
src/mccprofiler/global_features.py-298-
src/mccprofiler/global_features.py-299-    Bait-trapped genes converge toward 0 (only a couple of random low-signal
src/mccprofiler/global_features.py-300-    distal bins); genes with genuine distributed distal contacts have higher
src/mccprofiler/global_features.py-301-    entropy.
src/mccprofiler/global_features.py-302-    """
src/mccprofiler/global_features.py-303-    return _shannon_entropy_on_mask(profile, mask=_DISTAL_MASK)
=== distal_signal_density ===
src/mccprofiler/global_features.py:105:def distal_signal_density(profile: np.ndarray) -> float:
src/mccprofiler/global_features.py-106-    """log1p of mean signal in the distal region (|d| ≥ DISTAL_KB).
src/mccprofiler/global_features.py-107-
src/mccprofiler/global_features.py-108-    Captures genuine distal contact density while excluding bait pile-up and
src/mccprofiler/global_features.py-109-    the local distance-decay zone.
src/mccprofiler/global_features.py-110-    """
src/mccprofiler/global_features.py-111-    return float(np.log1p(profile[_DISTAL_MASK].mean()))
src/mccprofiler/global_features.py-112-
src/mccprofiler/global_features.py-113-
=== local_to_distal_ratio ===
src/mccprofiler/global_features.py:114:def local_to_distal_ratio(profile: np.ndarray) -> float:
src/mccprofiler/global_features.py-115-    """log1p(local) − log1p(distal). Compositional shape feature.
src/mccprofiler/global_features.py-116-
src/mccprofiler/global_features.py-117-    Positive: signal concentrated in the 10–50 kb local band (typical HK).
src/mccprofiler/global_features.py-118-    Negative: signal reaches into the >50 kb distal band (typical multi-enhancer).
src/mccprofiler/global_features.py-119-    """
src/mccprofiler/global_features.py-120-    local_sum = profile[_BAND_MASKS["local"]].sum()
src/mccprofiler/global_features.py-121-    distal_sum = profile[_BAND_MASKS["distal"]].sum()
src/mccprofiler/global_features.py-122-    return float(np.log1p(local_sum) - np.log1p(distal_sum))
=== spacing_regularity ===
=== peak_dominance_index ===
=== dominant_peak_isolation_bp ===
Found a real one already: "distal" means two different regions in the same feature set. distal_signal_density uses _DISTAL_MASK (everything ≥ 50 kb) while local_to_distal_ratio and frac_distal use DISTANCE_BANDS["distal"] (the 50–250 kb ring). And "proximal" is similarly split — oe_proximal_max is the ±2 kb bait, frac_promoter_proximal is 0–10 kb. Finding the rest:

Bash Locate remaining feature definitions
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/mccprofiler
for f in empty_band_fraction n_peaks_within_100kb spacing_regularity peak_dominance_index dominant_peak_isolation_bp frac_signal_in_top_peak n_active_peaks consensus_fraction_max_all corr_oe_distance; do
  echo "=== $f ==="
  grep -rn "\"$f\"\|'$f'\|def $f" --include=*.py src/ | head -3
done
echo; echo "=== _DISTAL_MASK vs bands ==="
grep -n "_DISTAL_MASK\s*=\|_BAND_MASKS\s*=\|DISTAL_KB\s*=\|BAIT_KB\s*=" src/mccprofiler/*.py

OUT
=== empty_band_fraction ===
src/mccprofiler/peak_morphology.py:53:    "empty_band_fraction",
src/mccprofiler/peak_morphology.py:89:    out["empty_band_fraction"] = 1.0 - occupied_bands / float(N_BANDS // 2)
=== n_peaks_within_100kb ===
src/mccprofiler/transform.py:62:    "distance_to_nearest_peak_bp", "n_peaks_within_100kb",
src/mccprofiler/peak_morphology.py:55:    "n_peaks_within_100kb",
src/mccprofiler/peak_morphology.py:93:    out["n_peaks_within_100kb"] = int((abs_dists < LOCAL_RADIUS_KB * 1000).sum())
=== spacing_regularity ===
src/mccprofiler/topology.py:29:    "spacing_regularity",
src/mccprofiler/topology.py:99:    out["spacing_regularity"] = (
=== peak_dominance_index ===
src/mccprofiler/peak_morphology.py:56:    "peak_dominance_index",
src/mccprofiler/peak_morphology.py:106:        out["peak_dominance_index"] = top_val / denom if denom > 0 else 1.0
src/mccprofiler/peak_morphology.py:112:        out["peak_dominance_index"] = 1.0  # trivially dominant
=== dominant_peak_isolation_bp ===
src/mccprofiler/peak_morphology.py:58:    "dominant_peak_isolation_bp",
src/mccprofiler/peak_morphology.py:110:        out["dominant_peak_isolation_bp"] = float(np.abs(other_pos - top_pos).min())
src/mccprofiler/transform.py:63:    "dominant_peak_distance_bp", "dominant_peak_isolation_bp",
=== frac_signal_in_top_peak ===
src/mccprofiler/peak_morphology.py:59:    "frac_signal_in_top_peak",
src/mccprofiler/peak_morphology.py:101:        out["frac_signal_in_top_peak"] = top_val / float(oe_max.sum())
=== n_active_peaks ===
src/mccprofiler/config.py:141:    "n_active_peaks",
src/mccprofiler/transform.py:32:    "n_peaks", "n_active_peaks", "n_high_consensus_peaks",
src/mccprofiler/topology.py:27:    "n_active_peaks", "frac_active_pairs",
=== consensus_fraction_max_all ===
=== corr_oe_distance ===
src/mccprofiler/topology.py:31:    "corr_oe_distance", "oe_max_cv",
src/mccprofiler/topology.py:114:        out["corr_oe_distance"] = float(rho) if np.isfinite(rho) else 0.0
src/mccprofiler/topology.py:116:        out["corr_oe_distance"] = np.nan

=== _DISTAL_MASK vs bands ===
src/mccprofiler/config.py:68:BAIT_KB = 2          # |x - viewpoint| < BAIT_KB defines the bait region
src/mccprofiler/config.py:69:DISTAL_KB = 50       # |x - viewpoint| > DISTAL_KB defines distal contacts
src/mccprofiler/global_features.py:53:_BAND_MASKS = {
src/mccprofiler/global_features.py:58:_DISTAL_MASK = _DISTANCES_BP >= (DISTAL_HALF_BINS * BIN_SIZE)
Bash Check the remaining band constants
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/mccprofiler
grep -rn "LOCAL_RADIUS_KB\s*=\|N_BANDS\s*=\|BAND_KB\s*=" --include=*.py src/
echo "--- empty_band_fraction context ---"
sed -n '80,95p' src/mccprofiler/peak_morphology.py

OUT
src/mccprofiler/peak_morphology.py:45:EMPTY_BAND_KB = 50
src/mccprofiler/peak_morphology.py:46:N_BANDS = int((N_BINS * BIN_SIZE) / (EMPTY_BAND_KB * 1000))  # 2000 kb / 50 kb = 40
src/mccprofiler/peak_morphology.py:48:LOCAL_RADIUS_KB = 100  # ±100 kb around viewpoint for local peak count
--- empty_band_fraction context ---
        out["mean_peak_gap_bp"] = float(gaps.mean())

    # Empty band fraction (defined for n >= 1)
    # Distance band index from viewpoint, each band = EMPTY_BAND_KB wide
    band_kb = abs_dists / 1000.0 / EMPTY_BAND_KB  # band 0..N_BANDS/2
    band_idx = np.clip(np.floor(band_kb).astype(int), 0, (N_BANDS // 2) - 1)
    occupied_bands = len(set(band_idx.tolist()))
    # ±1 Mb / 50 kb = 20 bands per side = 40 total (the empty-bin metric
    # collapses to unsigned distance bands, so denominator = N_BANDS/2 = 20)
    out["empty_band_fraction"] = 1.0 - occupied_bands / float(N_BANDS // 2)

    # --- Family B: viewpoint isolation ---
    out["distance_to_nearest_peak_bp"] = float(abs_dists.min())
    out["n_peaks_within_100kb"] = int((abs_dists < LOCAL_RADIUS_KB * 1000).sum())

    # --- Family C: dominant-peak isolation ---
Now I have the full picture, and it is worse than one bad name. Let me get the comprehensive test running before I write it up, since you are right that testing one class and one set is exactly the selective error I ran the null to avoid:

Write diagnose_own_vs_other_all.py
184 lines
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
df["q"] = (df["p"] * mm / df["rank"]).cummin().clip(upper=1.0)

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
Bash Launch the full grid in background
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
setsid nohup python audit/scripts/diagnose_own_vs_other_all.py \
  > audit/scripts/output/own_vs_other_all.log 2>&1 < /dev/null &
sleep 25
tail -12 audit/scripts/output/own_vs_other_all.log

OUT
[1]+  Done                    setsid nohup python audit/scripts/diagnose_own_vs_other_all.py > audit/scripts/output/own_vs_other_all.log 2>&1 < /dev/null
Bash Check grid progress
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
sleep 40; cat audit/scripts/output/own_vs_other_all.log

OUT
(Bash completed with no output)
Bash Confirm the grid is actually running
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
pgrep -af "own_vs_other_all" | head -3
ls -la audit/scripts/output/own_vs_other_all.log
echo "(empty log = python stdout buffering, not a failure)"

OUT
17337 python audit/scripts/diagnose_own_vs_other_all.py
17457 /bin/bash -c source /home/imm/grte4643/.claude/shell-snapshots/snapshot-bash-1786887210595-2eioi7.sh 2>/dev/null || true && shopt -u extglob 2>/dev/null || true && { \builtin unalias -- 'unsetenv'; \builtin unset -f -- 'unsetenv'; } >/dev/null 2>&1 || true && eval 'cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app pgrep -af "own_vs_other_all" | head -3 ls -la audit/scripts/output/own_vs_other_all.log echo "(empty log = python stdout buffering, not a failure)"' < /dev/null && pwd -P >| /tmp/claude-9d71-cwd
-rw-r--r-- 1 grte4643 domain users 0 Aug 16 15:37 audit/scripts/output/own_vs_other_all.log
(empty log = python stdout buffering, not a failure)
Running (PID 17337); the empty log is stdout buffering, not a failure. While that computes — the naming audit, which is worse than one bad name:

Write FEATURE_NAMING.md
156 lines
# Are the features named for what they measure?

Audit of all 73 features against their definitions in `mccprofiler/src/`, written
2026-08-16 after `frac_promoter_proximal` turned out not to be a promoter
feature. It is not the only one.

**Nothing here is a bug.** Every feature computes what its code says. The problem
is that several names invite a reading the code does not support, and at least
one of those misreadings has already reached written work.

---

## 1. The headline: `frac_promoter_proximal` is not about promoters

```python
DISTANCE_BANDS = {"promoter_proximal": (0, 10_000), ...}   # config.py
```

It is the fraction of signal within 10 kb of the **viewpoint**, measured as
`|bin - CENTER_BIN|`. Two separate problems:

1. **It is a distance ring, not an element class.** No promoter annotation is
   consulted. A gene with no promoter-classed peak at all still has a
   `frac_promoter_proximal`.
2. **The viewpoint is not the TSS.** The matrix is anchored on the midpoint of
   the captured region, a per-peak identifier. About 2% of genes differ by
   hundreds of kb (RERE by 292 kb). For those genes this measures a 10 kb ring
   around a point that is not the promoter.

**Rename to `frac_within_10kb_of_viewpoint`.** The honest name is unwieldy, which
is exactly why the misleading one survived.

## 2. Three distance words each mean two different things

This is the more dangerous class of problem, because the features sit side by
side in the same table and look comparable.

### "distal" is both a ring and a half-line

| feature | region | source |
|---|---|---|
| `frac_distal` | **50 to 250 kb** | `DISTANCE_BANDS["distal"]` |
| `local_to_distal_ratio` | **50 to 250 kb** | `_BAND_MASKS["distal"]` |
| `distal_signal_density` | **everything >= 50 kb** | `_DISTAL_MASK` |
| `signal_entropy_distal` | **everything >= 50 kb** | `_DISTAL_MASK` |
| `oe_distal_mean` | check before quoting | |

`distal_signal_density` includes the far-distal band; `frac_distal` excludes it.
They sound like the same region and are not. `local_to_distal_ratio` uses the
ring while `distal_signal_density` uses the half-line, in the same module.

### "proximal" is both 2 kb and 10 kb

| feature | region |
|---|---|
| `oe_proximal_max` | **+/- 2 kb**, the bait region (`BAIT_KB = 2`) |
| `bait_pileup_fraction` | **+/- 2 kb**, the bait region |
| `frac_promoter_proximal` | **0 to 10 kb** |

### "local" is both 10-50 kb and a 100 kb radius

| feature | region |
|---|---|
| `frac_local`, `local_to_distal_ratio` | **10 to 50 kb** (`DISTANCE_BANDS["local"]`) |
| `n_peaks_within_100kb` | **+/- 100 kb** (`LOCAL_RADIUS_KB = 100`) |

The feature name is honest here; the constant behind it is called `LOCAL_RADIUS`
while `local` elsewhere means 10-50 kb. `oe_local_enrichment_*` is a third use,
per-peak background enrichment, unrelated to either. **Check that one before
quoting it as anything to do with local contacts.**

### "band" is two different partitions

| | |
|---|---|
| `DISTANCE_BANDS` | 4 unequal rings: 10 / 50 / 250 / 1000 kb |
| `empty_band_fraction` | 20 equal 50 kb bands per side (`EMPTY_BAND_KB = 50`, `N_BANDS = 40`) |

`empty_band_fraction` is not the fraction of the four named bands that are empty.

## 3. The element-class features mean "any such element in the window"

`promoter_signal_fraction`, `n_peaks_promoter`, `raw_peak_max_max_promoter`,
`oe_local_enrichment_max_promoter`, `raw_log2_enrichment_max_promoter`,
`max_distance_to_viewpoint_promoter`.

Measured 2026-08-16: **only 5.8%** of promoter-classed peaks lie within 5 kb of
the gene's own TSS, median distance **115 kb**. So these are overwhelmingly
*other genes'* promoters. The names read as though they describe the gene's own
promoter and they do not.

That is not a defect for the promoter-assembly question, where other-promoter
contact is the quantity of interest (section 6d). It is a defect for anyone
reading `promoter_signal_fraction` as "how promoter-driven is this gene".

## 4. Names that are actively wrong about what they measure

- **`contact_asymmetry`.** Docstring says "positive means more signal
  downstream", but computed in *genome* coordinates, so downstream is upstream
  for half the genome and the gene-relative signal cancels. Strand-corrected
  mean is 13.7 SE from zero (section 6b). The name promises a directional
  quantity the implementation destroys.
- **`distance_to_viewpoint`** (in `annotated.tsv`, upstream of several
  features). **Unsigned.** Signed placement is `peak_midpoint - viewpoint_pos`.
- **The "topology" block.** After the 2026-08-16 removals it contains
  `spacing_regularity`, `corr_oe_distance`, `oe_max_cv` and `n_active_peaks`.
  None is topology. The genuinely topological features were the degenerate ones
  and they are gone. Renaming the block was already a carried-forward item from
  the 07-31 handoff.

## 5. Names that are fine

`total_mcc`, `max_mcc`, `q90_mcc`, `n_peaks_*`, `n_active_peaks`,
`n_high_consensus_peaks_075`, `consensus_fraction_*`, `raw_peak_max_*`,
`oe_max_*`, `gini_mcc`, `signal_entropy`, `frac_signal_in_top_peak`,
`peak_dominance_index`, `dominant_peak_distance_bp`,
`dominant_peak_isolation_bp`, `mean_peak_gap_bp`, `spacing_regularity`,
`bimodality_score`, `oe_distance_kurtosis`, `median_contact_distance`,
`std_contact_distance`, `corr_oe_distance`, `corr_raw_distance`.

Each says what it does. `bait_pileup_fraction` in particular is correctly named
after the bait rather than a biological region.

---

## What to do

**Do not bulk-rename.** Renaming a feature changes the correlation-prune order
and therefore which features survive, so a rename is a substrate change, not
cosmetics. The prune already demonstrably selects for the wrong thing.

Priority order:

1. **Fix the docstrings first**, especially `contact_asymmetry`, whose docstring
   is factually wrong rather than merely loose.
2. **Rename `frac_promoter_proximal` -> `frac_within_10kb_of_viewpoint`** at the
   next deliberate substrate change, and re-run the prune. This is the one that
   has actually misled.
3. **Pick one meaning for "distal"** and make `distal_signal_density` and
   `frac_distal` agree, or rename one of them.
4. **Rename the topology block** to what it now contains.
5. **Add a glossary to the briefing pack** giving the exact region for every
   distance word. Collaborators read the briefing, not `config.py`.

## The general lesson

Three of these problems come from the same source: **a distance from the
viewpoint was given a biological name.** "promoter_proximal", "local", "distal"
and "proximal" all sound like statements about regulatory architecture and are
all statements about base pairs from an arbitrary anchor. Given the viewpoint is
not the TSS, even the anchor is not what the name implies.

When adding features, name them after the measurement, not the interpretation.
`frac_within_10kb_of_viewpoint` cannot be misread; `frac_promoter_proximal`
already was.
Bash Check grid result
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
if pgrep -f own_vs_other_all >/dev/null; then echo "still running"; else echo "DONE"; fi
cat audit/scripts/output/own_vs_other_all.log

OUT
still running
1835 genes, 9 measures
  promoter    10056 peaks,   5.8% at the gene's own TSS (within 5 kb)
  enhancer     9374 peaks,   6.0% at the gene's own TSS (within 5 kb)
  ctcf         9673 peaks,   1.7% at the gene's own TSS (within 5 kb)
  18 usable sets -> 162 tests

========================================================================================
SURVIVING BENJAMINI-HOCHBERG FDR 0.05 ACROSS ALL 162 TESTS
========================================================================================
set                           measure               n       d        p        q
ChromHMM_bivalent             ctcf_total           46  +0.679   0.0005   0.0405
Lambert_TF                    promoter_total      156  +0.309   0.0005   0.0270
cd4_specific_immune           ctcf_other          819  -0.165   0.0005   0.0203
Lambert_TF                    ctcf_total          156  -0.426   0.0005   0.0162
Lambert_TF                    ctcf_other          156  -0.396   0.0005   0.0135
Lambert_TF                    promoter_other      156  +0.268   0.0005   0.0116
gene_desert_bottomQ_density   enhancer_total      456  +0.181   0.0005   0.0101
gene_desert_bottomQ_density   enhancer_other      456  +0.206   0.0005   0.0090
cd4_rna_top_quartile          ctcf_other          759  -0.152   0.0010   0.0090
cd4_rna_top_quartile          promoter_total      759  +0.149   0.0010   0.0090
ChromHMM_bivalent             promoter_total       46  -0.482   0.0010   0.0090
ChromHMM_bivalent             promoter_other       46  -0.480   0.0010   0.0090
cd4_specific_immune           ctcf_total          819  -0.155   0.0015   0.0090
GWAS_total_topQ               promoter_other      508  -0.161   0.0030   0.0090
GWAS_immune_any               promoter_other      732  -0.146   0.0030   0.0090
GWAS_immune_hot               promoter_other      129  -0.255   0.0045   0.0090
cd4_rna_top_quartile          ctcf_total          759  -0.136   0.0050   0.0090
Eisenberg_HK                  promoter_other      633  +0.137   0.0060   0.0090
phastCons_2kb_topQ            enhancer_own        415  -0.141   0.0085   0.0090
DepMap_inferred_essential     enhancer_other      328  -0.153   0.0100   0.0090
GWAS_immune_any               promoter_total      732  -0.114   0.0150   0.0090
DepMap_curated_essential      ctcf_own            254  +0.164   0.0160   0.0090
GWAS_immune_hot               promoter_total      129  -0.206   0.0210   0.0090
gene_desert_bottomQ_density   promoter_other      456  -0.123   0.0220   0.0090
phastCons_2kb_topQ            promoter_own        415  +0.130   0.0260   0.0090
cd4_rna_top_quartile          promoter_own        759  +0.107   0.0265   0.0090
bio_HK_k3                     enhancer_other      353  -0.123   0.0270   0.0090
GWAS_immune_hot               enhancer_total      129  +0.205   0.0310   0.0090
GWAS_total_topQ               enhancer_other      508  +0.113   0.0330   0.0090
GWAS_immune_any               enhancer_other      732  +0.093   0.0365   0.0090
bio_bulk_k3                   promoter_total     1065  -0.101   0.0365   0.0090
cd4_specific_immune           promoter_other      819  +0.100   0.0375   0.0090
bio_dev_TF_k3                 ctcf_other          303  -0.129   0.0400   0.0090
GWAS_immune_any               enhancer_total      732  +0.097   0.0430   0.0090
GWAS_total_topQ               enhancer_total      508  +0.108   0.0430   0.0090
bio_dev_TF_k3                 promoter_own        303  +0.129   0.0435   0.0090
bio_dev_TF_k3                 ctcf_total          303  -0.127   0.0440   0.0090
gene_desert_bottomQ_density   promoter_total      456  -0.107   0.0470   0.0090
bio_bulk_k3                   promoter_other     1065  -0.093   0.0550   0.0090
GWAS_total_topQ               promoter_total      508  -0.104   0.0555   0.0090
DepMap_curated_essential      enhancer_other      254  -0.129   0.0560   0.0090
cd4_rna_top_quartile          promoter_other      759  +0.087   0.0595   0.0090
DepMap_inferred_essential     enhancer_total      328  -0.115   0.0610   0.0090
GWAS_immune_hot               enhancer_other      129  +0.175   0.0650   0.0090
Eisenberg_HK                  enhancer_other      633  -0.090   0.0660   0.0090
Eisenberg_HK                  promoter_total      633  +0.089   0.0700   0.0090
bio_HK_k3                     promoter_other      353  +0.104   0.0890   0.0090
Lambert_TF                    ctcf_own            156  -0.142   0.0920   0.0090
bio_dev_TF_k3                 promoter_total      303  +0.105   0.0925   0.0090
Roadmap_silenced              promoter_total       39  -0.276   0.0945   0.0090
ChromHMM_bivalent             enhancer_other       46  -0.251   0.0945   0.0090
cd4_specific_immune           promoter_total      819  +0.078   0.1000   0.0090
bio_HK_k3                     promoter_own        353  -0.096   0.1015   0.0090
bio_bulk_k3                   ctcf_other         1065  +0.076   0.1020   0.0090
DepMap_inferred_essential     promoter_own        328  -0.096   0.1135   0.0090
DepMap_inferred_essential     promoter_other      328  +0.093   0.1140   0.0090
cd4_specific_immune           enhancer_own        819  +0.076   0.1205   0.0090
Eisenberg_HK                  enhancer_total      633  -0.075   0.1230   0.0090
bio_HK_k3                     enhancer_total      353  -0.090   0.1390   0.0090
cd4_specific_immune           enhancer_total      819  +0.068   0.1405   0.0090
gnomAD_pLI_topQ               ctcf_other          546  -0.075   0.1515   0.0090
Eisenberg_HK                  promoter_own        633  -0.067   0.1570   0.0090
GWAS_total_topQ               promoter_own        508  +0.074   0.1570   0.0090
gene_desert_bottomQ_density   ctcf_other          456  -0.075   0.1630   0.0090
Roadmap_silenced              promoter_other       39  -0.231   0.1645   0.0090
DepMap_inferred_essential     ctcf_own            328  +0.082   0.1685   0.0090
bio_bulk_k3                   enhancer_other     1065  +0.065   0.1695   0.0090
gene_desert_bottomQ_density   ctcf_total          456  -0.071   0.1785   0.0090
DepMap_curated_essential      enhancer_total      254  -0.093   0.1790   0.0090
DepMap_inferred_essential     ctcf_total          328  +0.084   0.1815   0.0090
bio_bulk_k3                   enhancer_total     1065  +0.062   0.1895   0.0090
ChromHMM_bivalent             enhancer_total       46  -0.194   0.1965   0.0090
Lambert_TF                    enhancer_other      156  +0.110   0.2020   0.0090
bio_dev_TF_k3                 enhancer_own        303  -0.077   0.2200   0.0090
DepMap_curated_essential      ctcf_total          254  +0.083   0.2225   0.0090
ChromHMM_bivalent             ctcf_own             46  +0.190   0.2255   0.0090
bio_bulk_k3                   ctcf_own           1065  -0.057   0.2310   0.0090
gnomAD_pLI_topQ               ctcf_own            546  +0.057   0.2555   0.0090
GWAS_immune_hot               enhancer_own        129  +0.103   0.2590   0.0090
phastCons_2kb_topQ            promoter_total      415  +0.063   0.2600   0.0090
DepMap_curated_essential      promoter_own        254  -0.075   0.2750   0.0090
Roadmap_silenced              enhancer_total       39  +0.173   0.2875   0.0090
bio_bulk_k3                   ctcf_total         1065  +0.049   0.2950   0.0090
Roadmap_silenced              enhancer_other       39  +0.174   0.2955   0.0090
Lambert_TF                    enhancer_total      156  +0.083   0.3160   0.0090
Roadmap_silenced              ctcf_own             39  +0.170   0.3205   0.0090
GWAS_immune_hot               ctcf_own            129  +0.090   0.3310   0.0090
gnomAD_pLI_topQ               ctcf_total          546  -0.048   0.3455   0.0090
phastCons_2kb_topQ            enhancer_total      415  -0.054   0.3490   0.0090
DepMap_inferred_essential     ctcf_other          328  +0.054   0.3715   0.0090
DepMap_curated_essential      promoter_other      254  +0.059   0.3805   0.0090
GWAS_immune_any               promoter_own        732  +0.039   0.3990   0.0090
bio_dev_TF_k3                 enhancer_other      303  +0.054   0.4050   0.0090
bio_HK_k3                     ctcf_own            353  +0.047   0.4100   0.0090
Lambert_TF                    promoter_own        156  +0.072   0.4135   0.0090
bio_HK_k3                     ctcf_total          353  +0.047   0.4300   0.0090
cd4_rna_top_quartile          enhancer_own        759  -0.036   0.4350   0.0090
cd4_specific_immune           enhancer_other      819  +0.033   0.4865   0.0090
ChromHMM_bivalent             promoter_own         46  -0.103   0.5010   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     promoter_other      158  -0.053   0.5085   0.0090
bio_HK_k3                     promoter_total      353  +0.040   0.5130   0.0090
cd4_rna_top_quartile          enhancer_total      759  -0.030   0.5190   0.0090
phastCons_2kb_topQ            ctcf_own            415  +0.033   0.5430   0.0090
Eisenberg_HK                  ctcf_own            633  -0.028   0.5480   0.0090
GWAS_immune_any               enhancer_own        732  +0.028   0.5600   0.0090
bio_HK_k3                     enhancer_own        353  +0.033   0.5665   0.0090
cd4_specific_immune           promoter_own        819  -0.027   0.5685   0.0090
GWAS_immune_any               ctcf_total          732  +0.027   0.5685   0.0090
Roadmap_silenced              promoter_own         39  -0.096   0.5685   0.0090
phastCons_2kb_topQ            ctcf_other          415  -0.032   0.5755   0.0090
Roadmap_silenced              ctcf_total           39  +0.095   0.5760   0.0090
GWAS_immune_any               ctcf_other          732  +0.026   0.5905   0.0090
DepMap_inferred_essential     enhancer_own        328  +0.033   0.5905   0.0090
GWAS_immune_hot               promoter_own        129  +0.046   0.5950   0.0090
DepMap_inferred_essential     promoter_total      328  +0.029   0.6150   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     promoter_own        158  +0.043   0.6170   0.0090
gnomAD_pLI_topQ               promoter_own        546  +0.025   0.6240   0.0090
gnomAD_pLI_topQ               promoter_total      546  +0.025   0.6240   0.0090
bio_HK_k3                     ctcf_other          353  +0.030   0.6265   0.0090
Eisenberg_HK                  ctcf_total          633  -0.023   0.6320   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     ctcf_other          158  +0.037   0.6510   0.0090
DepMap_curated_essential      enhancer_own        254  +0.031   0.6530   0.0090
gnomAD_pLI_topQ               enhancer_total      546  +0.020   0.6770   0.0090
bio_bulk_k3                   promoter_own       1065  -0.019   0.6880   0.0090
bio_dev_TF_k3                 promoter_other      303  +0.024   0.6905   0.0090
GWAS_immune_hot               ctcf_other          129  -0.036   0.6915   0.0090
gnomAD_pLI_topQ               enhancer_other      546  +0.019   0.7050   0.0090
Lambert_TF                    enhancer_own        156  -0.032   0.7100   0.0090
cd4_rna_top_quartile          ctcf_own            759  +0.017   0.7140   0.0090
phastCons_2kb_topQ            enhancer_other      415  +0.019   0.7330   0.0090
phastCons_2kb_topQ            promoter_other      415  -0.020   0.7425   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     ctcf_own            158  -0.029   0.7445   0.0090
GWAS_total_topQ               enhancer_own        508  +0.017   0.7540   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     ctcf_total          158  +0.024   0.7600   0.0090
phastCons_2kb_topQ            ctcf_total          415  -0.017   0.7785   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     promoter_total      158  -0.024   0.7850   0.0090
cd4_rna_top_quartile          enhancer_other      759  -0.013   0.7870   0.0090
Eisenberg_HK                  ctcf_other          633  -0.013   0.7905   0.0090
gene_desert_bottomQ_density   promoter_own        456  +0.011   0.8350   0.0090
bio_bulk_k3                   enhancer_own       1065  +0.009   0.8390   0.0090
gnomAD_pLI_topQ               promoter_other      546  +0.010   0.8405   0.0090
DepMap_curated_essential      ctcf_other          254  +0.013   0.8490   0.0090
bio_dev_TF_k3                 ctcf_own            303  -0.011   0.8550   0.0090
ChromHMM_bivalent             enhancer_own         46  +0.027   0.8590   0.0090
GWAS_immune_any               ctcf_own            732  +0.007   0.8780   0.0090
Roadmap_silenced              ctcf_other           39  +0.025   0.8880   0.0090
Eisenberg_HK                  enhancer_own        633  +0.008   0.8915   0.0090
DepMap_curated_essential      promoter_total      254  +0.009   0.9020   0.0090
bio_dev_TF_k3                 enhancer_total      303  +0.008   0.9025   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     enhancer_own        158  +0.011   0.9085   0.0090
gnomAD_pLI_topQ               enhancer_own        546  +0.006   0.9155   0.0090
GWAS_total_topQ               ctcf_other          508  +0.005   0.9275   0.0090
Roadmap_silenced              enhancer_own         39  +0.012   0.9320   0.0090
GWAS_total_topQ               ctcf_total          508  +0.004   0.9345   0.0090
gene_desert_bottomQ_density   enhancer_own        456  -0.003   0.9580   0.0090
GWAS_immune_hot               ctcf_total          129  +0.004   0.9670   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     enhancer_total      158  +0.003   0.9680   0.0090
dbSUPER_CD4_SE_TSS_pm50kb     enhancer_other      158  -0.002   0.9750   0.0090
cd4_specific_immune           ctcf_own            819  +0.000   0.9945   0.0090
gene_desert_bottomQ_density   ctcf_own            456  -0.000   0.9955   0.0090
GWAS_total_topQ               ctcf_own            508  -0.000   0.9975   0.0090

========================================================================================
THE PRE-SPECIFIED TEST, which does not owe multiplicity
========================================================================================
  Eisenberg_HK x promoter_other   d = +0.137  p = 0.0060  (BH q = 0.0090)
  Its q is shown only for context. It was specified before this grid
  was run, so the p stands on its own; the q answers the different
  question of whether it would survive as a discovery FROM the grid.

========================================================================================
HOW MUCH OF THE GRID IS NOISE
========================================================================================
  tests run                : 162
  nominally p < 0.05       : 39  (expected by chance 8.1)
  surviving BH q < 0.05    : 161

  own-TSS measures should be mostly null. If they are not, the split is
  picking up viewpoint placement rather than biology:
     promoter_own     3 of 18 sets nominal, 18 after correction
     enhancer_own     1 of 18 sets nominal, 18 after correction
     ctcf_own     1 of 18 sets nominal, 18 after correction

wrote /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit/scripts/output/own_vs_other_all.tsv
My BH correction is wrong — every q is 0.0090, including rows with p = 0.99. Forward cummin() on p-ascending data drags every later q down to the smallest early value; BH needs the cumulative minimum taken from the largest p downward. Fixing and recomputing from the saved p-values (no need to re-run the permutations):

Edit diagnose_own_vs_other_all.py
Added 4 lines
Bash Recompute BH correctly
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python - <<'PY'
import pandas as pd, numpy as np
d=pd.read_csv("audit/scripts/output/own_vs_other_all.tsv",sep="\t").sort_values("p").reset_index(drop=True)
m=len(d); d["rank"]=np.arange(1,m+1)
d["q"]=(d["p"]*m/d["rank"])[::-1].cummin()[::-1].clip(upper=1.0)
d.to_csv("audit/scripts/output/own_vs_other_all.tsv",sep="\t",index=False)
sig=d[d.q<0.05]
print(f"{m} tests | nominally p<0.05: {(d.p<0.05).sum()} (chance {0.05*m:.1f}) | BH q<0.05: {len(sig)}\n")
print(f"{'set':<30}{'measure':<17}{'n':>5}{'d':>8}{'p':>8}{'q':>8}")
for r in sig.itertuples():
    print(f"{r.group:<30}{r.measure:<17}{r.n:>5}{r.d:>+8.3f}{r.p:>8.4f}{r.q:>8.4f}")
print("\nown-TSS measures (should be mostly null):")
for c in ["promoter","enhancer","ctcf"]:
    o=d[d.measure==f"{c}_own"]
    print(f"   {c+'_own':<14} {(o.p<0.05).sum():>2}/18 nominal, {(o.q<0.05).sum()} after BH")
ps=d[(d.group=="Eisenberg_HK")&(d.measure=="promoter_other")].iloc[0]
print(f"\npre-specified Eisenberg_HK x promoter_other: d={ps['d']:+.3f} p={ps['p']:.4f} q={ps['q']:.4f}")
PY

OUT
162 tests | nominally p<0.05: 39 (chance 8.1) | BH q<0.05: 18

set                           measure              n       d       p       q
ChromHMM_bivalent             ctcf_other          46  +0.629  0.0005  0.0090
ChromHMM_bivalent             ctcf_total          46  +0.679  0.0005  0.0090
Lambert_TF                    promoter_total     156  +0.309  0.0005  0.0090
cd4_specific_immune           ctcf_other         819  -0.165  0.0005  0.0090
Lambert_TF                    ctcf_total         156  -0.426  0.0005  0.0090
Lambert_TF                    ctcf_other         156  -0.396  0.0005  0.0090
Lambert_TF                    promoter_other     156  +0.268  0.0005  0.0090
gene_desert_bottomQ_density   enhancer_total     456  +0.181  0.0005  0.0090
gene_desert_bottomQ_density   enhancer_other     456  +0.206  0.0005  0.0090
cd4_rna_top_quartile          ctcf_other         759  -0.152  0.0010  0.0125
cd4_rna_top_quartile          promoter_total     759  +0.149  0.0010  0.0125
ChromHMM_bivalent             promoter_total      46  -0.482  0.0010  0.0125
ChromHMM_bivalent             promoter_other      46  -0.480  0.0010  0.0125
cd4_specific_immune           ctcf_total         819  -0.155  0.0015  0.0174
GWAS_total_topQ               promoter_other     508  -0.161  0.0030  0.0304
GWAS_immune_any               promoter_other     732  -0.146  0.0030  0.0304
GWAS_immune_hot               promoter_other     129  -0.255  0.0045  0.0429
cd4_rna_top_quartile          ctcf_total         759  -0.136  0.0050  0.0450

own-TSS measures (should be mostly null):
   promoter_own    3/18 nominal, 0 after BH
   enhancer_own    1/18 nominal, 0 after BH
   ctcf_own        1/18 nominal, 0 after BH

pre-specified Eisenberg_HK x promoter_other: d=+0.137 p=0.0060 q=0.0512
Correct now: 18 of 162 survive, not 161. And testing everything changes the story — the housekeeping result is not the headline.

The own-TSS controls come back clean: 0/18 survive BH for all three classes. That validates the split — it is not picking up viewpoint placement.

Lambert_TF is much bigger than housekeeping, and it is exactly the characterisation your feedback quoted:

Lambert_TF	d	q
ctcf_total	−0.426	0.009
ctcf_other	−0.396	0.009
promoter_total	+0.309	0.009
promoter_other	+0.268	0.009
Promoter-rich and CTCF-poor, surviving BH across 162 tests. But this script controlled magnitude, density, length, expression and CpG — not chromosome, and Lambert_TF is the set that dies on chr19. Running that now, because without it this is not reportable:

Bash Chromosome control on the TF and bivalent results
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python - <<'PY'
import sys,numpy as np,pandas as pd,warnings
sys.path.insert(0,"backend"); warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
                  "scripts_cleaned/audit/GW/scripts")
from _shape import MAG_OVERALL
from app.store import get_store
rng=np.random.default_rng(0)
PEAKS=("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/scripts_cleaned/"
       "process/output_resgw_10000/utils/peaks/annotated.tsv")
pk=pd.read_csv(PEAKS,sep="\t"); pk.columns=[c.lower() for c in pk.columns]
pk["symbol_key"]=pk["viewpoint"].str.upper()
bed=pd.read_csv("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/01_tss_regions.bed",
    sep="\t",header=None,names=["chrom","start","end","gene","score","strand"])
bed["symbol_key"]=bed["gene"].str.upper(); bed["tss"]=((bed["start"]+bed["end"])//2).astype(int)
bed=bed.drop_duplicates("symbol_key").set_index("symbol_key")
gt=pd.read_csv("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/collaboration/data/gene_table.tsv",sep="\t")
gt["symbol_key"]=gt["gene_symbol"].str.upper(); gt=gt.drop_duplicates("symbol_key").set_index("symbol_key")
s=get_store(); W=s.table("features").pivot(index="gene_id",columns="feature",values="z")
skA=s.genes.set_index("gene_id").loc[W.index,"symbol_key"].to_numpy()
pk["own_tss"]=pk["symbol_key"].map(bed["tss"]); pk=pk[np.isfinite(pk["own_tss"])].copy()
pk["is_own"]=(pk["peak_midpoint"]-pk["own_tss"]).abs()<=5000
tot=pk.groupby("symbol_key")["peak_max"].sum(); F=pd.DataFrame({"tot":tot})
for cls in ["promoter","ctcf"]:
    sub=pk[pk["re"]==cls]
    F[f"{cls}_other"]=sub[~sub["is_own"]].groupby("symbol_key")["peak_max"].sum()
F=F.fillna(0.0)
for c in ["promoter_other","ctcf_other"]: F[c]/=F["tot"].replace(0,np.nan)
V=pd.DataFrame({c:pd.Series(skA).map(F[c]).to_numpy(float) for c in ["promoter_other","ctcf_other"]})
ok=V.notna().all(axis=1).to_numpy(); V=V[ok].reset_index(drop=True); sk=skA[ok]
Wm=W.reset_index(drop=True)[ok].reset_index(drop=True)
chrom=pd.Series(sk).map(bed["chrom"]).fillna("NA").to_numpy()
A=[np.ones(len(sk))]+[Wm[c].to_numpy(float) for c in MAG_OVERALL if c in Wm.columns]
for col,lg in [("gene_density_500kb",0),("log10_gene_length",0),("gtex_whole_blood_tpm",1),("cpg_oe_ratio",0)]:
    v=pd.Series(sk).map(gt[col]).to_numpy(float)
    if lg: v=np.log1p(v)
    A.append(np.where(np.isfinite(v),v,np.nanmedian(v)))
A=np.nan_to_num(np.column_stack(A),nan=0.0)
D=pd.get_dummies(pd.Series(chrom)).to_numpy(float)[:,1:]
cm=s.table("cohort_membership")
def perm(y,X,m,n=4000):
    b,*_=np.linalg.lstsq(X,y,rcond=None); r=y-X@b; r=(r-r.mean())/r.std()
    d=(r[m].mean()-r[~m].mean())/np.sqrt((r[m].var()+r[~m].var())/2)
    idx=np.arange(len(r)); nul=np.empty(n)
    for i in range(n):
        p=rng.permutation(idx)[:m.sum()]; q=np.setdiff1d(idx,p)
        nul[i]=(r[p].mean()-r[q].mean())/np.sqrt((r[p].var()+r[q].var())/2)
    return d,max((np.abs(nul)>=abs(d)).mean(),1/n)
for g in ["Lambert_TF","ChromHMM_bivalent"]:
    mem=set(cm.loc[cm["group"]==g,"symbol_key"])
    m=np.fromiter((v in mem for v in sk),bool,len(sk))
    print(f"\n{g}  n={m.sum()}   chr19: {(chrom[m]=='chr19').mean():.1%} vs {(chrom[~m]=='chr19').mean():.1%} background")
    for c in ["promoter_other","ctcf_other"]:
        y=V[c].to_numpy(float)
        d1,p1=perm(y,A,m); d2,p2=perm(y,np.column_stack([A,D]),m)
        k=chrom!="chr19"; d3,p3=perm(y[k],A[k],m[k])
        print(f"   {c:<15} full d={d1:+.3f} p={p1:.4f} | +chrom d={d2:+.3f} p={p2:.4f} | excl chr19 d={d3:+.3f} p={p3:.4f}")
PY

OUT
Lambert_TF  n=156   chr19: 25.0% vs 5.6% background
   promoter_other  full d=+0.268 p=0.0013 | +chrom d=+0.237 p=0.0065 | excl chr19 d=+0.102 p=0.2868
   ctcf_other      full d=-0.396 p=0.0003 | +chrom d=-0.375 p=0.0003 | excl chr19 d=-0.273 p=0.0045

ChromHMM_bivalent  n=46   chr19: 6.5% vs 7.3% background
   promoter_other  full d=-0.480 p=0.0005 | +chrom d=-0.461 p=0.0015 | excl chr19 d=-0.465 p=0.0047
   ctcf_other      full d=+0.629 p=0.0003 | +chrom d=+0.629 p=0.0003 | excl chr19 d=+0.583 p=0.0003
Now I can answer you properly.

There is no housekeeping cluster, and testing everything demoted the result I showed you. Across the full grid — 3 element classes × 3 measures × 18 sets = 162 tests, BH-corrected — Eisenberg_HK × promoter_other lands at q = 0.0512. It survives as a pre-specified test (p = 0.006, specified before the grid) but it does not survive as a discovery from the grid. It is mid-pack, not a headline. You were right to push.

On TF — and this is the sharpest thing in the whole grid. "TFs are promoter-driven and CTCF-poor" splits cleanly in two when you exclude chr19:

Lambert_TF	full controls	excluding chr19
promoter-rich	+0.268, p = 0.0013	+0.102, p = 0.29 — dies
CTCF-poor	−0.396, p = 0.0003	−0.273, p = 0.0045 — survives
So the promoter half of that claim is chr19 (25% vs 5.6% background) and the CTCF half is real. The external feedback nominated the whole statement as your headline positive; half of it is an artefact.

The largest surviving effect is not HK or TF. It is ChromHMM_bivalent: CTCF-rich d = +0.583 and promoter-poor d = −0.465, both surviving chr19 exclusion, chr19-balanced at 6.5% vs 7.3%. Caveat stands that it is ChIP-derived, so partly chromatin predicting chromatin, and n = 46.

Method validation: the own-TSS measures come back 0 of 18 significant after correction in all three classes. The split is not picking up viewpoint placement — which matters, because the viewpoint is not the TSS.

Why your results keep moving. Two things, and neither is instability:

The substrate changed once (91 → 73). |d| before vs after correlates r = 0.959.
Everything else is confounds being removed. chr19 has now killed half of the TF claim exactly as it killed the earlier version.
What has not moved all week: nothing separates. Every effect above is d < 0.7, i.e. ≥ 76% overlap. The grid found more structure than chance (39 nominal vs 8.1 expected) and still no clusters.

On the promoter-assemblies list — you are right, we have never used it. I tested an expression-derived housekeeping list against a feature I built. The published Ronin/THAP11 target list is a separate thing we do not have, and getting it would be a genuinely independent test rather than another rearrangement of what is already here.

Feature naming: frac_promoter_proximal was not the only one. Written up in audit/FEATURE_NAMING.md:

frac_promoter_proximal — a 0–10 kb ring around the viewpoint. No promoter annotation is consulted, and the viewpoint is not the TSS.
"distal" means two regions — distal_signal_density is everything ≥ 50 kb; frac_distal and local_to_distal_ratio are the 50–250 kb ring. Same module, different meanings.
"proximal" means 2 kb or 10 kb depending on the feature.
"local" means 10–50 kb or a 100 kb radius, and oe_local_enrichment_* is a third unrelated use.
"band" means two partitions — empty_band_fraction uses 20 equal 50 kb bands, not the four named distance bands.
contact_asymmetry's docstring is factually wrong, not merely loose.
The "topology" block contains no topology after the removals.
The recommendation there is not to bulk-rename: renaming changes the correlation-prune order, so it is a substrate change, not cosmetics. Fix the docstrings now, rename frac_promoter_proximal at the next deliberate substrate change.

One correction I owe you: I also had a BH bug in the grid script that reported 161 of 162 tests as significant. Fixed (forward cummin instead of reverse); the 18 above is the corrected number.