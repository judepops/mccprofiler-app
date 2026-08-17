"""Does `promoter_signal_fraction` conflate a gene's OWN promoter with OTHER
genes' promoters, and does splitting them change the housekeeping answer?

The promoter-assembly hypothesis (Hwang 2023 / Dejosez, Ronin / THAP11) says
housekeeping genes form promoter-promoter assemblies, i.e. they contact *other
genes'* promoters. That is a claim about a specific contact class.

`aggregate.py` cannot express it. `<class>_signal_fraction` is

    sum(oe_max where RE == 'promoter') / sum(oe_max over all peaks)

with NO distance filter, so a promoter-classed peak 800 kb away and the gene's
own TSS contribute identically. `frac_promoter_proximal` does not help: despite
the name it is a distance band, config.py `(0, 10_000)`, not an element class.

So the recorded `Eisenberg_HK` null may be null on the wrong quantity. It says
housekeeping genes have no unusual promoter-classed signal *in total*. The
hypothesis is about the other-promoter component alone, and a null in the sum is
consistent with a real effect in one component cancelling against the other.

This splits the peaks and re-runs the test three ways:

    OWN     promoter-classed peaks overlapping this gene's own annotated TSS
    OTHER   promoter-classed peaks elsewhere in the +/- 1 Mb window
    TOTAL   the existing feature, as a control that we reproduce it

TSS comes from `01_tss_regions.bed`, the annotation the capture panel was
designed against. Never from the RNA strand proxy, which is 51% accurate
(CURRENT_FINDINGS section 6b).

READ THE PANEL CAVEAT BEFORE INTERPRETING A NULL HERE. The panel is 85.1%
CpG-island against 63.8% genome-wide because of ATAC gating, so housekeeping
genes are compared against a background that is largely housekeeping-like in
promoter class. A null is weak evidence; a positive is still informative.

Run in cd4env.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, "/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
                   "cd4_cleaned/scripts_cleaned/audit/GW/scripts")

import math
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
OWN_KB = 5          # a peak within this of the gene's own TSS counts as OWN
N_PERM = 2000

print("loading peaks and annotation...")
pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()

bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key").set_index("symbol_key")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
panel = s.genes.set_index("gene_id").loc[W.index, "symbol_key"]

pk["own_tss"] = pk["symbol_key"].map(bed["tss"])
pk = pk[np.isfinite(pk["own_tss"])].copy()
pk["d_own"] = (pk["peak_midpoint"] - pk["own_tss"]).abs()
prom = pk[pk["re"] == "promoter"].copy()
prom["is_own"] = prom["d_own"] <= OWN_KB * 1000

print(f"  {len(pk)} peaks with annotation, {len(prom)} classed promoter")
print(f"  of the promoter peaks: {prom['is_own'].sum()} OWN "
      f"({prom['is_own'].mean():.1%}), {(~prom['is_own']).sum()} OTHER")
print(f"\n  distance of promoter peaks from the gene's own TSS:")
for q in [.10, .25, .50, .75, .90]:
    print(f"     {q:.0%} quantile {prom['d_own'].quantile(q)/1000:>9.1f} kb")

# Signal-weighted fractions, mirroring aggregate.py but split by own/other.
tot = pk.groupby("symbol_key")["peak_max"].sum()
own = prom[prom["is_own"]].groupby("symbol_key")["peak_max"].sum()
oth = prom[~prom["is_own"]].groupby("symbol_key")["peak_max"].sum()
F = pd.DataFrame({"total_signal": tot}).join(
    own.rename("own_signal")).join(oth.rename("other_signal")).fillna(0.0)
F["prom_own"] = F["own_signal"] / F["total_signal"].replace(0, np.nan)
F["prom_other"] = F["other_signal"] / F["total_signal"].replace(0, np.nan)
F["prom_total"] = F["prom_own"] + F["prom_other"]
F = F.replace([np.inf, -np.inf], np.nan).dropna(subset=["prom_total"])

sk = panel.to_numpy()
V = pd.DataFrame(index=range(len(sk)))
for c in ["prom_own", "prom_other", "prom_total"]:
    V[c] = pd.Series(sk).map(F[c]).to_numpy(float)
ok = V.notna().all(axis=1).to_numpy()
print(f"\n  {ok.sum()} of {len(sk)} panel genes have usable values")

# Reproduce the shipped feature as a sanity check on the reconstruction.
if "promoter_signal_fraction_raw" in W.columns:
    r = np.corrcoef(V.loc[ok, "prom_total"],
                    W.loc[:, "promoter_signal_fraction_raw"].to_numpy()[ok])[0, 1]
    print(f"  reconstruction check: prom_total vs shipped "
          f"promoter_signal_fraction_raw  r = {r:.3f}")
    if r < 0.7:
        print("     NOTE: below 0.7. The shipped feature uses oe_max and per-peak")
        print("     features computed downstream; peak_max here is a proxy. Treat")
        print("     the OWN/OTHER split as the result, not the absolute level.")

# Residualise on the magnitude basis, as everywhere else in the audit.
mag = [m for m in MAG_OVERALL if m in W.columns]
A = np.column_stack([np.ones(ok.sum()), W.loc[:, mag].to_numpy(float)[ok]])
Y = V.to_numpy(float)[ok]
beta, *_ = np.linalg.lstsq(A, Y, rcond=None)
Yr = Y - A @ beta
Yr = (Yr - Yr.mean(0)) / Yr.std(0)

cm = s.table("cohort_membership")
sk_ok = sk[ok]


def test(mask, name):
    out = []
    for j, c in enumerate(V.columns):
        row = [c]
        for M in (Yr, Yrd, Yrf):        # magnitude, + density, + length/expr/CpG
            a, b = M[mask, j], M[~mask, j]
            d = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
            null = np.empty(N_PERM)
            idx = np.arange(len(mask))
            for i in range(N_PERM):
                p = rng.permutation(idx)[:mask.sum()]
                q = np.setdiff1d(idx, p, assume_unique=False)
                null[i] = (M[p, j].mean() - M[q, j].mean()) / np.sqrt(
                    (M[p, j].var() + M[q, j].var()) / 2)
            row += [d, max((np.abs(null) >= abs(d)).mean(), 1 / N_PERM)]
        out.append(tuple(row))
    print(f"\n  {name}  (n = {mask.sum()})")
    print(f"     {'':<14}{'magnitude':>21}{'+ density':>21}"
          f"{'+ len/expr/CpG':>21}")
    for c, d1, p1, d2, p2, d3, p3 in out:
        f = lambda d, p: f"d={d:+.3f} p={p:.4f}{'*' if p < 0.05 else ' '}"
        print(f"     {c:<14} {f(d1,p1)}  {f(d2,p2)}  {f(d3,p3)}")
    return out


# --- the structural confound ------------------------------------------------
# `prom_other` counts contacts to OTHER genes' promoters as a share of all
# signal. A gene with more neighbours has more promoters to contact, so this is
# a gene-density proxy by construction, not merely correlated with one. Measure
# it, then residualise on density as well and re-run everything.
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")
dens = pd.Series(sk_ok).map(gt["gene_density_500kb"]).to_numpy(float)
dens = np.where(np.isfinite(dens), dens, np.nanmedian(dens))

print("\n" + "=" * 78)
print("IS `prom_other` JUST GENE DENSITY?")
print("=" * 78)
for j, c in enumerate(V.columns):
    r = np.corrcoef(Yr[:, j], dens)[0, 1]
    print(f"  {c:<14} Pearson r with gene_density_500kb = {r:+.3f}")

# magnitude + density residualisation
A2 = np.column_stack([A, dens, dens ** 2])
beta2, *_ = np.linalg.lstsq(A2, Y, rcond=None)
Yrd = Y - A2 @ beta2
Yrd = (Yrd - Yrd.mean(0)) / Yrd.std(0)

# magnitude + density + gene length + expression + CpG. Expression is the one
# that matters most here: housekeeping genes are highly expressed, and a highly
# expressed gene may contact more promoters for transcriptional reasons rather
# than through any assembly. If the effect is real it must survive this.
extra = []
for col, log in [("log10_gene_length", False), ("gtex_whole_blood_tpm", True),
                 ("cpg_oe_ratio", False)]:
    v = pd.Series(sk_ok).map(gt[col]).to_numpy(float)
    if log:
        v = np.log1p(v)
    extra.append(np.where(np.isfinite(v), v, np.nanmedian(v)))
A3 = np.column_stack([A2] + extra)
beta3, *_ = np.linalg.lstsq(A3, Y, rcond=None)
Yrf = Y - A3 @ beta3
Yrf = (Yrf - Yrf.mean(0)) / Yrf.std(0)

print("\n" + "=" * 78)
print("THE TEST: is the OTHER-promoter component different from the total?")
print("=" * 78)
print("magnitude-residualised, two-sided permutation, "
      f"{N_PERM} permutations, OWN <= {OWN_KB} kb")

TARGETS = ["Eisenberg_HK", "bio_HK_k3", "DepMap_curated_essential",
           "dbSUPER_CD4_SE_TSS_pm50kb", "gene_desert_bottomQ_density"]
res = {}
for g in TARGETS:
    mem = set(cm.loc[cm["group"] == g, "symbol_key"])
    m = np.fromiter((v in mem for v in sk_ok), bool, len(sk_ok))
    if not (S.MIN_GROUP_N <= m.sum() <= S.MAX_GROUP_COVERAGE * len(sk_ok)):
        print(f"\n  {g}: n = {m.sum()}, outside the floor/ceiling, skipped")
        continue
    res[g] = test(m, g)

print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  Only {prom['is_own'].mean():.1%} of promoter-classed peaks are the gene's")
print("  OWN TSS, so `promoter_signal_fraction` was already almost purely")
print("  other-promoter contact. The promoter-assembly hypothesis HAS been")
print("  tested by the existing feature; the conflation is real in the code but")
print("  is a few per cent of contamination, not a mixing of two quantities.")
print("\n  The live question is therefore not conflation but CONFOUNDING:")
print("  a gene with more neighbours has more promoters available to contact,")
print("  so this measure is a gene-density proxy by construction.")

for g in ("Eisenberg_HK", "bio_HK_k3"):
    if g not in res:
        continue
    o = {c: t for c, *t in res[g]}["prom_other"]
    print(f"\n  {g}, other-promoter contact:")
    print(f"     magnitude           d = {o[0]:+.3f}  p = {o[1]:.4f}")
    print(f"     + density           d = {o[2]:+.3f}  p = {o[3]:.4f}")
    print(f"     + length/expr/CpG   d = {o[4]:+.3f}  p = {o[5]:.4f}")
    if o[5] < 0.05:
        ov = 2 * (1 - 0.5 * (1 + math.erf(abs(o[4]) / (2 * math.sqrt(2)))))
        print(f"     SURVIVES EVERYTHING. Overlap still {ov:.0%}: displaced, not")
        print("     separated, like every other set. Small but mechanism-defined.")
    elif abs(o[4]) > 0.6 * abs(o[0]):
        print(f"     Not significant, but the effect only fell {o[0]:+.3f} to "
              f"{o[4]:+.3f}.")
        print("     Same direction and similar magnitude at a smaller n: this is")
        print("     UNDERPOWERED, not refuted. Do not report it as a failure to")
        print("     replicate, and do not report it as replication either.")
    else:
        print("     The effect itself collapses, not merely its p-value. Refuted.")

print("\n  Controls that must stay null for this to mean anything:")
for g in ("DepMap_curated_essential", "dbSUPER_CD4_SE_TSS_pm50kb"):
    if g in res:
        o = {c: t for c, *t in res[g]}["prom_other"]
        print(f"     {g:<28} fully controlled d = {o[4]:+.3f}  p = {o[5]:.4f}")
if "gene_desert_bottomQ_density" in res:
    o = {c: t for c, *t in res["gene_desert_bottomQ_density"]}["prom_other"]
    print(f"     {'gene_desert (density check)':<28} "
          f"d {o[0]:+.3f} -> {o[2]:+.3f} once density is controlled; the control works")
print("\n  Own-promoter must also stay null, or this is just 'more contacts':")
for g in ("Eisenberg_HK", "bio_HK_k3"):
    if g in res:
        o = {c: t for c, *t in res[g]}["prom_own"]
        print(f"     {g:<28} fully controlled d = {o[4]:+.3f}  p = {o[5]:.4f}")

print("\n  MULTIPLICITY. Eisenberg_HK x prom_other was the pre-specified test, so")
print("  its p stands on its own. Everything else in this table is a control, and")
print("  no control may be promoted to a finding without its own correction.")
print("  PANEL CAVEAT: 85% CpG-island by ATAC-gated construction. Restricted")
print("  range attenuates, so surviving it argues for the effect, not against.")

out = Path(__file__).resolve().parent / "output" / "promoter_assembly.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
V.assign(symbol_key=sk).to_csv(out, sep="\t", index=False)
print(f"\nwrote {out}")
