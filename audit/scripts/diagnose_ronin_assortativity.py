"""Do RONIN-target genes contact OTHER RONIN targets' promoters specifically?

This replaces `diagnose_ronin_targets.py`, which asked the wrong question.

THE ERROR IN THE PREVIOUS TEST. Dejosez et al. 2023 claim housekeeping promoters
cluster *with each other* into multi-promoter assemblies. That is a claim about
PARTNER IDENTITY. `promoter_other` measures the share of contact signal landing
in promoter-classed peaks and does not care whose promoter it is, so a gene
contacting random neighbouring promoters and a gene assembled with other RONIN
targets score identically. No per-gene summary feature can express assembly,
because assembly is pairwise. The whole feature set has this limitation.

THE TEST THAT DOES ADDRESS IT. For each gene, take its promoter-classed contact
peaks, ask which genes' TSSs they land on, and compute the share of those
partners that are RONIN targets. Compare against what the neighbourhood makes
available: the share of RONIN targets among ALL TSSs in the same +/- 1 Mb window.

    observed_g = share of CONTACTED promoters that are RONIN targets
    expected_g = share of AVAILABLE promoters in the window that are targets
    enrichment = observed - expected

This is a homophily/assortativity statistic. The `expected` term is what makes it
a real test: RONIN targets are common (THAP11 binds most CpG-island promoters)
and concentrated in gene-dense regions, so an unadjusted `observed` would be high
for reasons that have nothing to do with assembly.

WHY THE 70% COVERAGE CEILING DOES NOT APPLY HERE. That ceiling exists because a
displacement test comparing an 82% group against an 18% group is meaningless.
Assortativity asks about partner composition instead, so a common label costs
POWER (expected sits near observed, compressing the range) but not validity.
THAP11, excluded from the previous test, is testable here and is the factor that
actually is RONIN.

TWO CONTROLS.
  1. Non-target genes must show no preference for target partners. If they do,
     the statistic is measuring neighbourhood composition, not assembly.
  2. Label permutation stratified by gene-density decile, so the null preserves
     the fact that targets cluster in gene-dense regions.

STANDING CAVEATS, unchanged from the previous script: ENCODE binding is a proxy
for the published functional target definition, and none of HepG2 / GM12878 /
K562 is CD4.

Run in cd4env.
"""
import sys
import gzip
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import warnings
import numpy as np
import pandas as pd
from app.store import get_store

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)

BASE = Path(__file__).resolve().parents[1] / "data" / "dejosez2023" / "encode"
TSS_BED = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/"
           "01_tss_regions.bed")
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
PEAKS = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
         "scripts_cleaned/process/output_resgw_10000/utils/peaks/annotated.tsv")
TSS_WINDOW = 1000        # TSS-to-ChIP-peak distance defining a RONIN target
WINDOW_BP = 1_000_000    # the MCC capture window, +/- 1 Mb
OWN_KB = 5
N_PERM = 2000

# --- all genes, not just the panel: contacted promoters are mostly non-panel --
bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key").reset_index(drop=True)
print(f"{len(bed)} annotated TSSs genome-wide")


def bound(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                rows.append((f[0], int(f[1]), int(f[2])))
    pk = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    hit = set()
    for c, sub in bed.groupby("chrom"):
        p = pk[pk["chrom"] == c]
        if p.empty:
            continue
        st = np.sort(p["start"].to_numpy())
        en = np.sort(p["end"].to_numpy())
        for symk, t in zip(sub["symbol_key"], sub["tss"]):
            i = np.searchsorted(st, t + TSS_WINDOW)
            j = max(0, i - 50)
            if ((st[j:i] - TSS_WINDOW <= t) & (en[j:i] + TSS_WINDOW >= t)).any():
                hit.add(symk)
    return hit


TARGETS = {}
for f in sorted(BASE.glob("*.bed.gz")):
    TARGETS[f.name.replace(".bed.gz", "")] = bound(f)
TARGETS["RONIN_module_all3"] = (TARGETS["THAP11_HepG2"]
                                & TARGETS["HCFC1_GM12878"]
                                & TARGETS["ZNF143_K562"])
print("\ngenome-wide target rates (the denominator that matters here):")
for k, v in TARGETS.items():
    print(f"  {k:<20} {len(v):>6} of {len(bed)} genes = {len(v)/len(bed):>5.1%}")

# --- contacted promoters, per panel gene ------------------------------------
s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
panel_sk = set(s.genes.set_index("gene_id").loc[W.index, "symbol_key"])

pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()
pk = pk[(pk["re"] == "promoter") & pk["symbol_key"].isin(panel_sk)].copy()
tss_by_chrom = {c: sub.sort_values("tss").reset_index(drop=True)
                for c, sub in bed.groupby("chrom")}
own_tss = bed.set_index("symbol_key")["tss"]
own_chr = bed.set_index("symbol_key")["chrom"]

print(f"\n{len(pk)} promoter-classed peaks across {pk['symbol_key'].nunique()} panel genes")
print("mapping each peak to the TSSs it overlaps...")

rows = []
for symk, g in pk.groupby("symbol_key"):
    if symk not in own_chr.index:
        continue
    c, t0 = own_chr[symk], own_tss[symk]
    T = tss_by_chrom.get(c)
    if T is None:
        continue
    tv = T["tss"].to_numpy()
    names = T["symbol_key"].to_numpy()

    partners, weights = [], []
    for st, en, w in zip(g["start"], g["end"], g["peak_max"]):
        lo, hi = np.searchsorted(tv, [st, en])
        for k in range(lo, hi):
            if abs(tv[k] - t0) <= OWN_KB * 1000:      # skip the gene's own TSS
                continue
            partners.append(names[k])
            weights.append(float(w))
    if len(partners) < 3:                              # need a few to form a share
        continue

    # what the neighbourhood makes available
    lo, hi = np.searchsorted(tv, [t0 - WINDOW_BP, t0 + WINDOW_BP])
    avail = [n for n, v in zip(names[lo:hi], tv[lo:hi])
             if abs(v - t0) > OWN_KB * 1000]
    if len(avail) < 5:
        continue
    rows.append(dict(symbol_key=symk, partners=partners,
                     weights=np.asarray(weights), avail=avail))

print(f"{len(rows)} panel genes with >=3 contacted promoters and >=5 available\n")

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")
dens = np.array([gt["gene_density_500kb"].get(r["symbol_key"], np.nan) for r in rows],
                dtype=float)
dens = np.where(np.isfinite(dens), dens, np.nanmedian(dens))
dec = pd.qcut(dens, 10, labels=False, duplicates="drop")


def enrichment(tset):
    """observed minus expected share of target partners, per gene."""
    obs = np.array([np.mean([p in tset for p in r["partners"]]) for r in rows])
    exp = np.array([np.mean([a in tset for a in r["avail"]]) for r in rows])
    wobs = np.array([
        float(np.sum(r["weights"] * np.fromiter((p in tset for p in r["partners"]),
                                                bool, len(r["partners"])))
              / max(r["weights"].sum(), 1e-9))
        for r in rows])
    return obs, exp, wobs


print("=" * 90)
print("ASSORTATIVITY: do target genes contact target promoters more than the")
print("neighbourhood makes available?")
print("=" * 90)
print(f"{'factor':<20}{'n tgt':>7}{'obs':>8}{'exp':>8}{'diff':>8}"
      f"{'d':>8}{'p':>9}   control (non-targets)")

out = []
for name, tset in TARGETS.items():
    is_t = np.array([r["symbol_key"] in tset for r in rows])
    if is_t.sum() < 25 or (~is_t).sum() < 25:
        print(f"{name:<20}  n = {is_t.sum()}, too few, skipped")
        continue
    obs, exp, wobs = enrichment(tset)
    diff = obs - exp

    dt, dn = diff[is_t], diff[~is_t]
    d = (dt.mean() - dn.mean()) / np.sqrt((dt.var() + dn.var()) / 2)

    # density-stratified label permutation
    nul = np.empty(N_PERM)
    idx = np.arange(len(rows))
    for i in range(N_PERM):
        perm = np.zeros(len(rows), bool)
        for dd in np.unique(dec):
            m = idx[dec == dd]
            k = int(is_t[m].sum())
            if k:
                perm[rng.choice(m, k, replace=False)] = True
        a, b = diff[perm], diff[~perm]
        nul[i] = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
    p = max((np.abs(nul) >= abs(d)).mean(), 1 / N_PERM)

    out.append(dict(factor=name, n=int(is_t.sum()), obs=dt.mean() + exp[is_t].mean(),
                    exp=exp[is_t].mean(), diff=dt.mean(), d=d, p=p,
                    ctrl=dn.mean(), wdiff=(wobs - exp)[is_t].mean()))
    print(f"{name:<20}{is_t.sum():>7}{obs[is_t].mean():>8.3f}{exp[is_t].mean():>8.3f}"
          f"{dt.mean():>+8.3f}{d:>+8.3f}{p:>9.4f}   {dn.mean():>+7.3f}"
          f"{'  *' if p < .05 else ''}")

df = pd.DataFrame(out)
o = Path(__file__).resolve().parent / "output" / "ronin_assortativity.tsv"
o.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(o, sep="\t", index=False)

print("\n" + "=" * 90)
print("HOW TO READ THIS")
print("=" * 90)
print("  diff  = target genes' contacted-partner target share MINUS what their")
print("          own +/-1 Mb neighbourhood makes available. Assembly predicts > 0.")
print("  ctrl  = the same quantity for NON-target genes. Should be ~0. If it")
print("          tracks diff, the statistic is reading neighbourhood composition.")
print("  p     = density-decile-stratified label permutation, so the null keeps")
print("          targets concentrated in gene-dense regions.")
if len(df):
    best = df.loc[df["p"].idxmin()]
    print(f"\n  strongest: {best['factor']} diff {best['diff']:+.3f} vs control "
          f"{best['ctrl']:+.3f}, d = {best['d']:+.3f}, p = {best['p']:.4f}")
    print(f"  signal-weighted version of the same: {best['wdiff']:+.3f}")
print("\n  This is a PAIRWISE statistic. No per-gene feature in the substrate can")
print("  express it, which is why the earlier promoter_other test could not have")
print("  addressed the assembly claim whatever it returned.")
print(f"\nwrote {o}")
