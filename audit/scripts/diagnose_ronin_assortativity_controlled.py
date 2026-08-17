"""The two controls the raw assortativity result needs.

`diagnose_ronin_assortativity.py` found that target genes contact target
promoters +0.27 above neighbourhood availability, BUT that non-target genes do
the same at +0.21. About 80% of the apparent homophily is therefore not
homophily: MCC contact peaks preferentially land on active, accessible,
CpG-island promoters, and those are exactly the promoters THAP11 binds, so any
gene's contacts look RONIN-enriched. Only the ~+0.05 target-versus-non-target
excess was ever a candidate finding.

Two controls, attacking that residual from opposite directions.

CONTROL A, ACTIVITY-MATCHED EXPECTATION.
The per-gene `exp` term adjusts for how many targets are AVAILABLE nearby but not
for how ACTIVE they are. Recompute it by reweighting the available promoters to
match the expression profile of the promoters actually contacted:

    exp_matched = sum over deciles d of
                  P(contacted partner in d) * target rate among AVAILABLE in d

If the residual is "target genes sit in more active neighbourhoods", this removes
it. Expression is GTEx whole blood, available for 89% of annotated TSSs.

CONTROL B, WITHIN-CONTACT SIGNAL.
Sidesteps the peak-calling bias entirely by never comparing contacted against
uncontacted. Among the promoters a gene ACTUALLY contacts, do the target ones
receive more contact signal than the non-target ones? Every promoter in the
comparison already cleared the accessibility bar, so that bar cannot drive the
result. Signal is ranked within gene, because genes differ in overall depth.

Assembly predicts a positive effect in BOTH, and specifically a LARGER one for
target genes than for non-target genes. The non-target control is reported for
both, and remains the number to read first.

Run in cd4env.
"""
import sys
import gzip
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata
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
TSS_WINDOW, WINDOW_BP, OWN_KB, N_PERM, N_DEC = 1000, 1_000_000, 5, 2000, 10

bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key").reset_index(drop=True)

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")

expr = bed["symbol_key"].map(gt["gtex_whole_blood_tpm"]).to_numpy(float)
expr = np.log1p(np.where(np.isfinite(expr), expr, np.nanmedian(expr)))
bed["act_dec"] = pd.qcut(expr, N_DEC, labels=False, duplicates="drop")
ACT = dict(zip(bed["symbol_key"], bed["act_dec"]))
print(f"{len(bed)} TSSs, activity deciles from GTEx whole blood")


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
        st, en = np.sort(p["start"].to_numpy()), np.sort(p["end"].to_numpy())
        for symk, t in zip(sub["symbol_key"], sub["tss"]):
            i = np.searchsorted(st, t + TSS_WINDOW)
            j = max(0, i - 50)
            if ((st[j:i] - TSS_WINDOW <= t) & (en[j:i] + TSS_WINDOW >= t)).any():
                hit.add(symk)
    return hit


TARGETS = {f.name.replace(".bed.gz", ""): bound(f) for f in sorted(BASE.glob("*.bed.gz"))}
TARGETS["RONIN_module_all3"] = (TARGETS["THAP11_HepG2"] & TARGETS["HCFC1_GM12878"]
                                & TARGETS["ZNF143_K562"])

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
panel_sk = set(s.genes.set_index("gene_id").loc[W.index, "symbol_key"])
pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()
pk = pk[(pk["re"] == "promoter") & pk["symbol_key"].isin(panel_sk)].copy()

tss_by_chrom = {c: sub.sort_values("tss").reset_index(drop=True)
                for c, sub in bed.groupby("chrom")}
own_tss, own_chr = bed.set_index("symbol_key")["tss"], bed.set_index("symbol_key")["chrom"]

print("mapping peaks to partner TSSs...")
rows = []
for symk, g in pk.groupby("symbol_key"):
    if symk not in own_chr.index:
        continue
    c, t0 = own_chr[symk], own_tss[symk]
    T = tss_by_chrom.get(c)
    if T is None:
        continue
    tv, names = T["tss"].to_numpy(), T["symbol_key"].to_numpy()
    partners, sig = [], []
    for st, en, w in zip(g["start"], g["end"], g["peak_max"]):
        lo, hi = np.searchsorted(tv, [st, en])
        for k in range(lo, hi):
            if abs(tv[k] - t0) <= OWN_KB * 1000:
                continue
            partners.append(names[k])
            sig.append(float(w))
    if len(partners) < 4:
        continue
    lo, hi = np.searchsorted(tv, [t0 - WINDOW_BP, t0 + WINDOW_BP])
    avail = [n for n, v in zip(names[lo:hi], tv[lo:hi]) if abs(v - t0) > OWN_KB * 1000]
    if len(avail) < 5:
        continue
    rows.append(dict(sk=symk, partners=np.array(partners), sig=np.array(sig),
                     avail=np.array(avail)))
print(f"{len(rows)} panel genes usable\n")

dens = np.array([gt["gene_density_500kb"].get(r["sk"], np.nan) for r in rows], float)
dens = np.where(np.isfinite(dens), dens, np.nanmedian(dens))
dec = pd.qcut(dens, N_DEC, labels=False, duplicates="drop")
p_act = [np.array([ACT.get(p, 0) for p in r["partners"]]) for r in rows]
a_act = [np.array([ACT.get(a, 0) for a in r["avail"]]) for r in rows]


def controls(tset):
    """(A) activity-matched enrichment, (B) within-contact signal preference."""
    A, B = np.full(len(rows), np.nan), np.full(len(rows), np.nan)
    for i, r in enumerate(rows):
        pt = np.fromiter((p in tset for p in r["partners"]), bool, len(r["partners"]))
        at = np.fromiter((a in tset for a in r["avail"]), bool, len(r["avail"]))
        # A: reweight AVAILABLE to the activity profile of the CONTACTED
        exp_m, tot = 0.0, 0
        for d in np.unique(p_act[i]):
            sel = a_act[i] == d
            if sel.sum() == 0:
                continue
            w = (p_act[i] == d).sum()
            exp_m += w * at[sel].mean()
            tot += w
        if tot:
            A[i] = pt.mean() - exp_m / tot
        # B: within contacted promoters only, rank signal, target vs non-target
        if pt.sum() >= 2 and (~pt).sum() >= 2:
            rk = rankdata(r["sig"]) / len(r["sig"])
            B[i] = rk[pt].mean() - rk[~pt].mean()
    return A, B


def test(stat, is_t):
    ok = np.isfinite(stat)
    a, b = stat[ok & is_t], stat[ok & ~is_t]
    if len(a) < 25 or len(b) < 25:
        return None
    d = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
    idx = np.arange(len(rows))
    nul = np.empty(N_PERM)
    for i in range(N_PERM):
        perm = np.zeros(len(rows), bool)
        for dd in np.unique(dec):
            m = idx[dec == dd]
            k = int(is_t[m].sum())
            if k:
                perm[rng.choice(m, k, replace=False)] = True
        x, y = stat[ok & perm], stat[ok & ~perm]
        nul[i] = (x.mean() - y.mean()) / np.sqrt((x.var() + y.var()) / 2)
    return a.mean(), b.mean(), d, max((np.abs(nul) >= abs(d)).mean(), 1 / N_PERM), len(a)


print("=" * 94)
print("CONTROL A: activity-matched expectation   |   CONTROL B: within-contact signal")
print("=" * 94)
print(f"{'factor':<20}{'n':>6}| {'A tgt':>7}{'A ctrl':>8}{'A d':>7}{'A p':>8} | "
      f"{'B tgt':>7}{'B ctrl':>8}{'B d':>7}{'B p':>8}")
out = []
for name, tset in TARGETS.items():
    is_t = np.array([r["sk"] in tset for r in rows])
    A, B = controls(tset)
    ra, rb = test(A, is_t), test(B, is_t)
    if ra is None or rb is None:
        print(f"{name:<20} too few, skipped")
        continue
    out.append(dict(factor=name, n=ra[4], A_tgt=ra[0], A_ctrl=ra[1], A_d=ra[2],
                    A_p=ra[3], B_tgt=rb[0], B_ctrl=rb[1], B_d=rb[2], B_p=rb[3]))
    print(f"{name:<20}{ra[4]:>6}| {ra[0]:>+7.3f}{ra[1]:>+8.3f}{ra[2]:>+7.3f}"
          f"{ra[3]:>8.4f}{'*' if ra[3] < .05 else ' '}| {rb[0]:>+7.3f}{rb[1]:>+8.3f}"
          f"{rb[2]:>+7.3f}{rb[3]:>8.4f}{'*' if rb[3] < .05 else ' '}")

df = pd.DataFrame(out)
o = Path(__file__).resolve().parent / "output" / "ronin_assortativity_controlled.tsv"
o.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(o, sep="\t", index=False)

print("\n" + "=" * 94)
print("VERDICT")
print("=" * 94)
print("  Read the ctrl columns first. They are the non-target genes, and they")
print("  carried ~80% of the uncontrolled result.")
if len(df):
    nA = int((df.A_p < 0.05).sum())
    nB = int((df.B_p < 0.05).sum())
    print(f"\n  Control A: {nA} of {len(df)} factors significant")
    print(f"  Control B: {nB} of {len(df)} factors significant")
    both = df[(df.A_p < 0.05) & (df.B_p < 0.05)]
    print(f"  BOTH     : {len(both)} of {len(df)}"
          + (f"  -> {', '.join(both.factor)}" if len(both) else ""))
    if len(both) == 0:
        print("\n  The residual does not survive. The uncontrolled assortativity was")
        print("  contact-peak activity bias, and there is no evidence here that")
        print("  RONIN targets assemble with each other in CD4 contact space.")
    elif len(both) == len(df):
        print("\n  Survives both, in every factor. Note this is still an ENCODE")
        print("  binding proxy in the wrong cell types; it is not the published")
        print("  functional target list.")
    else:
        print("\n  Mixed. Report per factor, and treat the stringent all-3 module")
        print("  as the one that matters: it is the least contaminated definition.")
print(f"\nwrote {o}")
