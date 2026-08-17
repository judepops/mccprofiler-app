"""Do RONIN-module target genes show the promoter-assembly architecture?

The promoter-assembly claim is Dejosez et al., Cell Reports 2023,
doi 10.1016/j.celrep.2023.112505 (NOT Hwang: the wiki page
`2023-Hwang-HousekeepingGeneArchitecture.md` has the wrong authors, resolved
against Crossref on 2026-08-16). Housekeeping genes are said to lack distal
enhancer interactions and instead rely on RONIN/THAP11-assembled promoter
complexes.

Every previous test of this used `Eisenberg_HK`, an EXPRESSION-derived list. That
tests "are highly expressed constitutive genes architecturally distinct", which
is not the published claim. The published claim is about genes actually bound by
the RONIN module. This builds that set from ChIP-seq and tests it directly.

TARGETS are defined by TSS proximity to IDR-thresholded peaks (ENCODE, GRCh38):

    THAP11   HepG2     ENCFF054TEP    the core RONIN factor
    HCFC1    GM12878   ENCFF722QBB    complex partner, closest cell type available
    HCFC1    HepG2     ENCFF485SRU
    ZNF143   K562      ENCFF978YEJ    complex partner

CELL-TYPE MISMATCH IS THE MAIN THREAT AND MUST BE STATED. None of these is CD4.
A cell-type-mismatched external set already produced one null in this project
(Karpinska). The mitigation is that the THAP11/HCFC1/ZNF143 module binds
constitutive housekeeping promoters, which are far more cell-type invariant than
enhancers, and that GM12878 is lymphoid. The overlap with `Eisenberg_HK` is
reported as a transferability check: high overlap argues the binding transfers.

This is NOT the paper's own supplementary gene list. PMC serves a JS challenge on
binary downloads, so supplementary tables 2-5 of PMC10329844 must be fetched by
hand. When they are, add them here as a fifth set and compare.

Run in cd4env.
"""
import sys
import gzip
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

BASE = Path(__file__).resolve().parents[1] / "data" / "dejosez2023" / "encode"
TSS_BED = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/"
           "01_tss_regions.bed")
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
PEAKS = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/"
         "scripts_cleaned/process/output_resgw_10000/utils/peaks/annotated.tsv")
TSS_WINDOW = 1000        # RONIN binds promoters; keep the window tight
N_PERM = 4000

bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key")

s = get_store()
W = s.table("features").pivot(index="gene_id", columns="feature", values="z")
sk_all = s.genes.set_index("gene_id").loc[W.index, "symbol_key"].to_numpy()
panel = bed[bed["symbol_key"].isin(set(sk_all))]
print(f"{len(sk_all)} panel genes, {len(panel)} with a TSS annotation\n")


def bound(path):
    """Panel gene symbols whose TSS lies within TSS_WINDOW of a peak."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                rows.append((f[0], int(f[1]), int(f[2])))
    pk = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    hit = set()
    for c, sub in panel.groupby("chrom"):
        p = pk[pk["chrom"] == c]
        if p.empty:
            continue
        st, en = p["start"].to_numpy(), p["end"].to_numpy()
        order = np.argsort(st)
        st, en = st[order], en[order]
        for symk, t in zip(sub["symbol_key"], sub["tss"]):
            i = np.searchsorted(st, t + TSS_WINDOW)
            j = max(0, i - 50)
            if ((st[j:i] - TSS_WINDOW <= t) & (en[j:i] + TSS_WINDOW >= t)).any():
                hit.add(symk)
    return hit


sets = {}
for f in sorted(BASE.glob("*.bed.gz")):
    sets[f.name.replace(".bed.gz", "")] = bound(f)
    print(f"  {f.name.replace('.bed.gz',''):<18} {len(sets[f.name.replace('.bed.gz','')]):>5} "
          f"panel genes bound")

thap = sets.get("THAP11_HepG2", set())
core = thap & sets.get("HCFC1_GM12878", set()) & sets.get("ZNF143_K562", set())
sets["RONIN_module_all3"] = core
print(f"  {'RONIN_module_all3':<18} {len(core):>5} panel genes bound by all three")

cm = s.table("cohort_membership")
hk = set(cm.loc[cm["group"] == "Eisenberg_HK", "symbol_key"])
print(f"\ntransferability check, overlap with Eisenberg_HK ({len(hk)} genes):")
for k, v in sets.items():
    if v:
        print(f"  {k:<18} {len(v & hk):>4}/{len(v):<5} = {len(v & hk)/len(v):>5.1%} "
              f"(panel base rate {len(hk)/len(sk_all):.1%})")

# --- the architecture measures -------------------------------------------
pk = pd.read_csv(PEAKS, sep="\t")
pk.columns = [c.lower() for c in pk.columns]
pk["symbol_key"] = pk["viewpoint"].str.upper()
pk["own_tss"] = pk["symbol_key"].map(bed.set_index("symbol_key")["tss"])
pk = pk[np.isfinite(pk["own_tss"])].copy()
pk["is_own"] = (pk["peak_midpoint"] - pk["own_tss"]).abs() <= 5000
tot = pk.groupby("symbol_key")["peak_max"].sum()
F = pd.DataFrame({"tot": tot})
for cls in ["promoter", "enhancer", "ctcf"]:
    sub = pk[pk["re"] == cls]
    F[f"{cls}_other"] = sub[~sub["is_own"]].groupby("symbol_key")["peak_max"].sum()
F = F.fillna(0.0)
MEAS = ["promoter_other", "enhancer_other", "ctcf_other"]
for c in MEAS:
    F[c] /= F["tot"].replace(0, np.nan)

V = pd.DataFrame({c: pd.Series(sk_all).map(F[c]).to_numpy(float) for c in MEAS})
ok = V.notna().all(axis=1).to_numpy()
V, sk = V[ok].reset_index(drop=True), sk_all[ok]
Wm = W.reset_index(drop=True)[ok].reset_index(drop=True)

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")
chrom = pd.Series(sk).map(bed.set_index("symbol_key")["chrom"]).fillna("NA").to_numpy()

A = [np.ones(len(sk))] + [Wm[c].to_numpy(float) for c in MAG_OVERALL if c in Wm.columns]
for col, lg in [("gene_density_500kb", 0), ("log10_gene_length", 0),
                ("gtex_whole_blood_tpm", 1), ("cpg_oe_ratio", 0)]:
    v = pd.Series(sk).map(gt[col]).to_numpy(float)
    if lg:
        v = np.log1p(v)
    A.append(np.where(np.isfinite(v), v, np.nanmedian(v)))
A = np.nan_to_num(np.column_stack(A), nan=0.0)
D = pd.get_dummies(pd.Series(chrom)).to_numpy(float)[:, 1:]
AD = np.column_stack([A, D])


def perm(y, X, m):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    r = (r - r.mean()) / r.std()
    d = (r[m].mean() - r[~m].mean()) / np.sqrt((r[m].var() + r[~m].var()) / 2)
    idx = np.arange(len(r))
    nul = np.empty(N_PERM)
    for i in range(N_PERM):
        p = rng.permutation(idx)[:m.sum()]
        q = np.setdiff1d(idx, p)
        nul[i] = (r[p].mean() - r[q].mean()) / np.sqrt((r[p].var() + r[q].var()) / 2)
    return d, max((np.abs(nul) >= abs(d)).mean(), 1 / N_PERM)


print("\n" + "=" * 86)
print("THE PUBLISHED PREDICTION: RONIN targets should be promoter-driven and")
print("enhancer-poor. Fully controlled (magnitude, density, length, expression,")
print("CpG, chromosome).")
print("=" * 86)
rows = []
for name, genes in sets.items():
    m = np.fromiter((v in genes for v in sk), bool, len(sk))
    if not (S.MIN_GROUP_N <= m.sum() <= S.MAX_GROUP_COVERAGE * len(sk)):
        print(f"\n  {name}: n = {m.sum()}, outside floor/ceiling, skipped")
        continue
    print(f"\n  {name}  (n = {m.sum()})")
    for j, c in enumerate(MEAS):
        y = V[c].to_numpy(float)
        d, p = perm(y, AD, m)
        rows.append(dict(set=name, measure=c, n=int(m.sum()), d=d, p=p))
        pred = {"promoter_other": "predicted UP", "enhancer_other": "predicted DOWN",
                "ctcf_other": ""}[c]
        print(f"     {c:<16} d = {d:+.3f}  p = {p:.4f}  {'*' if p < .05 else ' '}  {pred}")

df = pd.DataFrame(rows)
m = len(df)
df = df.sort_values("p").reset_index(drop=True)
df["rank"] = np.arange(1, m + 1)
df["q"] = (df["p"] * m / df["rank"])[::-1].cummin()[::-1].clip(upper=1.0)
out = Path(__file__).resolve().parent / "output" / "ronin_targets.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)

print("\n" + "=" * 86)
print("VERDICT")
print("=" * 86)
sig = df[df["q"] < 0.05]
print(f"  {m} tests, {int((df.p < 0.05).sum())} nominal, {len(sig)} surviving BH q<0.05")
for r in sig.itertuples():
    print(f"     {r.set:<20}{r.measure:<16} d = {r.d:+.3f}  q = {r.q:.4f}")
if not len(sig):
    print("     nothing survives. Read the cell-type-mismatch caveat before")
    print("     concluding the published architecture is absent in CD4.")
print(f"\nwrote {out}")
