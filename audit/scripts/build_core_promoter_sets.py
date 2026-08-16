"""S1a: build core-promoter sequence classes, the last independent control.

Why this set and not another. Three of the four gene sets that separate in
contact-architecture space are ChIP-derived (ChromHMM active TSS, bivalent,
Roadmap silenced), and chromatin state and 3D contact are both downstream of the
same biology, so they are partly two assays measuring one thing. After
Lambert_TF was qualified for chr19 enrichment, `GWAS_immune_hot` is the only
positive control that touches no chromatin assay at all. One line of evidence is
not enough.

Core promoter elements are defined by DNA sequence and position relative to the
transcription start site. They require no measurement in any cell type, so they
are independent of both chromatin assays and expression. They also pose the
housekeeping question in its proper form: Haberle and Stark distinguish a
focused, TATA-containing, developmental class from a broad, CpG-island,
TATA-less housekeeping class. `Eisenberg_HK`, an expression list, has no
architectural signature once magnitude and density are controlled. Whether a
sequence-defined housekeeping class does is the open question.

POSITION IS WHAT MAKES THESE MEANINGFUL. A TATA-like string anywhere in a
promoter is noise; one at -31 to -24 relative to the TSS is a core promoter
element. Every motif here is scanned in its canonical window only.

STRAND. Taken from `Lab/Protocol_20k/Genes/Output/01_tss_regions.bed`, the
annotation used to design the capture panel. **Never** from
`cd4_rna_plus_promoter` versus `cd4_rna_minus_promoter`: checked against the
annotation that proxy is 51% accurate, i.e. chance, and it already invalidated
one recorded conclusion (see CURRENT_FINDINGS.md section 6b).

Writes `audit/scripts/output/core_promoter_sets.tsv`, one row per panel gene.
S1b then tests displacement of these sets on the final substrate.

Run in cd4env.
"""
import sys
import re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import warnings
import numpy as np
import pandas as pd
from pyfaidx import Fasta
from app.store import get_store
from app import store_schema as S

warnings.filterwarnings("ignore")

GENOME = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Input/"
          "GENOME/hg38_primary.fa")
TSS_BED = ("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/"
           "01_tss_regions.bed")
GENE_TABLE = ("/home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/"
              "collaboration/data/gene_table.tsv")
OUT = Path(__file__).resolve().parent / "output" / "core_promoter_sets.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)

FLANK = 120                      # bases either side of the TSS to extract

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "W": "[AT]", "S": "[CG]",
         "R": "[AG]", "Y": "[CT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "."}


def iupac_re(consensus: str) -> re.Pattern:
    return re.compile("".join(IUPAC[c] for c in consensus))


# (name, consensus, window relative to TSS at +1). Windows follow Haberle and
# Stark 2018; the TSS base itself is position +1 and there is no position 0.
MOTIFS = [
    ("TATA", "TATAWAWR", (-34, -18)),   # canonical -31 to -24, widened slightly
    ("Inr",  "YYANWYY",  (-4, +6)),     # the A sits at +1
    ("DPE",  "RGWYVT",   (+24, +36)),   # rare in human, reported for completeness
]


def rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def idx(pos: int) -> int:
    """TSS-relative coordinate to index in the extracted window.

    The window runs from -FLANK to +FLANK with no zero, so +1 is the TSS.
    """
    return FLANK + (pos - 1 if pos > 0 else pos)


print("loading panel, annotation and genome...")
s = get_store()
panel = s.genes[["gene_id", "symbol_key"]].copy()

bed = pd.read_csv(TSS_BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene", "score", "strand"])
bed["symbol_key"] = bed["gene"].str.upper()
bed["tss"] = ((bed["start"] + bed["end"]) // 2).astype(int)
bed = bed.drop_duplicates("symbol_key").set_index("symbol_key")

gt = pd.read_csv(GENE_TABLE, sep="\t")
gt["symbol_key"] = gt["gene_symbol"].str.upper()
gt = gt.drop_duplicates("symbol_key").set_index("symbol_key")

fa = Fasta(GENOME)
compiled = [(n, iupac_re(c), w) for n, c, w in MOTIFS]

rows = []
missing_ann = missing_seq = 0
for gid, sk in zip(panel["gene_id"], panel["symbol_key"]):
    if sk not in bed.index:
        missing_ann += 1
        continue
    b = bed.loc[sk]
    chrom, tss, strand = b["chrom"], int(b["tss"]), b["strand"]
    if chrom not in fa:
        missing_seq += 1
        continue
    seq = str(fa[chrom][max(0, tss - FLANK - 1):tss + FLANK]).upper()
    if len(seq) < 2 * FLANK:
        missing_seq += 1
        continue
    if strand == "-":
        seq = rc(seq)

    rec = {"gene_id": gid, "symbol_key": sk, "chrom": chrom,
           "tss": tss, "strand": strand}
    for name, pat, (lo, hi) in compiled:
        sub = seq[idx(lo):idx(hi)]
        rec[name] = bool(pat.search(sub))
    # GC and CpG observed/expected of the proximal promoter, the standard
    # CpG-island discriminant, computed on the same window as the motifs so the
    # two classes are defined on identical sequence.
    prox = seq[idx(-100):idx(+100)]
    g, c = prox.count("G"), prox.count("C")
    cg = prox.count("CG")
    rec["gc"] = (g + c) / max(len(prox), 1)
    rec["cpg_oe"] = (cg * len(prox)) / max(g * c, 1)
    rows.append(rec)

df = pd.DataFrame(rows)
print(f"  {len(df)} panel genes with sequence, "
      f"{missing_ann} without annotation, {missing_seq} without sequence")

# CpG island promoter by the standard Gardiner-Garden criterion, computed here
# rather than taken from gene_table so it is on the same window as the motifs.
df["CpG_island"] = (df["gc"] > 0.55) & (df["cpg_oe"] > 0.65)

print("\n" + "=" * 74)
print("MOTIF PREVALENCE, in the canonical window only")
print("=" * 74)
for name, cons, (lo, hi) in MOTIFS:
    n = int(df[name].sum())
    print(f"  {name:<6} {cons:<10} at {lo:+d} to {hi:+d}   "
          f"{n:>5} genes ({n / len(df):.0%})")
print(f"  {'CpG':<6} {'GC>0.55, o/e>0.65':<10} at -100 to +100   "
      f"{int(df['CpG_island'].sum()):>5} genes ({df['CpG_island'].mean():.0%})")

# The Haberle and Stark classes.
df["focused_TATA"] = df["TATA"] & ~df["CpG_island"]
df["broad_CpG_TATAless"] = df["CpG_island"] & ~df["TATA"]

SETS = {
    "coreprom_TATA": df["TATA"],
    "coreprom_Inr": df["Inr"],
    "coreprom_CpG_island": df["CpG_island"],
    "coreprom_focused_TATA_noCpG": df["focused_TATA"],
    "coreprom_broad_CpG_TATAless": df["broad_CpG_TATAless"],
}

print("\n" + "=" * 74)
print("CANDIDATE SETS, against the floor and the coverage ceiling")
print("=" * 74)
n_panel = len(df)
chrom_all = df["chrom"].to_numpy()
print(f"{'set':<34}{'n':>6}{'% panel':>9}{'chr19':>8}  usable")
keep = []
for name, mask in SETS.items():
    n = int(mask.sum())
    cov = n / n_panel
    ch19 = (chrom_all[mask.to_numpy()] == "chr19").mean() if n else 0.0
    ok = S.MIN_GROUP_N <= n <= S.MAX_GROUP_COVERAGE * n_panel
    flag = "yes" if ok else ("TOO SMALL" if n < S.MIN_GROUP_N else "OVER CEILING")
    warn = "  <-- chr19-enriched, check" if ch19 > 0.15 else ""
    print(f"{name:<34}{n:>6}{cov:>8.0%}{ch19:>7.0%}   {flag}{warn}")
    if ok:
        keep.append(name)

for name in SETS:
    df[name] = SETS[name]
df.to_csv(OUT, sep="\t", index=False)
print(f"\nwrote {OUT}")
print(f"{len(keep)} set(s) usable for the S1b displacement test: {', '.join(keep)}")
print("\nS1b: test these on the FINAL substrate, after S2.1, S2.1b and S2.1c.")
