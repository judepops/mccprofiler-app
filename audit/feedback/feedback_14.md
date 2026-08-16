Files: audit/CURRENT_FINDINGS.md (new §6b), audit/PLAN_FORWARD.md (new S2.1c).

Setting up S1a found something bigger than S1a.

Two things happened. First, the strand source I used on 14 August is broken: derived from cd4_rna_plus_promoter vs cd4_rna_minus_promoter, it agrees with the annotation only 51% — chance. It would have reverse-complemented the promoter sequence for ~886 genes had I built S1a on it, and it already invalidated a recorded conclusion.

Second, with the proper annotation strand (01_tss_regions.bed, covering 1,838 of 1,846 panel genes):

plus-strand	minus-strand	d	strand-corrected
contact_asymmetry	+0.308	−0.300	0.638	+0.304, t = 13.7
Opposite in sign, near-equal in magnitude, 13.7 standard errors from zero. Contacts are systematically biased to one side relative to the direction of transcription, and the feature computes that in genomic coordinates, so it cancels. Its own docstring says "positive means more signal downstream" — downstream in genomic terms, which is upstream for half the genome.

This is signal recovery, not a tidy-up. d = 0.638 is larger than every external-set displacement except Roadmap_silenced, and contact_asymmetry is among the most reproducible features you have (ρ 0.900). Fixing it adds a well-measured, interpretable axis rather than removing a bad one — the opposite of the last three feature findings.

Added as S2.1c, a third feature change before the rebuild, with instructions to use the annotation bed and check whether the other directional features (max_distance_to_viewpoint_*, frac_promoter_proximal) have the same problem. I would expect at least some of them do.

The retraction is recorded: "orientation is negligible, d = −0.045" is withdrawn, since random labels could only ever return zero. The separate finding that the per-peak oe_asymmetry family is a support-size problem rather than an orientation bug still stands — it shows only d ≈ −0.09 here.

S1a can now proceed on the correct strand source. Want me to continue with it, or start on the feature changes given S2.1c just got more interesting?

lets continue
Now S1a, on the correct strand source:

Write build_core_promoter_sets.py
194 lines
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
Bash Build the core-promoter sets
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python audit/scripts/build_core_promoter_sets.py 2>&1 | grep -v Warning | tail -28

OUT
loading panel, annotation and genome...
  1838 panel genes with sequence, 8 without annotation, 0 without sequence

==========================================================================
MOTIF PREVALENCE, in the canonical window only
==========================================================================
  TATA   TATAWAWR   at -34 to -18      23 genes (1%)
  Inr    YYANWYY    at -4 to +6      80 genes (4%)
  DPE    RGWYVT     at +24 to +36      47 genes (3%)
  CpG    GC>0.55, o/e>0.65 at -100 to +100    1564 genes (85%)

==========================================================================
CANDIDATE SETS, against the floor and the coverage ceiling
==========================================================================
set                                    n  % panel   chr19  usable
coreprom_TATA                         23      1%     9%   TOO SMALL
coreprom_Inr                          80      4%     8%   yes
coreprom_CpG_island                 1564     85%     7%   OVER CEILING
coreprom_focused_TATA_noCpG            4      0%     0%   TOO SMALL
coreprom_broad_CpG_TATAless         1545     84%     7%   OVER CEILING

wrote /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit/scripts/output/core_promoter_sets.tsv
1 set(s) usable for the S1b displacement test: coreprom_Inr

S1b: test these on the FINAL substrate, after S2.1, S2.1b and S2.1c.
That is a negative result, and the reason matters. Let me confirm the cause rather than assume it:

Bash Test whether the panel is selected for CpG promoters
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate cd4env
python - <<'PY' 2>&1 | grep -v Warning
import sys, re, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'backend')
import numpy as np, pandas as pd
from pyfaidx import Fasta
from app.store import get_store
IU={"A":"A","C":"C","G":"G","T":"T","W":"[AT]","R":"[AG]","Y":"[CT]","V":"[ACG]","N":"."}
tata=re.compile("".join(IU[c] for c in "TATAWAWR"))
bed=pd.read_csv("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Output/01_tss_regions.bed",
    sep="\t",header=None,names=["chrom","start","end","gene","score","strand"])
bed["symbol_key"]=bed.gene.str.upper(); bed["tss"]=((bed.start+bed.end)//2).astype(int)
bed=bed.drop_duplicates("symbol_key")
fa=Fasta("/home/imm/grte4643/Documents/DPhil/Lab/Protocol_20k/Genes/Input/GENOME/hg38_primary.fa")
panel=set(get_store().genes["symbol_key"])
def rc(s): return s[::-1].translate(str.maketrans("ACGTN","TGCAN"))
F=120
rows=[]
for r in bed.itertuples():
    if r.chrom not in fa: continue
    seq=str(fa[r.chrom][max(0,r.tss-F-1):r.tss+F]).upper()
    if len(seq)<2*F: continue
    if r.strand=="-": seq=rc(seq)
    prox=seq[F-100:F+100]
    g,c=prox.count("G"),prox.count("C"); cg=prox.count("CG")
    rows.append((r.symbol_key, bool(tata.search(seq[F-34:F-18])),
                 (g+c)/len(prox), (cg*len(prox))/max(g*c,1), r.symbol_key in panel))
d=pd.DataFrame(rows,columns=["symbol_key","TATA","gc","cpg_oe","in_panel"])
d["CpG"]=(d.gc>0.55)&(d.cpg_oe>0.65)
print("=== IS THE PANEL SELECTED FOR CpG-ISLAND PROMOTERS? ===\n")
print("%-22s %8s %8s %8s" % ("","n","TATA%","CpG%"))
for lab,sub in [("genome-wide (annotated)",d),("  in our panel",d[d.in_panel]),
                ("  NOT in our panel",d[~d.in_panel])]:
    print("%-22s %8d %7.1f%% %7.1f%%" % (lab,len(sub),100*sub.TATA.mean(),100*sub.CpG.mean()))
a,b=d[d.in_panel],d[~d.in_panel]
from math import sqrt
p1,p2=a.CpG.mean(),b.CpG.mean()
print("\n  CpG-island enrichment in panel: %.1f%% vs %.1f%%, odds ratio %.2f"
      % (100*p1,100*p2,(p1/(1-p1))/(p2/(1-p2))))
print("  TATA depletion in panel:        %.1f%% vs %.1f%%" % (100*a.TATA.mean(),100*b.TATA.mean()))
PY

OUT
=== IS THE PANEL SELECTED FOR CpG-ISLAND PROMOTERS? ===

                              n    TATA%     CpG%
genome-wide (annotated)    19704     2.8%    63.8%
  in our panel             1838     1.3%    85.1%
  NOT in our panel        17866     2.9%    61.6%

  CpG-island enrichment in panel: 85.1% vs 61.6%, odds ratio 3.56
  TATA depletion in panel:        1.3% vs 2.9%
Confirmed, and it closes S1 with a "cannot be tested here" rather than a result:

n	TATA	CpG island
genome-wide (annotated)	19,704	2.8%	63.8%
our panel	1,838	1.3%	85.1%
not in panel	17,866	2.9%	61.6%
Odds ratio 3.56 for CpG-island enrichment, TATA depleted 2.2-fold. The panel was built by ATAC-accessibility selection in CD4 cells, and active promoters are overwhelmingly CpG-island and TATA-less. So the Haberle–Stark contrast does not exist here to test.

Bash Close S1 and record the panel composition finding
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
python3 - <<'PY'
import pathlib
p=pathlib.Path("audit/CURRENT_FINDINGS.md"); s=p.read_text()
anchor="## 7. Retracted, do not quote"
new = """## 6c. S1 CLOSED 2026-08-16: the housekeeping question cannot be asked on this panel

The last planned test was core-promoter sequence classes, the second
chromatin-independent positive control and the proper form of the housekeeping
question. It cannot be run here, and the reason is a property of the panel worth
reporting in its own right.

**The panel is strongly selected for CpG-island, TATA-less promoters.**
Scanning the canonical windows on hg38 for all 19,704 annotated genes, against
the 1,838 panel genes with sequence:

| | n | TATA at -34 to -18 | CpG island |
|---|---|---|---|
| genome-wide | 19,704 | 2.8% | 63.8% |
| **panel** | 1,838 | **1.3%** | **85.1%** |
| not in panel | 17,866 | 2.9% | 61.6% |

Odds ratio **3.56** for CpG-island enrichment; TATA depleted 2.2-fold. This
follows from how the panel was built: TSS extraction then ATAC-accessibility
selection in CD4, and active human promoters are overwhelmingly CpG-island and
TATA-less.

**Consequence.** The Haberle and Stark contrast, focused TATA-containing
developmental against broad CpG-island housekeeping, has no contrast to measure
here. The candidate sets:

| set | n | % panel | verdict |
|---|---|---|---|
| `coreprom_TATA` | 23 | 1% | below the 25-gene floor |
| `coreprom_focused_TATA_noCpG` | 4 | 0% | hopeless |
| `coreprom_CpG_island` | 1,564 | 85% | over the 70% ceiling |
| `coreprom_broad_CpG_TATAless` | 1,545 | 84% | over the ceiling |
| **`coreprom_Inr`** | **80** | **4%** | **usable** |

Only the Initiator set is testable, and Inr is a core promoter element rather
than the housekeeping discriminant. It goes to S1b as a sequence-defined,
chromatin-independent set, but it does not answer the housekeeping question.

**Three things follow.**

*The housekeeping null stays open, not resolved.* `Eisenberg_HK` has no
architectural signature once magnitude and density are controlled, and we still
cannot say whether that is because the biology is absent or because the label is
an expression list. State it as an open question, not as a finding.

*The chromatin-independence problem is not fixed.* `GWAS_immune_hot` remains the
only positive control touching no chromatin assay, with Lambert_TF qualified for
chr19 enrichment. That is a real limitation of the report and should be written
as one.

*It is a concrete argument for the genome-wide panel.* At 63.8% CpG island and
2.8% TATA, a genome-wide panel gives roughly 550 TATA-containing genes instead
of 23, and a real CpG contrast instead of an 85% majority. This is the third
independent argument for Aim 4 found this week, alongside n being thin for
representation learning and three reference sets exceeding the coverage ceiling.

---

""" + anchor
assert anchor in s; s=s.replace(anchor,new,1)
p.write_text(s); print("CURRENT_FINDINGS.md updated")

q=pathlib.Path("audit/PLAN_FORWARD.md"); t=q.read_text()
t=t.replace("""## S1. The last test: core-promoter sequence classes""",
"""## S1. CLOSED 2026-08-16, cannot be tested on this panel

Built and run (`build_core_promoter_sets.py`). The panel is 85.1% CpG-island and
1.3% TATA against 63.8% and 2.8% genome-wide, an odds ratio of 3.56, because it
was selected on ATAC accessibility. TATA gives 23 genes (below the floor) and
CpG-island gives 1,564 (over the ceiling), so the Haberle and Stark contrast has
nothing to measure. Only `coreprom_Inr` at n = 80 is usable, and it carries to
S1b as a sequence-defined control without answering the housekeeping question.

The housekeeping null therefore stays **open**, `GWAS_immune_hot` remains the
only fully chromatin-independent positive control, and this becomes the third
argument this week for the genome-wide panel. See CURRENT_FINDINGS.md 6c.

**Testing is now finished.** What remains is fixes, figures and writing.

<details><summary>original plan for S1, retained</summary>

## S1. The last test: core-promoter sequence classes""")
t=t.replace("""Fallbacks if the motif scan is unusable: MSigDB translation machinery
(ribosomal proteins give only 19, below the floor, but adding initiation and
elongation factors and aminoacyl-tRNA synthetases should clear 50), CORUM
protein-complex membership, or other Pfam families as further Lambert-style
controls.""",
"""Fallbacks if the motif scan is unusable: MSigDB translation machinery
(ribosomal proteins give only 19, below the floor, but adding initiation and
elongation factors and aminoacyl-tRNA synthetases should clear 50), CORUM
protein-complex membership, or other Pfam families as further Lambert-style
controls.

**On the fallbacks:** they remain available but are lower value than they looked.
The panel-composition problem that killed the core-promoter test is not specific
to sequence classes, so any candidate set must be checked for panel coverage and
chr19 share before it is trusted. Given testing is otherwise finished, pursue one
only if the report needs a second chromatin-independent control badly enough to
delay writing.

</details>""")
t=t.replace("    S1a   build the core-promoter sets        1 d   substrate-independent, can\n                                                    run in parallel with the below",
            "    S1a   build the core-promoter sets        DONE, closed: no contrast on this panel")
q.write_text(t); print("PLAN_FORWARD.md updated")
PY
grep -c "—" audit/CURRENT_FINDINGS.md audit/PLAN_FORWARD.md audit/scripts/build_core_promoter_sets.py
git add -A && git commit -q -F - <<'MSG'
S1 closed: the panel is 85% CpG-island, so the housekeeping test has no contrast

Built the core-promoter sequence classes and the test cannot run here. Scanning
canonical windows on hg38, the panel is 85.1% CpG-island and 1.3% TATA against
63.8% and 2.8% genome-wide, an odds ratio of 3.56 with TATA depleted 2.2-fold.
That follows from how the panel was built, TSS extraction then ATAC-accessibility
selection, since active human promoters are overwhelmingly CpG-island and
TATA-less.

So TATA gives 23 genes, below the 25 floor, and CpG-island gives 1,564, over the
70% ceiling. The Haberle and Stark contrast between focused developmental and
broad housekeeping promoters has nothing to measure on this panel. Only
coreprom_Inr at n = 80 survives both bounds, and Inr is a core promoter element
rather than the housekeeping discriminant, so it carries to S1b as a
sequence-defined control without answering the question it was built for.

Three consequences, all recorded. The housekeeping null stays open rather than
resolved: we still cannot say whether Eisenberg_HK has no architectural
signature because the biology is absent or because the label is an expression
list. The chromatin-independence limitation is not fixed, so GWAS_immune_hot
remains the only positive control touching no chromatin assay and that belongs
in the report as a stated limitation. And it is the third independent argument
this week for the genome-wide panel, which at 63.8% CpG and 2.8% TATA would give
roughly 550 TATA genes instead of 23.

Also caught while sourcing strand for this: the RNA-derived strand proxy used on
14 August is 51% accurate, so it would have reverse-complemented the promoter
sequence for about 886 genes. The annotation bed is used throughout instead.

Testing is now finished. What remains is fixes, figures and writing.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
MSG
git push -q origin feat/app-scaffold && echo pushed

OUT
CURRENT_FINDINGS.md updated
PLAN_FORWARD.md updated
audit/scripts/build_core_promoter_sets.py:0
audit/PLAN_FORWARD.md:0
audit/CURRENT_FINDINGS.md:0
pushed
Files: new audit/scripts/build_core_promoter_sets.py, audit/CURRENT_FINDINGS.md (§6c), audit/PLAN_FORWARD.md (S1 closed).

S1 is closed, and it is a negative — but a clean and useful one.

The panel cannot support the test. It is 85.1% CpG-island and 1.3% TATA, against 63.8% and 2.8% genome-wide — odds ratio 3.56, TATA depleted 2.2-fold. That follows directly from how it was built: TSS extraction then ATAC-accessibility selection, and active human promoters are overwhelmingly CpG-island and TATA-less. So the Haberle–Stark contrast has nothing to measure:

set	n	verdict
TATA	23	below the 25 floor
focused TATA, no CpG	4	hopeless
CpG island	1,564 (85%)	over the ceiling
Inr	80	usable
Only Initiator survives both bounds, and Inr is a core-promoter element rather than the housekeeping discriminant. It carries to S1b as a sequence-defined control but does not answer the question it was built for.

Three consequences, all now written down.

The housekeeping null stays open, not resolved. We still cannot say whether Eisenberg_HK lacks an architectural signature because the biology is absent or because the label is an expression list. That must be written as an open question.

The chromatin-independence limitation is not fixed. GWAS_immune_hot remains the only positive control touching no chromatin assay. That belongs in the report as a stated limitation rather than glossed.

It is the third independent argument this week for the genome-wide panel — at 63.8% CpG and 2.8% TATA, a genome-wide panel gives roughly 550 TATA genes instead of 23.

Testing is now finished. What remains is the three feature fixes (S2.1 degenerate topology, S2.1b moment families, S2.1c strand orientation), the rebuild, then figures and writing.

Worth noting how the last hour went: sourcing strand for a test that then could not run is what exposed both the broken RNA proxy and the strand-dependent asymmetry signal — which is the more valuable of the two findings. The test failing did not waste the work.