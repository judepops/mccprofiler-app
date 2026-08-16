# Where the project stands, 2026-08-16

Written after two days of adversarial review that retracted eleven claims and
cancelled a research phase. Purpose: state what is true now, where it sits in
the DPhil, what the literature says about it, and what to do next.

Numbers here are quoted from `CURRENT_FINDINGS.md`, which is the authoritative
document. `PLAN_FORWARD.md` has the task order. This one is orientation.

---

## 1. What the project asked, and what it found

**The original question.** Do genes fall into distinct regulatory archetypes
based on their 3D contact architecture, and can those archetypes be discovered
without supervision?

**The answer is no, and it is answered rather than abandoned.** That distinction
is the whole value of the position. The landscape is continuous, structure is
dimensional rather than partitional, and external regulatory categories are
displaced within that space but never separated from it.

### The three results, in the order they should be presented

**The instrument** (formerly Aim 1, now Methods). MCCProfiler: 91
position-invariant features from base-pair-resolution MCC. Justified by what it
measures, not by what it predicts. Gene length, expression and CpG density
explain **2.6%** of the feature variance and **0.6%** of the retained
components; no feature exceeds R2 0.25 and no component exceeds 0.10. That
unexplained variance reproduces (median cross-capture rho **0.752**, 19
components above a parallel-analysis noise ceiling, kappa **0.72** for the
imposed partition). Orthogonal, reproducible, and non-random: all three together
are the argument, and any one alone is weak.

**Result 1: the landscape is continuous.** HDBSCAN returns zero clusters with
100% of genes unassigned across **sixteen conditions** (four substrates, whole
panel and amount tertiles, two definitions of amount). Every dip test unimodal,
gap statistic still rising at k=8 everywhere. Structure exists and is
substantial (19 components above noise, 78% of variance) but it is not
partitioned.

**Result 2: categories are directions, not regions.** 17 of 21 external
reference sets are displaced from the panel centroid; **none is separated**.
Largest effect |d| = 0.73 leaves 71% overlap; most sit above 90%. Super-enhancer
genes do not separate at all (p = 0.15 to 0.18 depending on treatment, 94%
overlap) at a set size where the test has 80% power to detect d = 0.35.

Controlled for **eight** things: overall magnitude, gene density, gene length,
expression, CpG density, chromosome, cell-state matching, and statistical power.

---

## 2. Where this sits in the DPhil

The stated aims were: build MCCProfiler (1), learn representations (2), unify
and validate (3), expand genome-wide (4), and challenge the super-enhancer
concept (5).

**Aim 5 is now the strongest thread, and it was originally the least
developed.** The super-enhancer argument was framed as a thesis position; it is
now a measurement. That is a better position to be in than the one planned.

**Aim 1 is complete but reframed.** The tool exists, is validated, and is
justified differently from how it was intended to be.

**Aim 2 (learned representations) is untouched** and is now Phase C, with a
pre-registered bar: beat the amount-only baseline on the targets that survive
confounding, and reproduce across the twice-captured genes.

**Aim 3 (consensus clustering) is largely moot.** Consensus clustering assumes
clusters. Sixteen conditions say there are none. The GNN work should be
redirected toward the peak-level analysis, which is where the untried
information is.

**Aim 4 (genome-scale) gained two independent arguments this week**: n = 1,846
is thin for representation learning, and three reference sets cover over 80% of
this panel, which would not happen genome-wide.

**Transfer of Status** target is HT Year 2 (Jan to Mar 2027), 5,000 words. Two
geometric results that have survived heavy scrutiny is a defensible report.

---

## 3. Literature context

Pulled from the project literature notebook (`1b49ca44`) on 2026-08-16.

### The continuum position is well precedented, and that is reassuring

This project is not inventing a heterodox position; it is applying an accepted
one to a new data type.

- **Crowley et al. (2026)**, Tabula Sapiens: gene expression in most cell types
  is constrained to **low-dimensional polytopes** rather than discrete clusters,
  and the type-versus-state distinction is a **continuum of specialised
  functions**.
- **Cano-Gamez et al. (2020)**: an **"effectorness gradient"** in CD4+ T cells,
  the same cell type used here, with activation state better modelled as a
  continuous trajectory than as fixed subtypes.
- **Shoval et al. (2012), Hart et al. (2015)**: Pareto front and archetypal
  analysis, where vertices are specialists and the interior is a continuous
  spectrum of intermediate strategies. This is the right formal framing for what
  we see, and it is a better fit than clustering.
- **Biggin (2011)**: animal transcription networks are "highly connected,
  quantitative continua" rather than discrete on/off switches.
- **Pott and Lieb (2015)**: super-enhancers are the extreme tail of a continuous
  distribution, argued from thresholding logic.
- **Loubiere et al. (2024)**: a "spectrum of regulatory strategies", extending
  the binary housekeeping-versus-developmental framework into a continuum.

**Consequence for the writeup.** Frame the result as *applying* archetypal or
polytope reasoning to contact architecture, citing Shoval and Crowley, rather
than as an isolated negative. Pott and Lieb argued the super-enhancer case
logically; this measures it in an assay they did not have.

### The complication, and it is the important one

**Sexton et al. (2024)** used k-means on a "sum of chromatin state by contact"
matrix to define 18 Chromatin Interaction Signatures, and reported that
**super-enhancers are significantly enriched in specific enhancer-contacting
clusters**. That is a published positive claim in the same territory as our
null, and it must be engaged with rather than ignored.

**The likely reconciliation is our own framing.** Enrichment is not separation.
A set can be reliably displaced toward a region while overlapping almost
completely with the rest, which is precisely "displaced but not separated". Our
super-enhancer set is displaced by |d| 0.16 to 0.24 with 94% overlap. Sexton
would call that an enrichment; we call it a non-separation. **Both can be true,
and saying so is stronger than pretending the disagreement does not exist.**

Two further differences to state: they cluster genomic regions where we cluster
genes, and their input is chromatin-state-by-contact, so it is partly two assays
measuring the same underlying biology. That is the same "chromatin predicts
chromatin" concern that demoted three of our own positive controls.

**Other positive claims to engage with:**
- **Promoter-promoter assemblies** for housekeeping genes (the Ronin/THAP11
  work; note the first-author discrepancy flagged in `CONTEXT.md`). Our
  Eisenberg-HK set has no architectural signature once magnitude and density are
  controlled, which is why the mechanism-defined housekeeping test (S1) matters:
  it is the difference between "the biology is absent" and "the label is wrong".
- **Karpinska et al. (2025)**: chromatin hubs as a defining feature of
  tissue-specific loci. Our `cd4_specific_immune` set sits at |d| 0.13 with 95%
  overlap. Same tension, same likely resolution.

### Two published negative results that support the approach

- **Standard PCA on raw binned contact profiles fails**: PC1 explains under 2%
  of variance and captures depth rather than shape. That is the direct
  justification for position-invariant features over bin-level comparison.
- **Goel et al. (2023)**: loop and compartment callers failed to detect
  microcompartments in gene-rich regions. Relevant to the peak-level phase.

### Confounders the field controls for, against what we did

| field-standard control | our status |
|---|---|
| O/E distance-decay normalisation | done, the `oe_*` feature family |
| open-chromatin ligation bias (~40% more junctions at ATAC peaks) | **partially**: ATAC gates membership and is never a feature, but the bias itself is not quantified in our substrate |
| positional heterogeneity | done by construction, position-invariant features |
| sequencing depth and enzyme choice | done, magnitude basis |
| gene length, expression, CpG | **not standard in the field, and it broke two of our claims** |
| gene density, chromosome | **not standard, and both mattered** |

The last two rows are worth stating explicitly in the report. We controlled
things the field does not routinely control, and two of them overturned
findings. That is a methods contribution in its own right.

---

## 4. Two warnings about the literature notebook

**The notebook contains our own draft.** Queries return "Transfer Report (Draft)"
as a source, and it was cited back quoting "56% vs 15% variance" for the feature
substrate. **Do not treat notebook answers about this project as external
literature.** They are partly our own prior claims reflected back, and some of
those claims are now retracted. Check any number attributed to the transfer
report against `CURRENT_FINDINGS.md`.

**Hennig (2015) is not in the notebook.** It has been cited repeatedly in this
project's working notes as the methodological reference for "what are the true
clusters". That citation comes from `CLAUDE.md`, not from a source in the
collection. Either add the paper or stop citing it.

---

## 5. Next steps

Full ordering in `PLAN_FORWARD.md`. In brief:

1. **S1a, build the core-promoter sequence classes** (TATA / Inr / DPE by JASPAR
   scan on hg38). Substrate-independent, runs in parallel with the feature work.
   **This is the last test.** It is the second chromatin-independent positive
   control, which the argument is short of, and it decides whether the
   housekeeping null is biology or a bad label.
2. **S2.1 and S2.1b, the two feature removals.** Drop the four degenerate
   topology features; re-support or retire the per-peak moment families at a
   pre-set threshold of rho 0.70. Fixes, not tests.
3. **S2.3, rebuild and re-run** `audit/scripts/run_all.sh`. Every number moves.
4. **S1b**, the core-promoter displacement test on the final substrate.
5. Remaining corrections, then figures, then write.

**Then stop testing.** The claim has survived magnitude, density, length,
expression, CpG, chromosome, cell state and power. Further checking has negative
expected return.

**After transfer**, the peak-level unit of analysis (Phase A) is the only
research bet, and it is the one with a measured asset behind it: summits
reproduce at 14 bp, 64x tighter than a jitter null, while the current features
use nothing below roughly 1 kb.
