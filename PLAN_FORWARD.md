# Forward plan, 2026-08-14

Written after the Aim 3 re-test (REVIEW.md Addendum 3). Transfer of Status
target is HT Year 2, roughly Jan to Mar 2027, 5,000 words.

## Where things actually stand

Three aims, and after today all three have an answer rather than a hope.

| aim | status |
|---|---|
| 1. MCCProfiler is worth building | **Yes, but narrowly.** 91 features beat the 11-feature magnitude basis on 9/9 targets, sign-test p = 0.002. But 6 of 9 margins sit inside one SD. Only LOEUF (4.1 SD) and immune GWAS (3.5 SD) are substantial, and they are the same two that survive amount correction. |
| 2. The landscape is a continuum | **Established.** HDBSCAN 0 clusters in 16 conditions across two amount definitions and four substrates; dip unimodal throughout; gap rising at k=8. 18-19 components above noise. |
| 3. External categories are displaced but not separated | **Survives, and now controlled twice.** 17 of 21 displaced after magnitude correction, largest d = 0.73 leaving 71% overlap. Positive controls retain 98-112% under gene-density stratification. Super-enhancers do not separate (p = 0.16 with density controlled). |

The claim that carries the thesis:

> On a substrate corrected for overall magnitude, and retained under gene-density
> stratification, chromatin-state-defined categories separate in contact
> architecture (ChromHMM bivalent d = -0.58, 112% retained; Roadmap silenced
> d = -0.73, 98% retained; both p = 0.0005). Super-enhancer genes do not separate
> (p = 0.16, 94% overlap), despite being defined as a distinct class of
> regulatory element.

That is a null against a positive control on the same test and substrate,
controlled for the two confounders that have broken every earlier version of
this analysis: overall magnitude and genomic context. Pott and Lieb argued it
from thresholding logic; this measures it.

---

## Phase 0. Corrections that change what goes in the report

Nothing should be written until these land. None is a research question.

**0.1 Drop the four degenerate topology features upstream.**
`mean_degree`, `mean_degree_raw`, `frac_active_pairs`, `n_isolates_raw`.
Confirmed algebraic at 1842/1842 genes. Re-derive PCA and dimension names.
PC1 becomes nameable at 18.41% variance: `signal_entropy` +0.79,
`oe_distal_mean` +0.75, `distal_signal_density` +0.72, `empty_band_fraction`
-0.69, i.e. **dispersed distal contact versus emptiness**.
*Half a day. Unblocks every figure.*

**0.2 Re-derive the archetype display names on the amount-corrected substrate.**
Current names come from unadjusted group means, and amount is now known to be
the dominant axis. `arch-HK` is already known not to be housekeeping; whether
"dispersed" survives correction is untested.
*Half a day. Blocks any naming in the report.*

**0.3 Fix the app's enrichment grid.** `/api/enrichment/grid` still computes on
raw components, so every displacement it shows is the uncorrected, overstated
one. Switch to the corrected substrate.
*Two hours.*

**0.4 Reconcile 116 vs 119 twice-captured genes.** The store has 116, the review
feedback has 119. The median rho 0.752 is quoted widely and needs one number.
*One hour.*

---

## Phase 1. Two open questions that could change the report's claims

**1.1 Is Eisenberg-HK failing because the label is wrong?**
It sits at d = 0.11, 96% overlap, the weakest meaningful set. The objection that
housekeeping genes form promoter assemblies and therefore *should* separate is
untested, because the mechanism-based check could not run: only 19 ribosomal
protein genes are in the panel, below the 25-gene floor.

Needs a larger mechanism-defined set. Options: Hwang 2023 promoter-assembly
genes, a broader ribosomal and translation-machinery set, or CpG-island-plus-
broad-promoter genes. If a mechanism-defined set separates where the expression
-defined list does not, that is a *positive* result and it strengthens Aim 1
considerably: it would show the features track mechanism rather than annotation.
*One day. **Highest-value item in Phase 1**, promoted above 1.2. It does two
jobs: it could add a positive finding, and it pre-empts the objection any
chromatin biologist will raise, that housekeeping genes form promoter assemblies
so they ought to separate. Addendum 4 makes it more urgent still: Eisenberg-HK
now fails entirely (p = 0.076) once density is controlled, so a mechanism-defined
set is the only remaining route to saying anything about housekeeping
architecture.*

**1.2 Re-support the per-peak moment families.**
Eight components rest on `oe_asymmetry` and `oe_tailedness` at rho 0.27 to 0.42.
Diagnosed as a support-size problem: median peak ~11 bins, 42% under 10 bins,
1-2 reads per non-zero bin. `contact_asymmetry`, the same statistic on the whole
profile, reproduces at 0.900. So compute the moments on aggregated support
(stacked peaks per gene, or per element class) rather than per peak.
*Two days. Deliberately second: it recovers variance in components that are not
load-bearing for any current claim.*

---

## Phase 2. Write the transfer report

Structure, following the framing that the work actually supports:

1. **Question.** Is regulatory architecture categorical? Not "we will find
   archetypes."
2. **Aim 1, the tool.** MCCProfiler, 91 features from bp-resolution MCC.
   Nested cross-validated baselines against `n_peaks` and against the
   11-feature magnitude basis. State the margins honestly: consistent wins,
   substantial only for constraint and immune GWAS.
3. **Aim 2, the answer is no.** Sixteen adversarial conditions. Structure is
   dimensional, not partitional. 18-19 components above a parallel-analysis
   ceiling. Include the caveat that only ~67% of that rests on reproducible
   features.
4. **Aim 3, what it is instead.** Categories are directions, not regions.
   Displaced but not separated. Lead with the super-enhancer null against the
   chromatin-state positive controls.
5. **Methods paragraph on the self-correction protocol.** Five overturned
   claims, each caught by a control that was built before the result was
   trusted. This is an asset; say so plainly and briefly.
6. **Forward programme.** Phase 3 below.

Two things to state explicitly rather than hide:
- We cannot exclude that a better representation would partition. We can exclude
  that these 91 features do, under sixteen conditions.
- `gene_desert` remains the largest displacement even after amount correction,
  so genomic-context confounding is real and separate from amount. The positive
  controls were separately shown to survive density stratification
  (Addendum 4).
- **That the analysis was redone.** The first pass ran on uncorrected components
  and overstated every displacement; the reported numbers come from a substrate
  corrected for magnitude, then checked against gene density. Saying so is
  stronger than presenting only final numbers, it is the same logic as the
  self-correction methods paragraph, and an examiner who finds the earlier
  version in a handoff will respect it.

*Three to four weeks of writing, after Phase 0, 1 and 2a.*

### Phase 2a. Figure production, scoped separately

Figures are not a by-product of writing and should not be assumed to fall out of
it. Every one needs the corrected substrate from Phase 0, so none can start
earlier, and seven composite figures is substantial work that will otherwise
compress against the deadline.

| # | figure | source |
|---|---|---|
| 1 | The assay and the pipeline | schematic |
| 2 | Feature substrate and nested baselines | `nested_baselines.tsv` |
| 3 | The continuum: gap, HDBSCAN, dip, silhouette vs null across 16 conditions | `experiment_cluster_search.py` |
| 4 | The dimensions: scree, noise ceiling, loadings, trust weighting | `diagnose_dimension_trust.py` |
| 5 | External sets displaced but not separated, corrected substrate | `diagnose_external_structure_corrected.py` |
| 6 | The super-enhancer null against its positive controls, density-stratified | `diagnose_external_density_stratified.py` |
| 7 | Resolution: 14 bp summits, sub-resolution collapse | `summit_precision.tsv` |

*Two weeks, in parallel with early writing but after Phase 0.*

---

## Phase 3. The research programme after transfer

Ranked by expected value, with the measured justification for each.

**3.1 Peak-level unit of analysis. The main bet.**
Contact summits reproduce at **14 bp median across independent captures, 64x
tighter than a jitter null**, 76.8% within 50 bp against 3.3% for the null. And
the current features use nothing below roughly 1 kb: at 5 kb binning, 18.7% of
peaks merge away and 90% of genes lose at least one.

So the assay resolves contacts two orders of magnitude finer than the
representation exploits. Changing the unit from gene to peak is the only route
in this project with a measured asset behind it rather than an argument. A
per-peak GNN with genes as subgraphs is the natural form, but it requires a real
adjacency first (0.1 removes the degenerate one; co-contact above an O/E
threshold is the obvious replacement).

**3.2 Learned representations on the raw profile.**
40,000 bins x 8 channels compressed to 91 position-invariant scalars. If
partitional structure exists and the summaries cannot hold it, this is the way
to find out. The bar is now explicit and should be pre-registered:
- beat the **amount-only** baseline on LOEUF and immune GWAS, not the `n_peaks`
  baseline
- reproduce across the twice-captured genes, held out

**3.3 Contrastive objective using technical replicates.**
The twice-captured genes are the same biological object measured twice. Training
a representation to place those pairs together learns invariance to capture
noise directly, which is the measured weakness (8 of 19 components on features
at rho < 0.5). Small n, so a regulariser on top of reconstruction, not a
standalone objective.

**3.4 Genome-scale expansion.** The 20k probe panel remains the route to n large
enough for representation learning to be well posed.

---

## Publication expectation, stated plainly

Specialist tier is realistic now: Genome Biology, Genome Research, NAR. A
validated tool, a rigorous negative, and a quantitative challenge to the
super-enhancer concept is a solid first paper.

High-impact venues want a mechanistic positive, which does not exist yet. The
candidate is 3.1: where contacts sit relative to motifs at 14 bp resolution, and
what that predicts. That is a genuine possibility rather than a consolation,
because the resolution asset is measured and unexploited.

For the transfer itself none of this matters. The transfer needs rigorous
science and a defensible programme, and both are in hand.

---

## Immediate next action

Phase 0.1. Drop the four degenerate features upstream in `mccprofiler`,
re-derive, and confirm PC1 lands where the app's experiment predicted
(18.41% variance, entropy-led). Everything else queues behind it.
