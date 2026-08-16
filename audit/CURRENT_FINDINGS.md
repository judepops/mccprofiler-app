# Current findings, as of 2026-08-14

**This is the document to write from.** Every number here has survived the
corrections and controls applied on 2026-08-14 and is the version to quote.
`REVIEW.md` is the provenance record: it shows how each number was arrived at
and, importantly, which earlier versions were wrong. Do not quote `REVIEW.md`
Parts 1 to 5 directly; several of its headline numbers were superseded the same
day, and Section 7 below lists exactly which.

Panel `gw_cd4_1`, human CD4+ T cells, 1,846 genes after QC and outlier removal,
91 features, store built 2026-08-03 at pipeline commit `96e13e4`.

---

## 1. The question

Is regulatory architecture, measured by base-pair-resolution contact profiles,
**categorical or continuous**?

Not "what are the archetypes". That framing presupposes the answer, and the
answer turns out to be no.

---

## 2. Aim 1: the feature substrate is worth building, narrowly

91 position-invariant features from viewpoint-anchored MCC profiles, nested and
cross-validated against two baselines.

**Against counting peaks:** wins on 9 of 9 targets, median gain 0.153.

**Against the 11-feature overall-magnitude basis** (`MAG_OVERALL`), which is the
fair comparison because it asks whether *shape* adds anything over *amount*:
wins on 9 of 9, **sign-test p = 0.002**.

But state the margins, because they are the honest part:

| target | 91 features | magnitude (11) | margin | margin / SD |
|---|---|---|---|---|
| gnomad_loeuf | 0.210 | 0.082 | +0.127 | **4.08** |
| GWAS_immune_hot | 0.746 | 0.695 | +0.051 | **3.47** |
| phastcons_2kb | 0.187 | 0.140 | +0.047 | 1.02 |
| blood TPM (log) | 0.380 | 0.359 | +0.021 | 0.94 |
| gtex_tau | 0.276 | 0.257 | +0.020 | 0.80 |
| cd4_rna_top_quartile | 0.707 | 0.697 | +0.009 | 0.44 |
| Lambert_TF | 0.661 | 0.647 | +0.014 | 0.36 |
| Eisenberg_HK | 0.637 | 0.628 | +0.009 | 0.24 |
| DepMap_curated_essential | 0.700 | 0.691 | +0.009 | 0.14 |

**Six of nine margins sit inside one standard deviation.** Only two clear 2 SD.
The defensible claim is therefore specific rather than sweeping:

> Contact architecture beyond overall magnitude predicts evolutionary constraint
> and immune-disease association, and little else.

Those two are also the only targets that survive amount correction in the
shape-corrected analysis, so two independent analyses converge on the same pair.

### The strongest single result in Aim 1: constraint

Do not let the honesty about the other seven flatten this one.

> **Evolutionary constraint is predicted almost entirely by contact shape, not by
> contact amount.** LOEUF: counting peaks gives r = 0.031, the 11-feature
> magnitude basis gives 0.082, the full 91 features give **0.210**, and the
> amount-corrected shape substrate alone still gives **0.191**. So 91% of the
> full-model performance survives removing magnitude entirely.

That is a 4.08 SD margin, by some distance the largest in the table, and it is
the one target where the architecture claim is doing real work rather than
riding on signal level. It is quotable as it stands.

---

## 3. Aim 2: the landscape is continuous, not categorical

Tested under **sixteen conditions**: four substrates (all 91 features; the 40
that reproduce at rho > 0.7; each with amount regressed out) across the whole
panel and within amount tertiles, under two different definitions of amount.

| test | result across all sixteen |
|---|---|
| HDBSCAN (min_cluster_size 25) | **0 clusters, 100% unassigned, every condition** |
| Hartigan dip test | **unimodal everywhere**, minimum p 0.42 |
| Gap statistic | **still rising at k=8 everywhere**, so k=1 never excluded |
| Silhouette | peaks at k=2, excess over permuted null +0.086 to +0.153 |

Cleaning the substrate did not help. Removing amount did not help. Conditioning
on amount did not help.

**Structure exists but is dimensional, not partitional.** 18 to 19 components
sit above a parallel-analysis noise ceiling, together 78% of variance, against 8
for permuted data. Effective rank 17.6 against 86.6 for noise. SigClust
z = -11.86 (trimmed).

**Caveat to state:** only **67%** of that retained structure rests on features
that reproduce across independent captures. Eight components, 14.9% of variance,
fall below 50% trusted and are all dominated by the `oe_asymmetry` and
`oe_tailedness` families, which reproduce at rho 0.27 to 0.42.

**The widest cut in the continuum**, if a binary description is wanted: k=2
splits 771 against 1075 on dispersed-and-long-range versus focal-and-local
(`signal_entropy` +1.38, `distal_signal_density` +1.12 against
`empty_band_fraction` -1.19, `frac_promoter_proximal` -1.12). Silhouette 0.096
under proper amount correction, so it is a cut through a continuum at its widest
point, not two clusters.

**RETRACTED, 2026-08-16: the k=2 cut is not a restatement of `arch-HK`.** An
earlier version of this section claimed it recovered 819 of `arch-HK`'s 844
genes. That number came from the superseded substrate (trusted features,
corrected against `total_mcc` alone). Recomputed on the correct
`MAG_OVERALL`-corrected substrate with the degenerate topology features dropped,
the cut captures **476 of 844 against 456 expected by chance, a 1.04x
enrichment**, i.e. essentially nothing. The claim is withdrawn.

What survives is only the weaker statement: the widest binary cut in the
continuum is dispersed-versus-focal, and it is amount-free. It is **not** the
same object as the imposed k=4 partition, and Sections 3 and 6 describe two
different things after all.

This was caught by drawing the figure rather than by a further review, which is
an argument for producing figures earlier than Phase 2a.

---

## 4. Aim 3: categories are directions, not regions

External reference sets, tested on a substrate corrected for overall magnitude
(11-feature `MAG_OVERALL` basis, basis features dropped, four degenerate
topology features dropped, 78 features, 18 components above a recomputed noise
ceiling), then checked against gene density.

**17 of 21 sets remain displaced** from the panel centroid after correction.
**None separates.** Largest effect d = 0.73 leaves **71% overlap**; most sets sit
above 90%.

### The headline claim

Lead with the matched-size comparison, because it forecloses the obvious
objection before it is raised.

> **Lambert transcription factors** (n = 156, defined by DNA-binding domain, so
> mechanism-defined rather than an expression list) separate in contact
> architecture at d = +0.46, p = 0.001, retaining 84% under gene-density
> stratification. They are displaced along an interpretable and well-measured
> axis: shape-PC1, 11.1% of corrected variance, 74% of it carried by features
> that reproduce across captures, loading **+promoter_signal_fraction against
> -CTCF signal, distance and consensus**. TFs are promoter-driven and CTCF-poor.
> **Super-enhancer genes** (n = 158, within two genes of the same set size, same
> substrate, same test) **do not separate**: d = -0.16, p = 0.18, 94% overlap.

**The super-enhancer list is state-matched.** It is derived by the primary
reference-set builder from dbSUPER `CD4p CD25- Il17- PMAstim Th` (867
super-enhancer regions, TSS +/- 50 kb), i.e. **stimulated** CD4, not naive. A
naive list (`CD4 Naive Primary 8pool`, 571 regions) exists in the same directory
but is used only by a comparison-only builder for a naive-versus-activated
sensitivity figure. So the null cannot be dismissed as a cell-state mismatch.
Running the naive list as a sensitivity check is still worth an hour
(`PLAN_FORWARD.md` Phase 0.5).

### Why super-enhancers are the case of interest

Not because their effect is smallest. It is not, and the earlier claim that it
was has been withdrawn twice. The reason is an asymmetry in what the literature
claims:

> **Nobody argues that Eisenberg housekeeping genes constitute a distinct
> three-dimensional architectural class. People do argue exactly that about
> super-enhancers.** A null for super-enhancers contradicts a live claim; a null
> for an expression-defined list contradicts nothing.

The methodological point sharpens it. ChromHMM and Roadmap fit a hidden Markov
model over many histone marks to assign a state, and both separate. ROSE ranks
stitched regions by a single mark and cuts at the inflection point of the ranked
curve, and that one does not. That is the precise procedure Pott and Lieb (2015)
argue has no biological warrant, and this is a direct measurement of the
consequence in a data type they did not have.

State it that way and the claim rests on the literature asymmetry and the
methodological contrast, not on an effect-size ranking that does not hold.

### Supporting observation: the kind of definition predicts displacement

Offered as context for the result above, not as the argument carrying it. An
earlier version of this section was a just-so story: a two-bin taxonomy drawn
over seven sets after seeing which ones separated. It has been redone properly
(`diagnose_definition_type.py`). All 21 sets are classified by the **data and
procedure their annotation was computed from**, which is a documented property
of the annotation, and the assignment is written into the script before any
effect size is consulted.

| category | rule | sets | median \|d\| |
|---|---|---|---|
| SEQUENCE | property of the gene product or DNA, no cell-type measurement | Lambert_TF 0.46, CpG_island 0.18, phastCons 0.16 | |
| GENETIC | association or constraint from population variation | GWAS_immune_hot 0.32, GWAS_total 0.24, GWAS_immune_any 0.20, pLI 0.17 | |
| STATE_MULTI | combinatorial state from an HMM over many marks | Roadmap_silenced 0.73, ChromHMM_bivalent 0.58, ChromHMM_active_TSS 0.29 | |
| **DIRECT (the three above)** | | **n = 10** | **0.27** |
| EXPRESSION | RNA abundance, specificity or stability | bio_dev_TF 0.24, DICE 0.23, bio_HK 0.16, cd4_rna 0.15, cd4_specific 0.14, bio_bulk 0.12, Eisenberg_HK 0.11 | |
| FITNESS | CRISPR dropout screen | DepMap_inferred 0.14, DepMap_curated 0.12 | |
| **OUTPUT (the two above)** | | **n = 9** | **0.14** |
| RANK_SINGLE | rank-order regions by ONE mark, cut at an inflection point | dbSUPER_CD4_SE 0.16 | n = 1 |
| CONTEXT | genomic neighbourhood, not the gene | gene_desert 0.71 | the density control |

**The pattern holds with the assignment fixed in advance.** Kruskal-Wallis
across the five testable categories H = 10.75, p = 0.030. DIRECT versus OUTPUT
by Mann-Whitney, one-sided, **U = 80, p = 0.0024**. Ten of ten DIRECT sets are
significant; six of nine OUTPUT sets are.

> Gene sets defined by what a gene or its regulation **is** (sequence, genetic
> consequence, multi-mark chromatin state) are displaced further in
> contact-architecture space than sets defined by what a gene **does**
> (expression, fitness). Median \|d\| 0.27 against 0.14, p = 0.0024.

**It is robust to how the sets are assigned.** The obvious objection is that
individual assignments are arguable and the result would move if they changed.
Tested four ways and it does not:

| check | p |
|---|---|
| baseline, DIRECT vs OUTPUT | 0.002 |
| leave-one-out across all 19 sets | worst case **0.003** |
| reassign each of the 4 most arguable sets in turn | 0.0006 to 0.0018 |
| reassign **all four simultaneously** | **0.0005** |
| 20,000 random 10/9 relabellings | 17 reach the observed p, i.e. **0.001** |

Moving the arguable sets *improves* the separation, which indicates the rule was
not drawn to maximise it. The four treated as arguable were
`CpG_island_promoter` (a promoter-class annotation as much as a sequence one),
`bio_dev_TF_k3` (about TFs, though the label is expression-derived),
`phastCons_2kb_topQ` (conservation is sequence but proxies function) and
`gnomAD_pLI_topQ` (population genetics, but reflects fitness).

**Three honest limits on this, all of which belong in the writeup.**

*It is not a preregistration.* The assignment is made from documented annotation
procedures rather than from the results, and the rule is written down so a
reader can check it was not drawn to fit. But the results had already been seen,
so this is a structured post-hoc test, not a blind one.

*It is a distribution shift, not a dichotomy.* Six of nine OUTPUT sets still
reach significance. Nothing here says output-defined sets have no architectural
signature; they have a weaker one.

*The super-enhancer category has n = 1 and cannot be tested.* RANK_SINGLE
contains one set. So no claim of the form "single-mark rank thresholds do not
correspond to architecture" is supported. Super-enhancers sit at the **38th
percentile** of \|d\| across all 21 sets: mid-pack, not last.

### The axis rank is a stronger statement of the null than the effect size

Effect size says how far a set is displaced. Axis rank says *where* it leans,
and for the super-enhancer null that is the more legible statement.

A set with real architectural structure should find its largest displacement in
a high-variance, well-measured component. A set with no architectural signature
has no preferred direction, so its largest lean lands wherever noise happens to
be biggest: a low-variance, poorly-reproducing component.

Measured across all 21 sets (`diagnose_external_axis_rank.py`):

| | sets clearing p < 0.05 (17) | sets that do not (4) |
|---|---|---|
| median axis rank | **3** | **9** |
| median trust of that axis | **74%** | **45%** |

Eleven of the seventeen significant sets displace along one of the top four
components. Spearman between -log10 p and axis rank is **-0.53**: the more
displaced a set is, the higher-variance the component it displaces along.

**Super-enhancers lean on sPC12** (2.7% of variance, 43% trusted). Only 3 of 21
sets lean on a lower-ranked component. So the finding is not "weakly displaced".
It is **no preferred direction anywhere in the well-measured part of the space**,
with the strongest lean landing exactly where noise would put it.

### Three independent positive controls, on three different well-measured axes

| set | defined by | n | axis | var | trust | d | p |
|---|---|---|---|---|---|---|---|
| Lambert_TF | DNA-binding domain (mechanism) | 156 | sPC1 | 11.1% | 74% | +0.46 | 0.0005 |
| GWAS_immune_hot | disease variants | 130 | sPC2 | 10.1% | 76% | +0.32 | 0.0005 |
| Roadmap_silenced | chromatin state | 41 | sPC3 | 7.8% | 74% | -0.73 | 0.0005 |
| **dbSUPER SE** | **H3K27ac signal threshold** | **158** | **sPC12** | **2.7%** | **43%** | **-0.16** | **0.18** |

Three different kinds of definition, three different top-4 components, all
well-measured, all significant. The super-enhancer set is defined by a threshold
on signal rather than by a mechanism, a variant, or a chromatin state, and it is
the one that finds no direction.

**One caveat on Lambert_TF.** It shares sPC1 with `gene_desert_bottomQ_density`,
in the opposite direction: TFs sit in gene-dense regions, deserts do not. Since
Lambert_TF retains 84% of its effect under gene-density stratification, the TF
signal is not merely density, but the shared axis should be stated rather than
discovered by a reader.

### Which positive control to qualify

The controls do not separate on equally trustworthy axes, and this matters more
than their effect sizes:

| set | separates on | variance | loading mass on reproducible features |
|---|---|---|---|
| Lambert_TF | shape-PC1 | 11.1% | **74%** |
| Roadmap_silenced | shape-PC3 | 7.8% | **74%** |
| ChromHMM_bivalent | shape-PC4 | 5.5% | **15%** |
| dbSUPER SE | shape-PC12 | 2.7% | 43% |

**Qualify ChromHMM_bivalent.** Its d = -0.58 sits on the `oe_asymmetry`
component, the feature family that reproduces at rho 0.27 to 0.42 (Section 3
caveat). It is a real displacement of a set, but it rests on the least
reproducible thing we measure and should not carry the argument alone.

Worth noting for completeness: the super-enhancer set's own non-significant
tendency lies on shape-PC12, a 2.7% component that is itself only 43% trusted.
Even the direction it leans is uninterpretable.

**The null is informative, not merely a failure to look.** Simulation on this
substrate gives the super-enhancer test **80% power to detect a displacement of
d = 0.35** at n = 158. The observed effect is 0.16. The positive controls show
0.46 to 0.73, all above that threshold. So an effect of the size seen elsewhere
in the table would have been found had it existed.

Note this **inverts** the natural worry that the positive controls are small
(n = 41, 46) while the null is large (n = 158). The super-enhancer test is the
**best** powered in the table; the chromatin-state sets need d > 0.50 for 80%
power and clear it only because their true effects are large.

That is a null measured against a mechanism-defined positive control of matched
size, on the same test and substrate, controlled for the two confounders that
broke every earlier version of this analysis: overall magnitude and genomic
context. Pott and Lieb (2015) argued the point from thresholding logic; this
measures it.

Lambert_TF is the positive control to foreground for a second reason: it is
mechanism-defined, so it establishes in advance that mechanism-defined sets can
separate in this space where expression-defined ones (Eisenberg-HK) do not.
That is the same distinction Phase 1.1 sets out to test directly.

### The full table, corrected substrate

| set | n | z | p | d | overlap |
|---|---|---|---|---|---|
| gene_desert_bottomQ_density | 456 | 20.8 | 0.0005 | -0.71 | 72% |
| ChromHMM_bivalent | 46 | 7.1 | 0.0005 | -0.58 | 77% |
| Lambert_TF | 156 | 5.4 | 0.0005 | +0.46 | 82% |
| Roadmap_silenced | 41 | 3.9 | 0.0005 | -0.73 | 71% |
| GWAS_immune_hot | 130 | 3.7 | 0.0005 | +0.32 | 87% |
| GWAS_immune_any | 733 | 3.6 | 0.0010 | +0.20 | 92% |
| Eisenberg_HK | 633 | 2.0 | 0.034 | +0.11 | 96% |
| DepMap_inferred_essential | 328 | 1.0 | 0.171 | -0.14 | 94% |
| dbSUPER_CD4_SE | 158 | 0.9 | 0.177 | -0.16 | 94% |
| DepMap_curated_essential | 254 | 0.5 | 0.295 | +0.12 | 95% |

### Three things that must travel with this table

**`gene_desert` is the density positive control, not a finding.** It is the
largest displacement in the table and it survives amount correction, so genomic
context is a real and separate confound. Its effect is by construction.

**Essentiality was amount.** Both DepMap sets fall to non-significance after
correction, having looked like real signal at d = 0.46 to 0.52 on raw
components. Nothing should be claimed about essentiality and architecture.

**Eisenberg housekeeping has no architectural signature** once magnitude and
density are both controlled (p = 0.034 -> 0.076). Whether that is because the
label is an expression list rather than an architecture list is **untested**: the
mechanism-based check could not run, only 19 ribosomal protein genes clear the
25-gene panel floor. This is the single most important open question.

---

## 5. Supporting results

**Reproducibility.** 116 genes captured in both panels, 63 shared features,
median rho **0.752**, 40 of 63 above 0.7, 1 below 0.3. Most reproducible are the
distance compositions (`frac_local` 0.984, `frac_far_distal` 0.977). Least are
the per-peak moments (`oe_asymmetry_mean_all` 0.274).
*(An external count of 119 genes exists and is unreconciled; the store has 116.)*

**Archetype label stability.** Cohen's kappa 0.72 at k=3-4 across independent
captures. State this as *the imposed partition is stable*, which licenses using
it as a descriptive device. It is **not** evidence that groups exist.

**Resolution.** Contact summits reproduce at **median 14.0 bp** across
independent captures against a null median of 895.3 bp, **64x tighter**; 76.8%
within 50 bp against 3.3% for the null. Sub-resolution collapse: at 5 kb
binning, 18.7% of peaks merge away and 90% of genes lose at least one; at 25 kb,
45.8% and 97%.

**The representation does not exploit that resolution.** The 91 features use
nothing below roughly 1 kb. The assay resolves contacts two orders of magnitude
finer than the features read them. This is the measured basis for the
peak-level research programme.

---

## 6. What the architecture groups are called, and the caveat

| canonical | display | top discriminating features |
|---|---|---|
| `arch-HK` | dispersed | signal_entropy +1.25, empty_band_fraction -1.17, frac_far_distal +1.09 |
| `arch-ME-constitutive` | promoter-local | promoter_signal_fraction_raw +1.13, n_peaks_promoter +0.93 |
| `arch-ME-effector` | enhancer-focal | raw_peak_max_max_enhancer +1.03, signal_entropy -0.97 |
| `arch-sparse` | sparse | n_high_consensus_peaks_075 -1.62, n_peaks_all -1.57 |
| `arch-off` | empty (QC) | all CTCF descriptors at floor |

Sizes 844 / 369 / 351 / 261 / 21. Core (posterior >= 0.8): 1,040 of 1,846; mean
max posterior 0.793; 44% of active genes are mixtures.

**`arch-HK` is not the housekeeping group**: Eisenberg-HK fraction 38.2% against
`arch-ME-constitutive`'s 40.2%, and it has the lowest median blood expression of
the three active groups (4.5 TPM vs 11.8 and 9.9).

**These names are provisional.** They were derived from *unadjusted* group means,
and amount is now known to dominate the geometry. They must be re-derived on the
amount-corrected substrate before use (Phase 0.2 in `PLAN_FORWARD.md`).

---

## 7. Retracted, do not quote

Listed explicitly because all of these appear in `REVIEW.md` Parts 1 to 5 and in
earlier handoffs.

| claim | status |
|---|---|
| "PC1 is not amount" (r = 0.033 with `total_mcc`) | **Retracted.** Against the `MAG_OVERALL` basis, PC1 correlates at **0.623**. PC1 is the amount axis. The error was using one feature as the amount proxy. |
| "Amount is not one quantity" | **Reframed.** Amount is multi-faceted, so it must be measured with the 11-feature basis, not with `total_mcc`. |
| "Regress `total_mcc` out before PCA" | **Superseded.** Use `MAG_OVERALL` via `_shape.corrected_shape`. |
| "Mechanism vs output taxonomy over 7 sets" | **Superseded.** Redone across all 21 with the assignment fixed by annotation procedure; the pattern holds at p = 0.0024, but RANK_SINGLE is n = 1 so nothing is claimable about that category. |
| "Super-enhancers are the weakest of 21" | **Withdrawn twice.** First as over-precise, then outright: on the corrected substrate the SE set has a LARGER effect (\|d\| 0.16) than Eisenberg-HK (0.11). Replaced by the kind-of-definition pattern in Section 4. |
| "The asymmetry family is a strand/orientation bug" | **Refuted.** `|value|` reproduces worse, and strand-relative asymmetry is d = -0.045. It is a support-size problem: median peak ~11 bins at 1-2 reads per bin. |
| "Most of the signal is amount, not shape" | **Softened.** Correct statement is *amount is sufficient for most targets*; residualisation can strip real architecture if amount is downstream of it. |
| Varimax rotation should be adopted | **Withdrawn.** It nearly doubles nameability for free, but its most concentrated factors are the `oe_asymmetry` families at rho 0.27-0.42. Do not adopt without weighting by reproducibility. |
| Displacement figures from the raw components (largest d = 1.61, 21 of 21) | **Superseded** by the corrected substrate in Section 4 (largest d = 0.73, 17 of 21). |
| The app's `/api/enrichment/grid` numbers | **Uncorrected.** Still computed on raw components. Phase 0.3. |
| "The k=2 cut recovers arch-HK almost exactly, 819 of 844" | **Retracted 2026-08-16.** On the correct `MAG_OVERALL` substrate it is 476 of 844 against 456 expected, 1.04x. The 819 figure came from the superseded `total_mcc` correction. |

---

## 8. Where to go next

See `PLAN_FORWARD.md`. In short: four corrections that change the report
(Phase 0), the mechanism-versus-annotation housekeeping test (Phase 1.1, now the
highest-value open question), figures as their own scoped phase (2a), then
writing. The research programme afterwards is led by the peak-level unit of
analysis, because the 14 bp summit reproducibility is a measured asset that the
current representation demonstrably does not use.

---

## Reproducing every number here

| section | script |
|---|---|
| 2 | `audit/continuous_methods/nested_baselines.tsv` |
| 3 | `backend/scripts/experiment_cluster_search.py` |
| 3 caveat | `backend/scripts/diagnose_dimension_trust.py` |
| 4 | `backend/scripts/diagnose_external_structure_corrected.py` |
| 4 density | `backend/scripts/diagnose_external_density_stratified.py` |
| 4 power | `backend/scripts/diagnose_external_power.py` |
| 4 axis rank | `backend/scripts/diagnose_external_axis_rank.py` |
| 5 resolution | `audit/continuous_methods/summit_precision.tsv`, `subresolution_collapse.tsv` |
| 7 (PC naming) | `backend/scripts/diagnose_pc_names.py`, `experiment_rotate_axes.py` |
| store invariants | `backend/scripts/verify_store.py`, 38 checks |
