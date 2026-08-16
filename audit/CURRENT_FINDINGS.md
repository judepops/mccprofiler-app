# Current findings

> **⚠ SUPERSEDED IN PART, 2026-08-16 (evening). The substrate changed.**
>
> The feature set went from 91 to 73 (four degenerate topology features and
> sixteen unreproducible per-peak moments removed, `n_active_peaks` protected)
> and the store was rebuilt. **Every component count, variance share and effect
> size below is from the old 91-feature substrate.** Section 0 states what
> changed and what survived. Numbers not yet re-checked against
> `audit/scripts/output_S2/` should not be quoted.

---

## 0b. Feature health audit, all 73 (`audit_feature_health.py`)

Four feature problems were found this week, each by accident and each a
different failure mode. This runs all four checks over every feature at once so
the next one is found deliberately. **It found more.**

### What was removed, and why it was right

| removed | n | why |
|---|---|---|
| `mean_degree`, `mean_degree_raw`, `frac_active_pairs`, `n_isolates_raw` | 4 | exact algebraic functions of two peak counts (1842/1842 genes), from a complete-graph adjacency with no contact criterion |
| `oe_asymmetry_*`, `oe_tailedness_*` | 16 | 3rd and 4th moments of ~11-bin peaks at 1-2 reads per bin; 14 of 15 measured below the pre-set rho 0.70 |
| **restored:** `n_active_peaks` | +1 | the interpretable primitive the four derivatives were functions of, wrongly pruned while they survived |

**Why it was positive, measured rather than asserted.** Retained structure
resting on reproducible features rose from **67% to 86%**; components below 50%
trusted went from **8 of 19 to 0 of 14**; PC1 rose 15.27% to 17.23% and its top
loading is now `n_active_peaks`. Cost: 0.3 points of variance. And the results
did not become unstable, they became sharper: effect sizes across all 21
external sets correlate **r = 0.959** before against after, with 16 of 21
increasing.

### What the audit found that we had not looked for

**Three feature families are compositional and sum to 1**, so one member of each
is redundant by construction, exactly as `mean_degree` was:

| family | members | row sum |
|---|---|---|
| distance bands: `frac_promoter_proximal`, `frac_local`, `frac_distal`, `frac_far_distal` | 4 | **1.0000** |
| element classes (O/E): `promoter_`, `enhancer_`, `ctcf_signal_fraction` | 3 | **1.0000** |
| element classes (raw): the same three, `_raw` | 3 | **1.0000** |

All four band fractions predict at **R2 = 1.000** from the others. Fifteen
features exceed R2 0.95 and thirty exceed 0.90. The correlation prune cannot
catch this: compositional redundancy is a linear dependency among three or four
variables, and a pairwise |r| filter only sees pairs.

**`total_mcc` is the weakest member of the magnitude basis.** R2 = **0.965**
predictable from the other 72 features, and cross-panel rho **0.65**, the
fourth-lowest in the substrate. It is the feature `MAG_OVERALL` leads with and
the one used as the amount proxy in the analysis that was already retracted for
being the wrong basis. The basis as a whole is fine; the habit of reaching for
`total_mcc` alone is not.

**`n_peaks_promoter` is 41% gene properties.** The only feature above R2 0.25
against length, expression, CpG, density and chr19, and it also predicts at 0.93
from the other features. It fails two checks.

**Eight of 45 measured features are below rho 0.70**, led by
`distance_to_nearest_peak_bp` at **0.32**, `oe_local_enrichment_max_all` 0.44,
`spacing_regularity` 0.55, `corr_raw_distance` 0.57.

**28 of 73 features have never been measured twice.** The reproducibility panel
covers 45. Nothing is known about the reliability of the other 28, and that is a
coverage gap rather than a clean bill of health.

### Should we be more careful: yes, and here is the specific answer

The four problems found this week were not bad luck. They are what happens when
features are added individually and audited never. Three rules follow, and the
script now enforces the first:

1. **Every new feature is checked for predictability from the existing set
   before it is adopted.** R2 > 0.95 means it is a derivative; R2 > 0.90 means
   it earns almost nothing.
2. **Compositional families must declare which member is dropped.** Fractions
   over an exhaustive partition always contain a redundancy, and no pairwise
   filter will find it.
3. **Reproducibility is a property of a feature, not of the substrate.** A
   feature with no cross-capture measurement should be labelled unknown rather
   than assumed sound; 28 of 73 currently are.

**What NOT to do.** Removing the redundant compositional members is not urgent
and may not be worth it. Redundancy is a mild inefficiency, whereas the four
removed features were degenerate or unmeasurable, which is a different problem.
The audit exists so the distinction is made deliberately rather than by whoever
happens to notice something in a loadings plot.

---

## 0. What the rebuild changed, 2026-08-16

Feature set 91 to 73. Store rebuilt, `verify_store.py` re-baselined and passing
38 of 38. Full re-run in `audit/scripts/output_S2/`.

### The substrate got substantially cleaner

| | before (91) | after (73) |
|---|---|---|
| components above the noise ceiling | 19 | **14** |
| variance they carry | 78.0% | 77.7% |
| PC1 variance | 15.27% | **17.23%** |
| **retained structure on reproducible features** | **67%** | **86%** |
| components below 50% trusted | 8 of 19 | **0 of 14** |

One percentage point of variance bought a 19-point rise in trust and removed
every weak component. PC1's top loading is now `n_active_peaks`, the primitive
protected from the prune, at +0.223.

### Aim 1 holds, essentially unchanged

Gene length, expression and CpG explain **2.9%** of the feature variance (was
2.6%) and **0.8%** of the retained components (was 0.6%). Zero of 18 components
exceed R2 0.10. The orthogonality claim is unaffected.

### Aim 2 holds

HDBSCAN returns **0 clusters, 100% unassigned**, every dip test unimodal,
silhouette peaking at k=2 at 0.123 to 0.159 against a null of 0.040. The
continuum result does not depend on the removed features.

### Aim 3 changed, and the super-enhancer claim must be restated

**The cleaner substrate detects more, not less.** Nineteen of 21 sets now clear
p < 0.05, where 17 did before. And the super-enhancer set is now among them:

| | before (91) | after (73) |
|---|---|---|
| dbSUPER SE | z 0.9, **p = 0.177**, d = -0.16, sPC12 | z 2.3, **p = 0.0145**, d = -0.25, sPC6 |

**"Super-enhancers do not separate at all" is withdrawn.** They are weakly
displaced, like almost every other set. What survives, and is now a cleaner
statement than the null was:

> Every external category tested is displaced in contact-architecture space and
> **none is separated**. Super-enhancers are displaced no more than an average
> set (d = 0.25, 90% overlap, 5 of 21 sets lean on a lower-ranked component) and
> far less than chromatin-state categories (ChromHMM bivalent d = -0.76 at 70%
> overlap, Roadmap silenced -0.65 at 74%), despite being defined as a distinct
> class of regulatory element.

The power argument still stands: at n = 158 the test detects d = 0.35 with 80%
power, and the observed 0.25 is below that, so this remains a weak effect rather
than a strong one measured precisely.

Other Aim 3 results on the new substrate:

- **Positive controls strengthen.** ChromHMM_bivalent -0.76 (was -0.58) and now
  on sPC3 at 87% trust rather than the old asymmetry-dominated component at 15%.
  That control is no longer compromised.
- **Definition-type taxonomy holds**: DIRECT median 0.29 against OUTPUT 0.16,
  Mann-Whitney p = 0.0024, Kruskal-Wallis p = 0.036. Unchanged.
- **Axis rank weakens as an argument**: SE now leans on sPC6 rather than sPC12,
  and the rank correlation falls from -0.53 to -0.30. Still directionally right,
  no longer a headline.
- `gene_desert` remains the largest at 0.74, still the density control.

### What must be redone

Figures (all six, drawn on the old substrate), the density and chromosome
re-tests quoted in Sections 4 and 6a, and every number in Sections 2 to 6 that
has not been checked against `output_S2/`.

---


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

### Aim 1, restated: the features measure something gene properties do not

The predictive justification is retracted (below). This replaces it, and it is a
better claim because it asks what the features **contain** rather than what they
predict. A readout can fail to beat simple covariates at predicting known labels
while still measuring a quantity nobody has measured. Those are different claims
and only the second is what a new assay readout is for.

`diagnose_feature_novelty.py`, cross-validated:

| | variance explained |
|---|---|
| gene length + expression + CpG, across the 91 features | **2.6%** |
| the same, across the 18 retained components | **0.6%** |
| for scale: the assay's own 11-feature magnitude basis | 27.3% |

Zero of 91 features exceed R2 0.25 and only three exceed 0.10, the largest being
`n_peaks_promoter` at 0.181. **Zero of 18 components exceed 0.10.** Several
features are at 0.000: `contact_asymmetry`, `max_mcc`, `oe_distance_kurtosis`,
`distance_to_nearest_peak_bp`.

> **Gene length, expression and CpG density explain 2.6% of the variance in the
> 91-feature space and 0.6% of its retained components. The substrate is not a
> re-description of gene properties.**

**The apparent tension with the retraction is the finding, not a problem.** The
features contain 97% variance that gene properties do not explain, yet they add
nothing to those properties for predicting existing annotations. Both are true,
and together they say something specific:

> The assay measures reproducible structure that existing gene annotations do
> not capture.

The reproducibility clause is what stops this being an elaborate way of saying
"noise". The unexplained variance reproduces: median cross-panel rho **0.752**
over 63 shared features, **19 components above a parallel-analysis noise
ceiling** carrying 78% of variance against 8 for permuted data, and Cohen's
kappa **0.72** for the imposed partition across independent captures. Noise does
not do that.

And it motivates Aims 2 and 3 directly rather than sitting beside them. If the
structure were categorical and matched existing categories, the annotations
would predict it. They do not, because the structure is continuous and the
categories are directions rather than regions.

### RETRACTED 2026-08-16: the predictive claim does not survive gene length

The nested baselines compare the 91 features against peak counting and against
the 11-feature magnitude basis. They never compared against **gene length,
expression level and CpG density**, and those three simple covariates predict
better than the features on every target tested
(`diagnose_loeuf_confounds.py`).

| target | confounders | shape | both | shape adds |
|---|---|---|---|---|
| gnomad_loeuf | **0.452** | 0.152 | 0.454 | **+0.002** |
| GWAS_immune_hot | 0.694 | 0.627 | 0.699 | +0.005 |
| DepMap_curated_essential | 0.576 | 0.520 | 0.582 | +0.006 |
| Eisenberg_HK | 0.680 | 0.533 | 0.673 | **-0.007** |
| gtex_tau | 0.218 | 0.080 | 0.230 | +0.012 |
| phastcons_2kb | 0.055 | 0.065 | 0.098 | +0.043 |
| Lambert_TF | 0.577 | 0.622 | 0.631 | +0.053 |

**Contact shape adds essentially nothing to three simple gene properties.** The
confound is specifically length: partialling length alone drops the LOEUF result
from 0.152 to 0.085, while expression and CpG cost nothing.

Why length. LOEUF is a depletion statistic, so it scales with coding sequence
length by construction, and length independently changes contact features
because a longer gene occupies more of the plus-or-minus 1 Mb window
(`r(length, n_peaks_promoter)` = 0.42, though the median across features is only
0.08).

**What this costs.** Aim 1's justification cannot be "the features predict
biology better than the baselines". They beat peak counting and they beat the
magnitude basis, both true and both now insufficient, because neither baseline
contained length. The honest version:

> The 91 features improve on peak counting and on overall magnitude, but their
> predictive advantage over gene length, expression and CpG density is
> negligible. Their value is not prediction of existing annotations.

**What survives.** Aims 2 and 3 are untouched, because they concern the geometry
of the feature space rather than prediction of external labels, and Aim 3 was
re-tested against length directly (below). And the tool's justification is not
lost, only relocated: see the section above, where gene properties explain 2.6%
of the feature space. It is justified by what it measures, not by what it
predicts.

### Aim 3 survives gene length: six confounders now controlled

| set | magnitude only | + length | + length, expression, CpG, density |
|---|---|---|---|
| ChromHMM_bivalent | 0.58 | 0.59 | 0.66, p 0.0005 |
| Roadmap_silenced | 0.73 | 0.74 | 0.47, p 0.0025 |
| Lambert_TF | 0.46 | 0.49 | 0.39, p 0.0005 |
| GWAS_immune_hot | 0.32 | 0.24 | 0.31, p 0.040 |
| **dbSUPER SE** | 0.16 | 0.19 | **0.24, p 0.15** |
| gene_desert (density control) | 0.71 | 0.57 | 0.20, p 0.0005 |

The super-enhancer null holds under all of them. `gene_desert` collapsing is
expected and is the check working: density is in the confounder set.

The claim is now controlled for **overall magnitude, gene density, gene length,
expression, CpG density, chromosome, cell-state matching, and statistical
power.**

### The former Aim 1 headline, retained for context only

Do not let the honesty about the other seven flatten this one.

> **Evolutionary constraint is predicted almost entirely by contact shape, not by
> contact amount.** LOEUF: counting peaks gives r = 0.031, the 11-feature
> magnitude basis gives 0.082, the full 91 features give **0.210**, and the
> amount-corrected shape substrate alone still gives **0.191**. So 91% of the
> full-model performance survives removing magnitude entirely.

That is a 4.08 SD margin, by some distance the largest in the table, and it is
the one target where the architecture claim is doing real work rather than
riding on signal level. It is quotable as it stands.

**And it is consistent with Aim 3 rather than in tension with it.** The
top-quartile constraint set (`gnomAD_pLI_topQ`) barely displaces, |d| = 0.17,
and conservation (`phastCons_2kb_topQ`) sits at 0.16. Architecture predicts
constraint well as a **continuous quantity** while the most-constrained genes do
**not** occupy a region of the space. That is the "directions, not regions"
claim stated twice from independent analyses, and it is worth a sentence in the
report because it makes an abstract claim concrete.

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
**None separates.** Largest effect |d| = 0.73 leaves **71% overlap**; most sets
sit above 90%.

**Report magnitudes, never signed d, when comparing sets.** Each set's effect is
measured on *its own* strongest component, PCA sign is arbitrary per component,
and different components are not a shared scale. A table sorted on signed d
implies a spectrum from most-positive to most-negative that does not exist. The
signs in the table below are retained only so each row can be traced back to its
component; they carry no cross-set meaning.

**A coverage ceiling applies, mirroring the existing size floor.** The cohort
convention already excludes sets below 25 genes as too small to interpret. The
symmetric rule was missing and is now `MAX_GROUP_COVERAGE = 0.70` in
`store_schema.py`: a set covering most of the panel cannot produce an
interpretable "set versus rest" contrast at any effect size, because the
comparison group is whatever is left over and is defined only by exclusion.

Three sets exceed it and are reported separately rather than mixed in:

| set | n | % of panel | "rest" |
|---|---|---|---|
| ChromHMM_active_TSS | 1,619 | 88% | 227 |
| DICE_top_TPM_quartile | 1,539 | 83% | 307 |
| CpG_island_promoter | 1,526 | 83% | 320 |

They are **flagged, not deleted**. Removing them invites the question of why a
set covering most of the panel is missing, and the honest answer is more useful
than a silent omission: the row exists, it is simply not evidence in either
direction. It also explains the otherwise surprising near-null for CpG-island
promoters, which are architecturally distinct in the literature but describe 83%
of this panel, so the contrast has little to work with.

**Nothing that matters depends on this.** The headline comparison sits well
inside the bound: super-enhancers cover 8.6% of the panel and Lambert_TF 8.4%.
The definition-type result is unchanged by dropping all three, p = 0.0019
against 0.0018, medians 0.28 against 0.14.

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

### Chromosome is a confound, and it costs us the Lambert_TF control

Found 2026-08-16 (`diagnose_chromosome_confound.py`). This is the fourth
confound class discovered in the external-set analysis and the only one that
removes a positive control.

Splitting Lambert TFs by DNA-binding domain showed the entire TF effect coming
from **C2H2 zinc fingers** (|d| 0.75) with non-C2H2 TFs at nothing (p = 0.52).
But C2H2-ZF genes are **38% chr19** in this panel against 6% elsewhere, because
of the KRAB-ZNF arrays, and **chr19 membership by itself displaces at |d| = 0.63,
p = 0.0005**. Every other chromosome sits at 0.23 to 0.30. So chr19 has
distinctive contact architecture and any set enriched there inherits it.

Two tests, answering different questions, both reported:

| set | all genes | excluding chr19 | chromosome controlled | chr19 share |
|---|---|---|---|---|
| **Lambert_TF** | 0.46, p 0.0005 | **0.20, p 0.16** | 0.39, p 0.0005 | **25%** |
| GWAS_immune_hot | 0.32, p 0.0005 | **0.33, p 0.0005** | 0.40, p 0.0005 | 7% |
| Roadmap_silenced | 0.73, p 0.0005 | 0.72, p 0.0020 | 0.77, p 0.0005 | 7% |
| ChromHMM_bivalent | 0.58, p 0.0005 | 0.58, p 0.0005 | 0.59, p 0.0005 | 7% |
| dbSUPER SE | 0.16, p 0.18 | 0.18, p 0.36 | 0.20, p 0.14 | 6% |
| gene_desert | 0.71, p 0.0005 | 0.68, p 0.0005 | 0.60, p 0.0005 | 0% |

**But chr19 is gene density wearing a chromosome label**, which partly reverses
this. chr19 is the most gene-dense chromosome in the genome and it is so in this
panel: median density 29.0 against 11.0 elsewhere, Cohen's d = +1.24. Its
displacement collapses under the density control already in use:

| chr19 membership | \|d\| |
|---|---|
| magnitude only | 0.63 |
| magnitude + density controlled | 0.28 |
| density-decile stratified | **0.10, retaining 16%** |

So chromosome is not a new confound class. It is the existing density confound
concentrated on one chromosome, and density is already controlled everywhere it
matters.

**Lambert_TF is therefore qualified rather than demoted.** With density
controlled it holds at **|d| 0.39, p = 0.0005**, meaning transcription factors
differ from equally-dense non-TFs. What it does not do is separate among
non-chr19 genes alone (0.20, p = 0.16). Those two facts are consistent: the TF
effect is real after density adjustment but is concentrated in the dense,
largely chr19 part of the set. Excluding chr19 removes a quarter of the set
including its densest members, which over-corrects. Report the density-adjusted
number and state the concentration.

**`GWAS_immune_hot` is nonetheless the cleanest positive control and should
lead.**
It is 7% chr19, unchanged by exclusion (0.32 to 0.33) and stronger with
chromosome controlled (0.40). It is variant-defined, so it touches no chromatin
assay, which is the property Lambert_TF was being relied on for. Set size 130
against the super-enhancer set's 158 is a looser match than 156 against 158, but
still comparable.

**The headline is unaffected.** Super-enhancers remain null under all three
treatments (0.16, 0.18, 0.20; p 0.18, 0.36, 0.14). The definition-type result
holds with chromosome controlled (p = 0.0078 either way, DIRECT median rising
from 0.28 to 0.36).

Two smaller casualties: `Eisenberg_HK` and `bio_bulk_k3` both fail the exclusion
test, so their marginal significance was partly locus-driven too.

**What this means for the chromatin-independence argument.** Three of the four
sets that separate are ChIP-derived (ChromHMM active TSS, bivalent, Roadmap
silenced), and chromatin state and 3D contact are both downstream of the same
biology, so those are partly two assays measuring one thing. With Lambert_TF
demoted, **`GWAS_immune_hot` is now the only genuinely chromatin-independent
positive control.** That is one line of evidence where there were two, and it is
the strongest argument yet for building a mechanism-defined set from outside
chromatin entirely (Phase 1.1).

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

## 6b. NEW 2026-08-16: contact asymmetry is strand-dependent, and the feature destroys it

Found while sourcing strand for the core-promoter test (S1a), and it reverses a
conclusion recorded on 14 August.

**The earlier test used a broken strand source.** Strand was derived from
`cd4_rna_plus_promoter` versus `cd4_rna_minus_promoter`. Checked against the
annotation in `Lab/Protocol_20k/Genes/Output/01_tss_regions.bed`, that proxy is
**51% accurate**, i.e. chance. It would have been wrong for about 886 of 1,807
panel genes. The conclusion drawn from it, that orientation is measurable but
negligible at d = -0.045, could only ever have returned zero and is withdrawn.

**With annotation strand (1,838 of 1,846 panel genes covered):**

| feature | plus-strand mean | minus-strand mean | d | strand-corrected mean |
|---|---|---|---|---|
| `contact_asymmetry` | **+0.308** | **-0.300** | **0.638** | **+0.304, t = 13.7** |
| `oe_asymmetry_mean_all` | -0.051 | +0.048 | -0.099 | -0.049, t = -2.1 |
| `oe_asymmetry_mean_promoter` | -0.036 | +0.042 | -0.079 | -0.039, t = -1.7 |

Opposite in sign, near-equal in magnitude, and a strand-corrected mean 13.7
standard errors from zero. **Contacts are systematically biased to one side
relative to the direction of transcription, and `contact_asymmetry` as computed
averages that away.** Its docstring says "positive means more signal
downstream", but downstream in *genomic* coordinates, which is upstream for half
the genome.

**Why this matters more than a tidy-up.** d = 0.638 is larger than every
external-set displacement in Section 4 except `Roadmap_silenced`. This is real
signal that the feature definition currently discards, and `contact_asymmetry`
is one of the most reproducible features in the set (rho 0.900). Recovering it
adds a well-measured, biologically interpretable axis rather than removing a bad
one.

The per-peak `oe_asymmetry` family shows only a weak strand effect, consistent
with those features being noise-dominated (rho 0.27 to 0.42) and independently
scheduled for re-supporting or retirement in S2.1b.

**Action: S2.1c.** Orient the profile by transcription direction before
computing asymmetry, using `01_tss_regions.bed` for strand, never the RNA
proxy. This is a third feature change and belongs with the other two, before
the rebuild.

---

## 6c. S1 CLOSED 2026-08-16: the housekeeping question cannot be asked on this panel

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

## 7. Retracted, do not quote

Listed explicitly because all of these appear in `REVIEW.md` Parts 1 to 5 and in
earlier handoffs.

| claim | status |
|---|---|
| "PC1 is not amount" (r = 0.033 with `total_mcc`) | **Retracted.** Against the `MAG_OVERALL` basis, PC1 correlates at **0.623**. PC1 is the amount axis. The error was using one feature as the amount proxy. |
| "Amount is not one quantity" | **Reframed.** Amount is multi-faceted, so it must be measured with the 11-feature basis, not with `total_mcc`. |
| "Regress `total_mcc` out before PCA" | **Superseded.** Use `MAG_OVERALL` via `_shape.corrected_shape`. |
| "Contact architecture predicts evolutionary constraint" | **Retracted 2026-08-16.** Shape adds +0.002 over gene length, expression and CpG. The confound is length, which LOEUF scales with by construction and which no baseline contained. |
| "The 91 features are justified because they predict biology" | **Retracted and replaced.** They add nothing over three gene properties for prediction, but those properties explain only 2.6% of the feature space. Justify by what the features contain, not by what they predict. |
| "Mechanism vs output taxonomy over 7 sets" | **Superseded.** Redone across all 21 with the assignment fixed by annotation procedure; the pattern holds at p = 0.0024, but RANK_SINGLE is n = 1 so nothing is claimable about that category. |
| "Super-enhancers are the weakest of 21" | **Withdrawn twice.** First as over-precise, then outright: on the corrected substrate the SE set has a LARGER effect (\|d\| 0.16) than Eisenberg-HK (0.11). Replaced by the kind-of-definition pattern in Section 4. |
| "Orientation is measurable but negligible, d = -0.045" | **Retracted 2026-08-16.** Computed with an RNA-derived strand proxy that is 51% accurate, i.e. random. With annotation strand, `contact_asymmetry` shows d = 0.638 and a strand-corrected mean at t = 13.7. See Section 6b. |
| "The PER-PEAK asymmetry family is a strand/orientation bug" | **Still refuted** for that family specifically. `|value|` reproduces worse, and strand-relative asymmetry is d = -0.045. It is a support-size problem: median peak ~11 bins at 1-2 reads per bin. |
| "Most of the signal is amount, not shape" | **Softened.** Correct statement is *amount is sufficient for most targets*; residualisation can strip real architecture if amount is downstream of it. |
| Varimax rotation should be adopted | **Withdrawn.** It nearly doubles nameability for free, but its most concentrated factors are the `oe_asymmetry` families at rho 0.27-0.42. Do not adopt without weighting by reproducibility. |
| Displacement figures from the raw components (largest d = 1.61, 21 of 21) | **Superseded** by the corrected substrate in Section 4 (largest d = 0.73, 17 of 21). |
| The app's `/api/enrichment/grid` numbers | **Uncorrected.** Still computed on raw components. Phase 0.3. |
| "Lambert_TF is the primary positive control, matched-n against SE" | **Qualified 2026-08-16, not withdrawn.** 25% of the set is chr19 and it does not separate among non-chr19 genes alone (\|d\| 0.20, p = 0.16), but chr19 is density and Lambert survives density control at 0.39, p = 0.0005. GWAS_immune_hot leads instead because it is unaffected by every treatment. |
| "The TF signature is promoter-driven and CTCF-poor" | **Narrowed.** The effect is C2H2 zinc fingers only (non-C2H2 TFs p = 0.52), and the C2H2 effect does not survive chr19 exclusion. |
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
