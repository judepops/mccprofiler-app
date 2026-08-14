# MCCProfiler: what the app shows, what the features found, and where the signal is not

Written 2026-08-14 for Jude, as input to a decision about direction. Every
number here is pulled from the built store or the audit outputs, not from
memory, and the scripts that produce them are named so each can be re-run.

The short version, stated first because it is the thing that matters:

> The continuum result is real and multiply confirmed, and no amount of extra
> features will turn it into clusters if the biology is continuous. But the
> current 91-feature representation is also genuinely underpowered in ways that
> are measurable and fixable, and most of its predictive power comes from
> **how much** signal a gene has rather than from the **shape** of its contacts.
> Those are two different problems and they need separating before choosing a
> direction.

---

## Part 1. What the app is, end to end

### 1.1 The chain

```
bigWigs (4 MCC experiments, 8 channels)
  └─ process/  ──────────────► results_combined_50bp_clean.pkl   (4.4 GB)
        signal extraction, QC, outlier removal        1,846 genes x 40,000 bins x 8 ch
  └─ mccprofiler/ ───────────► features_zscored.pkl
        91 position-invariant features                1,846 x 91
  └─ cluster/ ───────────────► archetype_labels.tsv   (k=4 KMeans)
  └─ audit/ ─────────────────► posteriors, dimension names, reproducibility,
                               nested baselines, structure vs noise
                                        │
                                        ▼
                    mccprofiler-app/backend/scripts/build_store.py
                                        │
                       profiles.h5 (1.6 GB) + parquet tables + manifest.json
                                        │
                                        ▼
                             FastAPI  ──►  React frontend
```

The app reads a **snapshot**. `build_store.py` runs in `cd4env`, opens the 4.4 GB
pickle once, and writes a self-contained store; the server never imports
`mccprofiler` and never opens the pickle. Current snapshot: panel `gw_cd4_1`,
built 2026-08-03, pipeline commit `96e13e4`.

### 1.2 What is read versus computed

This distinction matters when reconciling the app against a figure elsewhere.

| Layer | Origin |
|---|---|
| Profiles | read: `results_combined_50bp_clean.pkl`, resampled to 3 pyramid levels |
| 91 features | read: `features_zscored.pkl` |
| Peaks | read: `annotated.tsv` (Nick Denny's calls, ATAC-intersected) |
| Group labels | read: `archetype_labels.tsv` composed with `me_subtypes_labels.tsv` |
| Posteriors | read: `posterior_membership.tsv` |
| **PCA, loadings, scree, noise ceiling** | **computed in the builder** |
| **UMAP and its permuted null** | **computed in the builder**, seed 0 |
| gcPCA | read: `gcpca_axes.npz` |

The PCA is computed but *checked*: `verify_store.py` asserts the variance shares
reproduce `dimension_names.tsv` exactly (15.27 / 9.47 / 8.21 / 6.33 / 5.88),
which is what licenses calling PC1-PC5 the named dimensions. The UMAP
corresponds to no upstream file and will not match any deck unless parameters
and seed match.

### 1.3 The three pages

**Gene.** One gene: viewpoint-aligned profile (canvas, brush zoom, 8 channels,
raw or O/E), peak lane, position readout with mixture composition, the continuum
map with this gene marked, the 91 features by block with percentiles, and a
compare box that marks the second gene on the map.

**Panel.** Properties of the population and the coordinate system: the
plain-English finder, ranked lists, cohort comparison, the continuum map (any 2
of 38 axes, UMAP with its null, loadings biplot), the enrichment grid, the scree
and loadings panel, and the explanation layer.

**Lab.** Exploratory views whose caveats are too heavy to ship: currently the
cross-panel reproducibility view.

### 1.4 Design rules the app enforces

1. **Position first, the label is a readout.** The science forbids a category
   badge as primary output, so genes render as coordinates with a mixture
   composition and every label carries its posterior and the continuum caveat.
2. **Every figure is data-backed.** The API returns numbers; the frontend draws.
   No pre-rendered image is ever served.
3. **ATAC gates, it is never a feature.** `load_signal` reads only the MCC
   channel. ATAC determines which genes and peaks exist, so ATAC-derived
   validation is not independent.
4. **External sets are post-hoc only.** The finder cannot select on them: the
   `cohort` filter type is absent from the schema. Membership is annotation.
5. **Principal components are not selectable either**, for the reasons in 3.2.

---

## Part 2. The biological findings

### 2.1 The genes do not form clusters. This is the strongest result in the project.

Six independent lines, none of which depends on the others:

| test | result |
|---|---|
| Gap statistic | **k = 1**, gap declines monotonically |
| HDBSCAN | **one** cluster at min_cluster_size 25 and 50 |
| Dip test | **unimodal** on PC1-PC5, min p 0.88-0.96 |
| Silhouette at k=4 | **0.062** (permuted noise gives 0.009) |
| Mixture fraction | **44%** of active genes are blends; 1,040 of 1,846 exceed posterior 0.8 |
| Mean max posterior | **0.793** |

The k=4 labels used throughout the app are therefore **an imposed resolution
choice, not a discovery**. What licenses using them at all is that they are
*reproducible once imposed*: Cohen's kappa 0.72 at k=3-4 across independent
captures of the same 116 genes.

Group sizes: arch-HK 844, arch-ME-constitutive 369, arch-ME-effector 351,
arch-sparse 261, arch-off 21 (a QC class, near-empty signal).

**The structure that does exist is dimensional, not partitional.** SigClust
z = -11.86 (trimmed) says the cloud departs from a single Gaussian, and 19
components sit above a parallel-analysis noise ceiling against 8 for permuted
noise, together 78.0% of variance. Effective rank 17.6 against 86.6 for noise.
So there is a great deal of real structure. It is simply not partitioned.

### 2.2 The architecture groups, renamed from what they measure

The pipeline's inherited names describe biological categories the groups do not
track. Measured discriminating features (z, group vs rest):

| canonical | app display | top discriminating features |
|---|---|---|
| `arch-HK` | **dispersed** | signal_entropy +1.25, empty_band_fraction -1.17, frac_far_distal +1.09, mean_degree +1.01 |
| `arch-ME-constitutive` | **promoter-local** | promoter_signal_fraction_raw +1.13, n_peaks_promoter +0.93, frac_local +0.80 |
| `arch-ME-effector` | **enhancer-focal** | raw_peak_max_max_enhancer +1.03, signal_entropy -0.97, enhancer_signal_fraction_raw +0.86 |
| `arch-sparse` | **sparse** | n_high_consensus_peaks_075 -1.62, n_peaks_all -1.57 |

**`arch-HK` is not the housekeeping group.** Eisenberg-HK fraction 38.2% against
`arch-ME-constitutive`'s **40.2%**. It has the *lowest* median blood expression
of the three active groups (4.5 TPM vs 11.8 and 9.9) and is average on
promoter-drivenness. Its one real distinction is essentiality (17.9% vs 13.8%).
Meanwhile `arch-ME-constitutive` **is** the promoter-driven group, which
"multi-enhancer" directly contradicts.

**Caveat that must travel with these names:** they come from *unadjusted* group
means. Given 2.4 below, re-derive on the amount-corrected substrate before any
of this is written up.

### 2.3 The named dimensions are much looser than their names imply

`scripts/diagnose_pc_names.py`. Share of each axis's squared loading mass in the
concept its name refers to, concept features chosen from the concept rather than
from the loadings so the test is not circular:

| PC | name | mass in named concept | lift | r with `total_mcc` |
|---|---|---|---|---|
| 1 | amount, richness, reach | 23.3% | 1.93x | **+0.033** |
| 2 | local vs long-range | 26.1% | 1.69x | +0.150 |
| 3 | promoter vs enhancer/CTCF | 30.9% | 3.12x | **-0.621** |
| 4 | enhancer vs CTCF | 34.7% | 3.95x | -0.266 |
| 5 | concentrated vs dispersed | 22.0% | 2.85x | -0.006 |

Every name is well above chance, so none is invented, but **each accounts for
roughly a quarter to a third of its own axis.** PC1-PC5 carry 45.2% of variance;
19 components clear the noise ceiling at 78.0%, so about a third of the real
structure has no name at all.

**PC1 is not amount.** It correlates with `total_mcc` at +0.033. It tracks peak
counts and reach. Total signal lives on **PC3 (-0.621)** and PC6 (-0.433). And
amount is not one quantity: `r(total_mcc, n_high_consensus_peaks_075) = -0.127`,
so the genes with the most signal are not the genes with the most peaks.

The consequence is serious: the promoter-vs-enhancer contrast, the most quotable
dimensional result, sits on the axis most contaminated by how much signal a gene
has.

### 2.4 Most of the biological signal is amount, not shape

`audit/continuous_methods/nested_baselines.tsv`, nested and cross-validated,
linear models, AUC for binary targets and r for continuous:

| target | n_peaks only | amount (3) | all 91 | **shape-corrected (80)** |
|---|---|---|---|---|
| GWAS_immune_hot | 0.534 | 0.691 | **0.746** | 0.629 |
| cd4_rna_top_quartile | 0.606 | 0.647 | **0.707** | 0.576 |
| DepMap_curated_essential | 0.547 | 0.682 | **0.700** | 0.550 |
| Eisenberg_HK | 0.542 | 0.630 | **0.637** | 0.541 |
| Lambert_TF | 0.581 | 0.566 | **0.661** | 0.639 |
| blood TPM (log) | 0.222 | 0.242 | **0.380** | 0.132 |
| gtex_tau | 0.055 | 0.242 | **0.276** | 0.078 |
| gnomad_loeuf | 0.031 | 0.071 | **0.210** | 0.191 |
| phastcons_2kb | 0.128 | 0.127 | **0.187** | 0.138 |

**The 91 features beat counting peaks on 9 of 9 targets**, median gain 0.153.
That is the result that justifies the whole feature substrate existing.

**But look at the last column.** Remove amount and most of the gain disappears:
Eisenberg_HK falls to 0.541, back to the n_peaks baseline. `cd4_rna` 0.707 →
0.576. blood TPM 0.380 → 0.132, *below* n_peaks. Three targets survive amount
removal in a meaningful way: **gnomad_loeuf** (0.210 → 0.191, so constraint is
carried almost entirely by shape), **GWAS_immune_hot** (0.746 → 0.629), and
**Lambert_TF** (0.661 → 0.639).

This is the finding that should drive the next decision. The interesting claim,
that *architecture* carries information, is currently supported strongly for
constraint and immune GWAS, weakly for TFs, and barely at all for expression and
housekeeping status, where amount does the work.

### 2.5 External annotations are displaced but not separated

`scripts/diagnose_external_structure.py`. All 19 above-noise dimensions,
label-permutation nulls, BH correction over 399 tests.

- **21 of 21 sets** are displaced from the panel centroid, every one p <= 0.0045
- **396 of 399** set-by-dimension tests survive BH correction
- Largest effect anywhere: **d = 1.61**, which still leaves **42% overlap**
- Most sets sit at d 0.2-0.6, i.e. **76-92% overlap**

| set | strongest | d | overlap |
|---|---|---|---|
| Roadmap_silenced | PC1 | -1.61 | 42% |
| ChromHMM_bivalent | PC2 | -1.15 | 57% |
| GWAS_immune_hot | PC3 | -0.91 | 65% |
| Eisenberg_HK | PC3 | +0.37 | 85% |
| **dbSUPER super-enhancers** | PC2 | +0.29 | **88%** |
| cd4_specific_immune | PC4 | +0.13 | 95% |

Two readings, both true: there **is** universal positional structure, and it is
**invisible** because displacement is nowhere near separation. This is the
continuum result one level up: external categories map to *directions* in this
space, not to *regions* of it. That is why the enrichment grid looks like the
sets are everywhere. They are.

Two warnings attached. The **largest displacement of all is the density positive
control** (`gene_desert_bottomQ_density`, z = 26.7), so the strongest apparent
signal in that figure is a known confound. And **super-enhancers are the weakest
of the 21**, which is the sharpest version yet of the thesis argument about that
label: SE genes are not a class in contact-architecture space.

### 2.6 Reproducibility, and which parts of the space to trust

116 genes captured in both panels, 63 shared features. Median rho **0.752**,
40 of 63 above 0.7, 1 below 0.3.

**The spread is not uniform, and it maps onto the dimensions in a way that
matters.** Reproducibility-weighted loading mass, the share of each component
carried by features with rho > 0.7:

| PC | var % | trustworthy mass | top feature (rho) |
|---|---|---|---|
| 1 | 15.27 | 80% | mean_degree (0.75) |
| 2 | 9.47 | 82% | frac_far_distal (0.98) |
| 3 | 8.21 | 81% | promoter_signal_fraction (0.85) |
| 4 | 6.33 | 78% | enhancer_signal_fraction_raw (0.86) |
| 5 | 5.88 | 82% | mean_peak_gap_bp (0.82) |
| 6 | 4.29 | 56% | oe_tailedness_mean_all (0.66) |
| 7 | 3.59 | 51% | signal_profile_kurtosis (0.74) |
| 9 | 3.03 | **22%** | oe_asymmetry_mean_all (**0.27**) |

Most reproducible features: `frac_local` 0.984, `frac_far_distal` 0.977,
`local_to_distal_ratio` 0.976, `frac_distal` 0.964. **Distance composition is
close to perfectly reproducible.**

Least reproducible: `oe_asymmetry_mean_all` **0.274**,
`distance_to_nearest_peak_bp` 0.317, `oe_asymmetry_mean_promoter` 0.366,
`oe_asymmetry_mean_ctcf` 0.366.

Run in full (`scripts/diagnose_dimension_trust.py`), the pattern is systematic
rather than a couple of bad components. Across the 19 components above the noise
ceiling, **67% of the retained structure rests on reproducible features**. Eight
components fall below 50% trusted, together **14.9% of variance**, and every one
of them is dominated by the same two families:

| PC | var % | trusted | top feature (rho) |
|---|---|---|---|
| 9 | 3.03 | 22% | oe_asymmetry_mean_all (0.27) |
| 11 | 2.12 | 48% | oe_tailedness_mean_promoter (0.57) |
| 14 | 1.87 | 24% | oe_tailedness_mean_enhancer (0.65) |
| 15 | 1.71 | 34% | oe_asymmetry_mean_enhancer (0.42) |
| 16 | 1.69 | 32% | oe_tailedness_mean_ctcf (0.68) |
| 17 | 1.62 | 26% | oe_asymmetry_mean_promoter (0.37) |
| 18 | 1.49 | 40% | oe_asymmetry_mean_enhancer (0.42) |
| 19 | 1.38 | 25% | oe_asymmetry_mean_promoter (0.37) |

So the entire tail of the "19 real dimensions" is the asymmetry and tailedness
families, and those are the least reproducible things we measure. The headline
"19 dimensions, 78% of variance" should be qualified: about two thirds of it is
solid, and the last 15% is currently indistinguishable from measurement noise
that happens to be structured.

**This raises the value of fixing asymmetry (4.1) considerably.** If the
orientation problem is real, one bug is degrading eight components.

**It also retracts a recommendation I made earlier today.**
Varimax rotation nearly doubles nameability (top-5 loading mass 36.5% → 64.6%)
for free, since it is orthogonal and preserves gene distances to 4.3e-14. But
its *most concentrated* factors are precisely the asymmetry families, and those
sit at rho 0.27-0.42. Rotation would produce beautifully nameable axes built on
the least reproducible thing we measure. Do not adopt it without weighting by
reproducibility first.

---

## Part 3. Diagnosis, before choosing a direction

### 3.1 Two separate problems, currently conflated

**Problem A: the biology may be continuous.** Six independent tests say so. If
true, no representation will produce clusters, and looking harder for them is a
category error. The honest deliverable becomes a well-characterised continuum
with axes that predict orthogonal biology.

**Problem B: the representation is lossy and partly contaminated.** This is
demonstrably true regardless of A, and it is fixable. Evidence:

1. **Compression.** Each gene is 40,000 bins x 8 channels reduced to 91 scalars.
   That is roughly a 400-fold compression of the MCC channel alone, and every
   feature is position-invariant by design, so the actual *shape* of the profile
   is discarded. If contact architecture has structure that is not expressible
   as a summary statistic, the current features cannot see it.
2. **Amount dominates.** Section 2.4. Remove amount and most biological
   prediction collapses.
3. **A degenerate feature is driving PC1.** `topology.py:104` builds adjacency as
   the complete graph on active peaks, with no contact, distance or interaction
   criterion, so `mean_degree = n_active(n_active-1)/n_peaks` is an algebraic
   function of two counts (confirmed numerically, max |diff| 0.00e+00).
   `mean_degree` is **PC1's second-strongest correlate at 0.763**. A meaningful
   fraction of our largest component is an artefact of a feature that measures
   nothing.
4. **The asymmetry features are near-noise.** rho 0.27-0.42. Note that
   `distance_to_viewpoint` in `annotated.tsv` is **unsigned**, and there is a
   known strand-orientation issue. An asymmetry statistic computed without
   consistent orientation would be exactly this unreproducible. **This looks
   like a bug, not a biological limit, and it is worth an hour to check.**
5. **Element classes are coarse.** Three classes (enhancer, CTCF, promoter) over
   ATAC-intersected peaks. Any finer regulatory grammar is invisible.
6. **Effective dimensionality is 17.6** of 91 (participation ratio), and 32
   components are needed for 90% variance. The correlation filter did its job,
   0 pairs above |r| 0.8, but the space is still far smaller than its nominal
   size.

### 3.2 What "more signal" would have to mean

Adding features to the same summary-statistic paradigm is unlikely to help much:
the space is already 91 features collapsing to ~18 effective dimensions, and the
marginal feature will be another summary of the same profile. The three routes
that could add genuinely new information are:

- **Stop summarising.** Learn from the profile itself.
- **Fix what is broken.** Amount contamination, the degenerate topology block,
  the orientation problem.
- **Change the unit of analysis.** Currently one vector per gene. The peaks, or
  the peak-pairs, or sub-regions of the profile could be the unit instead.

---

## Part 4. Options, with my honest ranking

### Tier 1: do these first, they are cheap and they clean the substrate

**4.1 Fix the orientation / asymmetry problem. This is now the single highest
value item in the list.** rho 0.27 on a feature family that varimax says is the
cleanest axis in the space is a contradiction that almost certainly resolves as
a bug. Section 2.6 shows the damage is not local: **eight of the nineteen real
components, carrying 14.9% of variance, are dominated by the asymmetry and
tailedness families**, and all sit below 50% trusted. Check whether asymmetry is
computed on a strand-oriented profile; `distance_to_viewpoint` in
`annotated.tsv` is unsigned and there is a known strand-flip issue on record. If
orientation is the cause, one fix repairs eight dimensions. Cost: hours to
diagnose. Upside: larger than anything else here.

**4.2 Retire or rebuild the topology block.** `mean_degree` and
`frac_active_pairs` are algebraic functions of counts, and they are loading
heavily on PC1. Either delete them and re-derive the dimensions, or define
adjacency properly (co-contact above an O/E threshold, or distance-constrained)
so the graph means something. Cost: a day. Upside: PC1 becomes interpretable.

**4.3 Amount-correct and re-derive everything.** Regress `total_mcc` out before
PCA. Costs 0.3 points of variance and takes the worst amount correlation from
0.621 to exactly 0.000 (`scripts/experiment_rotate_axes.py`). Then re-derive the
archetype names on that substrate, which is already flagged as overdue. This
also makes section 2.4's shape-corrected column the *default* framing rather
than a robustness check, which is the more honest claim anyway.

### Tier 2: the representation-learning arm, now with a quantitative motivation

**4.4 Masked autoencoder on the raw profile.** This is Aim 2 and it already
exists but is undiagnosed. The argument for it is now numerical rather than
aspirational: 91 scalars from 40,000 bins, position-invariant by construction,
and the shape-corrected column shows the hand-crafted shape features carry
little biology. If there is structure in the profile that is not a summary
statistic, this is the only arm that can find it.

Two design notes that follow from this review:
- Evaluate it on the **shape-corrected** targets from 2.4, not the raw ones.
  Beating `n_peaks` on expression is easy and uninformative; beating the
  amount-only baseline on LOEUF and immune GWAS is the real test.
- Use the **116 twice-captured genes** as a held-out reproducibility check on
  the learned embedding, exactly as done for the features. An embedding that
  scores well but does not reproduce across captures is not usable.

**4.5 Contrastive learning with technical replicates as positive pairs.** The
116 genes captured in both panels are the same biological object measured
twice. Training a representation to place those pairs close together learns
invariance to capture noise directly, which is precisely the failure mode
section 2.6 identifies. This is a small, well-posed, and to my knowledge
unexplored idea in this project, and it addresses the measured weakness rather
than a hypothetical one. n = 116 pairs is small but sufficient for a contrastive
regulariser on top of a reconstruction objective.

### Tier 3: the GNN and the segmentation idea

**4.6 A GNN over peaks, not genes.** You raised segmenting promoters. The
obstacle is that the current graph is the complete graph, so there is no
structure for a GNN to exploit; 4.2 is a prerequisite. If adjacency were defined
by measured co-contact, a per-peak GNN with genes as subgraphs would change the
unit of analysis from gene to peak, which is one of the few ways to get genuinely
new information out of the same data. This is the most interesting of the
options and also the one with the most groundwork required.

### Tier 4: reframe the question

**4.7 Look for conditional structure rather than global clusters.** Every
clustering test so far has been global. Two variants have not been tried and
would be cheap:
- Cluster **within amount strata**. If amount dominates the geometry, clusters
  in shape may exist inside a stratum and be washed out globally.
- Cluster **in the shape-only space** after amount removal. Section 2.4 shows
  this space is nearly independent of the raw one for several targets.

If clusters do not appear there either, that is a much stronger negative result
than the current one, and it is worth having.

**4.8 Accept the continuum and make it the finding.** The strongest version of
the current work is not "we found archetypes". It is: *contact architecture is
continuous, external regulatory categories are displaced but not separated
within it (21 of 21 sets, all p <= 0.0045, all with 42-95% overlap), and the
super-enhancer label in particular has no territory at all.* That is a coherent,
defensible, and genuinely interesting negative result about a widely used label,
and it is exactly the argument the thesis set out to make.

---

## Part 5. What I would do

If it were my project, in order:

1. **4.1 and 4.2 this week.** Both are probable bugs contaminating headline
   dimensions. Neither is a research bet.
2. **4.3 immediately after**, and re-derive the archetype names on the corrected
   substrate. Nothing should be written up on the current names.
3. **4.7 as the decisive test of the clustering question.** One week. If
   conditional structure exists, everything changes; if not, 4.8 becomes the
   thesis argument with much stronger support than it has now.
4. **4.4 and 4.5 in parallel** as the Aim 2 arm, evaluated against the
   shape-corrected baselines and the 116-gene reproducibility check.
5. **4.6 only after 4.2**, because it is otherwise built on a degenerate graph.

The thing I would guard against is treating "find more signal" as the goal. The
measurements say the data has a great deal of real structure (19 dimensions, 78%
of variance, kappa 0.72 reproducibility) and that this structure is continuous
and amount-dominated. More features will not change continuity. Fixing the
contaminated features and removing amount will change *what the structure means*,
and learning from the raw profile is the only route to information the summaries
cannot hold.

---

## Appendix: scripts that produce these numbers

| claim | script |
|---|---|
| PC name accuracy, PC1 is not amount | `backend/scripts/diagnose_pc_names.py` |
| Rotation and amount-removal experiment | `backend/scripts/experiment_rotate_axes.py` |
| External set displacement, all 19 dims | `backend/scripts/diagnose_external_structure.py` |
| Store invariants, 38 checks | `backend/scripts/verify_store.py` |
| Nested baselines | `audit/continuous_methods/nested_baselines.tsv` |
| Structure vs noise | `audit/continuous_methods/structure_vs_noise.tsv` |
| Cross-panel reproducibility | `audit/continuous_methods/cross_panel_reproducibility.tsv` |

| Reproducibility-weighted trust per dimension | `backend/scripts/diagnose_dimension_trust.py` |

All five diagnostic scripts run in `mccapp` except `experiment_rotate_axes.py`,
which needs sklearn and therefore `cd4env`.
