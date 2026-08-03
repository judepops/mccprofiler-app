# Vocabulary the gene finder needs to understand

This file is loaded into the prompt of the translation step in `app/translate.py`.
It exists so that when someone types a term from the literature, the model knows
what the term means and, more importantly, what it does and does not correspond
to in this dataset.

It is a glossary, not a review. Every entry answers three questions: what does
the term mean, what does it map to here, and where does the mapping break. The
third question is the one that matters. Several widely used terms have no
counterpart in this store, and inventing a filter for them produces a confident
answer to a question the data cannot address.

Definitions are drawn from the project literature notebook (notebook
`1b49ca44`, pulled 2026-08-03). Citations are given so a claim can be traced.

---

## TWO RULES THAT OVERRIDE EVERY MAPPING BELOW

**1. Queries are built from features, not principal components.** There is no
`axis` filter type in the schema. A component is a mixture: measured in
`scripts/diagnose_pc_names.py`, the concept each PC is named after carries only
22 to 35 percent of that axis, so selecting on one selects on several things at
once. A feature is a single measured quantity and means what it says. Where an
entry below explains a contrast in terms of a component, that is background for
understanding the data, not a filter you can build.

**2. Selection is on contact architecture alone. External gene sets are never a
filter.**

Where an entry below says a term "maps to" a named set such as `Eisenberg_HK`,
`GWAS_immune_any` or `DepMap_curated_essential`, that identifies the set the
term corresponds to **for reading results afterwards**. It is not permission to
select on it. The `cohort` filter type is absent from the query schema, so there
is no vocabulary for it.

The reason is the whole argument of the project. The claim is that contact
architecture carries information that borrowed labels do not. Choosing genes by
those labels and then describing their architecture reverses the direction of
evidence and guarantees the answer. It is the same error as calling
super-enhancers a distinct class after defining them by a threshold on the
signal you then report.

So: if a question asks for something biological but not architectural (immune,
housekeeping, essential, conserved, expressed, disease-associated), that clause
goes in `unsupported`. Do not approximate it with an axis or a group either.
Every returned gene is annotated with the sets it belongs to, so the overlap is
still visible, as an observation about an architecture-selected list rather than
a property it was selected for.

---

## 0. What is actually measured here

Micro Capture-C contact profiles from human CD4+ T cells, one profile per gene,
anchored on an experimental viewpoint, plus or minus 1 Mb at 50 bp resolution.
From each profile, 91 position-invariant features are computed. Genes are placed
in a continuous coordinate system by PCA on those features, but **queries are
built from the 91 features, not from the components**. See the rule below.

Three consequences govern every mapping below.

**This measures contact architecture, not transcription.** No feature reads
expression, chromatin state, or sequence. When a user asks for something
transcriptional ("highly expressed", "silenced"), no feature answers it and the
clause is `unsupported`.

**The groups are a resolution choice, not a discovery.** Gap statistic returns
k=1, HDBSCAN returns one cluster, the dip test is unimodal on PC1 to PC5, and
44% of active genes are mixtures. The groups are reproducible once imposed
(Cohen's kappa 0.72 across independent captures) but they are not discoverable.
Treat every group label as a region of a continuum.

**Micro Capture-C resolves punctate contacts, not domains.** Hi-C sees a TAD as
a triangle of enriched signal. MCC sees the individual base-pair-resolution
contacts inside it. Any question phrased at domain scale is being asked of the
wrong assay.

---

## 1. Gene class terms

### housekeeping gene
Literature: a gene expressed ubiquitously and at stable level across tissues.
The standard human reference list is Eisenberg and Levanon (2013). Operational
definitions vary widely and disagree: the Human Protein Atlas annotates roughly
9,000 genes as housekeeping while stricter ubiquity criteria return far fewer,
so "housekeeping" is a family of overlapping lists rather than one set.

In the regulatory sense (Haberle and Stark 2018; Zabidi et al. 2015),
housekeeping genes have broad or dispersed core promoters that overlap CpG
islands and lack a TATA box, draw on TSS-proximal enhancers that are shared
across cell types, and combine those enhancers additively.

**Maps to:** cohort `Eisenberg_HK`. Related external sets: `bio_HK_k3`,
`CpG_island_promoter`.

**Does NOT map to:** the group `arch-HK`, despite the name. This is measured and
is the single most important trap in this file. `arch-HK` is 38.2% Eisenberg
housekeeping against `arch-ME-constitutive`'s 40.2%, so it is not even the most
housekeeping-enriched group. It also has the lowest median blood expression of
the three active groups (4.5 TPM against 11.8 and 9.9) and is average on
promoter-drivenness. Its real distinction is dispersed contact architecture.
Use the cohort for the biological question and the group only for the
architectural one.

### developmental gene, cell identity gene, master regulator
Literature: genes whose expression is restricted and stage- or lineage-specific.
Regulatory signature is the mirror of housekeeping: focused or sharp core
promoters with TATA, Inr and DPE elements, distal and cell-type-specific
enhancers, super-additive enhancer combination (Loubiere et al.), and
transcription factors with long intrinsically disordered regions. "Cell identity
gene" and "master regulator" in the Whyte and Hnisz sense are defined by
proximity to a super-enhancer, which makes them a super-enhancer derivative
rather than an independent category.

**Maps to:** `bio_dev_TF_k3`, `Lambert_TF` for transcription factors,
`cd4_specific_immune` for lineage-restricted immune genes.

**Careful:** because cell identity genes are defined via super-enhancers,
answering an identity question with the super-enhancer cohort is circular. Say
so rather than hiding it.

### essential gene, constrained gene
Literature: essentiality is measured by CRISPR-Cas9 dropout screens, DepMap
being the standard; core essential genes are roughly 10% of the genome.
Constraint is a population-genetic measure of depletion of loss-of-function
variation, reported as LOEUF or pLI from gnomAD. The two are related but not the
same: essentiality is cell-autonomous fitness in culture, constraint is
selection in the human population.

**Maps to:** `DepMap_curated_essential`, `DepMap_inferred_essential` for
essentiality; `gnomAD_pLI_topQ` for constraint; `phastCons_2kb_topQ` for
sequence conservation.

### tissue-specific gene, expression breadth
Literature: breadth is the number of tissues in which a gene is expressed;
specificity is often summarised by the tau statistic. DICE defines immune
cell-type specificity as TPM at least 1.0 with z at least 2 and fold change at
least 2.

**Maps to:** `DICE_top_TPM_quartile`, `cd4_rna_top_quartile`,
`cd4_specific_immune`. Silenced is `Roadmap_silenced`.

**Not a feature.** Expression breadth is not measured from contacts. If a user
asks for "broadly expressed genes with X architecture", the breadth half is
`unsupported` and only the X half becomes a filter.

### immune GWAS gene
Literature: a gene implicated by variants at genome-wide significance
(p < 5e-8) in immune-mediated disease, or an eGene with a significant cis-eQTL
in immune cells. Assigning the causal gene at a locus is the hard part, and one
published route uses the promoter with the highest number of unique Micro
Capture-C contacts, which is this assay.

**Maps to:** `GWAS_immune_any` (broad, 733 genes), `GWAS_immune_hot`
(restrictive), `GWAS_total_topQ` (GWAS burden regardless of trait).

### gene desert
Literature: a large region, often over 500 kb, without protein-coding genes.
Disease variants in deserts frequently act over long range on a distant target.

**Maps to:** `gene_desert_bottomQ_density`, genes in the lowest quartile of
local gene density.

**Note:** this cohort is the positive control for the density confound. Its
apparent architectural signal collapses under gene-density stratification while
insulation-based effects survive. If a user's question depends on it, the
stratified result is the one worth reporting.

---

## 2. Regulatory element terms

### enhancer, promoter, CTCF site
Here these are peak classes on the annotated peak set, assigned by chromatin
context. They are the element classes over which roughly half the 91 features
aggregate.

**Maps to:** element-class features, for example `n_peaks_enhancer`,
`enhancer_signal_fraction_raw`, `n_peaks_promoter`,
`promoter_signal_fraction_raw`, `n_peaks_ctcf`, `ctcf_signal_fraction`. Also
the element-class features directly. For a property OF one class, use the feature naming both, for example `max_distance_to_viewpoint_enhancer` for long-range enhancer contacts.

### super-enhancer
Literature: introduced by Whyte et al. (2013) and Hnisz et al. (2013) for large
clusters of enhancers densely occupied by master transcription factors and
coactivators. Operationally called by ROSE in three steps: call peaks on
H3K27ac or a master factor, stitch enhancers within 12.5 kb, rank by total
signal, and take everything above the inflection point where the slope of the
signal-versus-rank curve equals 1, typically under 3% of enhancers.

The critique matters as much as the definition. Pott and Lieb (2015) argue
super-enhancers are the extreme tail of a continuous signal distribution with no
biological rationale for either the 12.5 kb stitching distance or the rank
cutoff. High activity may be explained by additive contributions of components
rather than emergent synergy, and roughly 15% of super-enhancers are single
unclustered enhancers, so clustering is not even necessary for the label.

**Maps to:** cohort `dbSUPER_CD4_SE_TSS_pm50kb`, and **only as a comparison
set**. Super-enhancer labels never enter the feature path. Using them as input
would be circular, because a central argument of this project is that
"super-enhancer" is an over-trusted threshold on a continuum and that
architecture-derived descriptions are more informative. The continuum finding
here is the same shape of argument Pott and Lieb make about super-enhancers, so
the app should not repeat the error it criticises.

### stretch enhancer, locus control region, enhancer hub, facilitator
Literature, all distinct from super-enhancer and from each other. Stretch
enhancers (Parker et al. 2013) are single enhancers over 3 kb, and are roughly
ten times more numerous than super-enhancers. Locus control regions (Grosveld et
al. 1987) are the older concept for element sets controlling a gene cluster,
the beta-globin LCR being the type case; super-enhancers are essentially a
genome-wide generalisation of it. Enhancer hub or cluster emphasises 3D spatial
co-localisation rather than 1D stitching, so it is the closer concept to what
this assay measures. Facilitators (Blayney et al. 2023) are elements inside a
cluster with no intrinsic enhancer activity that amplify their neighbours.

**Maps to:** nothing directly. There is no stretch-enhancer, LCR or facilitator
annotation in this store. "Hub" is the one term with a real architectural
analogue: many distal enhancer contacts spread out rather than concentrated,
which is high `n_peaks_enhancer` with high `signal_entropy` and low
`frac_signal_in_top_peak`. Build that from features and say that is what you did,
rather than pretending a hub annotation exists.

### CpG island promoter, bivalent or poised promoter
Literature: CpG island promoters are the dispersed, TATA-less mammalian
housekeeping-type class. Bivalent promoters carry both H3K4me3 and H3K27me3 and
are the poised developmental class in stem cells.

**Maps to:** `CpG_island_promoter`, `ChromHMM_bivalent`, `ChromHMM_active_TSS`.

---

## 3. 3D contact terms

### local versus long-range, proximal versus distal
The primary architectural contrast in this dataset.

**Maps to:** `frac_local`, `frac_distal`,
`frac_far_distal`, `local_to_distal_ratio`, `distal_signal_density`,
`frac_promoter_proximal`, and the `max_distance_to_viewpoint_*` family.

PCA sign is arbitrary, so the pole must be read from the axis poles supplied in
the prompt and never guessed. This is the most common way to get a query exactly
backwards.

### chromatin loop, CTCF loop, convergent motif rule
Literature: a loop is a physical contact between two distal sites held by
protein complexes, spanning 1 kb to over 2 Mb. CTCF loops are anchored by
cohesin and CTCF, average around 360 kb, and over 90% have convergently oriented
CTCF motifs (Rao et al. 2014). Micro Capture-C localises the contact to the
central consensus motif rather than to a general anchor region.

**Maps to:** the CTCF element-class features. There is no loop-calling
step in this pipeline, so there is no per-loop object to filter on. A question
about "genes with strong CTCF loops" becomes a question about CTCF contact
features.

### loop extrusion, TAD, sub-TAD, compartment, insulation, boundary
Literature: loop extrusion is cohesin processively extruding DNA until stalled
at convergent CTCF sites, and generates TADs at 200 kb to 1 Mb. Sub-TADs are
nested, weaker, more cell-type specific. A and B compartments are the
megabase-scale segregation of active and inactive chromatin, read off Hi-C by
PCA. Insulation is the drop in contact frequency across a boundary.

**Maps to:** nothing. **These are all domain-scale Hi-C concepts and this store
has no domain calls, no insulation score and no compartment assignment.** The
window is plus or minus 1 Mb around a single viewpoint, which is a different
object from a genome-wide contact matrix. If a user asks about TADs,
compartments or insulation, say the store cannot answer it rather than
substituting a long-range contact filter, which measures something else.

### contact frequency, distance decay, observed over expected
Contact frequency falls steeply with genomic distance, so raw contact counts are
dominated by distance. Observed over expected divides by the panel-mean expected
curve at each distance, which is what makes contacts at different separations
comparable.

**Maps to:** the `oe_*` feature family, computed server-side from the stored
expected curve using the pipeline's own epsilon so the app agrees with
mccprofiler. `total_mcc` is the undivided amount.

### contact amount, signal depth
How much contact signal a gene has in total, before any question of shape.

**Maps to:** feature `total_mcc`. **Not `pc1`.**

Measured 2026-08-03: `r(PC1, total_mcc) = +0.033`, essentially zero. PC1 tracks
peak counts and reach, not signal. Total signal sits on **PC3** (`-0.621`) and
PC6 (`-0.433`). "Amount" is not even a single quantity here: total signal and
number of confident peaks correlate at `-0.127`, so the genes with the most
signal are not the genes with the most peaks.

**Two consequences.** Filtering on PC1 does not filter on amount. And the
promoter-versus-enhancer contrast on PC3 is entangled with amount, so a result
read off PC3 may be partly a statement about how much signal a gene has. Amount
is not quarantined on PC1; it never was.

### how far a name can be trusted

Every axis name is a summary of one contrast within a mixed axis, not a
definition. Measured share of each axis's squared loading mass that sits in the
concept its name refers to:

    PC1  23%   PC2  26%   PC3  31%   PC4  35%   PC5  22%

All are well above chance, so no name is invented, but each accounts for roughly
a quarter to a third of its axis. PC1 to PC5 carry 45% of total variance while
19 components clear the noise ceiling at 78%, so about a third of the real
structure has no name at all. When a question names something a **feature**
measures directly, prefer the feature over the axis.

---

## 4. Terms with no counterpart here

State plainly that the store cannot answer, rather than approximating. This list
is not exhaustive; the test is whether a real filter exists.

- TAD, sub-TAD, boundary, insulation score, A/B compartment, microcompartment
- called loops, loop anchors, CTCF motif orientation
- stretch enhancer, locus control region, facilitator annotations
- eQTL effect sizes, fine-mapped causal variants, individual SNPs
- transcription rate, nascent transcription, Pol II occupancy
- **membership of any external gene set, as a selection criterion.** The sets
  exist and are reported alongside results, but they cannot choose the results.
  See the rule at the top.
- allele-specific or single-cell contacts; every profile is a population average
- any tissue other than human CD4+ T cells
- anything about a gene not in the panel; the genome-wide panel is 1,846 genes
  after QC and outlier removal

---

## 5. Traps, in priority order

1. **Housekeeping is a cohort, never the `arch-HK` group.** Measured above.
   `arch-ME-constitutive` is the promoter-driven and slightly more
   housekeeping-enriched group, which is precisely what its own inherited name
   contradicts.
2. **PC1 is amount, not biology.**
3. **PC3's negative pole is not purely enhancer.** Its strongest negative
   loading is `total_mcc` (-0.227), with `n_peaks_enhancer` (-0.216) and
   `raw_peak_max_max_ctcf` (-0.193). So the negative pole is enhancer plus CTCF
   plus amount. This is one reason axis filters were removed: ask for
   enhancer-driven with `enhancer_signal_fraction_raw` high, which means only
   that.
4. **Axis poles must be read, not guessed.** PCA sign is arbitrary. Long-range
   is the low end of PC2.
5. **Super-enhancer labels are display-only** and are never a feature.
6. **ATAC gates membership and is not a feature.** ATAC determines which genes
   and peaks exist at all, so validation against ATAC-derived quantities is not
   independent. It is never a clustering input, and "drop the ATAC features" is
   a no-op because there are none.
7. **Prefer few filters.** With 1,846 genes, a conjunction of four or more
   narrow filters lands on a handful by chance. A small result is a warning, not
   a finding.
8. **Group labels carry a posterior.** 44% of active genes are mixtures. A
   group filter selects a region of a continuum, not a category.

---

## Sources

Haberle and Stark 2018, eukaryotic core promoters. Zabidi et al. 2015,
enhancer-core-promoter specificity. Loubiere et al., enhancer cooperativity.
Whyte et al. 2013 and Hnisz et al. 2013, super-enhancers. Pott and Lieb 2015,
the super-enhancer critique. Parker et al. 2013, stretch enhancers. Grosveld et
al. 1987, locus control regions. Blayney et al. 2023, facilitators. Rao et al.
2014, convergent CTCF loops. Eisenberg and Levanon 2013, housekeeping genes.
Hennig 2015 and Rousseeuw 1987, on what licenses naming clusters at all.

Pulled from the project literature notebook `1b49ca44` on 2026-08-03. Regenerate
by re-querying that notebook if the source set changes.
