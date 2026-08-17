# Replacing the archetypes with reach x composition

Written 2026-08-16 after the Leiden work. This is the plan for landing it and,
first, what it actually means.

---

## 1. What changes, and what it means

**Before.** Five flat groups (`arch-HK` 844, `arch-ME-constitutive` 369,
`arch-ME-effector` 351, `arch-sparse` 261, `arch-off` 21), fit by KMeans on the
**91-feature** substrate on 2026-07-21, named partly by overlap with external
gene sets (Eisenberg for HK, Lambert for ME-vs-sparse).

**After.** Two orthogonal levels:

| | far-reaching (1,074) | mid-range (772) |
|---|---|---|
| CTCF-dominated | 401 | 196 |
| enhancer-dominated | 356 | 340 |
| promoter-dominated | 317 | 236 |

**Level 1, reach.** How far contacts extend: the 50-250 kb band against
250 kb-1 Mb. This is sPC1, the largest axis in the substrate.

**Level 2, composition.** Which element class dominates the gene's contacts.
This is the sPC2/sPC3 plane.

### Why this is better, in four specific ways

1. **The levels are independent.** Level 2 carries 1-4% of its variance on sPC1.
   Composition is genuinely orthogonal to reach, which is why no flat method
   found it: reach dominates the distance metric and swamps composition, so
   KMeans spent k=2 on reach and k>=4 on 3-gene outlier pockets.
2. **It is not a dichotomised axis.** The k=2 KMeans partition put **97%** of its
   between-centroid variance on sPC1. It scored kappa 0.859 precisely because it
   was trivial: thresholding a strong continuous axis reproduces beautifully and
   tells you nothing. Selecting on kappa alone selects for uninformativeness.
3. **Larger, cleaner effects.** Discriminants run +/-1.0 to +/-1.36 against a
   previous maximum of 1.07, and are computed on the **amount-corrected**
   substrate, so no group can mean "these genes have more signal". That retires
   the caveat standing on the archetype names since 2026-08-03.
4. **No external sets anywhere.** Naming is by architecture only. The old
   procedure picked the HK cluster by best Eisenberg Fisher p, which makes any
   downstream "arch-HK is enriched for housekeeping genes" circular, and which
   misfired anyway: `arch-ME-constitutive` had the higher Eisenberg fraction
   (40.2% vs 38.2%).

### What does NOT change

**The continuum result stands, and is now stronger.** Leiden returns a **single
community** at every resolution below 0.4. That is a third independent method
agreeing with HDBSCAN (0 clusters, 16 conditions) and the gap statistic (still
rising at k=8). These six groups are **named regions of a continuum, not
discovered clusters**, and every caveat that travelled with the old labels
travels with these.

The old partition was not a coarser view of this structure, it cut across it:
`arch-HK`, the largest group, splits **48/52** on the reach axis.

---

## 2. The plan

### Phase 1 - generate and validate the labels (blocking, ~2h)

1. **`audit/scripts/build_taxonomy.py`**, writing a two-column labelling
   (`reach`, `composition`) plus the composed six-way id, for all 1,846 genes.
   Amount-corrected substrate, Leiden at resolution 0.6 within each reach half,
   `random_state` pinned, membership written with the substrate hash so it is
   reproducible.
2. **Re-measure reproducibility for THIS partition.** The kappa 0.72 figure
   belongs to the old KMeans k=3-4 on 91 features and does not transfer. Needed:
   agreement on the 116 twice-captured genes for the reach level and the
   composition level separately. **If composition does not reproduce, it does not
   ship** - that is the whole justification for imposing a partition at all.
3. **Check the ARI = 1.000 at resolution 0.4.** Perfect agreement is suspicious;
   confirm it is not an artefact of both panels collapsing to a trivial split.
4. **Stability under reseeding and subsampling.** Leiden is stochastic. Report
   the fraction of gene pairs co-assigned across 50 seeds.

### Phase 2 - store and schema (~2h)

5. Write `taxonomy_labels.tsv` beside `archetype_labels.tsv`. **Do not overwrite
   the old file**: it is referenced by the manifest, the audit outputs and five
   handoffs. Add it to `paths.py` as a new key.
6. `store_schema.py`: replace `ARCHETYPE_DISPLAY` with `REACH_DISPLAY` and
   `COMPOSITION_DISPLAY`, keep `CONTINUUM_CAVEAT`, drop `ARCHETYPE_PROVISIONAL`
   (it exists because the old labels were stale; these are not).
7. `build_store.py`: read the two-level file, keep `group` as the composed
   six-way id for backwards compatibility, add `reach` and `composition` columns.
   **Retire the `me_subtypes_labels.tsv` join** - it exists only to split the old
   `arch-ME` and has no meaning here.
8. `verify_store.py`: new expected counts (1074/772 and 401/356/317/196/340/236),
   and an assertion that the two levels are near-orthogonal (sPC1 share of the
   composition split stays under 10%).

### Phase 3 - app (~2h)

9. Readout shows both levels, not one badge. Six colours is too many for the
   map; colour by composition, shape or opacity by reach.
10. Cohort and enrichment views group by both levels.
11. Remove the `provisional` chip added this morning.

### Phase 4 - documents (~2h)

12. `CURRENT_FINDINGS.md` section 6, rewritten. Section 6 currently describes the
    five old groups and their caveat.
13. `CLAUDE.md` section 5 (the naming table) and the app `CLAUDE.md`.
14. **`PROJECT_STATUS.md` and the briefing pack**, which still say `arch-HK` is
    "housekeeping / essential / promoter-driven". This has been outstanding all
    day and gets worse with this change, not better. `PROJECT_STATUS.md` is
    auto-generated: fix the generator.
15. Regenerate the six figures; figure 6 is explicitly about the old partition.

### Phase 5 - post-hoc validation only (~1h)

16. Now, and only now, cross the six groups against the external sets. This is
    the **only** legitimate use of Eisenberg, Lambert and dbSUPER here, and it is
    a description, not a validation: the groups were not fit to them.

---

## 3. Risks, stated plainly

- **This is the fourth substrate-or-labelling change today.** Every number in
  `CURRENT_FINDINGS` has moved twice already. The evidence has been stable
  throughout (|d| correlation 0.904 across the last change, no conclusion
  altered), but the bookkeeping debt is real and compounding.
- **Composition may not reproduce.** Untested as of writing. Phase 1 step 2 is a
  genuine go/no-go: if the composition level does not hold across captures, ship
  the reach split alone or nothing.
- **Six groups on 1,846 genes** averages 308 each; the smallest is 196. Fine for
  description, thin for any per-group statistic against an external set.
- **"Reach" is partly a panel property.** The window is +/-1 Mb by construction,
  so "far-reaching" means "far within 1 Mb". Genome-wide capture would redefine
  the axis.

## 4. What this does NOT license

- Calling these clusters. Leiden finds one community below resolution 0.4.
- Claiming the taxonomy is discovered. It is a resolution choice, defended by
  reproducibility and interpretability, not by separation.
- Using the groups as anything but a description until Phase 1 step 2 passes.
