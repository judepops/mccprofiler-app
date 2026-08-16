Yes, you did — let me check the consolidated plan:

Bash Check which feature removals survive in the plan
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
echo "=== feature-removal items in the current plan ==="
grep -n "degenerate\|topology\|asymmetry\|tailedness\|re-support" audit/PLAN_FORWARD.md | head -12
echo
echo "=== is the moment re-supporting still an action item? ==="
grep -n "1.2\|moment famil" audit/PLAN_FORWARD.md | head -6

OUT
=== feature-removal items in the current plan ===
55:**S2.1 Drop the four degenerate topology features** (`mean_degree`,
113:   the literature asymmetry: a null here contradicts a live claim, where the
228:**0.1 Drop the four degenerate topology features upstream.**
324:Eight components rest on `oe_asymmetry` and `oe_tailedness` at rho 0.27 to 0.42.
326:1-2 reads per non-zero bin. `contact_asymmetry`, the same statistic on the whole
343:   and cell state. Carried by the literature asymmetry: a null here contradicts
449:adjacency first (0.1 removes the degenerate one; co-contact above an O/E
490:Phase 0.1. Drop the four degenerate features upstream in `mccprofiler`,

=== is the moment re-supporting still an action item? ===
203:| 3. External categories are displaced but not separated | **Survives, controlled three ways.** 17 of 21 displaced after magnitude correction, largest d = 0.73 leaving 71% overlap. Positive controls retain 84-112% under gene-density stratification. Super-enhancers do not separate (p = 0.18) at n = 158, where the test has 80% power to detect d = 0.35. |
272:confounding (Addendum 4, positive controls retain 84-112% under density
315:*One day. **Highest-value item in Phase 1**, promoted above 1.2. It does two
323:**1.2 Re-support the per-peak moment families.**
326:1-2 reads per non-zero bin. `contact_asymmetry`, the same statistic on the whole
You're right, and it is two items, not one. The consolidated plan kept the degenerate topology drop (S2.1) but lost the moment re-supporting — it survives only in the appendix as the old Phase 1.2.

That matters more than it looks: eight of nineteen components rest on oe_asymmetry/oe_tailedness at ρ 0.27–0.42, carrying 14.9% of variance, and ChromHMM_bivalent's separation sits on one of them at 15% trust. So it touches an Aim 3 control, not just tidiness.

Bash Restore the moment re-supporting item
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app
python3 - <<'PY'
import pathlib
p=pathlib.Path("audit/PLAN_FORWARD.md"); s=p.read_text()

s=s.replace("""## S1. The last test: core-promoter sequence classes""",
"""## S0. DONE 2026-08-16: the orthogonality measurement

Complete, and it is Aim 1's replacement claim. Gene length, expression and CpG
density explain **2.6%** of the 91-feature variance and **0.6%** of the retained
components; zero features exceed R2 0.25, zero components exceed 0.10.
`diagnose_feature_novelty.py`.

---

## S1. The last test: core-promoter sequence classes""")

s=s.replace("""**S2.2 Re-derive the archetype display names**""",
"""**S2.1b Re-support or retire the per-peak moment families.** This was dropped
when the plan was consolidated and should not have been. Eight of the nineteen
retained components rest on `oe_asymmetry` and `oe_tailedness`, which reproduce
at rho 0.27 to 0.42, and together they carry 14.9% of variance. Diagnosed as a
support-size problem: the median peak is about 11 bins at 1 to 2 reads per
non-zero bin, so a third or fourth central moment estimated there has enormous
sampling variance. `contact_asymmetry`, the same statistic computed over the
whole window, reproduces at 0.900, so the concept works and only the support is
wrong. Recompute on aggregated support (stacked peaks per gene, or per element
class) rather than per peak, or retire them.

**This is not cosmetic.** `ChromHMM_bivalent`, one of the Aim 3 corroborating
controls, takes its displacement on shape-PC4, which is 15% trustworthy and
asymmetry-dominated. Fixing or removing those families changes what that control
means. *Two days.*

**S2.2 Re-derive the archetype display names**""")

s=s.replace("""**S2.3 Rebuild the store, re-run the nine diagnostic scripts, regenerate the six
figures.** All are written; this is mechanical. *One hour.*""",
"""**S2.3 Rebuild the store, re-run the eleven diagnostic scripts, regenerate the
figures.** All are written; this is mechanical. **Note the consequence:** S2.1
and S2.1b both change the feature set, so the component count, the variance
shares, the trust weighting and every displacement number move. Aim 2's "19
components, 78% of variance" and Aim 3's effect sizes must be re-quoted from the
rebuilt store, not carried over. This is the reason figures and writing come
after the corrections rather than beside them. *One hour to run, then re-read
every number in CURRENT_FINDINGS.md.*""")

s=s.replace("""**A. Peak-level unit of analysis. Now the only research bet.**""",
"""**A. Peak-level unit of analysis. Now the only research bet.**

*Pre-registered before any of it is run:* peak-level analysis asks where contacts
land rather than what gene annotations they predict, so the length confound bites
differently, but it bites. **Any peak-level predictive claim gets tested against
peak-level equivalents of length, expression and local sequence composition
before it is believed.** Same rule as S0, one level down. No baseline in this
project has ever contained them, which is how a whole phase came to be built on
a length artefact.""")

s=s.replace("""## The rule for this phase

**One more test, then stop auditing and write.**""",
"""## Order of work

    S0 orthogonality  DONE
    S1 core promoters (the last test)
    S2.6 PROJECT_STATUS.md first, then the rest of S2, including BOTH feature
         removals (S2.1 degenerate topology, S2.1b the moment families)
    S3 figures
    S4 write

## The rule for this phase

**One more test, then stop auditing and write.**

Note what this rule does not cover. S2.1 and S2.1b are **fixes, not tests**:
they remove features that are provably degenerate or unmeasurable. Doing them is
not a continuation of the audit, and skipping them would leave known-bad features
in the substrate that the figures are drawn from.""")

s=s.replace("""2. **Aim 1, the tool.** MCCProfiler, 91 features from bp-resolution MCC.""",
"""2. **Methods, not a result: the instrument.** Aim 1 moves toward Methods. A
   description tool is justified by what it measures and how reliably, not by
   what it predicts, which is also what CellProfiler is: nobody justifies it by
   "it predicts cell type better than cell diameter". The report becomes two
   results rather than three, both geometric, both having survived everything.

   MCCProfiler, 91 features from bp-resolution MCC.""")
p.write_text(s); print("PLAN_FORWARD.md updated")
PY
grep -c "—" audit/PLAN_FORWARD.md
git add -A && git commit -q -F - <<'MSG'
Restore the second feature removal, which consolidation had dropped

Jude spotted it: the consolidated plan kept the degenerate topology drop but
lost the per-peak moment re-supporting, which survived only in the appendix as
the old Phase 1.2. It matters more than tidiness. Eight of the nineteen retained
components rest on oe_asymmetry and oe_tailedness at rho 0.27 to 0.42, carrying
14.9% of variance, and ChromHMM_bivalent takes its displacement on shape-PC4,
which is 15% trustworthy and asymmetry-dominated. Fixing or retiring those
families changes what one of the Aim 3 corroborating controls means.

Also recorded the consequence that both removals share and that nothing in the
plan had stated: changing the feature set moves the component count, the
variance shares, the trust weighting and every displacement number. Aim 2's "19
components, 78%" and Aim 3's effect sizes must be re-quoted from the rebuilt
store rather than carried over, which is the actual reason figures and writing
come after corrections.

Clarified that the stopping rule covers tests, not fixes. S2.1 and S2.1b remove
features that are provably degenerate or unmeasurable; skipping them would leave
known-bad features in the substrate the figures are drawn from.

Aim 1 moved toward Methods, as suggested. A description tool is justified by
what it measures and how reliably, not by what it predicts, which is what
CellProfiler is too. Two results instead of three, both geometric.

The orthogonality measurement is marked done, since it was run before the plan
asked for it: 2.6% of feature variance and 0.6% of the retained components.

And the peak-level phase now carries a pre-registered rule: any predictive claim
there gets tested against peak-level equivalents of length, expression and local
sequence composition first. That is how a phase came to be built on a length
artefact in the first place.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
MSG
git push -q origin feat/app-scaffold && echo pushed

OUT
PLAN_FORWARD.md updated
0
pushed
File updated: audit/PLAN_FORWARD.md.

On what you forgot: two feature removals, and only one survived consolidation.

S2.1 drop the four degenerate topology features — was there ✓
S2.1b re-support or retire the per-peak moment families — had fallen out, now restored
The second is not cosmetic. Eight of nineteen components rest on oe_asymmetry/oe_tailedness at ρ 0.27–0.42 (14.9% of variance), and ChromHMM_bivalent — an Aim 3 corroborating control — takes its displacement on shape-PC4, which is 15% trustworthy and asymmetry-dominated. Fixing those families changes what that control means. The diagnosis is support size: median peak ~11 bins at 1–2 reads, while contact_asymmetry computes the same statistic over the whole window at ρ 0.900. So the concept works and only the support is wrong.

On the plan itself — I agree with all four changes. Moving Aim 1 into Methods is right, and the CellProfiler analogy is the correct one: nobody justifies it by "it predicts cell type better than cell diameter". Two geometric results that survived everything beats three where one limps.

Three things I have added that the plan did not cover.

First, a consequence nobody stated: both feature removals change the numbers downstream. Component count, variance shares, trust weighting and every displacement value move. Aim 2's "19 components, 78%" and all of Aim 3's effect sizes must be re-quoted from the rebuilt store, not carried over. That is the real reason figures and writing come after corrections, and it was implicit rather than written down.

Second, a clarification to "now stop testing": the stopping rule covers tests, not fixes. S2.1 and S2.1b remove features that are provably degenerate or unmeasurable. Skipping them because of a testing moratorium would leave known-bad features in the substrate the figures are drawn from.

Third, the pre-registered rule for Phase A is now written into the plan rather than left as advice: any peak-level predictive claim gets tested against peak-level equivalents of length, expression and local sequence composition before it is believed. That is precisely how a phase came to be built on a length artefact.

Order of work is now explicit at the top: S0 done → S1 → S2.6 first, then the rest of S2 including both removals → figures → write.