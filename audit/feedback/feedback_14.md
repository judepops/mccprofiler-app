Finished, and it answers the question cleanly — no combination rescues separation, and searching harder actively makes things worse.

search level	real	null	gap
all 64 features, no selection	0.576	0.492	+0.084
best of 120 random subsets	0.591	0.545	+0.046
greedy forward selection	0.625	0.581	+0.044
Read the null column first. On labels that mean nothing, greedy selection climbs from 0.492 to 0.581 — the search manufactures +0.089 of AUC out of pure noise, which is twice the +0.044 that separates real labels from their own null at that level. Without the permutation you would have reported "optimal feature selection reaches AUC 0.625" and almost all of it would have been the search.

The direction is the important part: the gap is largest with no selection at all and shrinks as you search harder. Searching raises the real ceiling and the null ceiling together, but raises the null faster. So feature selection here degrades the signal-to-artefact ratio rather than improving it — the honest all-features baseline is the best discriminator you have, not the worst.

Super-enhancers: greedy 0.630 against a null of 0.620, gap +0.010. That is a considerably stronger result than the displacement test. It is no longer "we tested SE and found nothing" but "we let an adversarial search over 64 features pick the best possible combination specifically to separate SE, and it beat shuffled labels by one AUC point." That closes the obvious objection that you simply used the wrong features.

Two sets clear +0.10, and neither is a surprise:

gene_desert_bottomQ_density (+0.218, AUC 0.811) — the known positive control. Note it scores 0.806 with no selection and 0.811 with full greedy search: it needs no searching because it is genuinely strong. It is also the set that collapses under density stratification, so it is measuring gene density, which contact architecture can genuinely see.
ChromHMM_bivalent (+0.107, n=46) — ChIP-derived, so it is partly chromatin predicting chromatin, and at n=46 its own null is already 0.683, which is where small-set overfitting lives.
One thing to be careful of when reading the per-set rows: Roadmap_silenced looks impressive at greedy 0.732, but its null is 0.716, so the gap is +0.016. The raw column is the trap the null column exists to catch.

Writing it into the findings doc:

Bash Show findings doc structure
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit
grep -nE "^#{1,3} " CURRENT_FINDINGS.md | head -40

OUT
1:# Current findings
14:## 0b. Feature health audit, all 73 (`audit_feature_health.py`)
20:### What was removed, and why it was right
36:### What the audit found that we had not looked for
71:### Should we be more careful: yes, and here is the specific answer
95:## 0. What the rebuild changed, 2026-08-16
100:### The substrate got substantially cleaner
114:### Aim 1 holds, essentially unchanged
120:### Aim 2 holds
126:### Aim 3 changed, and the super-enhancer claim must be restated
162:### What must be redone
183:## 1. The question
193:## 2. Aim 1: the feature substrate is worth building, narrowly
227:### Aim 1, restated: the features measure something gene properties do not
272:### RETRACTED 2026-08-16: the predictive claim does not survive gene length
316:### Aim 3 survives gene length: six confounders now controlled
334:### The former Aim 1 headline, retained for context only
358:## 3. Aim 2: the landscape is continuous, not categorical
409:## 4. Aim 3: categories are directions, not regions
454:### The headline claim
478:### Chromosome is a confound, and it costs us the Lambert_TF control
551:### Why super-enhancers are the case of interest
572:### Supporting observation: the kind of definition predicts displacement
639:### The axis rank is a stronger statement of the null than the effect size
665:### Three independent positive controls, on three different well-measured axes
685:### Which positive control to qualify
728:### The full table, corrected substrate
743:### Three things that must travel with this table
761:## 5. Supporting results
786:## 6. What the architecture groups are called, and the caveat
809:## 6b. NEW 2026-08-16: contact asymmetry is strand-dependent, and the feature destroys it
854:## 6c. S1 CLOSED 2026-08-16: the housekeeping question cannot be asked on this panel
912:## 7. Retracted, do not quote
938:## 8. Where to go next
949:## Reproducing every number here
Read CURRENT_FINDINGS.md (lines 640-669)
Edit CURRENT_FINDINGS.md
Added 51 lines
Bash Inspect scripts folder and runner
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit/scripts
ls
echo "--- run_all.sh tail ---"
tail -25 run_all.sh 2>/dev/null

OUT
audit_feature_health.py
build_core_promoter_sets.py
diagnose_adhoc_checks.py
diagnose_chromosome_confound.py
diagnose_definition_type.py
diagnose_dimension_trust.py
diagnose_external_axis_rank.py
diagnose_external_density_stratified.py
diagnose_external_power.py
diagnose_external_structure_corrected.py
diagnose_external_structure.py
diagnose_feature_novelty.py
diagnose_loeuf_confounds.py
diagnose_pc_names.py
experiment_cluster_search.py
experiment_feature_search.py
experiment_rotate_axes.py
make_audit_figures.py
output
output_S2
README.md
run_all.sh
--- run_all.sh tail ---

echo; echo "--- 2. Aim 2: is the structure categorical? ---"
run experiment_cluster_search       "16 conditions, four cluster tests each"
run experiment_rotate_axes          "would rotation make the axes nameable"

echo; echo "--- 3. Aim 3: external sets, in confound order ---"
run diagnose_external_structure_corrected  "displacement, magnitude-corrected"
run diagnose_external_density_stratified   "+ gene density"
run diagnose_chromosome_confound           "+ chromosome"
run diagnose_external_power                "minimum detectable effect"
run diagnose_external_axis_rank            "where each set leans"
run diagnose_definition_type               "does definition type predict displacement"

echo; echo "--- 4. predictive claims, which need the length control ---"
run diagnose_loeuf_confounds        "length, expression, CpG vs the shape substrate"

echo; echo "--- 5. checks that were once ad hoc ---"
run diagnose_adhoc_checks           "coverage, taxonomy robustness, DBD split, chr19"

echo; echo "--- 6. figures ---"
run make_audit_figures              "six composites into audit/figures/"

echo; echo "done. outputs in $OUT"
echo "Now re-read audit/CURRENT_FINDINGS.md against these: every component count,"
echo "variance share and effect size in it will have moved."
Edit run_all.sh
Added 2 lines
Bash Find where to document the new script
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit/scripts
grep -nE "axis_rank|definition_type|feature_health" README.md | head

OUT
56:| `diagnose_external_axis_rank` | which component each set leans on |
57:| `diagnose_definition_type` | does the KIND of definition predict displacement |
Edit README.md
Added 1 line
Recorded in CURRENT_FINDINGS.md (new subsection in §4, ahead of the positive controls), added to run_all.sh and the scripts README.

Where that leaves the two things you asked for:

Feature search — done, and it strengthens the argument rather than weakening it. You asked which combination best separates the external sets; the answer is that the best combination is all of them, because every selection step buys less real signal than artefact. Worth keeping the standing caveat in mind: no effect size from this experiment can go in the report, since every subset was chosen on the outcome it was scored against. Its admissible use is purely negative — the search was run and did not rescue anything.
Figures — no work needed. All six were already on the 73-feature substrate; re-rendering produced byte-identical files, and fig4 confirms it (PC1 at 17.2%, 14 components above the ceiling). That item comes off PLAN_FORWARD.md.
Still outstanding from the plan:

S2.2	re-derive the archetype display names on the corrected substrate — the CLAUDE.md caveat that "dispersed" may partly mean "has more signal" is still unresolved
S2.4	app enrichment grid still runs on raw components, and MAX_GROUP_COVERAGE is not yet applied in the app
S2.5	reconcile 116 vs 119 twice-captured genes
—	re-check CURRENT_FINDINGS.md §2–6 against output_S2/ (§4 partly done today)
S2.2 is the one with thesis consequences, since those names travel into the writeup and the briefing pack still contradicts them. I'd take that next unless you want the app fixes first.