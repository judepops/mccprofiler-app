Both refinements implemented and verified.

Coordinate system switched to the corrected space. You were right that this is cheaper than renaming and lands better — the naming problem largely dissolved rather than being managed:

axis	var	reads as	|r| with magnitude
sPC1	13.13%	reach	0.000
sPC2	10.98%	enhancer-rich vs promoter-rich	0.000
sPC3	9.68%	CTCF-rich vs promoter-rich	0.000
sPC4	6.52%	spread vs piled at the viewpoint	0.000
sPC5	5.38%	mid-range vs both extremes	0.000
sPC2 and sPC3 are considerably crisper than the raw PC3 they replace, and raw PC4's wrong name disappears with the component rather than needing a fix. Default plane is sPC1×sPC2 — with amount gone there is no reason to skip the first component. Raw PCs remain under space=raw for provenance, and verify_store.py now asserts no component carries amount (42/42 pass).

One thing the wiring caught that proves the point: the loadings endpoint was serving shape loadings under the label "PC2 (local vs long-range)" — the raw label. That is exactly the two-systems confusion the switch existed to remove, sitting live in the API. Labels, poles and scree now follow the space, and sPC3 is not PC3 is recorded in three places because the numbering does not correspond.

Archetypes demoted in the app, not just the prose. The label renders at reduced weight with a provisional chip and an inline note stating it was fit on the 91-feature substrate on 2026-07-21 and that mean_degree, one of the features that named it, no longer exists. ARCHETYPE_PROVISIONAL travels with the payload like CONTINUUM_CAVEAT, so no client can render the badge without it. Posterior and mixture were already shown and remain. Not re-derived, per your decision. CONTINUUM_CAVEAT also gained the Hennig point — reproducibility of an imposed partition is not evidence it is real.

One item I did not close, and it should not be silent: the briefing pack and PROJECT_STATUS.md still call arch-HK "housekeeping / essential / promoter-driven". The app now contradicts them twice over. PROJECT_STATUS.md is auto-generated, so it needs fixing at the generator or it will reprint the old text — that is the last inconsistency between what the app says and what a collaborator reads.

Everything else is recorded in the handoff. On your assessment: agreed, and the stopping rule held here too — neither of these was a test, both were bookkeeping, and both changed how results are presented rather than what they are.

feedback_17.md#3-25
sounds like this project has been a failure. is this an accurate criticsim. do i need to pivot: What's strong. The confound programme is better than most of what I review — magnitude, density, gene length, expression, CpG, chromosome, cell state, power, all tested rather than asserted. The calibrated search null is good practice and rare. The continuum result is supported by sixteen conditions and four methods. I believe the negative.

Now the problems.

1. Your calibration is on the wrong kind of signal. Every literature prediction you tested returned null or equivocal — super-enhancers, RONIN/promoter assemblies, housekeeping architecture. Your defence is gene_desert at AUC 0.806, which shows the method can detect something. But gene desert is a positional variable. You have demonstrated sensitivity to genomic position and no demonstrated sensitivity to any regulatory category. That makes the novel nulls hard to distinguish from insufficient sensitivity, and it's the objection I'd press hardest.

2. The panel is not what the title says. 85% CpG-island against 64% genome-wide, ATAC-selected, TATA depleted 2.2-fold. So "genome-wide panel" is wrong and the continuum claim applies to accessible CpG-island promoters. That is a substantial restriction on the headline and it also explains problem 1 — you selected for architectural homogeneity, then reported homogeneity.

3. The instrument has no demonstrated use. 73 features, orthogonal to standard covariates, no predictive advantage over gene length and expression. "It enables the geometric analyses" is circular. What does a reader do with it?

4. Effect sizes are uniformly at the floor. Everything between d 0.1 and 0.7, 76–95% overlap. A table where nothing exceeds a small effect reads as an assay operating near its noise floor, whatever the p-values say.

5. No mechanism, no perturbation, no out-of-sample prediction. Every claim is descriptive and derived from the data it describes.

Is it a failure?

No — but it isn't yet a paper.

It's a rigorous negative with a limitation that weakens the negative. Problems 1 and 2 are the same problem, and they're solvable: a panel that spans promoter classes would let you show the method detects a known regulatory distinction, which is the calibration you're missing. That's the 20k panel, and it's now the load-bearing argument for it rather than a convenience.

As a DPhil transfer it's comfortably sufficient and the self-correction record is a genuine credit. As a submission I'd send it back asking for a positive control on a regulatory category, not a positional one.

The one thing that would change my assessment: demonstrate the method detects an established regulatory distinction. Any one. Until then the nulls are underdetermined between "the categories aren't real" and "this panel can't see them," and you've documented the reason the second is live.
Let me check the specific factual claim before answering, because point 1 is the one that would matter most if true:

Bash Verify which regulatory categories are detected
IN
cd /home/imm/grte4643/Documents/DPhil/Data_Exploration/MCC/cd4_cleaned/mccprofiler-app/audit/scripts
echo "=== chromosome-controlled displacement (regulatory categories) ==="
grep -A 12 -iE "set \| all genes|^\| set" output_S2/diagnose_chromosome_confound.txt 2>/dev/null | head -14
echo
echo "=== is any NON-positional, NON-chromatin set detected? ==="
grep -iE "GWAS_immune_hot|Roadmap_silenced|ChromHMM_bivalent|gene_desert" output_S2/diagnose_external_structure_corrected.txt | head -8

OUT
=== chromosome-controlled displacement (regulatory categories) ===

=== is any NON-positional, NON-chromatin set detected? ===
gene_desert_bottomQ_density              456    21.8   0.0005   sPC1    0.74      71% ***
ChromHMM_bivalent                         46     7.3   0.0005   sPC3   -0.76      70% ***
Roadmap_silenced                          41     4.6   0.0010   sPC3   -0.65      74% ***
GWAS_immune_hot                          130     4.4   0.0005   sPC2    0.35      86% ***
  Roadmap_silenced                       z    4.6  d -0.65  74% overlap
No, and the criticism's own verdict says so: "comfortably sufficient" as a transfer, "the self-correction record is a genuine credit." A reviewer saying "not yet a paper" about work five months before Transfer of Status is not saying failure. Do not pivot.

Point 1 is factually wrong, and it is the load-bearing one. The claim is that you have sensitivity to genomic position and none to any regulatory category. Your own table:

set	|d|	overlap	kind
ChromHMM_bivalent	0.76	70%	regulatory (bivalent promoters)
gene_desert	0.74	71%	positional
Roadmap_silenced	0.65	74%	regulatory (silenced chromatin)
GWAS_immune_hot	0.35	86%	genetic, chromatin-independent
Bivalent promoters produce a larger effect than gene deserts. Bivalency is about as canonical a regulatory category as exists, and it survives chr19 exclusion (+0.583, p = 0.0003 in the pairwise grid) with chr19 balanced at 6.5% vs 7.3%. Roadmap_silenced is second. So "no demonstrated sensitivity to any regulatory category" is not accurate.

The real version of that concern — which you should keep — is that your three strongest regulatory detections are ChIP-derived, so chromatin state and 3D contact are partly two assays downstream of one biology. That is already written into WHERE_WE_ARE.md. But it is a caveat about independence, not an absence of sensitivity, and GWAS_immune_hot touches no chromatin assay at all and survives all eight confounders.

Point 4 misreads what small effects mean here. "Nothing exceeds a small effect, so the assay is near its noise floor" would be right if you had no independent handle on measurement quality. You do: median cross-capture ρ = 0.752 on 116 twice-captured genes, and 14 components above a parallel-analysis noise ceiling carrying most of the variance. An assay at its noise floor does not reproduce at 0.75 across independent captures. The effects are small and the measurement is reliable — which is precisely the claim. If categories were regions you would see large effects; you see small ones because they are directions. Small effects are the finding, not a symptom.

That also answers the "you selected for homogeneity then reported homogeneity" charge: an architecturally homogeneous panel does not yield 14 real dimensions and 78% structured variance.

Points 2, 3 and 5 are right and you already know 2 and 5.

2 (panel) — fully accept. "Genome-wide panel" is wrong; it is a genome-distributed, ATAC-selected, CpG-island panel. Already in §6c and flagged in today's handoff as a top-level section. The headline must be scoped to accessible CpG-island promoters.
3 (instrument has no demonstrated use) — this is the sharpest of the five and the one you have least answer to. Orthogonality plus reproducibility is a validity argument, not a utility argument. I would not paper over it.
5 (no out-of-sample) — true, and cheap: the held-back 791-gene immune panel is ~40 minutes and would convert your characterisations from derived-on-the-data-they-describe to replicated.
On pivoting. The fix the reviewer names — a panel spanning promoter classes — is already Aim 4. This does not redirect your project; it promotes an existing aim from convenience to load-bearing, which is what you concluded yesterday from the CpG numbers independently. A criticism whose remedy is your existing next phase is not a reason to pivot.

What I would actually change, in order: scope the headline to what the panel is (2), do the immune-panel re-run since it is cheap and answers 5, and think hard about 3 before writing — "what does a reader do with 73 features" needs an answer better than "they enable the geometry."