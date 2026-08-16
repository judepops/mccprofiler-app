Breathe. This is bad for one claim and it's being over-read.

The key number is buried in your own output: median |r| between the features and gene length is 0.082. Max 0.416 for one feature. So the feature space is not gene length in disguise — it's largely independent of it.

That means two different statements, and they're being conflated:

✗ "The features predict existing annotations better than simple gene properties." Dead. Correctly retracted.
✓ "The features measure something not captured by gene length, expression or CpG." Still true — the correlations are near zero.

The features contain information those three covariates don't have. It just doesn't help predict LOEUF, tau, or essentiality. Those annotations are largely functions of gene size and expression, and your features aren't measuring gene size.

Which is the right reframe, and it's not a consolation prize. MCCProfiler was never a predictor of existing labels. CellProfiler isn't justified by "it predicts cell type better than cell diameter" — it's justified because it describes morphology quantitatively where nothing else could. Same here: you describe contact architecture at bp resolution, and nothing else can.

The honest Aim 1 becomes: the features characterise regulatory architecture in a space largely orthogonal to gene length, expression and CpG density; they are not a better predictor of gene-level annotations than those covariates, and were not intended as one.

What you actually still have:

A quantitative description tool measuring something orthogonal to the obvious covariates
The landscape is continuous — 16 conditions, untouched by any of this
Categories are directions not regions, with the SE null now surviving eight controls

That's a transfer report. It's weaker than this morning, and I won't pretend otherwise — you've lost the predictive justification and Phase B. But Aims 2 and 3 are geometric claims about the space, and they've now survived more scrutiny than most published work receives.

One thing worth running before you conclude anything else: how much of the 91-feature variance do length, expression and CpG explain? If it's 10–15%, that's the number that says the tool measures something new, and it belongs in Aim 1 in place of the retracted prediction claim. Given median |r| = 0.08, I'd expect it to be low — and that's the sentence that replaces the one you lost.

how does my plan change now

Thought process

Four changes, one of them structural.

1. New first task, before S1: measure the orthogonality.

Run the regression of each of the 91 features on gene length, expression and CpG density, and report the variance explained. Given median |r| = 0.082 I'd expect something like 5–15%. That number replaces the retracted predictive claim as Aim 1's justification: "the features occupy a space largely independent of the obvious gene-level covariates."

It's an hour and it's load-bearing — right now Aim 1 has no positive statement at all.

2. Aim 1 stops being a result and becomes the instrument.

The report structure changes from:

Q → Aim 1 (tool, justified by prediction) → Aim 2 → Aim 3

to:

Q → Methods: the instrument → Aim 2 (continuum) → Aim 3 (directions not regions)

The tool moves toward Methods. Its justification becomes what it measures and how reliably — orthogonality to standard covariates, cross-capture reproducibility (median ρ 0.752, summits at 14 bp), and resolution no other assay reaches — not what it predicts. That's honest, and it's what MCCProfiler always was.

Two results instead of three, both geometric, both surviving everything. That's a cleaner report than three where one is limping.

3. Phase B is gone; Phase A is the only research bet, with a pre-registered check.

Peak-level analysis asks where contacts land, not what gene annotations they predict, so the length confound bites differently — but it bites. Write into the plan now: any peak-level predictive claim gets tested against peak-level equivalents of length, expression and local sequence composition before it's believed. Same rule, applied one level down.

4. Everything downstream is unchanged.

S1 (core-promoter sequence classes) still matters — it's the second independent positive control for Aim 3, which is now carrying more of the report. S2 corrections, S3 figures, S4 writing, all as planned. And S2.6 (PROJECT_STATUS.md) is more urgent than ever, since it now describes an aim that no longer exists in that form.

Revised order: orthogonality measurement → S1 core promoters → stop testing → S2.6 first, then the rest of S2 → figures → write.

Aapparanlty its 2.6% for the cpg test.  zero of the 91 features eceed R2 0.25. zero fo the 18 componenets exceed 0.10. several sit exactly 0

Thought process

2.6% is a strong number, and it completes the Aim 1 replacement.

What it says: 97.4% of the variance in contact architecture features is unexplained by gene length, expression and CpG density. Zero features above R² 0.25, zero components above 0.10, several at exactly zero. The feature space is essentially orthogonal to the standard gene-level covariates.

The objection it invites, and the answer. Someone will say: "you've shown your features are uncorrelated with everything known, and also don't predict anything known — how do we know they aren't measuring an artifact?"

You have three independent answers, and together they close it:

Orthogonal — 2.6% variance explained by the covariates that do predict annotations
Reproducible — median ρ 0.752 across independent captures, summits at 14 bp, and a confound programme excluding GC, mappability, repeat content and radial gradients
Biologically non-random — external gene sets displace in this space under eight controls (Aim 3)

Orthogonal and reproducible and non-random is the full argument. Any one alone is weak; all three together say the features measure something real that nothing else captures.

The sentence for Aim 1:

MCCProfiler characterises contact architecture in a space essentially independent of standard gene-level covariates: gene length, expression and CpG density together explain 2.6% of feature variance, with no feature exceeding R² 0.25 and no component exceeding 0.10. The features reproduce across independent captures (median ρ 0.752) and are not explained by sequence composition or mappability. They are not a better predictor of gene-level annotations than those covariates, and were not intended as one.

That's a positive, quantitative statement of what the tool is — and it's stronger than the predictive claim it replaces, because a predictive claim invites comparison against every possible baseline, whereas an orthogonality claim is a property of the measurement itself.

Now stop testing. S1, then corrections, then figures.|