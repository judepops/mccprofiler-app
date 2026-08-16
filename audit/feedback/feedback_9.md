The left panel works. The right one is doing no work, and there's a better use for that space.

The right panel is a restatement, not a second view. For two normal distributions, overlap is a deterministic function of |d| — roughly 2Φ(−|d|/2). So the scatter is a curve by construction, and every point's y-position is derivable from its x-position. A reader who understands that will discount it; one who doesn't will think you've shown two independent things. Either way it costs you half the figure.

Replace it with the thing the claim is actually about: overlapping distributions. Pick three exemplars — Lambert TF (largest real separation), dbSUPER SE (the null), and one mid-pack set — and plot each set's density against the panel density on that set's strongest axis. Three small stacked panels. The reader sees two curves lying on top of each other, which is far more persuasive than a number saying 82% overlap. That's the whole argument in one image, and it's genuinely additional information rather than a transform of the left panel.

On the left panel, three fixes:

Colour by definition type (DIRECT vs OUTPUT) rather than singling out SE in orange. Right now the highlight implies SE is special, which is exactly the claim you just retracted. Colouring by category makes the taxonomy visible and lets SE sit mid-pack among the OUTPUT sets, which is the honest position.
Mark significance. Several of the small-|d| bars are non-significant and there's no way to tell. An asterisk or an open-vs-filled bar.
Show n. d = −0.73 at n=41 and d = +0.20 at n=733 are not equally reliable, and a reader can't see that. A small right-hand column of n values is enough.

Minor: the SE label in the right panel is truncated ("pm5") and floating away from its point. And the green gene-desert bar needs a legend note saying it's the positive control for genomic context, not a finding — otherwise it reads as the second-biggest result.

what do you tuink of the actual left panel distributions though

Thought process

Looking at the actual values, there's a problem more serious than the layout.

The signs aren't comparable across sets, so the sorting is misleading.

Each set's d is measured on its own strongest component — Lambert TF on shape-PC1, Roadmap silenced on shape-PC3, SE on shape-PC12. Principal component signs are arbitrary (sklearn's convention is essentially a coin flip per component), and different components aren't a shared scale anyway.

So the chart reads as a spectrum from +0.46 down to −0.73, implying the top sets are opposite in character to the bottom ones. They aren't. They're 21 unrelated magnitudes on 21 different axes with arbitrary sign. "Lambert TF is positive and Roadmap silenced is negative" means nothing.

Fix: plot |d| only, sorted descending. If you want direction, it has to be on a fixed axis with a stated sign convention (e.g. all sets projected onto shape-PC1), and then the comparison is meaningful.

On the biology, using magnitudes:

There's an internal tension worth noticing. phastCons 2kb topQ at 0.15 and gnomAD pLI topQ at 0.17 barely displace — yet LOEUF was your single strongest prediction target in Aim 1 (r = 0.210, 4.08 SD margin, 91% surviving amount correction).

That's not a contradiction, it's the continuum story stated twice: contact architecture predicts constraint well as a continuous quantity, but the top-quartile constraint genes don't occupy a region. Same data, both framings, and it's a nice piece of internal consistency between Aims 1 and 3. Worth a sentence — it makes the "directions not regions" claim concrete rather than abstract.

CpG island promoter at 0.19 is the one I'd double-check. CpG-island promoters are one of the most architecturally distinct classes in the literature, and a near-null there is either a real and interesting finding or a sign the set is badly defined in your panel (it's probably very large — check n). If most of the panel is CpG-island-promoter genes, the "set vs rest" contrast is weak by construction.