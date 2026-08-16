The reframe is better and the retraction was right. Two problems with the new taxonomy, both fixable, and one of them a reviewer will find immediately.

1. The categories are post-hoc, built from 7 of 21 sets.

"Mechanism / consequence / chromatin state" vs "output / threshold" was chosen after seeing which sets separated. With seven sets and two bins you can nearly always find a scheme that fits. As written it's a just-so story.

The fix is straightforward and makes it a real result: classify all 21 sets by definition type before looking at any effect size, write the assignment down, then test whether the categories differ (Mann-Whitney on d across categories, or just report the two distributions). If the pattern holds across all 21 with the assignment fixed in advance, it's a finding. If it only holds for the 7 you happened to tabulate, it isn't.

That's an hour and it's the difference between a claim and an observation.

2. The chromatin-state / signal-threshold boundary doesn't survive scrutiny.

dbSUPER super-enhancers are a rank-and-threshold on H3K27ac. ChromHMM bivalent is a threshold on H3K4me3 and H3K27me3. Roadmap silenced is a threshold on chromatin marks. All three are thresholds on ChIP signal — so why do two count as "chromatin state" and one as "threshold on signal"?

You need a principled distinction or the taxonomy collapses. The defensible one is combinatorial multi-mark state vs rank-cutoff on a single mark: ChromHMM and Roadmap use hidden Markov models over many marks to assign states; ROSE ranks regions by one mark and cuts at an inflection point. That's a real methodological difference and it's exactly the Pott & Lieb critique. But it has to be stated explicitly, because right now the categories look drawn to fit.

On the SE claim specifically — it survives the reframe and is still the interesting case, for a reason worth stating: nobody claims Eisenberg housekeeping genes constitute a distinct 3D architectural class. People do claim that about super-enhancers. So a null for SE contradicts a live claim in the literature; a null for Eisenberg contradicts nothing. That asymmetry is why SE stays the headline even though its effect size is marginally larger.

And the figure-as-control lesson is right. Drawing figure 6 falsified the arch-HK claim that four written review rounds missed. Move figure production before writing in the plan, not after.