Strong document — the diagnostics are real and it self-corrects (retracting varimax on reproducibility grounds is exactly right). But I tested two of its load-bearing claims and its top recommendation targets the wrong cause, and one headline is a semantic sleight.

Three things it gets right
The mean_degree degeneracy (§3.1.3) is a genuine bug and the most valuable find in the document. If adjacency is the complete graph, mean_degree = n_active(n_active−1)/n_peaks is algebra on two counts, and it's PC1's second-strongest correlate at 0.763. That's checkable and it should be fixed.

"Displaced but not separated" (§2.5) is the best articulation of the continuum result anyone has produced here. 21/21 sets displaced at p ≤ 0.0045, largest effect d = 1.61 still leaving 42% overlap. And super-enhancers weakest of the 21 at 88% overlap is the sharpest form of the thesis argument. Flagging the density positive control (z = 26.7) as the largest displacement, i.e. a known confound, is honest and easy to have omitted.

arch-HK is not the housekeeping group (38.2% vs arch-ME-constitutive's 40.2%, and the lowest median blood expression of the three active groups). Names that contradict measurements should go.

Where I'd push back
1. "PC1 is not amount" rests on a broken proxy
total_mcc has CV 0.216 — median 4,737, IQR 4,063–5,477. It barely varies, because the bigwigs are depth-normalised to ~20,000 per viewpoint. So total_mcc measures what fraction of normalised signal stayed inside ±1 Mb, not contact abundance. Testing PC1 against it is testing against a normalisation residue, and r = +0.033 is uninformative.

This also explains the review's own puzzle: r(total_mcc, n_high_consensus) = −0.16 isn't a deep insight about two kinds of amount, it's what you get when one variable is nearly constant by construction.

Consequence for §4.3: "regress total_mcc out before PCA" would remove a normalisation artifact, not amount, and could be actively harmful. The canonical amount basis in this project is the 11-feature MAG_OVERALL; raw pileup totals per viewpoint are also available (median 18,580, q05–q95 3,772–38,719) and are the real depth measure.

And PC1 is an amount axis — a count-amount axis. The review's own text says it "tracks peak counts and reach." "PC1 is not amount" only works if amount means total signal.

2. §4.1 — orientation is ruled out. It's sparsity.
The two hypotheses make opposite predictions: an orientation bug flips signs between panels, so |value| would reproduce much better than signed value. Sparsity produces variance, so neither reproduces. On the 119 shared genes:

feature	ρ(signed)	ρ(|value|)
oe_asymmetry_mean_all	0.276	0.151
oe_asymmetry_mean_ctcf	0.376	0.145
oe_asymmetry_mean_enhancer	0.420	0.208
contact_asymmetry (whole profile)	0.900	0.803
frac_far_distal (reference)	0.978	—
|value| is consistently worse, never better. Orientation is not the cause.

The real cause is visible in the last two rows: contact_asymmetry, the same concept computed over the whole ±1 Mb profile, reproduces at 0.900. The per-peak version sits at 0.28–0.42. Same statistic, different support — median peak is 11 bins, 41.7% span fewer than 10, 12.3% fewer than 6, and the median non-zero 50 bp bin holds 1–2 reads. A third or fourth central moment estimated from six bins carrying 0–2 counts has enormous sampling variance no matter how it's oriented.

So the fix isn't "check orientation" — it's that high moments are unmeasurable at peak scale at this depth. Either drop the per-peak moment families, or compute them on aggregated support (stacked peaks per gene, or the whole profile, where contact_asymmetry shows they work fine). That still repairs the eight weak components; it's just a different intervention, and it doesn't cost a week chasing a sign convention.

3. "Most of the signal is amount, not shape" over-reads a residualisation
corrected_shape residualises every feature against magnitude and drops the magnitude features. If amount is a consequence of architecture — more elements contacted → more contacts → more signal — then residualising it removes real architecture along with the confound. The drops in that column are consistent with over-correction, not only with "shape carries little."

The defensible statement is "amount is sufficient for most targets", not "shape carries little". The review also omits the 11-feature magnitude basis row, which is the canonical amount baseline here — 91 features beat it 9/9 (p = 7.5e-06), which is the fair version of the comparison.

What's missing
The review predates the resolution work and it matters for its own conclusions:

Contact summits reproduce at 14 bp median across the two panels (64× tighter than a jitter null), and CTCF-class summits sit 9 bp from the CTCF motif centre, 4.48× a centrality-matched null, with an enhancer negative control at 1.01×. That is the concrete, measured version of §3.1.1's "position-invariant by design, so the actual shape is discarded" — and it's the strongest argument for §4.6 (peak-level GNN), which the review ranks third.
The resolution ladder says the 91 features use nothing below ~1 kb. So the 400-fold compression is real but its cost is now bounded, not hypothetical.
amount_vs_shape_partition was corrected on 2026-08-03 (its old "shape" view kept 96% of features). If any app number derives from it, check the vintage.
Minor: I get 119 shared genes, not 116 — probably a QC-stage difference, worth reconciling.

Revised ranking
4.2 mean_degree — unchanged, real bug, tier 1.
4.1 restated — not orientation; retire or re-support the per-peak moment families. Cheaper than proposed and now diagnosed.
4.3 amended — amount-correct against MAG_OVERALL, not total_mcc.
4.6 promoted — peak-level unit of analysis is the one route with a measured asset behind it (14 bp summits, motif localisation) rather than an argument from compression.
4.7 stays the decisive clustering test; 4.8 remains the strongest honest framing.
The one line I'd change in the summary: "most of its predictive power comes from how much signal a gene has" — say how many elements it contacts, since total signal is nearly constant by construction here. Different claim, and it's the one the data supports.