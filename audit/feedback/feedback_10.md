All four landed cleanly, and the panel-coverage finding extends the point further than I did.

One symmetric fix worth adding: you have a minimum-n floor (25 genes) but no maximum-coverage rule. A set covering 88% of the panel can't produce an interpretable set-vs-rest contrast at any effect size — the comparison group is a 227-gene remainder defined only by not being in a very large set. Those rows shouldn't be read as evidence either way.

Adding a coverage ceiling (say, flag or exclude above 70%) mirrors the existing convention and stops a reader treating ChromHMM_active_TSS at |d| = 0.29 as comparable to Lambert TF at 0.46. It also pre-empts someone asking why a set covering most of your panel is in a "displaced vs not" table at all.

Worth noting for reassurance: this doesn't touch the headline. SE is 8.6% of the panel and Lambert TF 8.4% — both well inside any sensible bound, and the matched-n comparison is unaffected. And you've already shown the taxonomy survives dropping all three oversized sets (0.0019 vs 0.0018).

The density-overlay panel is the right call. That figure is now the single most persuasive thing in the report — a reader sees the SE curve sitting on the panel curve and doesn't need the statistics explained.