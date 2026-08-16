Files updated this run: audit/CURRENT_FINDINGS.md and audit/PLAN_FORWARD.md. No new script — the robustness checks ran ad hoc, and the table is now in the document.

To answer your question directly: three of its four points are correct, and the fourth is right in its recommendation but wrong in its reasoning.

Correct, and adopted:

The section ordering. Frame → super-enhancers → taxonomy as support. This is now written into both documents, and §4 of the findings has been physically reordered.
Making the 38th-percentile number prominent. Volunteering that super-enhancers sit mid-pack is what makes the surviving claim credible rather than less so.
Treating Pott and Lieb's multi-mark-versus-single-mark distinction as motivation for why super-enhancers were worth examining, not as something the data demonstrates. RANK_SINGLE has n = 1 and that is a hard limit.
Right recommendation, wrong reason: "don't make the taxonomy load-bearing."

It quotes p = 0.030, which is the omnibus Kruskal-Wallis across five categories. The test matching the hypothesis is the directional DIRECT-vs-OUTPUT contrast at ~0.002. And the fragility it asserts does not survive testing — worst-case leave-one-out 0.003, all four arguable sets reassigned simultaneously gives 0.0005, and only 17 of 20,000 random relabellings reach the observed value.

So I demoted the taxonomy anyway, because leading a section with its strongest argument is good writing and the super-enhancer claim is the one that contradicts something the field currently believes. But I put the robustness table in the document rather than accept the premise, because "we invented a scheme and got p = 0.03" is a fair thing for an examiner to probe, and the answer to that probe should be four robustness checks rather than a shrug.

One thing worth noting for when you write: the fact that moving the arguable sets improves the separation is itself the best defence against "you drew the categories to fit." It is worth one sentence.