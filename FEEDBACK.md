The plan is strong and the §1b pivot (axes as the product, archetypes demoted) is the right call. But comparing it against what we hypothesised earlier, three components are missing — and two of them map directly onto your headline results.

What we originally sketched: gene lookup, gene comparison, ranked lists/cohort query, and export for programmatic use. Only the first is in the plan.

Missing 1 — Cohort view (the important one). Your single strongest result is that functional gene classes have distinguishable, position-independent architecture (6 external groups, 6–36× above random floor, validated positive control). Nothing in the plan demonstrates it. A view where you paste a gene list — or pick a built-in set (Eisenberg HK, DepMap essential, GWAS-immune, Lambert TF) — and see that set's distribution on each named axis against the panel, is a live demo of your best finding. It's also the thing a collaborator would actually use: "here are my 40 hits, where do they sit?"

This should be P3 or P4, not an afterthought.

Missing 2 — Reproducibility demo. You have 116 genes measured twice in independent experiments. Showing both measurements for the same gene — two profiles overlaid, two positions on the continuum map — makes κ=0.72 and the technical noise floor (trap #7) visceral rather than a number in a caption. It's a small view and it pre-empts the obvious "how do I know this is real" question. Cheap, high value.

Missing 3 — Export. The CellProfiler analogy only works if the output is portable. A download button giving the gene × 91 features + gcPCA coordinates + axis percentiles + archetype proximity, for one gene, a gene list, or the whole panel. Minimum viable "this is a tool, not a demo."

Missing 4 — Gene comparison. Two genes side by side: profiles, axis positions, distance in feature space. Lower priority than the cohort view but it's the natural follow-up question after a lookup.

One thing to cut or defer: igv.js genomic context (P6). It's the view UCSC/IGV already do well, it's the only one needing bigWig range serving, and it's secondary by your own §5. Push it behind the cohort view.

Revised phase order I'd suggest: P0–P2 unchanged, then dimension profile → cohort view → continuum map → export → radar/P(s) → reproducibility demo → explanation layer → igv.js → SE panel.

The therapeutic framing from earlier — enhancer-driven vs promoter-driven implying different CRISPR strategies and redundancy expectations — should live in the explanation layer attached to PC3. That's the axis that carries it, and it's the sentence that makes a clinician care.