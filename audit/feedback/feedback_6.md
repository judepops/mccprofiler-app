Both documents are in good shape, and the power analysis inverting the small-n worry is the single most satisfying thing here — it forecloses "you weren't looking hard enough" completely. Three remaining things, one of which matters.

1. Which cell state is the dbSUPER super-enhancer list from?

This is the objection I'd raise as an examiner and it isn't addressed anywhere. Super-enhancers are defined as marking cell-identity genes, so they're among the most cell-state-specific annotations that exist. Your MCC is stimulated CD4. If dbSUPER's CD4 entry is naive CD4 — which is the more commonly deposited dataset — then a substantial fraction of your 158 SE genes may not be super-enhancer-associated in the cells you actually measured.

That doesn't rescue the SE concept, but it does give a reviewer an easy out: "your SE list is from the wrong cell state, so of course it doesn't separate."

Check the dbSUPER sample provenance and state it. If it's naive, either find a stimulated-CD4 SE call (Roadmap/ENCODE H3K27ac on activated CD4 run through ROSE) or state the mismatch explicitly as a limitation. Given the SE null is your headline, it's worth an hour.

2. State which axis Lambert_TF separates on.

d=+0.46 is your matched positive control, but the document doesn't say what direction that is. If TFs are displaced along an interpretable axis — promoter-driven vs enhancer-driven, say — the control becomes a mini-result in its own right and the comparison to SE is much more legible. If it's on an opaque high-order component, that's worth knowing too.

3. The k=2 split recovering arch-HK almost exactly (819/844) is underused.

That means your four-way partition is essentially a two-way split (dispersed vs focal) with one side subdivided. It's a cleaner description of what the geometry actually contains than "four archetypes, imposed," and it connects the widest-cut result in §3 to the group table in §6, which currently read as separate observations. One sentence linking them would tighten the document.

Everything else reads well. Section 7 and the "state that the analysis was redone" instruction in Phase 2 are the two things that will make an examiner trust the rest, and they're both in.