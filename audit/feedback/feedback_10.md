Correct — ChromHMM is a segmentation, not a list. It runs a hidden Markov model over histone ChIP tracks and assigns every genomic segment a state. To turn that into a gene set you have to take each gene's TSS, look up which state it falls in, and group accordingly. That's a derived categorisation from a measurement, not a curated annotation.

Contrast with Lambert TF: someone read the literature and identified proteins containing a DNA-binding domain. That's a property of the protein, independent of any chromatin assay.

Two consequences:

The DIRECT/OUTPUT taxonomy is shakier than it looks for the ChromHMM rows. Cluster Claude's distinction was multi-mark HMM vs single-mark rank-cutoff — a real methodological difference — but both are still thresholding operations on ChIP signal. ChromHMM sits closer to dbSUPER than the category boundary implies. Worth acknowledging rather than defending.

More importantly, it sharpens which controls are load-bearing. Chromatin state and 3D contact are both downstream of the same underlying chromatin biology, so "ChromHMM-derived sets separate in MCC space" is partly two assays measuring one thing. That's the "chromatin predicts chromatin" objection, and it applies to active_TSS, bivalent, and Roadmap_silenced alike.

Your genuinely independent controls are the two that touch no chromatin assay at all:

Lambert TF — protein sequence
GWAS_immune_hot — variant association

Those are the ones that show contact architecture corresponds to something outside chromatin measurement. Everything ChromHMM-derived should be presented as corroboration with that caveat stated, not as independent evidence.

Which makes Phase 1.1 more valuable again — a mechanism-defined housekeeping set would be a third control from outside chromatin, and you'd have three independent lines instead of two.

do we have a mechanism defined housekeeping set and anythign else like this for other genes?

Ran a command

No mechanism-defined housekeeping set has been built. But several are available, and the best one is free.

Strongest candidate, already on disk: Lambert DBD subfamilies.

The Lambert file carries DNA-binding-domain classes — C2H2 zinc finger, bHLH, IRF, Rel, forkhead, and so on. You're currently collapsing all 1,600+ human TFs into one set. Splitting by DBD gives you multiple mechanism-defined sets at finer grain, all sequence-derived and independent of any chromatin assay.

It also tests something specific: does the "promoter-driven, CTCF-poor" signature hold uniformly across TF families, or does it differ? That connects directly to your earlier finding that C2H2-ZF and IRF/Rel families split across archetypes. Zero fetch, and it could turn one positive control into three or four.

For housekeeping specifically, in order of independence:

Core promoter sequence classes — TATA-box, Initiator, DPE-containing promoters. Scan hg38 (on disk, same reference as your probes) with the JASPAR core-promoter motifs. TATA-containing is the classic sharp/developmental class and CpG-island-broad is the housekeeping class — this is the Haberle/Stark distinction, mechanism-defined by sequence, and completely independent of chromatin. This is the right test for the housekeeping question.
Translation machinery — ribosomal proteins alone give only 19, but RP + translation initiation and elongation factors + aminoacyl-tRNA synthetases should clear 50. MSigDB is on disk and unused; KEGG_RIBOSOME and the translation GO terms would build it.
Hwang 2023 promoter-assembly genes — the paper is in your folder; check whether it has a supplementary gene list.

Other mechanism-defined sets worth considering (all protein/sequence-defined, chromatin-independent):

Protein complex membership (CORUM) — genes whose products are in the same complex
Imprinted genes — defined by allele-specific regulation mechanism
Other Pfam domain families — kinases, GPCRs, as further Lambert-style controls

I'd do the Lambert DBD split first since it's free, then core promoter sequence classes since that's the one that answers the housekeeping question properly.