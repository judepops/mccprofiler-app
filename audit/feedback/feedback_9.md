Files updated this run: backend/scripts/make_audit_figures.py (fig 3 rebuilt) and audit/CURRENT_FINDINGS.md.

Your diagnosis of the right panel was exactly right — overlap is a deterministic function of |d| for two normals, so it was a curve by construction and every y was derivable from its x. It is now three density overlays: Lambert TF, a mid-pack set, and the super-enhancer null, each against the panel on its own component. Seeing the SE curve sitting on top of the panel curve makes the argument in a way "94% overlap" never could.

The sign point was a genuine error, not a presentation choice, and I have propagated the fix. The figure now plots |d| only, and the findings document carries the same warning: each set is measured on its own strongest component, PCA sign is arbitrary per component, and signed values are not comparable between sets. The signs remain in the table only so rows can be traced to their component.

Left panel also now: coloured by definition type (so SE sits visibly mid-pack among the output-defined sets rather than being singled out in orange, which implied exactly the specialness I had just retracted), significance marked, n printed per row, and the gene-desert bar labelled in the legend as the density control.

Your CpG instinct was right and the problem is broader than that one set. Three sets cover more than 80% of the panel:

set	n	% of panel	"rest"
ChromHMM_active_TSS	1,619	88%	227
DICE_top_TPM_quartile	1,539	83%	307
CpG_island_promoter	1,526	83%	320
Their "set versus rest" contrast is weak by construction — the comparison group is a small unusual remainder. That explains the CpG near-null without needing it to be biology. I checked whether the definition-type result depends on them: dropping all three gives p = 0.0019 against 0.0018, medians 0.28 and 0.14. Unchanged.

And your LOEUF observation is now in the document as internal consistency. Architecture predicts constraint strongly as a continuous quantity (4.08 SD, 91% surviving amount correction) while the top-quartile constraint set barely displaces (|d| = 0.17). Same claim reached twice from independent analyses, and it makes "directions, not regions" concrete rather than abstract. That is a better argument than either half alone.