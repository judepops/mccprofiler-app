/**
 * The coordinate space, chosen once for the whole page.
 *
 * Why this module exists. Every view used to hold its own axis state, so the
 * continuum map could be showing corrected axes while the ranked lists showed
 * raw PCs, with nothing on screen saying which was which. Two coordinate
 * systems in one page is the same confusion that produced the sPC/PC mix-ups
 * in the written analysis, and it had already reached the API: the loadings
 * endpoint was serving shape loadings under the label "PC2 (local vs
 * long-range)".
 *
 * So the space is a page-level decision and every view follows it.
 *
 *   corrected (sPC)  magnitude is projected out before PCA, so every component
 *                    is uncorrelated with overall signal BY CONSTRUCTION and a
 *                    position can never mean "this gene has more signal". This
 *                    is the space every Aim 3 result is computed in, and the
 *                    default.
 *
 *   raw (PC)         provenance only. On the current 73-feature substrate PC1
 *                    correlates 0.505 with the magnitude basis and PC3
 *                    correlates 0.519, so amount is smeared across two of the
 *                    top three components and neither is cleanly nameable.
 *
 * The numbering does NOT correspond between spaces. sPC3 is not PC3. Never
 * compare them by index, which is the same error as joining dimension_names.tsv
 * on PC index.
 */

export type Space = 'corrected' | 'raw'

export const SPACE_PREFIX: Record<Space, string> = {
  corrected: 'spc',
  raw: 'pc',
}

/** Default 2D plane per space.
 *
 * corrected: sPC1 x sPC2. With amount removed there is no reason to skip the
 * first component.
 *
 * raw: PC2 x PC3, the historical default. It was chosen because PC1 was the
 * amount axis and worth skipping; that justification no longer holds, which is
 * part of why this space is provenance-only now.
 */
export const DEFAULT_PLANE: Record<Space, [string, string]> = {
  corrected: ['spc1', 'spc2'],
  raw: ['pc2', 'pc3'],
}

/** Component n in the given space: axisKey('corrected', 2) -> 'spc2'. */
export function axisKey(space: Space, n: number): string {
  return `${SPACE_PREFIX[space]}${n}`
}

/** Which space an axis key belongs to. UMAP axes belong to neither. */
export function spaceOf(axis: string): Space | null {
  if (axis.startsWith('spc')) return 'corrected'
  if (axis.startsWith('umap')) return null
  if (axis.startsWith('pc')) return 'raw'
  return null
}

/** Does this axis belong to the active space? UMAP passes in both, since it is
 *  a separate embedding and is labelled as such wherever it is offered. */
export function inSpace(axis: string, space: Space): boolean {
  const s = spaceOf(axis)
  return s === null ? true : s === space
}

/** Move a plane into another space, keeping the component numbers where the
 *  target has them. Used when the toggle flips: staying on sPC1 x sPC2 while
 *  the page claims to be showing raw PCs would be exactly the mix-up this
 *  module prevents. UMAP planes are left alone. */
export function planeInSpace(
  plane: [string, string],
  space: Space,
): [string, string] {
  const moved = plane.map((a) => {
    if (a.startsWith('umap')) return a
    const n = a.replace(/^s?pc/, '')
    return /^\d+$/.test(n) ? axisKey(space, Number(n)) : a
  }) as [string, string]
  // If the source plane was UMAP, or anything failed to map, fall back to the
  // target space's default rather than leaving a mixed pair on screen.
  return moved.every((a) => inSpace(a, space)) ? moved : DEFAULT_PLANE[space]
}

export const SPACE_LABEL: Record<Space, string> = {
  corrected: 'amount-corrected (sPC)',
  raw: 'raw (PC)',
}

export const SPACE_NOTE: Record<Space, string> = {
  corrected:
    'Magnitude is projected out before PCA, so no axis carries overall signal ' +
    'and a position cannot mean "this gene has more signal". This is the space ' +
    'the analysis uses.',
  raw:
    'Provenance only. PC1 and PC3 both carry magnitude (r 0.505 and 0.519) on ' +
    'the current substrate, so position here is partly a statement about how ' +
    'much signal a gene has, not only about its architecture.',
}
