/**
 * Scree + loadings, the evidence for the axis names.
 *
 * The continuum map relabels PC1-PC5 with prose ("promoter-driven vs
 * enhancer-driven"). Those readings come from the audit's 14_name_dimensions.py,
 * but a name is an *interpretation of loadings*, and showing the name without
 * them asks the reader to take it on trust. This panel is where the
 * interpretation can be checked and, if it is wrong, seen to be wrong.
 *
 * Panel-level, not gene-level, it describes the coordinate system, not any
 * particular gene.
 */

import { useEffect, useState } from 'react'
import { api, type Loadings, type Scree } from './api'
import { type Space } from './space'

export function DimensionPanel({ space = 'corrected' }: { space?: Space }) {
  const [scree, setScree] = useState<Scree | null>(null)
  const [pc, setPc] = useState(1)
  const [load, setLoad] = useState<Loadings | null>(null)

  // Both requests carry the space. They used to disagree: the scree came from
  // the raw PCs while the loadings came from the corrected ones, so the panel
  // showed "PC1, 17.23%" above sPC1's feature list.
  useEffect(() => {
    api.scree(space).then(setScree).catch(() => setScree(null))
  }, [space])
  useEffect(() => {
    api.loadings(pc, 15, space).then(setLoad).catch(() => setLoad(null))
  }, [pc, space])

  const maxVar = scree ? Math.max(...scree.rows.map((r) => r.variance_pct)) : 1
  const maxLoad = load ? Math.max(...load.loadings.map((l) => Math.abs(l.loading))) : 1

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Coordinate system, components and what they are made of
      </h2>
      <p className="mb-3 text-[11px] leading-relaxed text-ink-500">
        Panel-level, not gene-specific. Axis names are readings of the loadings below;
        if a name does not match its loadings, the name is wrong.
      </p>

      {/* Measured, not asserted. Reproduce with
          backend/scripts/diagnose_pc_names.py. This sits at the top of the
          panel because it governs how every name below should be read. */}
      <div className="mb-4 rounded border border-element-enhancer/30 bg-element-enhancer/5 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
        <strong>How much of a name is true.</strong> Each axis name summarises one
        contrast inside a mixed axis. The share of an axis's squared loading mass
        that sits in the concept its name refers to is{' '}
        <span className="font-mono">23%</span> for PC1,{' '}
        <span className="font-mono">26%</span> PC2,{' '}
        <span className="font-mono">31%</span> PC3,{' '}
        <span className="font-mono">35%</span> PC4,{' '}
        <span className="font-mono">22%</span> PC5. All are well above chance, so no
        name is invented, but a name accounts for roughly a quarter to a third of
        its axis. Two specifics worth carrying:{' '}
        <strong>PC1 is not signal amount</strong> (it correlates with{' '}
        <span className="font-mono">total_mcc</span> at{' '}
        <span className="font-mono">+0.03</span>; it tracks peak counts and reach),
        and <strong>amount lives on PC3</strong> (
        <span className="font-mono">-0.62</span>), so the promoter-versus-enhancer
        contrast there is partly a statement about how much signal a gene has.
        PC1-5 cover 45% of variance while 19 components clear the noise ceiling at
        78%, so about a third of the real structure is unnamed.
      </div>

      <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
        {/* scree ------------------------------------------------------------ */}
        <div>
          <h3 className="mb-2 text-[11px] font-medium text-ink-700">
            Variance per component
            {scree && (
              <span className="ml-2 font-normal text-ink-500">
                {scree.n_above_noise} above the noise ceiling
              </span>
            )}
          </h3>

          <div className="space-y-0.5">
            {scree?.rows.slice(0, 20).map((r) => {
              const active = r.pc === pc
              return (
                <button
                  key={r.pc}
                  onClick={() => setPc(r.pc)}
                  className={`flex w-full items-center gap-2 rounded px-1 py-0.5 text-left hover:bg-ink-50 ${
                    active ? 'bg-ink-50' : ''
                  }`}
                >
                  <span
                    className={`w-10 shrink-0 font-mono text-[10px] ${
                      active ? 'font-semibold text-ink-900' : 'text-ink-500'
                    }`}
                  >
                    PC{r.pc}
                  </span>
                  <div className="relative h-3 flex-1 rounded-sm bg-ink-50">
                    <div
                      className="absolute inset-y-0 left-0 rounded-sm"
                      style={{
                        width: `${(r.variance_pct / maxVar) * 100}%`,
                        background: r.above_noise ? '#2b5070' : '#c3d6e4',
                      }}
                    />
                    {/* noise ceiling marker */}
                    <div
                      className="absolute inset-y-0 w-px bg-element-enhancer"
                      style={{ left: `${(r.noise_pct / maxVar) * 100}%` }}
                    />
                  </div>
                  <span className="w-10 shrink-0 text-right font-mono text-[10px] text-ink-500">
                    {r.variance_pct.toFixed(1)}%
                  </span>
                </button>
              )
            })}
          </div>

          <p className="mt-2 text-[10px] leading-relaxed text-ink-500">
            <span className="inline-block h-2 w-2 rounded-sm bg-ink-600" /> above noise ·{' '}
            <span className="inline-block h-2 w-2 rounded-sm bg-ink-200" /> below ·{' '}
            <span className="inline-block h-2 w-px bg-element-enhancer align-middle" /> ceiling.
            {scree && ` ${scree.note}`}
          </p>
        </div>

        {/* loadings --------------------------------------------------------- */}
        <div>
          <h3 className="mb-2 text-[11px] font-medium text-ink-700">
            {load?.label ?? `PC${pc}`}
            {load?.variance_pct != null && (
              <span className="ml-2 font-normal text-ink-500">
                {load.variance_pct.toFixed(2)}% of variance
              </span>
            )}
          </h3>

          {load && (
            <div className="space-y-1">
              {load.loadings.map((l) => {
                const frac = Math.abs(l.loading) / maxLoad
                const pos = l.loading >= 0
                return (
                  <div key={l.feature} className="flex items-center gap-2">
                    <span
                      className="w-52 shrink-0 truncate font-mono text-[10px] text-ink-600"
                      title={l.feature}
                    >
                      {l.feature}
                    </span>
                    <div className="relative h-3 flex-1">
                      <div className="absolute inset-y-0 left-1/2 w-px bg-ink-200" />
                      <div
                        className="absolute inset-y-0 rounded-sm"
                        style={{
                          left: pos ? '50%' : `${50 - frac * 50}%`,
                          width: `${frac * 50}%`,
                          background: pos ? '#c2703d' : '#2b5070',
                          opacity: 0.8,
                        }}
                      />
                    </div>
                    <span className="w-12 shrink-0 text-right font-mono text-[10px] text-ink-500">
                      {l.loading >= 0 ? '+' : ''}
                      {l.loading.toFixed(3)}
                    </span>
                  </div>
                )
              })}
            </div>
          )}

          <p className="mt-2 text-[10px] leading-relaxed text-ink-500">
            Top {load?.loadings.length ?? 0} features by absolute loading, of{' '}
            {load?.n_features ?? 0}. Sign is arbitrary in PCA, only the *contrast*
            between the two ends is meaningful, not which end is positive.
          </p>
        </div>
      </div>
    </div>
  )
}
