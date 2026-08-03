/**
 * Position-first archetype readout.
 *
 * The label is rendered as a consequence of position, never as the headline.
 * A gene below the core threshold reads "mixture" and shows its composition,
 * because 44% of active genes are blends and rounding them to a corner is the
 * exact dichotomisation this project argues against.
 */

import type { Archetype } from './api'

const BAR_COLOR: Record<string, string> = {
  'arch-HK': '#3d6b91',
  'ME-constitutive': '#2b5070',
  'ME-effector': '#c2703d',
  'arch-sparse': '#8badc9',
  'arch-off': '#c3d6e4',
}

const DISPLAY: Record<string, string> = {
  'arch-HK': 'dispersed',
  'ME-constitutive': 'promoter-local',
  'ME-effector': 'enhancer-focal',
  'arch-sparse': 'sparse',
  'arch-off': 'empty (QC)',
}

export function ArchetypeReadout({ a }: { a: Archetype }) {
  const mix = Object.entries(a.mixture ?? {})
    .filter(([, v]) => v > 0.001)
    .sort((x, y) => y[1] - x[1])

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex items-baseline justify-between">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Position in the continuum
        </h2>
        {a.is_qc_class && (
          <span className="rounded bg-ink-100 px-2 py-0.5 text-[11px] text-ink-600">
            QC class, not biology
          </span>
        )}
      </div>

      <div className="mb-1 flex items-baseline gap-3">
        <span className="text-2xl font-semibold text-ink-900">
          {a.display ?? a.group ?? '—'}
        </span>
        {a.is_mixture ? (
          <span className="rounded bg-element-enhancer/10 px-2 py-0.5 text-xs font-medium text-element-enhancer">
            mixture
          </span>
        ) : (
          <span className="rounded bg-ink-100 px-2 py-0.5 text-xs font-medium text-ink-700">
            core
          </span>
        )}
      </div>

      <p className="mb-3 font-mono text-[11px] text-ink-400">
        {a.group}
        {a.max_posterior != null && <> · p {a.max_posterior.toFixed(2)}</>}
        {a.entropy != null && <> · entropy {a.entropy.toFixed(2)}</>}
      </p>

      {a.architecture && (
        <p className="mb-4 text-[13px] leading-relaxed text-ink-700">{a.architecture}</p>
      )}

      {/* Composition, always shown — including for core genes, so the display
          is the same object in both cases rather than a badge that sometimes
          hides its own uncertainty. */}
      <div className="space-y-1.5">
        {mix.map(([k, v]) => (
          <div key={k} className="flex items-center gap-2">
            <span className="w-32 shrink-0 text-[11px] text-ink-600">
              {DISPLAY[k] ?? k}
            </span>
            <div className="h-2.5 flex-1 overflow-hidden rounded-sm bg-ink-50">
              <div
                className="h-full rounded-sm"
                style={{ width: `${v * 100}%`, background: BAR_COLOR[k] ?? '#5b89ae' }}
              />
            </div>
            <span className="w-10 shrink-0 text-right font-mono text-[11px] text-ink-500">
              {v.toFixed(2)}
            </span>
          </div>
        ))}
      </div>

      {a.caveat && (
        <p className="mt-4 border-t border-ink-100 pt-3 text-[11px] leading-relaxed text-ink-500">
          {a.caveat}
        </p>
      )}
    </div>
  )
}
