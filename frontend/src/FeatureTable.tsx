/**
 * The 91 MCCProfiler features for one gene, grouped by block.
 *
 * Replaces the feature table at the foot of the J1 case-study figures. Two
 * changes of substance, not just formatting:
 *
 *  - **Percentile leads, z is secondary.** A z of +1.74 means nothing without
 *    knowing the distribution; "95th percentile of the panel" is directly
 *    readable. The original figure showed rank as `#45/791`, which requires
 *    mental arithmetic to interpret.
 *  - **One colour scale, stated.** The original coloured each cell against a
 *    different reference and said so in small print. Here a single diverging
 *    scale runs from bottom to top of the panel distribution, and the legend
 *    is on screen rather than in a caption.
 */

import { useState } from 'react'
import type { FeatureSet } from './api'
import { FeatureExplainer } from './FeatureExplainer'

/** Diverging low→high, colourblind-safe, avoiding red/green opposition. */
function scaleColor(pct: number): string {
  const t = Math.max(0, Math.min(100, pct)) / 100
  if (t < 0.5) {
    const k = t / 0.5 // 0 = extreme low
    return `rgba(43, 80, 112, ${0.42 * (1 - k) + 0.03})`
  }
  const k = (t - 0.5) / 0.5 // 1 = extreme high
  return `rgba(194, 112, 61, ${0.42 * k + 0.03})`
}

export function FeatureTable({
  data,
  onPick,
}: {
  data: FeatureSet
  onPick?: (s: string) => void
}) {
  const [open, setOpen] = useState<Set<string>>(new Set(['G2', 'P1']))
  const [onlyExtreme, setOnlyExtreme] = useState(false)
  // Clicking a feature name opens "what does this measure?", the gap the
  // percentile bars leave.
  const [explain, setExplain] = useState<string | null>(null)

  function toggle(block: string) {
    const next = new Set(open)
    next.has(block) ? next.delete(block) : next.add(block)
    setOpen(next)
  }

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex items-baseline justify-between">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          MCCProfiler features · {data.n_features}
        </h2>
        <label className="flex items-center gap-1.5 text-[11px] text-ink-600">
          <input
            type="checkbox"
            checked={onlyExtreme}
            onChange={(e) => setOnlyExtreme(e.target.checked)}
          />
          only top/bottom 10%
        </label>
      </div>

      {/* colour legend, on screen rather than in a caption */}
      <div className="mb-4 flex items-center gap-2 text-[10px] text-ink-500">
        <span>bottom of panel</span>
        <div className="flex h-2 flex-1 overflow-hidden rounded-sm">
          {Array.from({ length: 40 }, (_, i) => (
            <div key={i} className="flex-1" style={{ background: scaleColor((i / 39) * 100) }} />
          ))}
        </div>
        <span>top of panel</span>
      </div>

      <div className="space-y-2">
        {data.blocks.map((b) => {
          const shown = onlyExtreme
            ? b.features.filter((f) => f.percentile >= 90 || f.percentile <= 10)
            : b.features
          if (!shown.length) return null
          const isOpen = open.has(b.block)
          return (
            <div key={b.block} className="rounded border border-ink-100">
              <button
                onClick={() => toggle(b.block)}
                className="flex w-full items-baseline justify-between px-3 py-2 text-left hover:bg-ink-50"
              >
                <span className="text-[13px]">
                  <span className="font-mono font-semibold text-ink-700">{b.block}</span>
                  <span className="ml-2 text-ink-600">{b.description}</span>
                </span>
                <span className="font-mono text-[11px] text-ink-400">
                  {shown.length} {isOpen ? '−' : '+'}
                </span>
              </button>

              {isOpen && (
                <table className="w-full border-t border-ink-100 text-[12px]">
                  <tbody>
                    {shown.map((f) => (
                      <tr key={f.name} className="border-b border-ink-50 last:border-0">
                        <td className="py-1 pl-3 pr-2">
                          <button
                            onClick={() => setExplain(f.name)}
                            className={`font-mono text-[11px] hover:underline ${
                              explain === f.name
                                ? 'font-semibold text-ink-900'
                                : 'text-ink-700'
                            }`}
                            title="What does this measure?"
                          >
                            {f.name}
                          </button>
                        </td>
                        <td className="w-40 py-1 pr-2">
                          <div className="h-3 w-full overflow-hidden rounded-sm bg-ink-50">
                            <div
                              className="h-full"
                              style={{
                                width: `${f.percentile}%`,
                                background: scaleColor(f.percentile),
                              }}
                            />
                          </div>
                        </td>
                        <td className="w-14 py-1 pr-3 text-right font-mono text-[11px] text-ink-600">
                          {f.percentile.toFixed(0)}%
                        </td>
                        <td className="w-16 py-1 pr-3 text-right font-mono text-[11px] text-ink-400">
                          {f.z >= 0 ? '+' : ''}
                          {f.z.toFixed(2)}σ
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
            </div>
          )
        })}
      </div>

      {explain && (
        <div className="mt-4">
          <FeatureExplainer
            feature={explain}
            onClose={() => setExplain(null)}
            onPick={onPick}
          />
        </div>
      )}

      <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
        Click a feature name to see what it measures, drawn on the panel's most and
        least extreme genes. Percentile is against the panel of 1,846 genes, computed per feature. Bars
        show percentile; σ is the z-score on the same scaled substrate the clustering uses.
        Features are correlation-pruned at |r| &gt; 0.8, so the set is not exhaustive of
        everything computed.
      </p>
    </div>
  )
}
