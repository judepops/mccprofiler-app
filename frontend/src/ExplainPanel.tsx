/**
 * The argument the app is making, with its numbers attached.
 *
 * Leads with the nested-baseline result, because that is the justification for
 * the 91-feature substrate existing at all — without it, counting peaks would
 * do the job. What the work is NOT sits at the same level as what it is, not in
 * a footnote: the effect sizes are small, most of the association is
 * position-mediated, and half a typical between-gene distance is noise.
 */

import { useEffect, useState } from 'react'
import { api, type Explain } from './api'

const ORDER: [keyof Explain, string][] = [
  ['why_these_features', 'Why these features'],
  ['why_not_clusters', 'Why not clusters'],
  ['why_name_regions_at_all', 'Why name regions at all'],
  ['what_it_is_not', 'What this is not'],
  ['noise_floor', 'The noise floor'],
  ['atac_is_not_a_feature', 'What ATAC does here'],
]

export function ExplainPanel() {
  const [data, setData] = useState<Explain | null>(null)
  const [open, setOpen] = useState<string | null>('why_these_features')

  useEffect(() => {
    api.explain().then(setData).catch(() => setData(null))
  }, [])

  if (!data) return null

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
        How this works, and what it does not show
      </h2>
      <p className="mb-4 text-[11px] text-ink-500">
        {data.provenance.panel} · {data.provenance.n_genes.toLocaleString()} genes ·
        built {data.provenance.built.slice(0, 10)} · pipeline commit{' '}
        <span className="font-mono">{data.provenance.scripts_cleaned_commit}</span>
      </p>

      <div className="space-y-2">
        {ORDER.map(([key, heading]) => {
          const section = data[key] as
            | { claim: string; detail: string; citations?: string[] }
            | undefined
          if (!section?.claim) return null
          const isOpen = open === key
          const isCaveat = key === 'what_it_is_not' || key === 'noise_floor'
          return (
            <div
              key={String(key)}
              className={`rounded border ${
                isCaveat ? 'border-element-enhancer/30' : 'border-ink-100'
              }`}
            >
              <button
                onClick={() => setOpen(isOpen ? null : String(key))}
                className="flex w-full items-baseline justify-between gap-3 px-3 py-2 text-left hover:bg-ink-50"
              >
                <span>
                  <span
                    className={`mr-2 text-[10px] font-semibold uppercase tracking-wide ${
                      isCaveat ? 'text-element-enhancer' : 'text-ink-400'
                    }`}
                  >
                    {heading}
                  </span>
                  <span className="text-[13px] text-ink-800">{section.claim}</span>
                </span>
                <span className="shrink-0 font-mono text-[11px] text-ink-400">
                  {isOpen ? '−' : '+'}
                </span>
              </button>

              {isOpen && (
                <div className="border-t border-ink-100 px-3 py-2">
                  <p className="text-[12px] leading-relaxed text-ink-700">
                    {section.detail}
                  </p>
                  {section.citations && (
                    <p className="mt-2 text-[11px] text-ink-500">
                      {section.citations.join(' · ')}
                    </p>
                  )}

                  {/* the nested-baseline table, where it belongs */}
                  {key === 'why_these_features' && data.why_these_features.table && (
                    <table className="mt-3 w-full text-[11px]">
                      <thead>
                        <tr className="border-b border-ink-100 text-ink-500">
                          {Object.keys(data.why_these_features.table[0] ?? {})
                            .slice(0, 5)
                            .map((c) => (
                              <th key={c} className="py-1 pr-3 text-left font-normal">
                                {c}
                              </th>
                            ))}
                        </tr>
                      </thead>
                      <tbody>
                        {data.why_these_features.table.slice(0, 8).map((r, i) => (
                          <tr key={i} className="border-b border-ink-50 last:border-0">
                            {Object.values(r)
                              .slice(0, 5)
                              .map((v, j) => (
                                <td key={j} className="py-1 pr-3 font-mono text-ink-700">
                                  {typeof v === 'number' ? v.toFixed(3) : String(v)}
                                </td>
                              ))}
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  )}
                </div>
              )}
            </div>
          )
        })}
      </div>
    </div>
  )
}
