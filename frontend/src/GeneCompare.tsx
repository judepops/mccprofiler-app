/**
 * Two genes side by side.
 *
 * The natural follow-up to a lookup: "how does this one differ from that one?".
 * Shows where they diverge most on the 91 features, and their positions on the
 * named axes, rather than two plots the reader has to diff by eye.
 */

import { useEffect, useState } from 'react'
import { api, type FeatureSet, type Gene } from './api'

interface Row {
  name: string
  a: number
  b: number
  aPct: number
  bPct: number
  gap: number
}

export function GeneCompare({
  primary,
  primaryFeatures,
}: {
  primary: Gene
  primaryFeatures: FeatureSet | null
}) {
  const [query, setQuery] = useState('')
  const [other, setOther] = useState<Gene | null>(null)
  const [otherFeatures, setOtherFeatures] = useState<FeatureSet | null>(null)
  const [err, setErr] = useState<string | null>(null)

  useEffect(() => {
    // Clear when the primary gene changes — a stale comparison is worse than none.
    setOther(null)
    setOtherFeatures(null)
    setQuery('')
  }, [primary.gene_id])

  async function load(symbol: string) {
    setErr(null)
    try {
      const [g, f] = await Promise.all([api.gene(symbol), api.features(symbol)])
      setOther(g)
      setOtherFeatures(f)
    } catch (e) {
      setErr(String((e as Error).message))
    }
  }

  let rows: Row[] = []
  if (primaryFeatures && otherFeatures) {
    const bIdx = new Map<string, { z: number; percentile: number }>()
    for (const blk of otherFeatures.blocks) {
      for (const f of blk.features) bIdx.set(f.name, { z: f.z, percentile: f.percentile })
    }
    for (const blk of primaryFeatures.blocks) {
      for (const f of blk.features) {
        const m = bIdx.get(f.name)
        if (!m) continue
        rows.push({
          name: f.name,
          a: f.z,
          b: m.z,
          aPct: f.percentile,
          bPct: m.percentile,
          gap: Math.abs(f.z - m.z),
        })
      }
    }
    rows.sort((x, y) => y.gap - x.gap)
    rows = rows.slice(0, 12)
  }

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-3 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Compare with another gene
      </h2>

      <div className="mb-4 flex gap-2">
        <input
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          onKeyDown={(e) => e.key === 'Enter' && query.trim() && load(query.trim())}
          placeholder={`Compare ${primary.gene_symbol} with…`}
          className="flex-1 rounded border border-ink-200 px-3 py-1.5 text-[13px] outline-none focus:border-ink-400"
        />
        <button
          onClick={() => query.trim() && load(query.trim())}
          className="rounded bg-ink-600 px-3 py-1.5 text-xs text-white hover:bg-ink-700"
        >
          Compare
        </button>
      </div>

      {err && <p className="text-[12px] text-element-enhancer">{err}</p>}

      {other && (
        <>
          <div className="mb-4 grid grid-cols-2 gap-4">
            {[
              { g: primary, tag: 'A', color: '#2b5070' },
              { g: other, tag: 'B', color: '#c2703d' },
            ].map(({ g, tag, color }) => (
              <div key={g.gene_id} className="rounded border border-ink-100 p-2.5">
                <div className="flex items-baseline gap-2">
                  <span
                    className="inline-block h-2.5 w-2.5 rounded-sm"
                    style={{ background: color }}
                  />
                  <span className="text-sm font-semibold">{g.gene_symbol}</span>
                  <span className="font-mono text-[10px] text-ink-400">{tag}</span>
                </div>
                <p className="mt-1 text-[11px] text-ink-600">
                  {g.archetype.display}
                  {g.archetype.is_mixture && (
                    <span className="ml-1 text-element-enhancer">· mixture</span>
                  )}
                </p>
                <p className="font-mono text-[10px] text-ink-400">
                  p {g.archetype.max_posterior?.toFixed(2)} ·{' '}
                  {g.viewpoint.chrom}:{g.viewpoint.pos?.toLocaleString()}
                </p>
              </div>
            ))}
          </div>

          <h3 className="mb-2 text-[11px] font-medium text-ink-700">
            Where they differ most
          </h3>
          <div className="space-y-1.5">
            {rows.map((r) => (
              <div key={r.name} className="flex items-center gap-2">
                <span
                  className="w-52 shrink-0 truncate font-mono text-[10px] text-ink-600"
                  title={r.name}
                >
                  {r.name}
                </span>
                {/* percentile positions on a shared 0-100 track */}
                <div className="relative h-4 flex-1 rounded-sm bg-ink-50">
                  <div
                    className="absolute top-1/2 h-px bg-ink-300"
                    style={{
                      left: `${Math.min(r.aPct, r.bPct)}%`,
                      width: `${Math.abs(r.aPct - r.bPct)}%`,
                    }}
                  />
                  <div
                    className="absolute top-1/2 h-2.5 w-2.5 -translate-x-1/2 -translate-y-1/2 rounded-sm"
                    style={{ left: `${r.aPct}%`, background: '#2b5070' }}
                    title={`A ${r.aPct.toFixed(0)}%`}
                  />
                  <div
                    className="absolute top-1/2 h-2.5 w-2.5 -translate-x-1/2 -translate-y-1/2 rounded-sm"
                    style={{ left: `${r.bPct}%`, background: '#c2703d' }}
                    title={`B ${r.bPct.toFixed(0)}%`}
                  />
                </div>
                <span className="w-12 shrink-0 text-right font-mono text-[10px] text-ink-500">
                  {r.gap.toFixed(2)}σ
                </span>
              </div>
            ))}
          </div>

          <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            Top 12 features by absolute difference in z. Squares are panel percentiles, so
            the track is comparable across features; the number is the gap in standard
            deviations. Remember roughly half a typical between-gene distance is technical
            noise (gene-to-itself 5.36 vs gene-to-other 10.49 across repeat captures), so
            small gaps are not meaningful.
          </p>
        </>
      )}
    </div>
  )
}
