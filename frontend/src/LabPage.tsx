/**
 * Lab — exploratory views, deliberately separate from the main app.
 *
 * Things worth looking at that are not part of the product: methodological
 * checks, diagnostics, and views whose caveats are too heavy for a page a
 * collaborator might read quickly. Nothing here should be quoted without the
 * caveat attached to it.
 */

import { useEffect, useState } from 'react'
import { api, type Repro } from './api'

function Scatter({
  points,
  feature,
}: {
  points: { symbol_key: string; gw: number; immune: number }[]
  feature: string
}) {
  const W = 340
  const H = 340
  const PAD = 34
  const all = points.flatMap((p) => [p.gw, p.immune])
  const lo = Math.min(...all)
  const hi = Math.max(...all)
  const sx = (v: number) => PAD + ((v - lo) / (hi - lo || 1)) * (W - PAD * 1.4)
  const sy = (v: number) => H - PAD - ((v - lo) / (hi - lo || 1)) * (H - PAD * 1.4)

  return (
    <svg width={W} height={H} className="max-w-full">
      {/* identity line: perfect reproducibility would put every point on it */}
      <line
        x1={sx(lo)} y1={sy(lo)} x2={sx(hi)} y2={sy(hi)}
        stroke="#c3d6e4" strokeDasharray="3 3"
      />
      <line x1={PAD} y1={PAD * 0.4} x2={PAD} y2={H - PAD} stroke="#c3d6e4" />
      <line x1={PAD} y1={H - PAD} x2={W - PAD * 0.4} y2={H - PAD} stroke="#c3d6e4" />
      {points.map((p) => (
        <circle
          key={p.symbol_key}
          cx={sx(p.gw)}
          cy={sy(p.immune)}
          r={2.6}
          fill="#2b5070"
          fillOpacity={0.55}
        >
          <title>{`${p.symbol_key}\nGW ${p.gw.toFixed(2)} · immune ${p.immune.toFixed(2)}`}</title>
        </circle>
      ))}
      <text x={W / 2} y={H - 6} textAnchor="middle" className="fill-ink-500 text-[10px]">
        {feature} — genome-wide capture
      </text>
      <text
        x={11} y={H / 2} textAnchor="middle" className="fill-ink-500 text-[10px]"
        transform={`rotate(-90 11 ${H / 2})`}
      >
        immune capture
      </text>
    </svg>
  )
}

export function LabPage() {
  const [data, setData] = useState<Repro | null>(null)
  const [feature, setFeature] = useState('frac_local')

  useEffect(() => {
    api.labReproducibility(feature).then(setData).catch(() => setData(null))
  }, [feature])

  return (
    <div className="space-y-5">
      <div className="rounded-lg border border-element-enhancer/40 bg-element-enhancer/5 px-4 py-3">
        <h2 className="text-sm font-semibold text-ink-900">Lab</h2>
        <p className="mt-1 text-[12px] leading-relaxed text-ink-700">
          Exploratory views, not part of the app. Methodological checks and diagnostics
          whose caveats are too heavy for a page someone might read quickly. Nothing here
          should be quoted without its caveat.
        </p>
      </div>

      <div className="rounded-lg border border-ink-200 bg-white p-4">
        <h3 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
          Technical reproducibility — the same gene, captured twice
        </h3>
        <p className="mb-3 text-[11px] leading-relaxed text-ink-500">
          {data?.n_genes} genes appear in both the genome-wide and immune panels. Two
          independent captures, same pipeline, so the spread is technical noise rather
          than biology. This is what the κ = 0.72 label-reproducibility and the 0.511
          noise floor are made of.
        </p>

        {data && (
          <div className="mb-4 rounded border border-ink-200 bg-ink-50 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
            <strong>Caveat.</strong> {data.caveat}
          </div>
        )}

        {data && (
          <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
            {/* per-feature reproducibility */}
            <div>
              <div className="mb-2 flex items-baseline gap-3 text-[11px]">
                <span className="text-ink-700">
                  median ρ <span className="font-mono">{data.median_rho.toFixed(3)}</span>
                </span>
                <span className="text-ink-500">
                  {data.n_above_0_7}/{data.n_features} above 0.7
                </span>
                <span className="text-ink-500">{data.n_below_0_3} below 0.3</span>
              </div>

              <div className="max-h-[380px] space-y-0.5 overflow-y-auto pr-1">
                {data.per_feature.map((f) => {
                  const bad = f.spearman_rho < 0.3
                  return (
                    <button
                      key={f.feature}
                      onClick={() => setFeature(f.feature)}
                      className={`flex w-full items-center gap-2 rounded px-1 py-0.5 text-left hover:bg-ink-50 ${
                        f.feature === feature ? 'bg-ink-50' : ''
                      }`}
                    >
                      <span
                        className="w-52 shrink-0 truncate font-mono text-[10px] text-ink-600"
                        title={f.feature}
                      >
                        {f.feature}
                      </span>
                      <div className="h-2.5 flex-1 rounded-sm bg-ink-50">
                        <div
                          className="h-full rounded-sm"
                          style={{
                            width: `${Math.max(f.spearman_rho, 0) * 100}%`,
                            background: bad ? '#c2703d' : '#2b5070',
                          }}
                        />
                      </div>
                      <span className="w-10 shrink-0 text-right font-mono text-[10px] text-ink-500">
                        {f.spearman_rho.toFixed(2)}
                      </span>
                    </button>
                  )
                })}
              </div>
              <p className="mt-2 text-[10px] leading-relaxed text-ink-500">
                {data.note} Click a feature to see its paired measurements.
              </p>
            </div>

            {/* paired scatter */}
            <div>
              {data.pairs ? (
                <>
                  <Scatter points={data.pairs.points} feature={data.pairs.feature} />
                  <p className="mt-1 text-[10px] leading-relaxed text-ink-500">
                    Each point is one gene measured twice. The dashed line is perfect
                    agreement — spread away from it is measurement noise, and it bounds
                    how well anything downstream can possibly do.
                  </p>
                </>
              ) : (
                <p className="text-[11px] text-ink-400">Select a feature.</p>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
