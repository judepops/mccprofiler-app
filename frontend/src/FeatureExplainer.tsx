/**
 * What a feature measures, shown on two real traces.
 *
 * Extends the spatial-zones schematic in report.tex from a diagram of the
 * cutoffs to an overlay on actual signal: the panel's highest and lowest gene
 * on the feature, side by side, with the region the feature reads shaded on
 * both. The contrast is what makes it legible, a single median profile shows
 * nothing.
 *
 * Features with no faithful geometry say so rather than getting a
 * plausible-looking overlay. Topology is the loud case: its adjacency is the
 * complete graph on active peaks, so a network drawing would imply structure
 * the data does not measure.
 */

import { useEffect, useRef, useState } from 'react'
import { api, formatBp, type FeatureExplain } from './api'

const PAD = { top: 8, right: 8, bottom: 18, left: 8 }

function MiniProfile({
  values,
  startBp,
  endBp,
  geom,
  color,
}: {
  values: number[]
  startBp: number
  endBp: number
  geom: FeatureExplain
  color: string
}) {
  const ref = useRef<HTMLCanvasElement>(null)
  const wrapRef = useRef<HTMLDivElement>(null)
  const [w, setW] = useState(340)
  const h = 120

  useEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver(([e]) => setW(e.contentRect.width))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  useEffect(() => {
    const c = ref.current
    if (!c) return
    const dpr = window.devicePixelRatio || 1
    c.width = w * dpr
    c.height = h * dpr
    const ctx = c.getContext('2d')
    if (!ctx) return
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, w, h)

    const pw = w - PAD.left - PAD.right
    const ph = h - PAD.top - PAD.bottom
    const span = endBp - startBp
    const x = (bp: number) => PAD.left + ((bp - startBp) / span) * pw
    const yMax = Math.max(...values, 1)

    // ---- the region the feature actually reads --------------------------
    // Bands are defined on |distance|, so each one is drawn twice.
    if (geom.regions) {
      for (const r of geom.regions) {
        const sides: [number, number][] = r.mirrored
          ? [
              [r.start_bp, r.end_bp],
              [-r.end_bp, -r.start_bp],
            ]
          : [[r.start_bp, r.end_bp]]
        for (const [a, b] of sides) {
          const x0 = Math.max(x(a), PAD.left)
          const x1 = Math.min(x(b), PAD.left + pw)
          if (x1 <= x0) continue
          ctx.fillStyle =
            r.role === 'denominator' ? 'rgba(139,173,201,0.30)' : 'rgba(194,112,61,0.22)'
          ctx.fillRect(x0, PAD.top, x1 - x0, ph)
        }
      }
    }

    // ---- the trace -------------------------------------------------------
    const colMax = new Float64Array(Math.ceil(pw))
    for (let i = 0; i < values.length; i++) {
      const col = Math.min(Math.floor((i / values.length) * pw), colMax.length - 1)
      if (values[i] > colMax[col]) colMax[col] = values[i]
    }
    ctx.fillStyle = color
    for (let col = 0; col < colMax.length; col++) {
      const bh = (colMax[col] / yMax) * ph
      if (bh > 0) ctx.fillRect(PAD.left + col, PAD.top + ph - bh, 1, bh)
    }

    // ---- threshold features: a line at the quantile -----------------------
    if (geom.kind === 'threshold') {
      const sorted = [...values].sort((a, b) => a - b)
      const q = geom.quantile ?? 1
      const v = sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))]
      const y = PAD.top + ph - (v / yMax) * ph
      ctx.strokeStyle = '#c2703d'
      ctx.setLineDash([4, 3])
      ctx.beginPath()
      ctx.moveTo(PAD.left, y)
      ctx.lineTo(PAD.left + pw, y)
      ctx.stroke()
      ctx.setLineDash([])
    }

    // ---- viewpoint --------------------------------------------------------
    ctx.strokeStyle = 'rgba(15,30,46,0.35)'
    ctx.beginPath()
    ctx.moveTo(x(0), PAD.top)
    ctx.lineTo(x(0), PAD.top + ph)
    ctx.stroke()

    ctx.font = '9px ui-monospace, monospace'
    ctx.fillStyle = '#8badc9'
    ctx.textAlign = 'left'
    ctx.fillText(formatBp(startBp), PAD.left, h - 5)
    ctx.textAlign = 'right'
    ctx.fillText(formatBp(endBp), PAD.left + pw, h - 5)
    ctx.textAlign = 'center'
    ctx.fillText('VP', x(0), h - 5)
  }, [values, startBp, endBp, geom, color, w])

  return (
    <div ref={wrapRef} className="w-full">
      <canvas ref={ref} style={{ width: '100%', height: h }} />
    </div>
  )
}

export function FeatureExplainer({
  feature,
  onClose,
  onPick,
}: {
  feature: string
  onClose: () => void
  onPick?: (s: string) => void
}) {
  const [data, setData] = useState<FeatureExplain | null>(null)

  useEffect(() => {
    setData(null)
    api.feature(feature).then(setData).catch(() => setData(null))
  }, [feature])

  if (!data) return null

  return (
    <div className="rounded-lg border border-ink-300 bg-white p-4">
      <div className="mb-2 flex items-baseline justify-between">
        <h2 className="font-mono text-[13px] font-semibold text-ink-900">
          {data.feature}
        </h2>
        <button onClick={onClose} className="text-[11px] text-ink-500 hover:text-ink-800">
          close
        </button>
      </div>

      {data.measures && (
        <p className="mb-3 text-[12px] leading-relaxed text-ink-700">{data.measures}</p>
      )}

      {!data.drawable ? (
        <div className="rounded border border-element-enhancer/40 bg-element-enhancer/5 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
          <strong>Not drawn.</strong> {data.why_not_drawable}
        </div>
      ) : (
        <>
          {data.reads && (
            <p className="mb-2 text-[11px] text-ink-500">
              The plot reads: <span className="text-ink-700">{data.reads}</span>
            </p>
          )}

          <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
            {data.examples.map((ex) => (
              <div key={ex.gene_id}>
                <div className="mb-1 flex items-baseline justify-between">
                  <button
                    onClick={() => onPick?.(ex.gene_symbol)}
                    className="text-[12px] font-medium hover:underline"
                  >
                    {ex.gene_symbol}
                  </button>
                  <span className="font-mono text-[10px] text-ink-500">
                    {ex.role === 'high' ? 'highest' : 'lowest'} · {ex.percentile.toFixed(0)}%
                  </span>
                </div>
                <MiniProfile
                  values={ex.profile.values}
                  startBp={ex.profile.start_bp}
                  endBp={ex.profile.end_bp}
                  geom={data}
                  color={ex.role === 'high' ? '#2b5070' : '#8badc9'}
                />
              </div>
            ))}
          </div>

          <p className="mt-2 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            {data.caveat}
          </p>
        </>
      )}
    </div>
  )
}
