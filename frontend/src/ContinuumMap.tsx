/**
 * All 1,846 genes as a scatter, with the selected gene highlighted.
 *
 * PCA is the default because it is linear and preserves distances, so a
 * continuum renders as a continuum. UMAP is offered but paired with a
 * per-feature-permuted null: if the null looks equally grouped, the grouping is
 * the method rather than the biology. That comparison is the point of including
 * UMAP at all.
 */

import { useEffect, useRef, useState } from 'react'
import type { Embedding } from './api'

const GROUP_COLOR: Record<string, string> = {
  'arch-HK': '#3d6b91',
  'arch-ME-constitutive': '#2b5070',
  'arch-ME-effector': '#c2703d',
  'arch-sparse': '#8badc9',
  'arch-off': '#c3d6e4',
}

const GROUP_DISPLAY: Record<string, string> = {
  'arch-HK': 'dispersed',
  'arch-ME-constitutive': 'promoter-local',
  'arch-ME-effector': 'enhancer-focal',
  'arch-sparse': 'sparse',
  'arch-off': 'empty (QC)',
}

const PAD = { top: 12, right: 12, bottom: 46, left: 62 }

interface Props {
  data: Embedding
  /** The REQUESTED axes, not the loaded ones. The selects must reflect intent:
   *  reading them back off `data` meant that while a fetch was in flight the
   *  other select still held the previous response's axis, so a second change
   *  silently reverted the first. */
  axes: [string, string]
  onAxisChange: (x: string, y: string) => void
  onPick?: (symbol: string) => void
  loading?: boolean
  error?: string | null
  height?: number
}

export function ContinuumMap({
  data,
  axes,
  onAxisChange,
  onPick,
  loading,
  error,
  height = 360,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const wrapRef = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(600)
  const [hover, setHover] = useState<{ x: number; y: number; label: string } | null>(null)
  const [dimByConfidence, setDimByConfidence] = useState(true)

  useEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver(([e]) => setWidth(e.contentRect.width))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  const plotW = Math.max(width - PAD.left - PAD.right, 10)
  const plotH = height - PAD.top - PAD.bottom

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const dpr = window.devicePixelRatio || 1
    canvas.width = width * dpr
    canvas.height = height * dpr
    const ctx = canvas.getContext('2d')
    if (!ctx) return
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, width, height)

    const pts = data.points.filter((p) => Number.isFinite(p.x) && Number.isFinite(p.y))
    if (!pts.length) return
    const xs = pts.map((p) => p.x)
    const ys = pts.map((p) => p.y)
    const x0 = Math.min(...xs)
    const x1 = Math.max(...xs)
    const y0 = Math.min(...ys)
    const y1 = Math.max(...ys)
    const sx = (v: number) => PAD.left + ((v - x0) / (x1 - x0 || 1)) * plotW
    const sy = (v: number) => PAD.top + plotH - ((v - y0) / (y1 - y0 || 1)) * plotH

    ctx.strokeStyle = '#e4edf4'
    ctx.beginPath()
    ctx.moveTo(PAD.left, PAD.top)
    ctx.lineTo(PAD.left, PAD.top + plotH)
    ctx.lineTo(PAD.left + plotW, PAD.top + plotH)
    ctx.stroke()

    for (const p of pts) {
      // Opacity carries assignment confidence, so the 44% mixture body reads as
      // the diffuse interior it is instead of looking like crisp membership.
      const a = dimByConfidence ? 0.15 + 0.5 * (p.max_posterior ?? 0.5) : 0.5
      ctx.beginPath()
      ctx.arc(sx(p.x), sy(p.y), 2.1, 0, Math.PI * 2)
      ctx.fillStyle = (GROUP_COLOR[p.group ?? ''] ?? '#8badc9') + Math.round(a * 255).toString(16).padStart(2, '0')
      ctx.fill()
    }

    const hlPoint = pts.find((p) => p.gene_id === data.highlight)
    if (hlPoint) {
      const px = sx(hlPoint.x)
      const py = sy(hlPoint.y)
      ctx.beginPath()
      ctx.arc(px, py, 7, 0, Math.PI * 2)
      ctx.strokeStyle = '#0f1e2e'
      ctx.lineWidth = 2
      ctx.stroke()
      ctx.beginPath()
      ctx.arc(px, py, 3.5, 0, Math.PI * 2)
      ctx.fillStyle = '#0f1e2e'
      ctx.fill()
      ctx.font = '11px ui-sans-serif, system-ui'
      ctx.fillStyle = '#0f1e2e'
      ctx.textAlign = px > PAD.left + plotW * 0.75 ? 'right' : 'left'
      ctx.fillText(hlPoint.gene_symbol, px + (ctx.textAlign === 'right' ? -11 : 11), py + 4)
    }

    // ---- axis labels, with the poles named -------------------------------
    // The label gives the contrast; the poles say which end you are looking at.
    ctx.font = '10px ui-monospace, monospace'
    ctx.fillStyle = '#5b89ae'
    ctx.textAlign = 'center'
    ctx.fillText(data.x_axis.label, PAD.left + plotW / 2, height - 6)

    const xp = data.x_axis.poles
    if (xp) {
      ctx.fillStyle = '#8badc9'
      ctx.textAlign = 'left'
      ctx.fillText(`\u2190 ${xp.neg}`, PAD.left + 2, PAD.top + plotH + 15)
      ctx.textAlign = 'right'
      ctx.fillText(`${xp.pos} \u2192`, PAD.left + plotW - 2, PAD.top + plotH + 15)
    }

    ctx.save()
    ctx.translate(11, PAD.top + plotH / 2)
    ctx.rotate(-Math.PI / 2)
    ctx.fillStyle = '#5b89ae'
    ctx.textAlign = 'center'
    ctx.fillText(data.y_axis.label, 0, 0)
    ctx.restore()

    const yp = data.y_axis.poles
    if (yp) {
      ctx.fillStyle = '#8badc9'
      ctx.save()
      ctx.translate(PAD.left - 6, PAD.top + plotH)
      ctx.rotate(-Math.PI / 2)
      ctx.textAlign = 'left'
      ctx.fillText(`\u2190 ${yp.neg}`, 0, 0)
      ctx.restore()
      ctx.save()
      ctx.translate(PAD.left - 6, PAD.top)
      ctx.rotate(-Math.PI / 2)
      ctx.textAlign = 'right'
      ctx.fillText(`${yp.pos} \u2192`, 0, 0)
      ctx.restore()
    }
  }, [data, width, height, plotW, plotH, dimByConfidence])

  function onMove(e: React.MouseEvent) {
    const rect = canvasRef.current!.getBoundingClientRect()
    const mx = e.clientX - rect.left
    const my = e.clientY - rect.top
    const pts = data.points
    const xs = pts.map((p) => p.x)
    const ys = pts.map((p) => p.y)
    const x0 = Math.min(...xs)
    const x1 = Math.max(...xs)
    const y0 = Math.min(...ys)
    const y1 = Math.max(...ys)
    let best: (typeof pts)[number] | null = null
    let bd = 10
    for (const p of pts) {
      const px = PAD.left + ((p.x - x0) / (x1 - x0 || 1)) * plotW
      const py = PAD.top + plotH - ((p.y - y0) / (y1 - y0 || 1)) * plotH
      const d = Math.hypot(px - mx, py - my)
      if (d < bd) {
        bd = d
        best = p
      }
    }
    setHover(best ? { x: mx, y: my, label: best.gene_symbol } : null)
  }

  const axisOptions = data.axes_available

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex flex-wrap items-baseline justify-between gap-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Continuum map · {data.n.toLocaleString()} genes
          {loading && <span className="ml-2 font-normal text-ink-400">updating…</span>}
          {error && (
            <span className="ml-2 font-normal text-element-enhancer">{error}</span>
          )}
        </h2>
        <div className="flex items-center gap-2 text-[11px]">
          <select
            value={axes[0]}
            onChange={(e) => onAxisChange(e.target.value, axes[1])}
            className="rounded border border-ink-200 px-1.5 py-0.5"
          >
            {axisOptions.map((a) => (
              <option key={a.key} value={a.key}>{a.label}</option>
            ))}
          </select>
          <span className="text-ink-400">×</span>
          <select
            value={axes[1]}
            onChange={(e) => onAxisChange(axes[0], e.target.value)}
            className="rounded border border-ink-200 px-1.5 py-0.5"
          >
            {axisOptions.map((a) => (
              <option key={a.key} value={a.key}>{a.label}</option>
            ))}
          </select>
        </div>
      </div>

      {data.is_umap && (
        <div
          className={`mb-3 rounded border px-3 py-2 text-[11px] leading-relaxed ${
            data.is_null
              ? 'border-ink-300 bg-ink-50 text-ink-700'
              : 'border-element-enhancer/40 bg-element-enhancer/5 text-ink-700'
          }`}
        >
          {data.is_null ? (
            <>
              <strong>This is the permuted null.</strong> Every feature keeps its
              marginal distribution; all joint structure is destroyed. There is provably
              no structure here, so any grouping you see is what UMAP does to
              structureless data.
            </>
          ) : (
            data.umap_caveat
          )}
          {!data.is_null && (
            <button
              onClick={() => onAxisChange('umap_null_1', 'umap_null_2')}
              className="ml-1 underline underline-offset-2"
            >
              Show the null →
            </button>
          )}
          {data.is_null && (
            <button
              onClick={() => onAxisChange('umap_1', 'umap_2')}
              className="ml-1 underline underline-offset-2"
            >
              ← Back to the real data
            </button>
          )}
        </div>
      )}

      <div ref={wrapRef} className="relative w-full">
        <canvas
          ref={canvasRef}
          style={{ width: '100%', height }}
          className="cursor-crosshair"
          onMouseMove={onMove}
          onMouseLeave={() => setHover(null)}
          onClick={() => hover && onPick?.(hover.label)}
        />
        {hover && (
          <div
            className="pointer-events-none absolute rounded border border-ink-200 bg-white/95 px-2 py-0.5 text-[11px] shadow-sm"
            style={{ left: Math.min(hover.x + 8, width - 90), top: hover.y - 22 }}
          >
            {hover.label}
          </div>
        )}
      </div>

      <div className="mt-2 flex flex-wrap items-center justify-between gap-2 text-[11px] text-ink-500">
        <span className="flex flex-wrap items-center gap-3">
          {Object.entries(GROUP_DISPLAY).map(([k, label]) => (
            <span key={k} className="flex items-center gap-1">
              <span
                className="inline-block h-2 w-2 rounded-full"
                style={{ background: GROUP_COLOR[k] }}
              />
              {label}
            </span>
          ))}
        </span>
        <label className="flex items-center gap-1.5">
          <input
            type="checkbox"
            checked={dimByConfidence}
            onChange={(e) => setDimByConfidence(e.target.checked)}
          />
          opacity = confidence
        </label>
      </div>
    </div>
  )
}
