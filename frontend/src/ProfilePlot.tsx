/**
 * Viewpoint-aligned contact profile, the app's primary view.
 *
 * This is NOT a genome browser. The x-axis is signed distance from the
 * experimental viewpoint, which sits at the window centre, and that alignment
 * is exactly what makes genes comparable to each other. UCSC and IGV are
 * coordinate-anchored and structurally cannot draw this.
 *
 * Canvas rather than SVG: at level 0 a window holds 40,000 points, which is far
 * past the point where one DOM node per point stays responsive.
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { formatBp, type Peak, type Profile } from './api'

const BAND_FILL: Record<string, string> = {
  viewpoint_proximal: 'rgba(43, 80, 112, 0.10)',
  local: 'rgba(43, 80, 112, 0.06)',
  distal: 'rgba(43, 80, 112, 0.035)',
  far_distal: 'rgba(43, 80, 112, 0.015)',
}

const BAND_LABEL: Record<string, string> = {
  viewpoint_proximal: '0–10 kb',
  local: '10–50 kb',
  distal: '50–250 kb',
  far_distal: '250 kb–1 Mb',
}

const PAD = { top: 14, right: 16, bottom: 30, left: 56 }
/** Dedicated lane for peak markers, below the signal area. */
const PEAK_LANE = 16

interface Props {
  profile: Profile
  peaks?: Peak[]
  onPeakClick?: (p: Peak) => void
  /** Called with the new window when the user brushes; null resets to full. */
  onZoom: (range: { start: number; end: number } | null) => void
  loading?: boolean
  height?: number
}

const ELEMENT_COLOR: Record<string, string> = {
  enhancer: '#c2703d',
  ctcf: '#4a7c59',
  promoter: '#2b5070',
}

export function ProfilePlot({
  profile,
  peaks = [],
  onPeakClick,
  onZoom,
  loading,
  height = 300,
}: Props) {
  const totalHeight = height + (peaks.length ? PEAK_LANE : 0)
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const wrapRef = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(900)
  const [drag, setDrag] = useState<{ x0: number; x1: number } | null>(null)
  const [hover, setHover] = useState<{ x: number; bp: number; value: number } | null>(null)

  // Track container width so the plot is responsive without a layout library.
  useEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver(([entry]) => setWidth(entry.contentRect.width))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  const plotW = Math.max(width - PAD.left - PAD.right, 10)
  const plotH = height - PAD.top - PAD.bottom

  const { values, start_bp, end_bp } = profile
  const spanBp = end_bp - start_bp

  const yMax = useMemo(() => {
    let m = 0
    for (const v of values) if (v > m) m = v
    return m > 0 ? m : 1
  }, [values])

  const bpToX = useCallback(
    (bp: number) => PAD.left + ((bp - start_bp) / spanBp) * plotW,
    [start_bp, spanBp, plotW],
  )
  const xToBp = useCallback(
    (x: number) => start_bp + ((x - PAD.left) / plotW) * spanBp,
    [start_bp, spanBp, plotW],
  )

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const dpr = window.devicePixelRatio || 1
    canvas.width = width * dpr
    canvas.height = totalHeight * dpr
    const ctx = canvas.getContext('2d')
    if (!ctx) return
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, width, totalHeight)

    // ---- distance bands, mirrored either side of the viewpoint -------------
    for (const band of profile.bands) {
      ctx.fillStyle = BAND_FILL[band.name] ?? 'rgba(0,0,0,0.02)'
      for (const sign of [1, -1]) {
        const a = sign > 0 ? band.start_bp : -band.end_bp
        const b = sign > 0 ? band.end_bp : -band.start_bp
        const x0 = Math.max(bpToX(a), PAD.left)
        const x1 = Math.min(bpToX(b), PAD.left + plotW)
        if (x1 > x0) ctx.fillRect(x0, PAD.top, x1 - x0, plotH)
      }
    }

    // ---- axes --------------------------------------------------------------
    ctx.strokeStyle = '#c3d6e4'
    ctx.lineWidth = 1
    ctx.beginPath()
    ctx.moveTo(PAD.left, PAD.top)
    ctx.lineTo(PAD.left, PAD.top + plotH)
    ctx.lineTo(PAD.left + plotW, PAD.top + plotH)
    ctx.stroke()
    if (peaks.length) {
      // left edge extended down past the peak lane, so the lane is enclosed
      ctx.beginPath()
      ctx.moveTo(PAD.left, PAD.top + plotH)
      ctx.lineTo(PAD.left, PAD.top + plotH + PEAK_LANE)
      ctx.stroke()
    }

    // y ticks
    ctx.fillStyle = '#5b89ae'
    ctx.font = '10px ui-monospace, monospace'
    ctx.textAlign = 'right'
    for (let i = 0; i <= 2; i++) {
      const frac = i / 2
      const y = PAD.top + plotH - frac * plotH
      const v = frac * yMax
      ctx.fillText(v >= 100 ? v.toFixed(0) : v.toFixed(v >= 10 ? 1 : 2), PAD.left - 6, y + 3)
      if (i > 0) {
        ctx.strokeStyle = '#e4edf4'
        ctx.beginPath()
        ctx.moveTo(PAD.left, y)
        ctx.lineTo(PAD.left + plotW, y)
        ctx.stroke()
      }
    }

    // ---- the profile -------------------------------------------------------
    // One vertical extent per pixel column: with up to 40k points over ~900px,
    // drawing every point would alias badly and hide real peaks.
    const colMax = new Float64Array(Math.ceil(plotW))
    const n = values.length
    for (let i = 0; i < n; i++) {
      const col = Math.min(Math.floor((i / n) * plotW), colMax.length - 1)
      if (values[i] > colMax[col]) colMax[col] = values[i]
    }

    ctx.fillStyle = profile.mode === 'oe' ? '#2b5070' : '#3d6b91'
    for (let col = 0; col < colMax.length; col++) {
      const h = (colMax[col] / yMax) * plotH
      if (h > 0) ctx.fillRect(PAD.left + col, PAD.top + plotH - h, 1, h)
    }

    // ---- peaks, in their own lane -------------------------------------------
    // Drawn under the signal rather than over it. Overlaying 20-40 markers on
    // the trace hid the data they were annotating, which is the same mistake
    // the J1 figures make.
    if (peaks.length) {
      const laneY = PAD.top + plotH + PEAK_LANE / 2
      ctx.strokeStyle = '#e4edf4'
      ctx.beginPath()
      ctx.moveTo(PAD.left, laneY)
      ctx.lineTo(PAD.left + plotW, laneY)
      ctx.stroke()

      const pmax = Math.max(...peaks.map((p) => p.peak_max), 1)
      for (const pk of peaks) {
        const x = bpToX(pk.offset_bp)
        if (x < PAD.left - 4 || x > PAD.left + plotW + 4) continue
        const r = 1.8 + 3.2 * Math.sqrt(pk.peak_max / pmax)
        // A faint stem ties each marker to the position it annotates.
        ctx.strokeStyle = (ELEMENT_COLOR[pk.re] ?? '#888') + '44'
        ctx.beginPath()
        ctx.moveTo(x, PAD.top + plotH)
        ctx.lineTo(x, laneY - r)
        ctx.stroke()
        ctx.beginPath()
        ctx.arc(x, laneY, r, 0, Math.PI * 2)
        ctx.fillStyle = (ELEMENT_COLOR[pk.re] ?? '#888') + 'dd'
        ctx.fill()
      }
    }

    // ---- viewpoint marker --------------------------------------------------
    const vpX = bpToX(0)
    if (vpX >= PAD.left && vpX <= PAD.left + plotW) {
      ctx.strokeStyle = '#c2703d'
      ctx.setLineDash([3, 3])
      ctx.beginPath()
      ctx.moveTo(vpX, PAD.top)
      ctx.lineTo(vpX, PAD.top + plotH)
      ctx.stroke()
      ctx.setLineDash([])
    }

    // ---- x ticks -----------------------------------------------------------
    ctx.fillStyle = '#5b89ae'
    ctx.textAlign = 'center'
    ctx.font = '10px ui-monospace, monospace'
    const nTicks = 7
    for (let i = 0; i <= nTicks; i++) {
      const bp = start_bp + (i / nTicks) * spanBp
      ctx.fillText(
        formatBp(Math.round(bp)),
        bpToX(bp),
        PAD.top + plotH + (peaks.length ? PEAK_LANE : 0) + 14,
      )
    }

    // ---- brush selection ---------------------------------------------------
    if (drag) {
      const x0 = Math.min(drag.x0, drag.x1)
      const x1 = Math.max(drag.x0, drag.x1)
      ctx.fillStyle = 'rgba(43, 80, 112, 0.15)'
      ctx.fillRect(x0, PAD.top, x1 - x0, plotH)
      ctx.strokeStyle = '#2b5070'
      ctx.strokeRect(x0, PAD.top, x1 - x0, plotH)
    }

    // ---- hover crosshair ---------------------------------------------------
    if (hover && !drag) {
      ctx.strokeStyle = 'rgba(15, 30, 46, 0.35)'
      ctx.beginPath()
      ctx.moveTo(hover.x, PAD.top)
      ctx.lineTo(hover.x, PAD.top + plotH)
      ctx.stroke()
    }
  }, [values, yMax, width, height, totalHeight, plotW, plotH, bpToX, drag, hover, profile, peaks, start_bp, spanBp])

  function localX(e: React.MouseEvent) {
    const rect = canvasRef.current!.getBoundingClientRect()
    return Math.min(Math.max(e.clientX - rect.left, PAD.left), PAD.left + plotW)
  }

  function onMove(e: React.MouseEvent) {
    const x = localX(e)
    if (drag) {
      setDrag({ ...drag, x1: x })
      return
    }
    const bp = xToBp(x)
    const idx = Math.min(
      values.length - 1,
      Math.max(0, Math.floor(((bp - start_bp) / spanBp) * values.length)),
    )
    setHover({ x, bp, value: values[idx] })
  }

  function nearestPeak(x: number): Peak | null {
    if (!peaks.length) return null
    let best: Peak | null = null
    let bestDist = 12 // px tolerance
    for (const pk of peaks) {
      const d = Math.abs(bpToX(pk.offset_bp) - x)
      if (d < bestDist) {
        bestDist = d
        best = pk
      }
    }
    return best
  }

  function onUp() {
    if (drag) {
      const a = xToBp(Math.min(drag.x0, drag.x1))
      const b = xToBp(Math.max(drag.x0, drag.x1))
      // Ignore a click-without-drag, and refuse to zoom below one native bin.
      if (Math.abs(b - a) > 500) {
        onZoom({ start: Math.round(a), end: Math.round(b) })
      } else if (onPeakClick) {
        // A click without a drag selects the nearest peak, if one is close.
        const pk = nearestPeak(drag.x1)
        if (pk) onPeakClick(pk)
      }
      setDrag(null)
    }
  }

  const zoomed = spanBp < 2_000_000

  return (
    <div ref={wrapRef} className="relative w-full">
      <canvas
        ref={canvasRef}
        style={{ width: '100%', height: totalHeight }}
        className={`select-none ${drag ? 'cursor-col-resize' : 'cursor-crosshair'}`}
        onMouseDown={(e) => {
          const x = localX(e)
          setDrag({ x0: x, x1: x })
        }}
        onMouseMove={onMove}
        onMouseUp={onUp}
        onMouseLeave={() => {
          setHover(null)
          setDrag(null)
        }}
      />

      {hover && !drag && (
        <div
          className="pointer-events-none absolute rounded border border-ink-200 bg-white/95 px-2 py-1 font-mono text-[11px] text-ink-800 shadow-sm"
          style={{ left: Math.min(hover.x + 8, width - 130), top: PAD.top + 4 }}
        >
          {formatBp(Math.round(hover.bp))} · {hover.value.toFixed(2)}
        </div>
      )}

      <div className="mt-1 flex items-center justify-between text-[11px] text-ink-500">
        <span className="flex items-center gap-3">
          <span>
            {profile.n.toLocaleString()} bins @ {profile.bp_per_bin} bp
            <span className="ml-2 text-ink-400">L{profile.level}</span>
            {loading && <span className="ml-2 text-ink-400">loading…</span>}
          </span>
          {peaks.length > 0 && (
            <span className="flex items-center gap-2 border-l border-ink-200 pl-3">
              {(['enhancer', 'ctcf', 'promoter'] as const).map((cls) => {
                const n = peaks.filter((p) => p.re === cls).length
                if (!n) return null
                return (
                  <span key={cls} className="flex items-center gap-1">
                    <span
                      className="inline-block h-2 w-2 rounded-full"
                      style={{ background: ELEMENT_COLOR[cls] }}
                    />
                    {cls} {n}
                  </span>
                )
              })}
            </span>
          )}
        </span>
        <span className="flex items-center gap-3">
          {profile.bands.map((b) => (
            <span key={b.name} className="flex items-center gap-1">
              <span
                className="inline-block h-2 w-2 rounded-sm"
                style={{ background: BAND_FILL[b.name]?.replace(/[\d.]+\)$/, '0.5)') }}
              />
              {BAND_LABEL[b.name] ?? b.name}
            </span>
          ))}
          {zoomed && (
            <button
              onClick={() => onZoom(null)}
              className="rounded border border-ink-300 px-2 py-0.5 text-ink-700 hover:bg-ink-50"
            >
              reset zoom
            </button>
          )}
        </span>
      </div>
    </div>
  )
}
