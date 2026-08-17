import { inSpace, type Space } from './space'
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
import { api, type Embedding, type PlaneLoadings } from './api'

// Colour is the CURRENT taxonomy (region), not the superseded archetypes.
// Far-reaching regions take the saturated tone of their element class and
// mid-range the muted one, so the two levels are legible at once: hue is
// composition, lightness is reach.
const GROUP_COLOR: Record<string, string> = {
  'extended-ctcf': '#4a7c59',
  'extended-enhancer': '#c2703d',
  'extended-promoter': '#2b5070',
  'contained-ctcf': '#7fa98b',
  'contained-enhancer': '#dba57f',
  'contained-promoter': '#6d94b8',
}

// Legend for the CURRENT taxonomy. The old entries (dispersed / promoter-local
// / enhancer-focal / sparse / empty) named the superseded archetypes and would
// have kept a retired vocabulary on screen under new colours.
const GROUP_DISPLAY: Record<string, string> = {
  'contained-promoter': 'contained, promoter',
  'contained-enhancer': 'contained, enhancer',
  'extended-promoter': 'extended, promoter',
  'extended-enhancer': 'extended, enhancer',
  'extended-ctcf': 'extended, CTCF',
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
  /** Page-level coordinate space; the axis menus are filtered to it. */
  space?: Space
  onPick?: (symbol: string) => void
  loading?: boolean
  error?: string | null
  height?: number
  /** Extra genes to mark alongside `data.highlight`, for the comparison view.
   *  Matched on gene_symbol because that is what the compare UI carries. */
  compare?: string[]
  /** Shown above the plot when the map is embedded in the gene page, where the
   *  question is "where does THIS gene sit" rather than "what is the panel". */
  subtitle?: string
}

/** Comparison marker colours. Distinct from the group palette on purpose: a
 *  compared gene is not a category, and reusing a group colour would read as
 *  one. */
const COMPARE_COLORS = ['#c2703d', '#4a7c59', '#7a5c9e', '#a8452f']

export function ContinuumMap({
  data,
  axes,
  onAxisChange,
  space = 'corrected',
  onPick,
  loading,
  error,
  height = 360,
  compare = [],
  subtitle,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const wrapRef = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(600)
  const [hover, setHover] = useState<{ x: number; y: number; label: string } | null>(null)
  const [dimByConfidence, setDimByConfidence] = useState(true)
  const [showLoadings, setShowLoadings] = useState(false)
  const [vecs, setVecs] = useState<PlaneLoadings | null>(null)

  // Fetched per plane rather than once, because the vectors are a property of
  // the pair of axes on screen, not of the dataset.
  useEffect(() => {
    if (!showLoadings) return
    let live = true
    api
      .planeLoadings(axes[0], axes[1], 8)
      .then((r) => live && setVecs(r))
      .catch(() => live && setVecs(null))
    return () => {
      live = false
    }
  }, [showLoadings, axes[0], axes[1]])

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
      ctx.fillStyle = (GROUP_COLOR[p.region ?? p.group ?? ''] ?? '#8badc9') + Math.round(a * 255).toString(16).padStart(2, '0')
      ctx.fill()
    }

    // ---- comparison genes -------------------------------------------------
    // Drawn before the primary highlight so the searched gene stays on top.
    // A connecting line makes the separation between compared genes readable
    // as a distance, which is the whole point of comparing on a map.
    const cmp = compare
      .map((sym, i) => ({
        p: pts.find((q) => q.gene_symbol === sym),
        color: COMPARE_COLORS[i % COMPARE_COLORS.length],
      }))
      .filter((c): c is { p: (typeof pts)[number]; color: string } => !!c.p)

    const hlPoint = pts.find((p) => p.gene_id === data.highlight)

    if (hlPoint && cmp.length) {
      ctx.strokeStyle = '#8badc9'
      ctx.lineWidth = 1
      ctx.setLineDash([3, 3])
      for (const c of cmp) {
        ctx.beginPath()
        ctx.moveTo(sx(hlPoint.x), sy(hlPoint.y))
        ctx.lineTo(sx(c.p.x), sy(c.p.y))
        ctx.stroke()
      }
      ctx.setLineDash([])
    }

    for (const c of cmp) {
      const px = sx(c.p.x)
      const py = sy(c.p.y)
      ctx.beginPath()
      ctx.arc(px, py, 6, 0, Math.PI * 2)
      ctx.strokeStyle = c.color
      ctx.lineWidth = 2
      ctx.stroke()
      ctx.beginPath()
      ctx.arc(px, py, 3, 0, Math.PI * 2)
      ctx.fillStyle = c.color
      ctx.fill()
      ctx.font = '11px ui-sans-serif, system-ui'
      ctx.fillStyle = c.color
      ctx.textAlign = px > PAD.left + plotW * 0.75 ? 'right' : 'left'
      ctx.fillText(c.p.gene_symbol, px + (ctx.textAlign === 'right' ? -10 : 10), py + 4)
    }

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

    // ---- loading vectors --------------------------------------------------
    // Drawn from the centre of the plot and scaled to a fixed fraction of it.
    // Loadings and scores live in different units, so arrow length here is
    // relative between arrows and carries no absolute meaning against the
    // point cloud. Said in the caption rather than implied by the drawing.
    if (showLoadings && vecs?.available && vecs.vectors?.length) {
      const cx = PAD.left + plotW / 2
      const cy = PAD.top + plotH / 2
      const maxLen = Math.max(...vecs.vectors.map((v) => v.length))
      const reach = Math.min(plotW, plotH) * 0.42

      for (const v of vecs.vectors) {
        // Screen y is inverted relative to data y, so the vertical component
        // is negated. Without this every arrow points at the wrong quadrant.
        const ex = cx + (v.x / maxLen) * reach
        const ey = cy - (v.y / maxLen) * reach

        ctx.strokeStyle = '#0f1e2e'
        ctx.globalAlpha = 0.55
        ctx.lineWidth = 1.4
        ctx.beginPath()
        ctx.moveTo(cx, cy)
        ctx.lineTo(ex, ey)
        ctx.stroke()

        const ang = Math.atan2(ey - cy, ex - cx)
        ctx.beginPath()
        ctx.moveTo(ex, ey)
        ctx.lineTo(ex - 6 * Math.cos(ang - 0.4), ey - 6 * Math.sin(ang - 0.4))
        ctx.lineTo(ex - 6 * Math.cos(ang + 0.4), ey - 6 * Math.sin(ang + 0.4))
        ctx.closePath()
        ctx.fillStyle = '#0f1e2e'
        ctx.fill()

        ctx.globalAlpha = 0.9
        ctx.font = '9px ui-monospace, monospace'
        ctx.fillStyle = '#2b5070'
        ctx.textAlign = ex >= cx ? 'left' : 'right'
        // Label sits just beyond the head, pushed clear of the arrow line.
        ctx.fillText(v.feature, ex + (ex >= cx ? 4 : -4), ey + (ey >= cy ? 10 : -4))
        ctx.globalAlpha = 1
      }
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
  }, [data, width, height, plotW, plotH, dimByConfidence, compare, showLoadings, vecs])

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

  // Only offer axes from the active space. Listing all 68 (both spaces) would
  // let a viewer cross sPC1 with PC3, which is meaningless: the numbering does
  // not correspond between spaces and the two are different decompositions.
  // UMAP axes survive the filter because they belong to neither and are
  // labelled separately.
  // Only axes from the active space, and GROUPED, because the raw list is 30
  // components per space and a flat menu of 30 is unusable. Only the first five
  // carry interpretations; past that a component is a number, and offering them
  // with equal prominence invites reading meaning into an axis that has none.
  const inThisSpace = data.axes_available.filter((a: { key: string }) =>
    inSpace(a.key, space),
  )
  const num = (k: string) => Number(k.replace(/^s?pc/, ''))
  const isComp = (k: string) => /^s?pc\d+$/.test(k)
  const AXIS_GROUPS: { label: string; opts: typeof inThisSpace }[] = [
    {
      label: 'Named axes',
      opts: inThisSpace.filter((a) => isComp(a.key) && num(a.key) <= 5),
    },
    {
      label: 'Higher components (unnamed)',
      opts: inThisSpace.filter((a) => isComp(a.key) && num(a.key) > 5),
    },
    {
      label: 'Nonlinear (see caveat)',
      opts: inThisSpace.filter((a) => a.key.startsWith('umap')),
    },
    {
      label: 'Other',
      opts: inThisSpace.filter((a) => !isComp(a.key) && !a.key.startsWith('umap')),
    },
  ].filter((g) => g.opts.length > 0)

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex flex-wrap items-baseline justify-between gap-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Continuum map · {data.n.toLocaleString()} genes
          {loading && <span className="ml-2 font-normal text-ink-400">updating…</span>}
          {error && (
            <span className="ml-2 font-normal text-element-enhancer">{error}</span>
          )}
          {subtitle && (
            <span className="ml-2 font-normal normal-case tracking-normal text-ink-400">
              {subtitle}
            </span>
          )}
        </h2>
        <div className="flex items-center gap-2 text-[11px]">
          <label
            className="flex items-center gap-1.5"
            title="Overlay the features that drive variation in this plane"
          >
            <input
              type="checkbox"
              checked={showLoadings}
              onChange={(e) => setShowLoadings(e.target.checked)}
            />
            loadings
          </label>
          <select
            value={axes[0]}
            onChange={(e) => {
              // Swap rather than duplicate; x == y is rejected by the API.
              const v = e.target.value
              onAxisChange(v, v === axes[1] ? axes[0] : axes[1])
            }}
            className="rounded border border-ink-200 px-1.5 py-0.5"
          >
            {AXIS_GROUPS.map((g) => (
              <optgroup key={g.label} label={g.label}>
                {g.opts.map((a) => (
                  <option key={a.key} value={a.key}>{a.label}</option>
                ))}
              </optgroup>
            ))}
          </select>
          <span className="text-ink-400">×</span>
          <select
            value={axes[1]}
            onChange={(e) => {
              const v = e.target.value
              onAxisChange(v === axes[0] ? axes[1] : axes[0], v)
            }}
            className="rounded border border-ink-200 px-1.5 py-0.5"
          >
            {AXIS_GROUPS.map((g) => (
              <optgroup key={g.label} label={g.label}>
                {g.opts.map((a) => (
                  <option key={a.key} value={a.key}>{a.label}</option>
                ))}
              </optgroup>
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

      {showLoadings && (
        <div className="mt-2 rounded border border-ink-100 bg-ink-50/60 px-3 py-2 text-[11px] leading-relaxed text-ink-600">
          {vecs?.available ? (
            <>
              <strong className="text-ink-800">
                Top {vecs.vectors?.length} of {vecs.n_features} features
              </strong>{' '}
              by loading length in this plane, so a feature that loads hard on one
              axis and not the other ranks below one that drives both. Each arrow
              points the way that feature increases. Arrow lengths are comparable
              with each other only: loadings and gene scores are in different
              units, so the arrows carry no scale against the point cloud. PCA
              sign is arbitrary, read the contrast between opposite arrows rather
              than the absolute orientation.
            </>
          ) : (
            <>{vecs?.reason ?? 'Loading vectors…'}</>
          )}
        </div>
      )}

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
          {compare.map((sym, i) => (
            <span key={sym} className="flex items-center gap-1 font-medium">
              <span
                className="inline-block h-2 w-2 rounded-full"
                style={{ background: COMPARE_COLORS[i % COMPARE_COLORS.length] }}
              />
              {sym}
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
