/**
 * Where each external reference set sits on the map, as small multiples.
 *
 * Legitimate because it is post-hoc: the embedding is built from MCC features
 * alone and these sets are painted on afterwards, so the overlap is an
 * observation about an architecture-derived map rather than something the map
 * was built to show. That is the same licence the super-enhancer comparison
 * has.
 *
 * Two readings, deliberately shown together, because they are different claims
 * and a panel can be strong on one and silent on the other. A PATCH is a cell
 * exceeding a label-permutation null, drawn outlined. A GRADIENT is a shift of
 * the whole set along an axis, printed as a number. Cell-wise testing cannot
 * see a smooth gradient, so without the second number several strongly
 * positioned sets read as nothing at all.
 *
 * Defaults to PCA, not UMAP. UMAP renders continuous data as apparent clusters,
 * so enrichments painted on it acquire crisp territories that are partly the
 * method. The UMAP and its permuted null are both selectable, which is the
 * point: if a territory survives on the null, it was never real.
 */

import { useEffect, useRef, useState } from 'react'
import { api, type EnrichmentGrid as Grid, type EnrichmentPanel } from './api'

/** Diverging, house palette. Steel for depleted, ochre for enriched. */
function colour(v: number | null, scale: number): string {
  if (v == null) return '#f4f7fa'
  const t = Math.max(-1, Math.min(1, v / (scale || 1)))
  if (t >= 0) {
    const a = 0.12 + 0.88 * t
    return `rgba(194, 112, 61, ${a.toFixed(3)})`
  }
  const a = 0.12 + 0.88 * -t
  return `rgba(43, 80, 112, ${a.toFixed(3)})`
}

function Panel({
  p,
  nx,
  ny,
  xKey,
  yKey,
  size,
  points,
  extent,
  showHeat,
}: {
  p: EnrichmentPanel
  nx: number
  ny: number
  xKey: string
  yKey: string
  size: number
  points: [number, number][]
  extent: { x0: number; x1: number; y0: number; y1: number }
  showHeat: boolean
}) {
  const ref = useRef<HTMLCanvasElement>(null)

  useEffect(() => {
    const c = ref.current
    if (!c) return
    const dpr = window.devicePixelRatio || 1
    c.width = size * dpr
    c.height = size * dpr
    const ctx = c.getContext('2d')
    if (!ctx) return
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, size, size)

    const { x0, x1, y0, y1 } = extent
    const sx = (v: number) => ((v - x0) / (x1 - x0 || 1)) * size
    // Screen y is inverted relative to data y.
    const sy = (v: number) => size - ((v - y0) / (y1 - y0 || 1)) * size

    const w = size / nx
    const h = size / ny
    // Scaled to this set's own null threshold, so colour means "relative to
    // what chance produces for a set of THIS size" and panels stay comparable.
    const scale = p.null_threshold || p.max_abs || 1

    if (showHeat) {
      for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
          const v = p.cells[i * ny + j]
          if (v == null) continue
          ctx.fillStyle = colour(v, scale)
          ctx.fillRect(i * w, size - (j + 1) * h, Math.ceil(w), Math.ceil(h))
        }
      }
    }

    // The map itself. Every gene faint, this set's members on top, so the tile
    // is recognisably the same cloud as the continuum map rather than an
    // abstraction of it.
    ctx.fillStyle = 'rgba(150, 170, 190, 0.30)'
    for (const [px, py] of points) {
      ctx.beginPath()
      ctx.arc(sx(px), sy(py), 0.9, 0, Math.PI * 2)
      ctx.fill()
    }
    ctx.fillStyle = '#c2703d'
    for (const idx of p.members) {
      const q = points[idx]
      if (!q) continue
      ctx.beginPath()
      ctx.arc(sx(q[0]), sy(q[1]), 1.5, 0, Math.PI * 2)
      ctx.fill()
    }

    // Cells beating the permutation null, outlined over the points. The
    // statistics annotate the picture; they do not replace it.
    if (p.null_threshold != null) {
      ctx.strokeStyle = '#0f1e2e'
      ctx.lineWidth = 1.25
      for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
          const v = p.cells[i * ny + j]
          if (v == null || Math.abs(v) <= p.null_threshold) continue
          ctx.strokeRect(i * w + 0.5, size - (j + 1) * h + 0.5, w - 1, h - 1)
        }
      }
    }
  }, [p, nx, ny, size, points, extent, showHeat])

  const sx = p.axis_shift?.[xKey] ?? 0
  const sy = p.axis_shift?.[yKey] ?? 0
  const strongest = Math.abs(sx) > Math.abs(sy) ? sx : sy
  const strongestAxis = Math.abs(sx) > Math.abs(sy) ? xKey : yKey
  const notable = Math.abs(strongest) >= 0.4 || (p.n_above_null ?? 0) > 0

  return (
    <div
      className={`rounded border p-2 ${
        p.is_super_enhancer
          ? 'border-element-enhancer/40'
          : p.is_positive_control
            ? 'border-ink-400 border-dashed'
            : 'border-ink-100'
      }`}
    >
      <div className="mb-1 flex items-baseline justify-between gap-1">
        <span className="truncate font-mono text-[10px] text-ink-700" title={p.group}>
          {p.group}
        </span>
        <span className="shrink-0 font-mono text-[9px] text-ink-400">{p.n}</span>
      </div>

      <canvas ref={ref} style={{ width: size, height: size }} className="block" />

      <div className="mt-1 space-y-0.5 text-[9px] leading-tight">
        {/* The headline number. Displacement across ALL above-noise dimensions,
            paired with the overlap it implies, because every set here clears
            significance and only the overlap says whether it is visible. */}
        {p.displacement && (
          <div className="text-ink-700">
            <span className="font-mono">
              d {p.displacement.cohens_d >= 0 ? '+' : ''}
              {p.displacement.cohens_d.toFixed(2)}
            </span>{' '}
            on {p.displacement.strongest_dim}
            {', '}
            <span className={p.displacement.overlap_pct > 80 ? 'text-ink-400' : 'text-ink-700'}>
              {p.displacement.overlap_pct.toFixed(0)}% overlap
            </span>
          </div>
        )}
        <div className={notable ? 'text-ink-500' : 'text-ink-400'}>
          {p.n_above_null ?? 0} patch{(p.n_above_null ?? 0) === 1 ? '' : 'es'}
          {' · '}
          <span className="font-mono">
            {strongest >= 0 ? '+' : ''}
            {strongest.toFixed(2)}
          </span>{' '}
          on {strongestAxis}
        </div>
        {p.stratification && p.stratification.pct_retained < 0.5 && (
          <div className="text-element-enhancer">
            {Math.round(p.stratification.pct_retained * 100)}% survives{' '}
            {p.stratification.stratifier}
          </div>
        )}
        {p.is_super_enhancer && (
          <div className="text-element-enhancer">comparison only, never an input</div>
        )}
        {p.is_positive_control && (
          <div className="text-ink-500">density positive control</div>
        )}
      </div>
    </div>
  )
}

export function EnrichmentGrid() {
  const [axes, setAxes] = useState<[string, string]>(['pc2', 'pc3'])
  const [data, setData] = useState<Grid | null>(null)
  const [busy, setBusy] = useState(false)
  const [err, setErr] = useState<string | null>(null)
  const [open, setOpen] = useState(false)
  // Off by default. The points are what makes a tile recognisable as the map;
  // the heat shading is a summary of them and, drawn underneath, it competes
  // with the thing it is summarising.
  const [showHeat, setShowHeat] = useState(false)
  const req = useRef(0)

  useEffect(() => {
    if (!open) return
    const token = ++req.current
    setBusy(true)
    setErr(null)
    api
      .enrichmentGrid(axes[0], axes[1], 12, 200)
      .then((d) => token === req.current && setData(d))
      .catch((e) => token === req.current && setErr(String((e as Error).message)))
      .finally(() => token === req.current && setBusy(false))
  }, [axes, open])

  const CANVASES: [string, string, string][] = [
    ['pc2', 'pc3', 'PCA (default)'],
    ['pc1', 'pc2', 'PCA, PC1 x PC2'],
    ['umap_1', 'umap_2', 'UMAP'],
    ['umap_null_1', 'umap_null_2', 'UMAP permuted null'],
  ]

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="flex flex-wrap items-baseline justify-between gap-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          External sets across the map
          {busy && <span className="ml-2 font-normal text-ink-400">computing…</span>}
          {err && <span className="ml-2 font-normal text-element-enhancer">{err}</span>}
        </h2>
        <button
          onClick={() => setOpen(!open)}
          className="rounded border border-ink-200 px-2 py-0.5 text-[11px] text-ink-600 hover:border-ink-400"
        >
          {open ? 'hide' : 'show'}
        </button>
      </div>

      <p className="mt-1 text-[11px] leading-relaxed text-ink-500">
        Post-hoc. The map is built from MCC features alone; these sets are painted on
        afterwards and never entered the coordinates. Runs a label-permutation null,
        so it takes a few seconds.
      </p>

      {open && (
        <>
          <div className="mt-3 flex flex-wrap items-center gap-2 text-[11px]">
            <span className="text-ink-500">canvas</span>
            {CANVASES.map(([x, y, label]) => {
              const active = axes[0] === x && axes[1] === y
              return (
                <button
                  key={label}
                  onClick={() => setAxes([x, y])}
                  className={`rounded border px-2 py-0.5 ${
                    active
                      ? 'border-ink-600 bg-ink-600 text-white'
                      : 'border-ink-200 text-ink-600 hover:border-ink-400'
                  }`}
                >
                  {label}
                </button>
              )
            })}
          </div>

          {data?.is_null && (
            <div className="mt-3 rounded border border-ink-300 bg-ink-50 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
              <strong>This is the permuted null.</strong> There is provably no joint
              structure here, so any territory you can see is what the embedding does
              to structureless data. Anything that survives here is not a finding.
            </div>
          )}
          {data?.is_umap && !data.is_null && (
            <div className="mt-3 rounded border border-element-enhancer/40 bg-element-enhancer/5 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
              UMAP renders continuous data as apparent clusters, so enrichments here
              acquire crisp territories that are partly the method rather than the
              biology. Compare against the null before believing a shape.
            </div>
          )}

          {data && (
            <>
              <div className="mt-3 flex flex-wrap items-baseline justify-between gap-2">
                <p className="max-w-3xl text-[11px] leading-relaxed text-ink-600">
                  <strong>Every tile is the same {data.x_axis.label} by{' '}
                  {data.y_axis.label} scatter as the map above</strong>, all{' '}
                  {data.n_genes.toLocaleString()} genes in pale grey, with that set's
                  members in ochre on top. Boxes mark cells where the local rate beats
                  the label-permutation null. Sorted by strength, patch or gradient,
                  whichever is larger.
                </p>
                <label className="flex shrink-0 items-center gap-1.5 text-[11px] text-ink-600">
                  <input
                    type="checkbox"
                    checked={showHeat}
                    onChange={(e) => setShowHeat(e.target.checked)}
                  />
                  shade cells by enrichment
                </label>
              </div>
              {showHeat && (
                <p className="mt-1 text-[10px] leading-relaxed text-ink-500">
                  Shading is log2 of the observed rate over the panel base rate,
                  scaled to each set's own permutation threshold so tiles stay
                  comparable. Ochre enriched, steel depleted, unshaded means fewer
                  than {data.min_n} genes in that cell.
                </p>
              )}

              <div className="mt-3 grid grid-cols-2 gap-2 sm:grid-cols-3 lg:grid-cols-5">
                {data.panels.map((p) => (
                  <Panel
                    key={p.group}
                    p={p}
                    nx={data.nx}
                    ny={data.ny}
                    xKey={data.x_axis.key}
                    yKey={data.y_axis.key}
                    size={112}
                    points={data.points}
                    extent={data.extent}
                    showHeat={showHeat}
                  />
                ))}
              </div>

              <div className="mt-3 space-y-1.5 border-t border-ink-100 pt-2 text-[10px] leading-relaxed text-ink-500">
                <p>
                  <strong className="text-ink-700">Two different claims.</strong>{' '}
                  {data.caveats.gradient}
                </p>
                {data.caveats.robustness && (
                  <p>
                    <strong className="text-ink-700">Which number to trust.</strong>{' '}
                    {data.caveats.robustness}
                  </p>
                )}
                <p>{data.caveats.null}</p>
                <p>
                  <strong className="text-ink-700">Not independent.</strong>{' '}
                  {data.caveats.atac} {data.caveats.density}
                </p>
              </div>
            </>
          )}
        </>
      )}
    </div>
  )
}
