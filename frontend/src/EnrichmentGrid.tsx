import { DEFAULT_PLANE, axisKey, planeInSpace, type Space } from './space'
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
import { api, type EnrichmentGrid as Grid, type EnrichmentPanel,
         type Loadings } from './api'

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

/** Member vs rest density along one axis.
 *
 * Drawn rather than summarised because the finding is that the two curves are
 * SHIFTED but almost entirely OVERLAPPING, and a reader given "d = -0.21" has
 * to be told that means 91% overlap, whereas a reader given the curves has
 * already seen it. Filled grey is the rest of the panel, outlined is the set;
 * solid and dashed rules mark the two means. Densities, not counts: the groups
 * differ by up to tenfold in n.
 */
/** The shared axis frame every tile is drawn in.
 *
 * Without this a tile is an unoriented cloud, so "shift +0.61 on sPC1" has
 * nothing visual to attach to. One labelled frame above the grid is enough,
 * because every tile uses identical extents.
 *
 * UMAP gets a warning instead of labels. UMAP axes have no interpretation:
 * there is no linear map back to features, distances between distant points are
 * not meaningful, and a pole label would be invented. That is the same reason
 * /api/embedding/loadings refuses to draw arrows over it.
 */
function AxisFrame({
  x,
  y,
  isUmap,
}: {
  x: Grid['x_axis']
  y: Grid['y_axis']
  isUmap: boolean
}) {
  // What actually DRIVES each axis. A pole label is an interpretation; these
  // are the features it was read off, so a viewer can check the name instead of
  // taking it on trust. Same reasoning as the loadings panel on the main map.
  // PER-AXIS loadings, not plane loadings. /api/embedding/loadings ranks by
  // length in the PLANE and returns only the top N, so a feature that dominates
  // sPC2 but is near zero on sPC1 gets cut before it is ever seen. That is how
  // this frame came to report n_peaks_within_100kb as sPC2's positive driver
  // when the true leaders are promoter_signal_fraction_raw +0.325 and
  // promoter_signal_fraction +0.318, and how the pole labels ended up inverted.
  const [load, setLoad] = useState<Record<string, Loadings | null>>({})
  useEffect(() => {
    if (isUmap) return
    let live = true
    for (const k of [x.key, y.key]) {
      const n = Number(k.replace(/^s?pc/, ''))
      if (!Number.isFinite(n)) continue
      api
        .loadings(n, 30, k.startsWith('spc') ? 'corrected' : 'raw')
        .then((d) => live && setLoad((s) => ({ ...s, [k]: d })))
        .catch(() => live && setLoad((s) => ({ ...s, [k]: null })))
    }
    return () => {
      live = false
    }
  }, [x.key, y.key, isUmap])

  const drivers = (axis: 'x' | 'y', sign: 1 | -1) => {
    const d = load[axis === 'x' ? x.key : y.key]
    return (d?.loadings ?? [])
      .filter((l) => Math.sign(l.loading) === sign)
      .sort((a, b) => Math.abs(b.loading) - Math.abs(a.loading))
      .slice(0, 3)
      .map((l) => l.feature)
  }

  const Drivers = ({ names }: { names: string[] }) =>
    names.length ? (
      <div className="font-mono text-[8px] leading-tight text-ink-400">
        {names.join(', ')}
      </div>
    ) : null

  if (isUmap) {
    return (
      <div className="mt-2 rounded border border-ink-200 bg-ink-50 p-2 text-[10px] text-ink-600">
        <strong className="text-ink-700">UMAP axes have no interpretation.</strong>{' '}
        There is no linear map from a UMAP coordinate back to the features, so
        neither axis can be labelled and a shift along one means nothing. Switch
        to the sPC axes to read displacement.
      </div>
    )
  }
  return (
    <div className="mt-2 flex items-stretch gap-2 rounded border border-ink-200 p-2">
      <div className="flex-1">
        <div className="font-mono text-[10px] text-ink-700">{y.label}</div>
        <div className="flex justify-between gap-4">
          <div className="min-w-0">
            <div className="text-[9px] text-ink-500">&darr; {y.poles?.neg}</div>
            <Drivers names={drivers('y', -1)} />
          </div>
          <div className="min-w-0 text-right">
            <div className="text-[9px] text-ink-500">{y.poles?.pos} &uarr;</div>
            <Drivers names={drivers('y', 1)} />
          </div>
        </div>

        <div className="my-1.5 h-px bg-ink-100" />

        <div className="font-mono text-[10px] text-ink-700">{x.label}</div>
        <div className="flex justify-between gap-4">
          <div className="min-w-0">
            <div className="text-[9px] text-ink-500">&larr; {x.poles?.neg}</div>
            <Drivers names={drivers('x', -1)} />
          </div>
          <div className="min-w-0 text-right">
            <div className="text-[9px] text-ink-500">{x.poles?.pos} &rarr;</div>
            <Drivers names={drivers('x', 1)} />
          </div>
        </div>
      </div>
      <div className="w-40 shrink-0 border-l border-ink-100 pl-2 text-[9px] leading-tight text-ink-500">
        A shift on the horizontal axis moves a set left or right in every tile;
        a shift on the vertical moves it up or down. The shifts are small, so
        the movement is a drift of the ochre points, never a separate cluster.
      </div>
    </div>
  )
}

function DistCurve({
  dist,
  width,
  height = 26,
}: {
  dist: NonNullable<EnrichmentPanel['distributions']>[string]
  width: number
  height?: number
}) {
  const n = dist.members.length
  const peak = Math.max(...dist.members, ...dist.rest) || 1
  const path = (vals: number[]) =>
    vals
      .map((v, i) => {
        const x = (i / (n - 1)) * width
        const y = height - (v / peak) * (height - 2) - 1
        return `${i === 0 ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`
      })
      .join(' ')
  const meanX = (m: number) => {
    const lo = dist.edges[0]
    const hi = dist.edges[dist.edges.length - 1]
    return ((m - lo) / (hi - lo || 1)) * width
  }
  return (
    <svg width={width} height={height} className="block">
      <path
        d={`${path(dist.rest)} L${width},${height} L0,${height} Z`}
        fill="#cbd5df"
        opacity={0.55}
      />
      <path d={path(dist.members)} fill="none" stroke="#2b5070" strokeWidth={1.3} />
      <line
        x1={meanX(dist.rest_mean)} x2={meanX(dist.rest_mean)}
        y1={0} y2={height} stroke="#8fa3b5" strokeWidth={1}
      />
      <line
        x1={meanX(dist.member_mean)} x2={meanX(dist.member_mean)}
        y1={0} y2={height} stroke="#2b5070" strokeWidth={1} strokeDasharray="2 2"
      />
    </svg>
  )
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

      {/* THE SHIFT IS THE RESULT, so it leads. A set displaced a fifth of a
          standard deviation across the whole map produces no hot cell anywhere,
          because the shift is spread over every cell -- which is why the
          heatmap below is empty for 20 of 21 sets and this bar is not. */}
      <div className="mb-1">
        <div className="relative h-3 rounded-sm bg-ink-50" style={{ width: size }}>
          <div className="absolute inset-y-0 left-1/2 w-px bg-ink-300" />
          <div
            className={`absolute inset-y-0 ${
              strongest >= 0 ? 'left-1/2' : ''
            } ${notable ? 'bg-ink-600' : 'bg-ink-300'}`}
            style={{
              width: `${Math.min(Math.abs(strongest) / 0.8, 1) * 50}%`,
              ...(strongest < 0
                ? { right: '50%' }
                : {}),
            }}
          />
        </div>
        <div className="font-mono text-[9px] leading-tight text-ink-500">
          shift {strongest >= 0 ? '+' : ''}{strongest.toFixed(2)} on {strongestAxis}
        </div>
      </div>

      {/* ALWAYS shown. Three layers live here and only one of them is empty:
          the log2 heat cells (behind `showHeat`, off by default, because 20 of
          21 sets have no cell above the null), and the SCATTER, which is the
          evidence for the whole result. Members drawn over the full cloud is
          what makes "displaced but not separated" visible; the shift bar above
          is the same fact with the proof removed. Hiding this was a mistake. */}
      <canvas ref={ref} style={{ width: size, height: size }} className="block" />

      {/* The distribution along x, directly under the map, so the shift and the
          overlap are read together rather than one being inferred from a
          number. */}
      {/* Draw the axis the SHIFT is on, not always x. Showing a shift measured
          on sPC2 against a distribution on sPC1 made every set look unmoved,
          which is a display bug masquerading as a result. */}
      {p.distributions?.[strongestAxis] && (
        <div className="mt-1">
          <DistCurve dist={p.distributions[strongestAxis]} width={size} />
          <div className="font-mono text-[8px] leading-none text-ink-400">
            {strongestAxis} alone ·{' '}
            {p.distributions[strongestAxis].overlap_pct.toFixed(0)}% overlap
          </div>
        </div>
      )}

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
            across all {' '}
            <span title="Displacement of the centroid across every component above the noise ceiling, not just the two on screen. This is a DIFFERENT quantity from the single-axis overlap above the line, and the two will not agree.">
              18 dims
            </span>
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

export function EnrichmentGrid({ space = 'corrected' }: { space?: Space }) {
  const [axes, setAxes] = useState<[string, string]>(DEFAULT_PLANE[space])

  // Follow the page-level space. Without this the grid could sit on corrected
  // axes while the rest of the page showed raw PCs.
  useEffect(() => {
    setAxes((cur) => planeInSpace(cur, space))
  }, [space])
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

  // Amount-corrected axes first and default. The raw PCs are kept for
  // provenance only: on the current substrate PC1 and PC3 both carry magnitude
  // (r 0.505 and 0.519), so a position on that plane partly means "more signal".
  // Planes are built from the page-level space rather than hardcoded, so this
  // panel cannot end up showing raw PCs while the rest of the page is on the
  // corrected axes. Labels differ per space because the axes mean different
  // things: sPC2 is enhancer-vs-promoter, PC2 is local-vs-long-range.
  const CANVASES: [string, string, string][] =
    space === 'corrected'
      ? [
          [axisKey(space, 1), axisKey(space, 2), 'reach x enhancer/promoter'],
          [axisKey(space, 2), axisKey(space, 3), 'enhancer/promoter x CTCF/promoter'],
          ['umap_1', 'umap_2', 'UMAP'],
          ['umap_null_1', 'umap_null_2', 'UMAP permuted null'],
        ]
      : [
          [axisKey(space, 2), axisKey(space, 3), 'PC2 x PC3 (historical default)'],
          [axisKey(space, 1), axisKey(space, 2), 'PC1 x PC2 (PC1 carries amount)'],
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
              <AxisFrame
                x={data.x_axis}
                y={data.y_axis}
                isUmap={data.is_umap}
              />

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

              {/* Sorted by how far the set has moved, so the panel reads as a
                  ranking. Unsorted, 21 near-identical tiles in arbitrary order
                  invite the eye to hunt for a hot spot that does not exist. */}
              <div className="mt-2 flex items-center gap-3 text-[10px] text-ink-500">
                <span>
                  sorted by displacement ·{' '}
                  {data.panels.filter((p) => (p.n_above_null ?? 0) > 0).length} of{' '}
                  {data.panels.length} sets have any cell above the permutation null
                </span>
              </div>

              <div className="mt-2 grid grid-cols-2 gap-2 sm:grid-cols-3 lg:grid-cols-5">
                {[...data.panels]
                  .sort((a, b) => {
                    const m = (q: typeof a) =>
                      Math.max(...Object.values(q.axis_shift ?? {}).map(Math.abs), 0)
                    return m(b) - m(a)
                  })
                  .map((p) => (
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
