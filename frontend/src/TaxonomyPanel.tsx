/**
 * The reach x composition taxonomy: how it was derived, and how the regions
 * differ.
 *
 * Two things this panel must do that a bare list of group names would not.
 *
 * It shows the DERIVATION. The old archetypes were five names with no visible
 * provenance, and it took until 2026-08-16 to notice they had been fit on a
 * substrate that no longer existed and named partly by overlap with external
 * gene sets. A taxonomy that cannot be interrogated gets trusted by default, so
 * the split-by-split trail and the k evidence are on screen, including the k
 * values that were REJECTED and why.
 *
 * It shows the regions as OVERLAPPING PROFILES, not as a partition. Every gene
 * in the panel has a top membership below 0.5, so these are named areas of a
 * continuum. The radar deliberately draws all regions on one set of axes so the
 * overlap is visible rather than implied.
 */

import { useEffect, useState } from 'react'
import { useRef } from 'react'
import { api, type Taxonomy, type TaxonomyMap } from './api'
import { DEFAULT_PLANE, type Space } from './space'

const COLOURS: Record<string, string> = {
  'extended-ctcf': '#4a7c59',
  'extended-enhancer': '#c2703d',
  'extended-promoter': '#2b5070',
  'contained-ctcf': '#7fa98b',
  'contained-enhancer': '#dba57f',
  'contained-promoter': '#6d94b8',
}

function Radar({ data, size = 300 }: { data: Taxonomy; size?: number }) {
  const axes = data.radar
  const regions = data.regions.map((r) => r.region)
  const [hover, setHover] = useState<string | null>(null)
  if (!axes.length) return null

  const cx = size / 2
  const cy = size / 2
  const R = size / 2 - 52
  const pt = (i: number, v: number) => {
    const a = (i / axes.length) * 2 * Math.PI - Math.PI / 2
    const r = (v / 100) * R
    return [cx + r * Math.cos(a), cy + r * Math.sin(a)]
  }

  return (
    <svg width={size} height={size} className="shrink-0">
      {[25, 50, 75, 100].map((g) => (
        <polygon
          key={g}
          points={axes.map((_, i) => pt(i, g).join(',')).join(' ')}
          fill="none"
          stroke={g === 50 ? '#b8c6d3' : '#e6ecf1'}
          strokeDasharray={g === 50 ? '3 3' : undefined}
        />
      ))}
      {axes.map((a, i) => {
        // Anchor by quadrant instead of always centring, and push the top and
        // bottom labels clear of the ring. With a fixed 112% radius and centred
        // anchoring, the top label sat exactly where the median caption is
        // drawn and the two overlapped.
        const [x, y] = pt(i, 114)
        const ang = (i / axes.length) * 2 * Math.PI - Math.PI / 2
        const cosA = Math.cos(ang)
        const anchor = cosA > 0.3 ? 'start' : cosA < -0.3 ? 'end' : 'middle'
        const nudge = Math.sin(ang) < -0.7 ? -6 : Math.sin(ang) > 0.7 ? 6 : 0
        return (
          <text
            key={a.axis}
            x={x}
            y={y + nudge}
            textAnchor={anchor}
            dominantBaseline="middle"
            className="fill-ink-500 text-[8px]"
          >
            {a.axis}
          </text>
        )
      })}
      {regions.map((r) => {
        const pts = axes.map((a, i) => pt(i, a.values[r] ?? 50).join(',')).join(' ')
        const dim = hover !== null && hover !== r
        return (
          <polygon
            key={r}
            points={pts}
            fill={COLOURS[r] ?? '#888'}
            fillOpacity={dim ? 0.03 : hover === r ? 0.28 : 0.1}
            stroke={COLOURS[r] ?? '#888'}
            strokeWidth={hover === r ? 2.2 : 1.3}
            strokeOpacity={dim ? 0.2 : 1}
            onMouseEnter={() => setHover(r)}
            onMouseLeave={() => setHover(null)}
            style={{ cursor: 'pointer' }}
          />
        )
      })}
      {/* 50 is the panel median: the ring, not the centre, is "average".
          Drawn at the very bottom of the viewBox, clear of the top axis label
          it used to sit on top of. */}
      <text x={cx} y={size - 3} textAnchor="middle" className="fill-ink-400 text-[7px]">
        dashed ring = panel median
      </text>
    </svg>
  )
}


/** The regions drawn on an actual embedding.
 *
 * Shown on TWO planes because the regions were defined in the 77-dimensional
 * corrected space, not on any projection. If they only look coherent on the
 * plane they happen to be plotted in, that is an artefact of the plane; seeing
 * them hold on a linear (sPC) and a nonlinear (UMAP) projection is a check.
 *
 * Opacity is the gene's top membership, so the faintness of most of the panel
 * is the continuum result rendered directly: 99.7% of genes sit below 0.5.
 */
function RegionScatter({
  initialX,
  initialY,
  colours,
  options,
  size = 300,
  onPick,
}: {
  initialX: string
  initialY: string
  colours: Record<string, string>
  /** Axis keys offered here. Kept short on purpose: the point of this pair of
   *  plots is to check the regions against a LINEAR and a NONLINEAR projection,
   *  not to browse all 60 components. */
  options: { key: string; label: string }[]
  size?: number
  onPick?: (s: string) => void
}) {
  const [x, setX] = useState(initialX)
  const [y, setY] = useState(initialY)
  const [d, setD] = useState<TaxonomyMap | null>(null)
  const ref = useRef<HTMLCanvasElement>(null)

  useEffect(() => {
    let live = true
    api.taxonomyMap(x, y).then((v) => live && setD(v)).catch(() => live && setD(null))
    return () => {
      live = false
    }
  }, [x, y])

  useEffect(() => {
    const c = ref.current
    if (!c || !d) return
    const dpr = window.devicePixelRatio || 1
    c.width = size * dpr
    c.height = size * dpr
    const ctx = c.getContext('2d')
    if (!ctx) return
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, size, size)
    const x0 = Math.min(...d.xs)
    const x1 = Math.max(...d.xs)
    const y0 = Math.min(...d.ys)
    const y1 = Math.max(...d.ys)
    const px = (v: number) => 4 + ((v - x0) / (x1 - x0 || 1)) * (size - 8)
    const py = (v: number) => size - 4 - ((v - y0) / (y1 - y0 || 1)) * (size - 8)
    for (let i = 0; i < d.n; i++) {
      const reg = d.regions[d.r[i]]
      ctx.globalAlpha = 0.18 + Math.min(d.w[i] / 0.5, 1) * 0.55
      ctx.fillStyle = colours[reg] ?? '#888'
      ctx.beginPath()
      ctx.arc(px(d.xs[i]), py(d.ys[i]), 1.6, 0, Math.PI * 2)
      ctx.fill()
    }
    ctx.globalAlpha = 1
  }, [d, size, colours])

  if (!d) return null
  return (
    <div>
      <canvas
        ref={ref}
        style={{ width: size, height: size }}
        className="block rounded border border-ink-100"
        onClick={(e) => {
          if (!onPick || !d) return
          const rect = e.currentTarget.getBoundingClientRect()
          const mx = e.clientX - rect.left
          const my = e.clientY - rect.top
          const x0 = Math.min(...d.xs), x1 = Math.max(...d.xs)
          const y0 = Math.min(...d.ys), y1 = Math.max(...d.ys)
          let best = -1, bd = 64
          for (let i = 0; i < d.n; i++) {
            const ex = 4 + ((d.xs[i] - x0) / (x1 - x0 || 1)) * (size - 8)
            const ey = size - 4 - ((d.ys[i] - y0) / (y1 - y0 || 1)) * (size - 8)
            const dd = (ex - mx) ** 2 + (ey - my) ** 2
            if (dd < bd) { bd = dd; best = i }
          }
          if (best >= 0) onPick(d.symbols[best])
        }}
      />
      <div className="mt-1 flex items-center gap-1 text-[10px]">
        <select
          value={x}
          onChange={(e) => {
            // Swap rather than duplicate. Picking the axis already on the other
            // dropdown used to send x == y, which the API answered with a 500.
            const v = e.target.value
            if (v === y) setY(x)
            setX(v)
          }}
          className="rounded border border-ink-200 px-1 py-0.5 font-mono text-[10px]"
        >
          {options.map((o) => (
            <option key={o.key} value={o.key}>{o.key}</option>
          ))}
        </select>
        <span className="text-ink-400">x</span>
        <select
          value={y}
          onChange={(e) => {
            const v = e.target.value
            if (v === x) setX(y)
            setY(v)
          }}
          className="rounded border border-ink-200 px-1 py-0.5 font-mono text-[10px]"
        >
          {options.map((o) => (
            <option key={o.key} value={o.key}>{o.key}</option>
          ))}
        </select>
      </div>
      <div className="mt-0.5 max-w-[19rem] text-[9px] leading-tight text-ink-400">
        {d.x_axis.poles && (
          <>{d.x_axis.poles.neg} &rarr; {d.x_axis.poles.pos}</>
        )}
      </div>
    </div>
  )
}

/** Component axes offered for the linear plot. Only the five named ones: past
 *  sPC5 a component has no interpretation, so browsing them here would invite
 *  reading meaning into an axis that has none. */
const LINEAR_AXES = (space: Space) =>
  [1, 2, 3, 4, 5].map((n) => ({
    key: `${space === 'corrected' ? 'spc' : 'pc'}${n}`,
    label: `${space === 'corrected' ? 'sPC' : 'PC'}${n}`,
  }))

const NONLINEAR_AXES = [
  { key: 'umap_1', label: 'UMAP 1' },
  { key: 'umap_2', label: 'UMAP 2' },
  { key: 'umap_null_1', label: 'UMAP null 1' },
  { key: 'umap_null_2', label: 'UMAP null 2' },
]

export function TaxonomyPanel({
  space = 'corrected',
  onPick,
}: {
  space?: Space
  onPick?: (s: string) => void
}) {
  const [d, setD] = useState<Taxonomy | null>(null)
  const [err, setErr] = useState<string | null>(null)
  const [showDeriv, setShowDeriv] = useState(false)

  useEffect(() => {
    api
      .taxonomy()
      .then(setD)
      .catch((e) => setErr(String(e)))
  }, [])

  if (err) return null
  if (!d) return null

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Contact-architecture taxonomy · {d.regions.length} regions
      </h2>
      <p className="mb-3 max-w-3xl text-[11px] leading-relaxed text-ink-600">
        Two orthogonal levels: how far contacts <strong>reach</strong>, and which
        element class <strong>dominates</strong> them. These are named areas of a
        continuum, not clusters, and the mixture statistic says so directly:{' '}
        <strong>{d.mixture.pct_below_half}% of genes have a top membership below
        0.5</strong> (median {d.mixture.median_top_weight}, against{' '}
        {d.mixture.uniform} for an even split across {d.regions.length} regions).
        A gene is a blend, never a member.
      </p>

      {/* Two rows. Radar and table first because they answer "how do the
          regions differ"; the projections come second because they answer the
          narrower "does that hold on an embedding". Previously all four sat in
          one wrapping flex row, so on a normal window the table was pushed
          below the scatters and read as a caption to them. */}
      <div className="flex flex-wrap items-start gap-6">
        <Radar data={d} />

        <div className="min-w-[15rem] flex-1">
          <table className="w-full text-[11px]">
            <thead>
              <tr className="border-b border-ink-200 text-left text-ink-500">
                <th className="pb-1 font-medium">region</th>
                <th className="pb-1 text-right font-medium">n</th>
                <th className="pb-1 text-right font-medium">mean top weight</th>
              </tr>
            </thead>
            <tbody>
              {d.regions.map((r) => (
                <tr key={r.region} className="border-b border-ink-50">
                  <td className="py-1">
                    <span
                      className="mr-1.5 inline-block h-2 w-2 rounded-sm align-middle"
                      style={{ background: COLOURS[r.region] ?? '#888' }}
                    />
                    <span className="font-mono">{r.region}</span>
                  </td>
                  <td className="py-1 text-right font-mono">{r.n}</td>
                  <td className="py-1 text-right font-mono text-ink-500">
                    {r.mean_top_weight?.toFixed(2)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          <p className="mt-2 text-[10px] leading-relaxed text-ink-400">
            {d.radar_note}
          </p>
        </div>
      </div>

      <div className="mt-4 flex flex-wrap items-start gap-6 border-t border-ink-100 pt-4">
        <RegionScatter
          initialX={DEFAULT_PLANE[space][0]}
          initialY={DEFAULT_PLANE[space][1]}
          colours={COLOURS}
          options={LINEAR_AXES(space)}
          onPick={onPick}
        />
        <RegionScatter
          initialX="umap_1"
          initialY="umap_2"
          colours={COLOURS}
          options={NONLINEAR_AXES}
          onPick={onPick}
        />
        <p className="max-w-[16rem] text-[10px] leading-relaxed text-ink-500">
          The regions were defined in 77 dimensions, not on either of these
          planes. Seeing them hold on a linear projection and a nonlinear one is
          a check: if they only cohere on the plane they are drawn in, that is
          the plane's doing. Point opacity is the gene's top membership, so the
          general faintness is the continuum result rendered directly.
        </p>
      </div>

      <button
        onClick={() => setShowDeriv((v) => !v)}
        className="mt-3 text-[11px] text-ink-600 underline decoration-ink-300 underline-offset-2"
      >
        {showDeriv ? 'hide' : 'show'} how these regions were derived
      </button>

      {showDeriv && (
        <div className="mt-2 space-y-3 border-t border-ink-100 pt-3">
          <ol className="space-y-2 text-[11px] leading-relaxed text-ink-600">
            {d.derivation.map((s) => (
              <li key={s.step} className="flex gap-2">
                <span className="shrink-0 font-mono text-ink-400">{s.step}</span>
                <span>
                  <strong className="text-ink-800">{s.name}.</strong> {s.detail}
                  {s.result && (
                    <span className="ml-1 font-mono text-ink-500">-&gt; {s.result}</span>
                  )}
                  {s.evidence && (
                    <span className="mt-0.5 block text-[10px] text-ink-400">
                      {s.evidence}
                    </span>
                  )}
                </span>
              </li>
            ))}
          </ol>

          <div>
            <p className="mb-1 text-[11px] font-medium text-ink-700">
              Why k differs between the two halves
            </p>
            <p className="mb-1.5 text-[10px] leading-relaxed text-ink-500">
              k was chosen inside each reach half rather than forced symmetric.
              A CTCF-dominated group only exists among far-reaching genes, which
              is what CTCF biology predicts: CTCF loops are long-range structural
              contacts. Forcing a third mid-range group manufactures four genes.
            </p>
            <table className="w-full max-w-lg text-[11px]">
              <thead>
                <tr className="border-b border-ink-200 text-left text-ink-500">
                  <th className="pb-1 font-medium">half</th>
                  <th className="pb-1 text-right font-medium">k</th>
                  <th className="pb-1 text-right font-medium">smallest group</th>
                  <th className="pb-1 text-right font-medium">seed stability</th>
                  <th className="pb-1 pl-2 font-medium">used</th>
                </tr>
              </thead>
              <tbody>
                {d.k_evidence.map((r, i) => (
                  <tr
                    key={i}
                    className={`border-b border-ink-50 ${
                      r.chosen ? 'text-ink-800' : 'text-ink-400'
                    }`}
                  >
                    <td className="py-1">{r.half}</td>
                    <td className="py-1 text-right font-mono">{r.k}</td>
                    <td className="py-1 text-right font-mono">
                      {r.smallest}
                      {r.smallest < 0.05 * r.n && (
                        <span className="ml-1 text-[9px]">pocket</span>
                      )}
                    </td>
                    <td className="py-1 text-right font-mono">
                      {r.stability.toFixed(3)}
                    </td>
                    <td className="py-1 pl-2">{r.chosen ? 'yes' : ''}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <p className="text-[10px] leading-relaxed text-ink-500">{d.caveat}</p>
        </div>
      )}
    </div>
  )
}
