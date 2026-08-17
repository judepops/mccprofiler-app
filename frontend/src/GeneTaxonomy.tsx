/**
 * One gene's contact architecture, in three sections that answer three
 * different questions:
 *
 *   REACH PROFILE   where this gene's contacts sit, against the panel median
 *   GROUP SIMILARITY  how close it is to each of the five clustered groups
 *   CONTACT TABLE   the full band x class breakdown with panel references
 *
 * Two rules learned the hard way today.
 *
 * 1. THE TABLE IS A CROSS-TABULATION, NOT A MEMBERSHIP. An earlier version
 *    listed the nine cells as "prox-promoter 49%, mid-promoter 15%, ..." with
 *    colour swatches in ranked order -- the same visual language as the
 *    clustered regions -- so bins read as cluster memberships. Every gene has
 *    all nine cells.
 *
 * 2. GROUP SIMILARITY IS NOT A CONTACT BREAKDOWN. The softmax over centroid
 *    distances says how close a gene sits to each group across all 77 features.
 *    Presenting it as composition produced "24% far-enhancer" beside "no peaks
 *    of this kind". The two are allowed to differ, and the beyond-50 kb line
 *    explains why they do.
 *
 * Prose lives in collapsible notes rather than inline paragraphs: the panel
 * carries a lot of caveat and inline it made the numbers hard to find.
 */

import { useEffect, useState } from 'react'
import { api, type GeneTaxonomy as GT } from './api'

/** Clustered-region palette, identical to the panel views: hue is composition,
 *  lightness is reach. */
const REGION_COLOUR: Record<string, string> = {
  'extended-ctcf': '#4a7c59',
  'extended-enhancer': '#c2703d',
  'extended-promoter': '#2b5070',
  'contained-ctcf': '#7fa98b',
  'contained-enhancer': '#dba57f',
  'contained-promoter': '#6d94b8',
}

const BAND_RANGE: Record<string, string> = {
  prox: '<50 kb',
  mid: '50-250 kb',
  far: '>250 kb',
}
const BAND_PHRASE: Record<string, string> = {
  prox: 'within 50 kb',
  mid: 'between 50 and 250 kb',
  far: 'beyond 250 kb',
}
const BAND_COLOUR: Record<string, string> = {
  prox: '#b8c6d3',
  mid: '#7b95ab',
  far: '#3d5f7d',
}

/** hue = element class (the dimension where the clustering and the evidence
 *  agree, 63% against 33% by chance); lightness = distance band. */
const CLASS_HUE: Record<string, [number, number, number]> = {
  promoter: [43, 80, 112],
  enhancer: [194, 112, 61],
  ctcf: [74, 124, 89],
}
const BAND_LIGHTEN: Record<string, number> = { prox: 0.55, mid: 0.25, far: 0 }

function cellColour(band: string, cls: string, pct: number): string {
  const hue = CLASS_HUE[cls] ?? [120, 130, 140]
  const l = BAND_LIGHTEN[band] ?? 0
  const [r, g, b] = hue.map((c) => Math.round(c + (255 - c) * l))
  const a = pct <= 0 ? 0 : 0.18 + Math.min(pct / 45, 1) * 0.72
  return `rgba(${r},${g},${b},${a})`
}

/** A collapsible note. Keeps explanation available without letting it crowd
 *  the numbers. */
function Note({ label = 'why', children }: { label?: string; children: React.ReactNode }) {
  const [open, setOpen] = useState(false)
  return (
    <>
      <button
        onClick={() => setOpen((v) => !v)}
        className="text-[10px] text-ink-400 underline decoration-ink-200 underline-offset-2 hover:text-ink-600"
      >
        {open ? 'hide' : label}
      </button>
      {open && (
        <div className="mt-1.5 space-y-1.5 rounded border border-ink-100 bg-ink-50/60 p-2 text-[10px] leading-relaxed text-ink-600">
          {children}
        </div>
      )}
    </>
  )
}

/**
 * The whole membership in one bar.
 *
 * Deliberately NOT a stacked bar of the five weights. That version was removed
 * because a share cannot express the uniform floor: with K=5 a weight of 20% is
 * zero evidence, so a stacked bar drew a fifth of the width for a region the
 * gene has no lean toward at all, and a gene with no preference looked identical
 * to a five-way blend.
 *
 * What is stacked instead is the ABOVE-FLOOR EXCESS. Only regions the gene leans
 * toward appear, sized by how far above the floor they sit. Regions at or below
 * the floor contribute nothing, which is the correct width for no evidence.
 *
 * The filled fraction of the track is then TVD exactly, because the weights sum
 * to 100 so the positive and negative deviations are equal in magnitude and the
 * total excess is TVD x 100. So one bar carries both questions at once: WHICH
 * regions from the segments, and HOW MUCH LEAN AT ALL from the total width. An
 * empty bar is a gene with no preference, which is the common case.
 */
function MembershipSummary({
  steps,
  floor,
  flatness,
}: {
  steps: { region: string; weight: number }[]
  floor: number
  flatness: GT['flatness'] | undefined
}) {
  const total = steps.reduce((s, x) => s + x.weight, 0) || 1
  const excess = steps
    .map((s) => ({
      region: s.region,
      pp: (s.weight / total) * 100 - floor,
    }))
    .filter((x) => x.pp > 0)
    .sort((a, b) => b.pp - a.pp)
  const sum = excess.reduce((s, x) => s + x.pp, 0)

  // Scale to just above the panel maximum so bars are comparable BETWEEN genes,
  // not normalised per gene (which would make every gene look equally decided).
  const full = (flatness?.panel_max_tvd ?? 0.5) * 100
  const medTick = ((flatness?.panel_median_tvd ?? 0.119) * 100 * 100) / full
  const filled = Math.min((sum / full) * 100, 100)
  const inflated = flatness?.distance_inflated

  return (
    <div>
      <div className="relative h-6 w-full overflow-hidden rounded bg-ink-100">
        {excess.map((x) => (
          <div
            key={x.region}
            className="absolute inset-y-0"
            style={{
              left: `${(excess.slice(0, excess.indexOf(x)).reduce((s, y) => s + y.pp, 0) / full) * 100}%`,
              width: `${(x.pp / full) * 100}%`,
              background: REGION_COLOUR[x.region] ?? '#8badc9',
              // Striped when the confidence is a distance artefact, so the one
              // glance summary cannot silently endorse it.
              backgroundImage: inflated
                ? 'repeating-linear-gradient(45deg, rgba(255,255,255,.45) 0 3px, transparent 3px 6px)'
                : undefined,
            }}
            title={`${x.region}: +${x.pp.toFixed(1)} pp above the ${floor.toFixed(0)}% floor`}
          />
        ))}
        {/* Panel median lean, so the width means something. */}
        <div
          className="absolute inset-y-0 w-px bg-ink-500"
          style={{ left: `${medTick}%` }}
          title={`panel median lean (TVD ${flatness?.panel_median_tvd})`}
        />
      </div>
      <div className="mt-1 flex flex-wrap items-baseline gap-x-3 gap-y-0.5 text-[10px] text-ink-500">
        <span>
          total lean{' '}
          <strong className="text-ink-700">{sum.toFixed(0)} pp</strong> above the{' '}
          {floor.toFixed(0)}% floor
        </span>
        <span className="text-ink-400">
          panel median {((flatness?.panel_median_tvd ?? 0) * 100).toFixed(0)} pp ·
          max {((flatness?.panel_max_tvd ?? 0) * 100).toFixed(0)} pp
        </span>
        {filled < 1 && (
          <span className="text-ink-400">
            empty bar = no lean toward any region
          </span>
        )}
        {inflated && (
          <span className="text-element-enhancer">
            striped: lean is a distance artefact
          </span>
        )}
      </div>
      <div className="mt-1.5 flex flex-wrap gap-x-3 gap-y-1 text-[10px]">
        {excess.map((x) => (
          <span key={x.region} className="flex items-center gap-1">
            <span
              className="inline-block h-2 w-2 rounded-sm"
              style={{ background: REGION_COLOUR[x.region] ?? '#8badc9' }}
            />
            <span className="font-mono text-ink-600">{x.region}</span>
            <strong className="text-ink-700">+{x.pp.toFixed(0)}</strong>
          </span>
        ))}
        {!excess.length && (
          <span className="text-ink-500">
            no region sits above the floor for this gene
          </span>
        )}
      </div>
    </div>
  )
}

function SectionHead({ title, note }: { title: string; note?: React.ReactNode }) {
  return (
    <div className="mb-1.5 flex items-baseline gap-2">
      <h3 className="text-[11px] font-medium text-ink-800">{title}</h3>
      {note && <Note>{note}</Note>}
    </div>
  )
}

export function GeneTaxonomy({ gene }: { gene: string }) {
  const [d, setD] = useState<GT | null>(null)
  const [showDeriv, setShowDeriv] = useState(false)

  useEffect(() => {
    let live = true
    setD(null)
    api.geneTaxonomy(gene).then((v) => live && setD(v)).catch(() => live && setD(null))
    return () => {
      live = false
    }
  }, [gene])

  if (!d) return null

  const ev = d.evidence
  const pm = ev?.panel_median
  const rows = ev?.rows ?? []
  const classes = ev?.classes ?? []
  const steps = d.derivation?.steps ?? []
  const a = d.assignment

  const mine = rows.map((r) => ({
    band: r.band,
    pct: (r.promoter ?? 0) + (r.enhancer ?? 0) + (r.ctcf ?? 0),
  }))
  const medBands = pm?.bands ?? {}
  const cells = rows.flatMap((r) =>
    classes.map((c) => ({
      phrase: BAND_PHRASE[r.band] ?? r.band,
      cls: c,
      pct: ((r as unknown as Record<string, number | null>)[c] ?? 0) as number,
    })),
  )
  const strongest = cells.length
    ? cells.reduce((x, y) => (y.pct > x.pct ? y : x))
    : null
  const simTotal = steps.reduce((s, x) => s + x.weight, 0) || 1
  // The uniform floor, 100/K. Served rather than hardcoded so it follows K if
  // the taxonomy ever changes shape.
  const floor = d.flatness?.uniform_pct ?? (steps.length ? 100 / steps.length : null)

  if (ev?.no_peaks_called) {
    return (
      <div className="rounded-lg border border-ink-200 bg-white p-4">
        <h2 className="mb-2 text-xs font-semibold uppercase tracking-wide text-ink-500">
          Contact architecture
        </h2>
        <p className="text-[11px] leading-relaxed text-ink-700">
          <strong>No peaks were called for this gene</strong>, so there is nothing
          to describe. Its MCC signal is normal; the annotation is what is
          missing.
        </p>
        {d.under_evidenced && (
          <p className="mt-2 text-[10px] leading-relaxed text-ink-500">
            {d.under_evidenced_note}
          </p>
        )}
      </div>
    )
  }

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-5">
      <h2 className="mb-3 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Contact architecture
      </h2>

      {/* Genes with 1 or 2 peaks still have a contact table to draw, so without
          this banner the panel would show an architecture breakdown with no
          membership beside it and no reason given. */}
      {d.under_evidenced && (
        <div className="mb-4 rounded border border-element-enhancer/40 bg-element-enhancer/5 p-2.5">
          <p className="text-[11px] leading-relaxed text-ink-700">
            <strong className="text-element-enhancer">
              Excluded from the taxonomy.
            </strong>{' '}
            Only {ev?.n_peaks} called peak{ev?.n_peaks === 1 ? '' : 's'}, so this
            gene has no region and no membership. Anything below is computed from
            those peaks alone and should not be read as an architecture.
          </p>
          <div className="mt-1">
            <Note label="why it was excluded">
              <p>{d.under_evidenced_note}</p>
            </Note>
          </div>
        </div>
      )}

      <div className="space-y-6">
        {/* ---------------- reach profile ---------------- */}
        <section>
          <SectionHead
            title="Reach profile"
            note={
              <>
                <p>
                  Where this gene's peak signal sits, against the panel median.
                  The bands are fixed cutoffs in the feature pipeline, not a
                  result: every distance feature uses the same ones.
                </p>
                <p>
                  Contact frequency falls steeply with distance, so a large
                  proximal share is normal for every gene (panel median about
                  55% within 50 kb). What distinguishes genes is how much reaches
                  beyond that.
                </p>
              </>
            }
          />
          {/* Identical colours in both rows. An earlier version faded the median
              row to 40% opacity, which changed every hue and made the two rows
              incomparable by eye -- exactly what the pairing exists for. The
              reference row is distinguished by height instead. */}
          <div className="space-y-1.5">
            {[
              { label: 'this gene', data: mine, total: 100, h: 'h-4' },
              {
                label: 'panel median',
                data: ['prox', 'mid', 'far'].map((b) => ({
                  band: b,
                  pct: medBands[b] ?? 0,
                })),
                total:
                  Object.values(medBands).reduce((x, y) => x + y, 0) || 1,
                h: 'h-2',
              },
            ].map((row) => (
              <div key={row.label} className="flex items-center gap-2">
                <span className="w-24 shrink-0 text-[10px] text-ink-500">
                  {row.label}
                </span>
                <div
                  className={`flex ${row.h} flex-1 overflow-hidden rounded`}
                >
                  {row.data.map((b) => (
                    <div
                      key={b.band}
                      style={{
                        width: `${(b.pct / row.total) * 100}%`,
                        background: BAND_COLOUR[b.band],
                      }}
                      title={`${BAND_RANGE[b.band]}: ${b.pct.toFixed(0)}%`}
                    />
                  ))}
                </div>
              </div>
            ))}
          </div>
          <div className="mt-1.5 flex flex-wrap gap-x-4 gap-y-1 pl-26 text-[10px] text-ink-500">
            {mine.map((b) => (
              <span key={b.band} className="flex items-center gap-1">
                <span
                  className="inline-block h-2 w-2 rounded-sm"
                  style={{ background: BAND_COLOUR[b.band] }}
                />
                {BAND_RANGE[b.band]}{' '}
                <strong className="text-ink-700">{b.pct.toFixed(0)}%</strong>
                <span className="text-ink-400">
                  vs {(medBands[b.band] ?? 0).toFixed(0)}%
                </span>
              </span>
            ))}
          </div>
        </section>

        {/* ---------------- group similarity ---------------- */}
        {steps.length > 0 && (
          <section>
            <SectionHead
              title="Similarity to each group"
              note={
                <>
                  <p>
                    How close this gene sits to each of the five clustered groups
                    across all 77 features. This is a similarity profile, NOT a
                    breakdown of the gene's contacts, so it is allowed to differ
                    from the table below.
                  </p>
                  <p>
                    <strong>Read against the {floor?.toFixed(0)}% floor, not as
                    a share.</strong>{' '}
                    A softmax over five groups cannot go below{' '}
                    {floor?.toFixed(0)}% for a gene that resembles nothing in
                    particular, so {floor?.toFixed(0)}% means no evidence either
                    way and below it means evidence against. Bars therefore
                    diverge from the floor. The room to move is small by
                    construction: nearest-centroid distance is 7.29 against a
                    furthest of 9.67 at the panel median, because in 77
                    dimensions every centroid is far from every point.
                    Consequently{' '}
                    {d.flatness?.panel_pct_near_uniform ?? 63}% of genes have
                    their entire profile within 10 pp of the floor and only{' '}
                    {d.flatness?.panel_pct_above_40 ?? 3}% reach any weight of
                    40%.
                  </p>
                  <p>
                    The groups reproduce across independent captures (reach ARI
                    0.741, composition 0.616, seed stability 0.977) but there are
                    no clusters in the strict sense: HDBSCAN finds none across
                    sixteen conditions and Leiden returns a single community below
                    resolution 0.4.
                  </p>
                  <p>
                    "Extended" is relative: it means more extended than the panel
                    average, not that the bulk of contacts sits beyond 250 kb.
                  </p>
                </>
              }
            />
            {/* DIVERGING AROUND THE UNIFORM FLOOR, not a stacked share.
                A softmax over 5 regions cannot fall below 20% for a gene that
                resembles nothing in particular, so a stacked bar reads 20% as a
                fifth of the architecture when it is in fact zero evidence. That
                misreading was reported directly: a gene 79% within 50 kb and 0%
                beyond 250 kb showed 20% extended-enhancer, which looked like a
                contradiction and was the floor. Below the floor is evidence
                AGAINST, which a share cannot express at all.

                The floor is not a presentation choice, it is forced: measured on
                the panel, nearest-centroid distance is 7.29 against a furthest
                of 9.67, so the softmax has almost no room to move and 63% of
                genes never leave the floor. */}
            {floor != null && (
              <>
                <div className="mb-1 flex items-center gap-2 text-[9px] text-ink-400">
                  <span className="w-40 shrink-0" />
                  <span className="flex-1 text-center">
                    less similar than chance &lt; {floor.toFixed(0)}% &gt; more
                    similar
                  </span>
                  <span className="w-20 shrink-0" />
                </div>
                <div className="space-y-1">
                  {steps.map((s) => {
                    const pct = (s.weight / simTotal) * 100
                    const dev = pct - floor
                    // Full width = the largest deviation seen anywhere in the
                    // panel, so bar length is comparable between genes.
                    const w = Math.min(Math.abs(dev) / 12, 1) * 50
                    return (
                      <div
                        key={s.region}
                        className="flex items-center gap-2 text-[11px]"
                      >
                        <span className="flex w-40 shrink-0 items-center gap-1.5">
                          <span
                            className="inline-block h-2 w-2 shrink-0 rounded-sm"
                            style={{
                              background: REGION_COLOUR[s.region] ?? '#8badc9',
                            }}
                          />
                          <span className="truncate font-mono text-ink-700">
                            {s.region}
                          </span>
                        </span>
                        <span className="relative h-4 flex-1">
                          <span className="absolute inset-y-0 left-1/2 w-px bg-ink-300" />
                          <span
                            className="absolute inset-y-0 rounded-sm"
                            style={{
                              left: dev >= 0 ? '50%' : `${50 - w}%`,
                              width: `${w}%`,
                              background: REGION_COLOUR[s.region] ?? '#8badc9',
                              opacity: dev >= 0 ? 0.85 : 0.35,
                            }}
                            title={`${s.region}: ${pct.toFixed(1)}%, ${
                              dev >= 0 ? '+' : ''
                            }${dev.toFixed(1)} pp vs the ${floor.toFixed(0)}% floor`}
                          />
                        </span>
                        <span className="w-20 shrink-0 text-right font-mono text-[10px]">
                          <span className="text-ink-600">
                            {dev >= 0 ? '+' : ''}
                            {dev.toFixed(0)} pp
                          </span>
                          <span className="ml-1 text-ink-400">
                            d {s.distance.toFixed(1)}
                          </span>
                        </span>
                      </div>
                    )
                  })}
                </div>
              </>
            )}

            {d.flatness?.near_uniform && (
              <p className="mt-2 rounded border border-ink-200 bg-ink-50 p-2 text-[11px] leading-relaxed text-ink-700">
                <strong>This profile is flat.</strong> No region is more than{' '}
                {d.flatness.max_deviation_pp?.toFixed(0)} pp from the{' '}
                {floor?.toFixed(0)}% floor, so the ranking below carries very
                little information about this gene. That is the common case, not a
                fault: {d.flatness.panel_pct_near_uniform}% of the panel sits
                within 10 pp of the floor and only{' '}
                {d.flatness.panel_pct_above_40}% of genes reach any weight of
                40%.
              </p>
            )}

            {/* The opposite warning, and the more dangerous one because it looks
                like a strong result rather than a weak one. */}
            {d.flatness?.distance_inflated && (
              <div className="mt-2 rounded border border-element-enhancer/40 bg-element-enhancer/5 p-2">
                <p className="text-[11px] leading-relaxed text-ink-700">
                  <strong className="text-element-enhancer">
                    Confidence here is inflated by distance.
                  </strong>{' '}
                  This gene sits{' '}
                  {d.flatness.nearest_distance?.toFixed(1)} from its NEAREST
                  group centroid, the{' '}
                  {d.flatness.nearest_distance_pct?.toFixed(0)}th percentile of
                  the panel against a median of{' '}
                  {d.flatness.panel_median_nearest_distance}. It is far from all
                  five groups rather than close to one, so a large percentage
                  above means "unlike everything" and not "clearly this group".
                </p>
                <div className="mt-1">
                  <Note label="why distance does this">
                    <p>
                      The weight is exp(-d&sup2;/tau), so the distance is
                      squared. A 10% gap between the nearest and furthest
                      centroid becomes an enormous gap once squared, which the
                      softmax reads as certainty. {d.flatness.distance_note}
                    </p>
                  </Note>
                </div>
              </div>
            )}

            {a && (
              <div className="mt-2.5 rounded border border-ink-100 bg-ink-50/50 p-2">
                <p className="text-[11px] leading-relaxed text-ink-700">
                  <span
                    className="mr-1.5 inline-block h-2.5 w-2.5 rounded-sm align-middle"
                    style={{ background: REGION_COLOUR[a.group] ?? '#8badc9' }}
                  />
                  <strong>Assigned: {a.group}</strong>
                  {a.runner_up && a.margin_pct != null && (
                    <>
                      {' '}
                      &mdash; {a.margin_pct.toFixed(0)}% closer than {a.runner_up}
                    </>
                  )}
                  {a.close_call && (
                    <span className="ml-1.5 rounded bg-element-enhancer/15 px-1.5 py-0.5 text-[9px] font-medium uppercase text-element-enhancer">
                      close call
                    </span>
                  )}
                </p>
                <div className="mt-1">
                  <Note label="how decisive is this?">
                    <p>
                      Panel median margin over the runner-up is 8.1%, and 59% of
                      genes have a runner-up within 10%. For those the single
                      label is close to arbitrary and should be read loosely.
                    </p>
                    {a.beyond_50kb?.dominant && a.beyond_50kb.total_pct > 0 && (
                      <p>
                        The grouping keys on reach and composition together, so
                        what drives it is the{' '}
                        {a.beyond_50kb.total_pct.toFixed(0)}% of contacts beyond
                        50 kb: those are{' '}
                        {a.beyond_50kb.shares[a.beyond_50kb.dominant]?.toFixed(0)}%{' '}
                        {a.beyond_50kb.dominant}. That is why the assigned group
                        can differ from the largest cell in the table, which is
                        usually proximal for every gene.
                      </p>
                    )}
                  </Note>
                </div>
              </div>
            )}
          </section>
        )}

        {/* ---------------- membership summary ---------------- */}
        {steps.length > 0 && floor != null && (
          <section>
            <SectionHead
              title="Membership, in one bar"
              note={
                <>
                  <p>
                    The five bars above, condensed. Segments are each region's{' '}
                    <strong>excess above the {floor.toFixed(0)}% floor</strong>,
                    so only regions the gene actually leans toward appear and the
                    rest of the track stays empty.
                  </p>
                  <p>
                    Deliberately not a stacked bar of the five percentages. A
                    share cannot express the floor: a stacked bar gives a fifth of
                    its width to a region carrying no evidence, which makes a gene
                    with no preference look identical to a genuine five-way blend.
                  </p>
                  <p>
                    The filled fraction is therefore the total lean, and it is
                    comparable between genes because the scale is fixed to the
                    panel maximum rather than normalised per gene. The tick marks
                    the panel median lean of{' '}
                    {((d.flatness?.panel_median_tvd ?? 0) * 100).toFixed(0)} pp.
                    An empty or near-empty bar is the common case, not a failure:
                    {' '}{d.flatness?.panel_pct_near_uniform ?? 63}% of the panel
                    sits within 10 pp of the floor.
                  </p>
                </>
              }
            />
            <MembershipSummary
              steps={steps}
              floor={floor}
              flatness={d.flatness}
            />
          </section>
        )}

        {/* ---------------- contact table ---------------- */}
        <section>
          <SectionHead
            title={`Contact breakdown (${ev?.n_peaks} peaks)`}
            note={
              <>
                <p>
                  Share of this gene's peak signal. Rows are distance from the
                  viewpoint, columns are element class. These are bins, not
                  groups: every gene has all nine, and the row totals give the
                  reach profile above.
                </p>
                <p>
                  Small grey numbers are the panel <strong>median</strong> and
                  this gene's <strong>percentile</strong>. A tilde marks cells
                  where over a quarter of genes sit at zero, so the percentile is
                  inflated by ties: in the far row 41-63% of genes have no
                  contact of that class at all, which is itself the useful fact.
                </p>
                <p>
                  <strong>These shares are NOT distance-corrected, so read the
                  median and percentile, not the raw number.</strong>{' '}
                  Cells weight each peak by its height (`peak_max`), which falls
                  with distance from the viewpoint: median 13.70 within 50 kb
                  against 6.17 mid and 4.73 far, a 2.9-fold drop, rho -0.555
                  with distance. So a large proximal number is mostly contact
                  decay, which is why the panel median for proximal signal is
                  about 55% for every gene. Comparing against the median or the
                  percentile cancels the decay, because the reference carries the
                  same geometry; the absolute share does not.
                </p>
                <p>
                  Within a row the comparison is clean, since those peaks sit at
                  similar distance. <strong>Down a column it is not.</strong>{' '}
                  Element class is confounded with distance in this panel:
                  median offset is 57.8 kb for enhancer peaks, 92.8 kb for
                  promoter and 112.5 kb for CTCF, and class composition shifts
                  hard across the bands (enhancers are 47% of proximal peaks but
                  21% of far ones; CTCF 21% rising to 37%). Enhancers are both
                  closer and rewarded for being closer, so a high overall
                  enhancer share partly reports that a gene's peaks are near the
                  viewpoint.
                </p>
                <p>
                  The proximal row is not one of the clustered reach levels; those
                  distinguish relative extension.
                </p>
              </>
            }
          />
          {strongest && (
            <p className="mb-2 text-[12px] leading-relaxed text-ink-700">
              Mostly <strong>{strongest.cls}</strong>-classed, concentrated{' '}
              {strongest.phrase} ({strongest.pct.toFixed(0)}% of signal).
            </p>
          )}
          <table className="w-full max-w-lg text-[11px]">
            <thead>
              <tr className="text-left text-ink-500">
                <th className="pb-1.5 pr-4 font-medium">distance</th>
                {classes.map((c) => (
                  <th key={c} className="pb-1.5 pr-4 text-right font-medium">
                    {c}
                    <div className="font-normal text-[9px] text-ink-400">
                      med · pct
                    </div>
                  </th>
                ))}
                <th className="pb-1.5 text-right font-medium">row</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((r) => {
                const rowTotal =
                  (r.promoter ?? 0) + (r.enhancer ?? 0) + (r.ctcf ?? 0)
                return (
                  <tr key={r.band} className="border-t border-ink-50">
                    <td className="py-1.5 pr-4 font-mono text-ink-700">
                      {BAND_RANGE[r.band] ?? r.band}
                    </td>
                    {classes.map((c) => {
                      const key = `${r.band}-${c}`
                      const v = ((r as unknown as Record<string, number | null>)[
                        c
                      ] ?? 0) as number
                      const tie = (pm?.zero_pct?.[key] ?? 0) > 25
                      return (
                        <td key={c} className="py-1.5 pr-4 text-right">
                          <div
                            className="inline-block rounded px-1.5 py-0.5 font-mono"
                            style={{
                              background: cellColour(r.band, c, v),
                              color: v > 28 ? '#ffffff' : '#3d4a57',
                            }}
                          >
                            {v === 0 ? '-' : `${v.toFixed(0)}%`}
                          </div>
                          <div className="font-mono text-[9px] text-ink-400">
                            {(pm?.cells?.[key] ?? 0).toFixed(0)} ·{' '}
                            {tie ? '~' : ''}
                            {(pm?.percentile?.[key] ?? 0).toFixed(0)}
                          </div>
                        </td>
                      )
                    })}
                    <td className="py-1.5 text-right font-mono text-ink-500">
                      {rowTotal.toFixed(0)}%
                    </td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </section>

        {/* ---------------- derivation ---------------- */}
        {steps.length > 0 && (
          <section className="border-t border-ink-100 pt-3">
            <button
              onClick={() => setShowDeriv((v) => !v)}
              className="text-[11px] text-ink-600 underline decoration-ink-300 underline-offset-2"
            >
              {showDeriv ? 'hide' : 'show'} how the assigned group was chosen
            </button>
            {showDeriv && (
              <div className="mt-2 space-y-2">
                <p className="font-mono text-[10px] text-ink-500">
                  {d.derivation.formula} · tau = {d.derivation.tau}
                </p>
                <table className="max-w-md text-[11px]">
                  <thead>
                    <tr className="text-left text-ink-500">
                      <th className="pb-1 pr-4 font-medium">group</th>
                      <th className="pb-1 pr-4 text-right font-medium">distance</th>
                      <th className="pb-1 text-right font-medium">
                        further than nearest
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {steps.map((s) => (
                      <tr key={s.region} className="border-t border-ink-50">
                        <td className="py-1 pr-4 font-mono text-ink-700">
                          {s.region}
                        </td>
                        <td className="py-1 pr-4 text-right font-mono">
                          {s.distance.toFixed(2)}
                        </td>
                        <td className="py-1 text-right font-mono text-ink-500">
                          {s.excess_over_nearest === 0
                            ? '--'
                            : `+${s.excess_over_nearest.toFixed(2)}`}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                <p className="max-w-2xl text-[10px] leading-relaxed text-ink-500">
                  {d.derivation.reading}
                </p>
              </div>
            )}
          </section>
        )}
      </div>
    </div>
  )
}
