/**
 * Three pages, with one rule: a page is about a gene, or about the panel, and
 * never both.
 *
 *   Gene   this gene's own data. Profile, position, peaks, features, compare.
 *   Panel  properties of the coordinate system and the population. Nothing
 *          here is about any particular gene, and the finding tools live here
 *          because you use them BEFORE you have one.
 *   Lab    exploratory and methodological views, not part of the product.
 *
 * Previously the panel-level views were repeated underneath the gene detail,
 * which made it ambiguous whether a number described the gene you had searched
 * or the panel it sits in. That ambiguity is the thing this split removes.
 */

import { useCallback, useEffect, useRef, useState } from 'react'
import {
  api,
  levelForSpan,
  type Embedding,
  type FeatureSet,
  type Gene,
  type GeneSummary,
  type Health,
  type Peak,
  type PeakSet,
  type Profile,
} from './api'
import { ProfilePlot } from './ProfilePlot'
import { ArchetypeReadout } from './ArchetypeReadout'
import { FeatureTable } from './FeatureTable'
import { ContinuumMap } from './ContinuumMap'
import { CohortView } from './CohortView'
import { DimensionPanel } from './DimensionPanel'
import { LabPage } from './LabPage'
import { ExplainPanel } from './ExplainPanel'
import { RankedLists } from './RankedLists'
import { AskPanel } from './AskPanel'
import { PeakDetail } from './PeakDetail'
import { GeneCompare } from './GeneCompare'

type Page = 'gene' | 'panel' | 'lab'

const PAGE_LABEL: Record<Page, string> = {
  gene: 'Gene',
  panel: 'Panel',
  lab: 'Lab',
}

const CHANNEL_LABEL: Record<string, string> = {
  mcc: 'MCC',
  atac: 'ATAC',
  ctcf: 'CTCF',
  rna_plus: 'RNA +',
  rna_minus: 'RNA −',
  h3k4me1: 'H3K4me1',
  h3k4me3: 'H3K4me3',
  h3k27ac: 'H3K27ac',
}

export default function App() {
  const [page, setPage] = useState<Page>('gene')
  const [health, setHealth] = useState<Health | null>(null)
  const [err, setErr] = useState<string | null>(null)

  const [query, setQuery] = useState('')
  const [hits, setHits] = useState<GeneSummary[]>([])
  const [gene, setGene] = useState<Gene | null>(null)

  const [channel, setChannel] = useState('mcc')
  const [mode, setMode] = useState<'raw' | 'oe'>('raw')
  const [range, setRange] = useState<{ start: number; end: number } | null>(null)
  const [profile, setProfile] = useState<Profile | null>(null)
  const [peaks, setPeaks] = useState<PeakSet | null>(null)
  const [features, setFeatures] = useState<FeatureSet | null>(null)
  const [peakLimit, setPeakLimit] = useState(25)
  const [selectedPeak, setSelectedPeak] = useState<Peak | null>(null)
  const [loading, setLoading] = useState(false)

  const [embedding, setEmbedding] = useState<Embedding | null>(null)
  const [axes, setAxes] = useState<[string, string]>(['pc2', 'pc3'])
  const [embLoading, setEmbLoading] = useState(false)
  const [embErr, setEmbErr] = useState<string | null>(null)
  const embReq = useRef(0)
  // Lifted out of GeneCompare so the map can mark the compared gene too.
  const [compareGenes, setCompareGenes] = useState<string[]>([])
  // Latched on first visit to the Panel page. Before that the panel is not
  // rendered at all, so an initial Gene-page load does not pay for cohort and
  // ranked-list fetches nobody asked for.
  const panelVisited = useRef(false)
  if (page === 'panel') panelVisited.current = true

  const plotWidth = useRef(900)

  useEffect(() => {
    api.health().then(setHealth).catch((e) => setErr(String(e.message ?? e)))
  }, [])

  useEffect(() => {
    if (query.trim().length < 1) {
      setHits([])
      return
    }
    const t = setTimeout(() => {
      api
        .searchGenes(query.trim(), 12)
        .then((r) => setHits(r.genes))
        .catch(() => setHits([]))
    }, 150)
    return () => clearTimeout(t)
  }, [query])

  const loadProfile = useCallback(
    async (symbol: string, r: { start: number; end: number } | null) => {
      setLoading(true)
      try {
        const span = r ? r.end - r.start : 2_000_000
        setProfile(
          await api.profile(symbol, {
            channel,
            mode,
            level: levelForSpan(span, plotWidth.current),
            start_bp: r?.start,
            end_bp: r?.end,
          }),
        )
      } catch (e) {
        setErr(String((e as Error).message))
      } finally {
        setLoading(false)
      }
    },
    [channel, mode],
  )

  /** Selecting a gene anywhere always lands you on the Gene page. */
  async function select(symbol: string) {
    setQuery('')
    setHits([])
    setRange(null)
    setErr(null)
    setSelectedPeak(null)
    setPage('gene')
    try {
      const [g, pk, ft] = await Promise.all([
        api.gene(symbol),
        api.peaks(symbol).catch(() => null),
        api.features(symbol).catch(() => null),
      ])
      setGene(g)
      setPeaks(pk)
      setFeatures(ft)
      await loadProfile(symbol, null)
    } catch (e) {
      setErr(String((e as Error).message))
    }
  }

  useEffect(() => {
    if (gene) loadProfile(gene.gene_symbol, range)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [channel, mode, range])

  // Fetched for both the Gene and Panel pages. The same view answers two
  // different questions: on the Panel page it is the population, on the Gene
  // page it is where this one gene sits within it. The Lab page does not need
  // it, so it is still not fetched there.
  useEffect(() => {
    if (page === 'lab') return
    const token = ++embReq.current
    setEmbLoading(true)
    setEmbErr(null)
    api
      .embedding(axes[0], axes[1], gene?.gene_symbol)
      .then((e) => token === embReq.current && setEmbedding(e))
      .catch((e) => token === embReq.current && setEmbErr(String((e as Error).message)))
      .finally(() => token === embReq.current && setEmbLoading(false))
  }, [axes, gene, page])

  return (
    <div className="min-h-screen bg-ink-50 font-sans text-ink-900">
      <header className="border-b border-ink-200 bg-white">
        <div className="mx-auto flex max-w-6xl items-center justify-between gap-4 px-6 py-3">
          <div>
            <h1 className="text-sm font-semibold tracking-tight">mccprofiler</h1>
            <p className="text-[11px] text-ink-500">
              Gene position in the MCC regulatory continuum · CD4+ T cells
            </p>
          </div>

          <nav className="flex gap-1 text-[11px]">
            {(Object.keys(PAGE_LABEL) as Page[]).map((p) => (
              <button
                key={p}
                onClick={() => setPage(p)}
                className={`rounded px-2.5 py-1 ${
                  page === p ? 'bg-ink-600 text-white' : 'text-ink-600 hover:bg-ink-50'
                }`}
              >
                {PAGE_LABEL[p]}
              </button>
            ))}
          </nav>

          <div className="text-right">
            <a
              href="http://localhost:8000/api/export/panel.csv"
              className="text-[11px] text-ink-500 underline-offset-2 hover:underline"
            >
              export panel
            </a>
            {health && (
              <p className="font-mono text-[10px] leading-4 text-ink-400">
                {health.panel} · {health.n_genes.toLocaleString()} genes
                <br />
                store v{health.store_version} · {health.scripts_cleaned_commit}
              </p>
            )}
          </div>
        </div>
      </header>

      <main className="mx-auto max-w-6xl px-6 py-6">
        {err && (
          <div className="mb-4 rounded border border-element-enhancer/30 bg-element-enhancer/5 px-3 py-2 text-[13px] text-element-enhancer">
            {err}
            {err.includes('no store') && (
              <span className="ml-1 text-ink-600">
                run <code className="font-mono">build_store.py</code> first.
              </span>
            )}
          </div>
        )}

        {/* ---------------------------------------------------------------- */}
        {/* GENE: this gene's own data, and nothing else                      */}
        {/* ---------------------------------------------------------------- */}
        {page === 'gene' && (
          <>
            <div className="relative mb-6">
              <input
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === 'Enter' && hits.length) select(hits[0].gene_symbol)
                }}
                placeholder="Search a gene, try IL7R, CTCF, EEF1A1"
                className="w-full rounded-lg border border-ink-200 bg-white px-4 py-2.5 text-sm outline-none placeholder:text-ink-300 focus:border-ink-400"
              />
              {hits.length > 0 && (
                <ul className="absolute z-10 mt-1 w-full overflow-hidden rounded-lg border border-ink-200 bg-white shadow-lg">
                  {hits.map((h) => (
                    <li key={h.gene_id}>
                      <button
                        onClick={() => select(h.gene_symbol)}
                        className="flex w-full items-baseline justify-between px-4 py-2 text-left text-sm hover:bg-ink-50"
                      >
                        <span className="font-medium">{h.gene_symbol}</span>
                        <span className="font-mono text-[11px] text-ink-400">
                          {h.viewpoint_chrom}:{h.viewpoint_pos?.toLocaleString()}
                        </span>
                      </button>
                    </li>
                  ))}
                </ul>
              )}
            </div>

            {!gene && !err && (
              <div className="py-16 text-center">
                <p className="text-sm text-ink-400">
                  Search a gene to see its contact architecture.
                </p>
                <button
                  onClick={() => setPage('panel')}
                  className="mt-2 text-[12px] text-ink-500 underline underline-offset-2 hover:text-ink-800"
                >
                  or find one on the Panel page
                </button>
              </div>
            )}

            {gene && (
              <div className="grid grid-cols-1 gap-5 lg:grid-cols-3">
                <section className="lg:col-span-2">
                  <div className="rounded-lg border border-ink-200 bg-white p-4">
                    <div className="mb-3 flex items-baseline justify-between gap-3">
                      <div>
                        <h2 className="text-lg font-semibold">{gene.gene_symbol}</h2>
                        <p className="font-mono text-[11px] text-ink-500">
                          viewpoint {gene.viewpoint.chrom}:
                          {gene.viewpoint.pos?.toLocaleString()}
                        </p>
                      </div>

                      <div className="flex flex-wrap items-center gap-2">
                        <select
                          value={channel}
                          onChange={(e) => setChannel(e.target.value)}
                          className="rounded border border-ink-200 px-2 py-1 text-xs"
                        >
                          {(health?.channels ?? ['mcc']).map((c) => (
                            <option key={c} value={c}>
                              {CHANNEL_LABEL[c] ?? c}
                            </option>
                          ))}
                        </select>

                        {channel === 'mcc' && (peaks?.n ?? 0) > 0 && (
                          <label className="flex items-center gap-1 text-[11px] text-ink-600">
                            top
                            <select
                              value={peakLimit}
                              onChange={(e) => setPeakLimit(Number(e.target.value))}
                              className="rounded border border-ink-200 px-1 py-0.5"
                            >
                              {[10, 25, 50, 200].map((n) => (
                                <option key={n} value={n}>
                                  {n >= (peaks?.n ?? 0) ? `all ${peaks?.n}` : n}
                                </option>
                              ))}
                            </select>
                            peaks
                          </label>
                        )}

                        <div className="flex overflow-hidden rounded border border-ink-200 text-xs">
                          {(['raw', 'oe'] as const).map((m) => (
                            <button
                              key={m}
                              onClick={() => setMode(m)}
                              className={`px-2.5 py-1 ${
                                mode === m ? 'bg-ink-600 text-white' : 'bg-white text-ink-600'
                              }`}
                            >
                              {m === 'raw' ? 'raw' : 'O/E'}
                            </button>
                          ))}
                        </div>

                        <a
                          href={`http://localhost:8000/api/export/${encodeURIComponent(
                            gene.gene_symbol,
                          )}.csv`}
                          className="rounded border border-ink-200 px-2 py-1 text-[11px] text-ink-600 hover:bg-ink-50"
                        >
                          export
                        </a>
                      </div>
                    </div>

                    {profile && (
                      <ProfilePlot
                        profile={profile}
                        peaks={
                          channel === 'mcc' ? (peaks?.peaks ?? []).slice(0, peakLimit) : []
                        }
                        loading={loading}
                        onZoom={setRange}
                        onPeakClick={setSelectedPeak}
                      />
                    )}

                    <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
                      Distance is signed bp from the <strong>experimental viewpoint</strong>,
                      not the canonical TSS. For about 2% of genes these differ by hundreds
                      of kb. Aligning every gene on its viewpoint is what makes profiles
                      comparable, and is why this is not a genome-browser view. Drag to zoom.
                      {channel !== 'mcc' && profile?.channel_role && (
                        <>
                          {' '}
                          <strong>{CHANNEL_LABEL[channel]}</strong> is {profile.channel_role}.
                        </>
                      )}
                    </p>
                  </div>

                  {gene.ps?.alpha_both != null && (
                    <div className="mt-5 rounded-lg border border-ink-200 bg-white p-4">
                      <h2 className="mb-2 text-xs font-semibold uppercase tracking-wide text-ink-500">
                        Contact decay P(s) ~ s<sup>−α</sup>
                      </h2>
                      <div className="flex flex-wrap gap-8 text-sm">
                        {(
                          [
                            ['α overall', gene.ps.alpha_both],
                            ['α near (10 to 100 kb)', gene.ps.alpha_near],
                            ['α far (0.1 to 1 Mb)', gene.ps.alpha_far],
                            ['fit R²', gene.ps.fit_r2_both],
                          ] as [string, number | undefined][]
                        ).map(([label, v]) =>
                          v == null ? null : (
                            <div key={label}>
                              <p className="text-[11px] text-ink-500">{label}</p>
                              <p className="font-mono text-base">{Number(v).toFixed(3)}</p>
                            </div>
                          ),
                        )}
                      </div>
                    </div>
                  )}
                </section>

                <aside className="space-y-5">
                  <ArchetypeReadout a={gene.archetype} />
                  {selectedPeak && (
                    <PeakDetail
                      peak={selectedPeak}
                      onClose={() => setSelectedPeak(null)}
                    />
                  )}
                </aside>

                {/* Where this gene sits among the rest. Placed above the
                    feature table because position is the headline and the
                    features are the detail behind it, and above the compare
                    box so that adding a gene there visibly moves the map. */}
                <section className="lg:col-span-3">
                  {embedding && (
                    <ContinuumMap
                      data={embedding}
                      axes={axes}
                      loading={embLoading}
                      error={embErr}
                      onAxisChange={(x, y) => setAxes([x, y])}
                      onPick={select}
                      compare={compareGenes}
                      subtitle={`${gene.gene_symbol} in black${
                        compareGenes.length ? `, compared with ${compareGenes.join(', ')}` : ''
                      }`}
                    />
                  )}
                </section>

                <section className="lg:col-span-3">
                  {features && <FeatureTable data={features} onPick={select} />}
                </section>

                <section className="lg:col-span-3">
                  <GeneCompare
                    primary={gene}
                    primaryFeatures={features}
                    onCompareChange={setCompareGenes}
                  />
                </section>
              </div>
            )}
          </>
        )}

        {/* ---------------------------------------------------------------- */}
        {/* PANEL: the coordinate system and the population. No single gene.  */}
        {/* ---------------------------------------------------------------- */}
        {/* Mounted on first visit and kept mounted, hidden with CSS rather than
            unmounted. Picking a gene here switches to the Gene page, and
            unmounting would throw away the search results that led you to that
            gene, so coming back would mean re-running the question. React keeps
            child state only while the component stays mounted, so `hidden` is
            doing real work here and is not a style choice. */}
        {panelVisited && (
          <div className="space-y-5" hidden={page !== 'panel'}>
            <p className="text-[11px] leading-relaxed text-ink-500">
              These views describe the panel and the coordinate system, not any one
              gene. {gene && (
                <>
                  <strong>{gene.gene_symbol}</strong> is highlighted on the map for
                  reference.{' '}
                </>
              )}
              Selecting a gene anywhere here opens it on the Gene page, and what you
              found here is kept.
            </p>

            <AskPanel onPick={select} />
            <RankedLists onPick={select} />
            <CohortView onPick={select} />
            {embedding && (
              <ContinuumMap
                data={embedding}
                axes={axes}
                loading={embLoading}
                error={embErr}
                onAxisChange={(x, y) => setAxes([x, y])}
                onPick={select}
              />
            )}
            <DimensionPanel />
            <ExplainPanel />
          </div>
        )}

        {page === 'lab' && <LabPage />}
      </main>
    </div>
  )
}
