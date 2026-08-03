import { useCallback, useEffect, useRef, useState } from 'react'
import {
  api,
  levelForSpan,
  type Gene,
  type GeneSummary,
  type Embedding,
  type FeatureSet,
  type Health,
  type PeakSet,
  type Profile,
} from './api'
import { ProfilePlot } from './ProfilePlot'
import { ArchetypeReadout } from './ArchetypeReadout'
import { FeatureTable } from './FeatureTable'
import { ContinuumMap } from './ContinuumMap'
import { CohortView } from './CohortView'
import { DimensionPanel } from './DimensionPanel'

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
  const [embedding, setEmbedding] = useState<Embedding | null>(null)
  const [axes, setAxes] = useState<[string, string]>(['pc2', 'pc3'])
  const [peakLimit, setPeakLimit] = useState(25)
  const [loading, setLoading] = useState(false)

  const plotWidth = useRef(900)

  useEffect(() => {
    api.health().then(setHealth).catch((e) => setErr(String(e.message ?? e)))
  }, [])

  // Debounced search.
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
        const p = await api.profile(symbol, {
          channel,
          mode,
          level: levelForSpan(span, plotWidth.current),
          start_bp: r?.start,
          end_bp: r?.end,
        })
        setProfile(p)
      } catch (e) {
        setErr(String((e as Error).message))
      } finally {
        setLoading(false)
      }
    },
    [channel, mode],
  )

  async function select(symbol: string) {
    setQuery('')
    setHits([])
    setRange(null)
    setErr(null)
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
    api
      .embedding(axes[0], axes[1], gene?.gene_symbol)
      .then(setEmbedding)
      .catch(() => setEmbedding(null))
  }, [axes, gene])

  // Re-fetch when channel, mode or zoom changes.
  useEffect(() => {
    if (gene) loadProfile(gene.gene_symbol, range)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [channel, mode, range])

  return (
    <div className="min-h-screen bg-ink-50 font-sans text-ink-900">
      <header className="border-b border-ink-200 bg-white">
        <div className="mx-auto flex max-w-6xl items-center justify-between px-6 py-3">
          <div>
            <h1 className="text-sm font-semibold tracking-tight">mccprofiler</h1>
            <p className="text-[11px] text-ink-500">
              Gene position in the MCC regulatory continuum · CD4+ T cells
            </p>
          </div>
          {health && (
            <p className="text-right font-mono text-[10px] leading-4 text-ink-400">
              {health.panel} · {health.n_genes.toLocaleString()} genes
              <br />
              store v{health.store_version} · {health.scripts_cleaned_commit}
            </p>
          )}
        </div>
      </header>

      <main className="mx-auto max-w-6xl px-6 py-6">
        {err && (
          <div className="mb-4 rounded border border-element-enhancer/30 bg-element-enhancer/5 px-3 py-2 text-[13px] text-element-enhancer">
            {err}
            {err.includes('no store') && (
              <span className="ml-1 text-ink-600">
                — run <code className="font-mono">build_store.py</code> first.
              </span>
            )}
          </div>
        )}

        {/* search ---------------------------------------------------------- */}
        <div className="relative mb-6">
          <input
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter' && hits.length) select(hits[0].gene_symbol)
            }}
            placeholder="Search a gene — try IL7R, CTCF, EEF1A1…"
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
                      {h.max_posterior != null && (
                        <span className="ml-2 text-ink-300">
                          p {h.max_posterior.toFixed(2)}
                        </span>
                      )}
                    </span>
                  </button>
                </li>
              ))}
            </ul>
          )}
        </div>

        {!gene && !err && (
          <div className="space-y-6">
            <p className="pt-8 text-center text-sm text-ink-400">
              Search a gene to see its contact architecture — or start from a gene set below.
            </p>
            <CohortView onPick={select} />
            <DimensionPanel />
          </div>
        )}

        {gene && (
          <div className="grid grid-cols-1 gap-5 lg:grid-cols-3">
            <section className="lg:col-span-2">
              <div className="rounded-lg border border-ink-200 bg-white p-4">
                <div className="mb-3 flex items-baseline justify-between">
                  <div>
                    <h2 className="text-lg font-semibold">{gene.gene_symbol}</h2>
                    <p className="font-mono text-[11px] text-ink-500">
                      viewpoint {gene.viewpoint.chrom}:
                      {gene.viewpoint.pos?.toLocaleString()}
                    </p>
                  </div>

                  <div className="flex gap-3">
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
                  />
                )}

                <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
                  Distance is signed bp from the <strong>experimental viewpoint</strong>, not
                  the canonical TSS — for ~2% of genes these differ by hundreds of kb.
                  Aligning every gene on its viewpoint is what makes profiles comparable, and
                  is why this is not a genome-browser view. Drag to zoom.
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
                  <div className="flex gap-8 text-sm">
                    {(
                      [
                        ['α overall', gene.ps.alpha_both],
                        ['α near (10–100 kb)', gene.ps.alpha_near],
                        ['α far (0.1–1 Mb)', gene.ps.alpha_far],
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
            </aside>

            {/* Gene-specific detail. */}
            <section className="lg:col-span-3">
              {features && <FeatureTable data={features} />}
            </section>

            {/* Panel-level below: these describe the coordinate system and the
                population, not this gene. */}
            <section className="lg:col-span-3 border-t border-ink-200 pt-6">
              <p className="mb-4 text-[11px] uppercase tracking-wide text-ink-400">
                Panel-level — the space {gene.gene_symbol} sits in
              </p>
              <div className="space-y-5">
                {embedding && (
                  <ContinuumMap
                    data={embedding}
                    onAxisChange={(x, y) => setAxes([x, y])}
                    onPick={select}
                  />
                )}
                <DimensionPanel />
                <CohortView onPick={select} />
              </div>
            </section>
          </div>
        )}
      </main>
    </div>
  )
}
