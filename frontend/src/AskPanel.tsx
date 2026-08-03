/**
 * Ask a question in words; get genes back.
 *
 * The model translates the sentence into a query and never sees the data. The
 * query renders BEFORE the results, as the filter chain with per-step counts,
 * so a misreading is visible as a wrong query rather than arriving disguised as
 * a gene list. Every filter is editable, and the whole panel works with no API
 * key, the ask box is the only part that needs one.
 */

import { useEffect, useState } from 'react'
import { api, type QueryResult, type Vocabulary } from './api'

const EXAMPLES = [
  'a gene with long-range contacts that is enhancer-driven and immune',
  'promoter-driven housekeeping genes',
  'genes with the most dispersed contacts',
]

export function AskPanel({ onPick }: { onPick?: (s: string) => void }) {
  const [vocab, setVocab] = useState<Vocabulary | null>(null)
  const [question, setQuestion] = useState('')
  const [query, setQuery] = useState<Record<string, unknown> | null>(null)
  const [result, setResult] = useState<QueryResult | null>(null)
  const [interpretation, setInterpretation] = useState<string | null>(null)
  const [unsupported, setUnsupported] = useState<string | null>(null)
  const [reasoning, setReasoning] = useState<string[] | null>(null)
  // Which stage the request is at, and how long it has been there. The
  // translation is a network round trip to a model and takes seconds; a button
  // that only greys out gives no way to tell "thinking" from "hung".
  const [stage, setStage] = useState<'translating' | 'filtering' | null>(null)
  const [elapsed, setElapsed] = useState(0)

  useEffect(() => {
    if (!stage) return
    const t0 = Date.now()
    const id = setInterval(() => setElapsed((Date.now() - t0) / 1000), 100)
    return () => clearInterval(id)
  }, [stage])
  const [err, setErr] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)

  useEffect(() => {
    api.vocabulary().then(setVocab).catch(() => setVocab(null))
  }, [])

  async function ask() {
    if (!question.trim()) return
    setBusy(true)
    setErr(null)
    setStage('translating')
    setElapsed(0)
    // Cleared up front. A previous answer left on screen under a new spinner
    // reads as the answer to the new question.
    setQuery(null)
    setResult(null)
    setInterpretation(null)
    setUnsupported(null)
    setReasoning(null)
    try {
      const r = await api.ask(question.trim())
      setStage('filtering')
      setQuery(r.query ?? null)
      setInterpretation(r.interpretation ?? null)
      setUnsupported(r.unsupported ?? null)
      setReasoning(r.reasoning ?? null)
      setResult(r.error ? null : r)
      if (r.error) setErr(r.error)
    } catch (e) {
      setErr(String((e as Error).message))
    } finally {
      setBusy(false)
      setStage(null)
    }
  }

  async function rerun(q: Record<string, unknown>) {
    setBusy(true)
    setErr(null)
    try {
      setResult(await api.query(q))
    } catch (e) {
      setErr(String((e as Error).message))
      setResult(null)
    } finally {
      setBusy(false)
    }
  }

  function dropFilter(i: number) {
    if (!query) return
    const filters = [...((query.filters as unknown[]) ?? [])]
    filters.splice(i, 1)
    const next = { ...query, filters }
    setQuery(next)
    rerun(next)
  }

  const keyMissing = err?.includes('credentials') || err?.includes('API_KEY')

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Find a gene by describing it
      </h2>
      <p className="mb-3 text-[11px] leading-relaxed text-ink-500">
        The model turns your sentence into a query and never sees the data. The query
        is shown below the box, check it before trusting the genes, and edit it if it
        misread you.
      </p>

      <div className="flex gap-2">
        <input
          value={question}
          onChange={(e) => setQuestion(e.target.value)}
          onKeyDown={(e) => e.key === 'Enter' && ask()}
          placeholder="e.g. a gene with long-range contacts that is enhancer-driven and immune"
          className="flex-1 rounded border border-ink-200 px-3 py-2 text-[13px] outline-none focus:border-ink-400"
        />
        <button
          onClick={ask}
          disabled={busy}
          className="rounded bg-ink-600 px-4 py-2 text-xs text-white hover:bg-ink-700 disabled:opacity-50"
        >
          {busy ? '…' : 'Ask'}
        </button>
      </div>

      <div className="mt-1.5 flex flex-wrap gap-1.5">
        {EXAMPLES.map((ex) => (
          <button
            key={ex}
            onClick={() => setQuestion(ex)}
            className="rounded border border-ink-100 px-1.5 py-0.5 text-[10px] text-ink-500 hover:border-ink-300"
          >
            {ex}
          </button>
        ))}
      </div>

      {/* Two stages, both named, with a running clock. The stages are real:
          translating is a call to the model, filtering is local pandas and is
          effectively instant, which is itself worth showing because it makes
          clear where the time actually goes. */}
      {stage && (
        <div className="mt-3 rounded border border-ink-200 bg-ink-50 px-3 py-2">
          <div className="flex items-center gap-2 text-[11px] text-ink-700">
            <span className="inline-block h-3 w-3 animate-spin rounded-full border-2 border-ink-300 border-t-ink-700" />
            <span className="font-medium">
              {stage === 'translating'
                ? 'Translating your question into a query'
                : 'Running the filters'}
            </span>
            <span className="ml-auto font-mono text-ink-400">
              {elapsed.toFixed(1)}s
            </span>
          </div>
          <div className="mt-1.5 space-y-0.5 text-[10px] text-ink-500">
            <div className={stage === 'translating' ? 'text-ink-700' : ''}>
              1. The model reads your wording and picks axes, features and
              directions. It never sees the data.
            </div>
            <div className={stage === 'filtering' ? 'text-ink-700' : ''}>
              2. The query runs locally over {vocab ? '1,846' : 'all'} genes. No
              model involved, so the same query always returns the same genes.
            </div>
          </div>
        </div>
      )}

      {err && (
        <div
          className={`mt-3 rounded border px-3 py-2 text-[11px] leading-relaxed ${
            keyMissing
              ? 'border-ink-300 bg-ink-50 text-ink-700'
              : 'border-element-enhancer/40 bg-element-enhancer/5 text-ink-700'
          }`}
        >
          {err}
          {keyMissing && vocab && (
            <span className="mt-1 block text-ink-500">
              {vocab.axes.length} axes, {vocab.cohorts.length} reference sets and{' '}
              {vocab.features.length} features are queryable directly. Only the
              plain-English box needs a key.
            </span>
          )}
        </div>
      )}

      {/* the translation, before the results */}
      {query && (
        <div className="mt-4 rounded border border-ink-200 bg-ink-50 p-3">
          <h3 className="mb-2 text-[11px] font-medium text-ink-700">
            How your question was read
          </h3>
          {interpretation && (
            <p className="mb-2 text-[12px] italic text-ink-700">“{interpretation}”</p>
          )}
          {/* Shown above the filter chain, not below the results. A gap in what
              the data can answer changes how you read the genes, so it has to
              arrive before them. */}
          {unsupported && (
            <p className="mb-2 rounded border border-element-enhancer/40 bg-element-enhancer/5 px-2 py-1.5 text-[11px] leading-relaxed text-ink-700">
              <strong className="text-element-enhancer">Not used for selection.</strong>{' '}
              {unsupported}
            </p>
          )}

          {/* Called "working", not "thinking". This is a rationale the model
              wrote alongside its answer, not a record of how it produced one,
              and the difference matters for how much weight to give it. */}
          {reasoning && reasoning.length > 0 && (
            <details className="mb-2 rounded border border-ink-100 bg-white px-2 py-1.5">
              <summary className="cursor-pointer text-[11px] font-medium text-ink-600">
                Its working, {reasoning.length} steps
              </summary>
              <ol className="mt-1.5 space-y-1 pl-4 text-[11px] leading-relaxed text-ink-600">
                {reasoning.map((r, i) => (
                  <li key={i} className="list-decimal">{r}</li>
                ))}
              </ol>
              <p className="mt-1.5 border-t border-ink-100 pt-1.5 text-[10px] text-ink-400">
                This is the model's stated rationale, written alongside the query
                rather than a record of how it arrived at one. The query below is
                the part that actually ran.
              </p>
            </details>
          )}
          <div className="space-y-1">
            {result?.steps?.map((s, i) => (
              <div key={i} className="flex items-center gap-2 text-[11px]">
                <span className="flex-1 text-ink-700">{s.reads_as}</span>
                <span className="font-mono text-ink-500">
                  {s.before.toLocaleString()} → {s.after.toLocaleString()}
                </span>
                <button
                  onClick={() => dropFilter(i)}
                  className="rounded border border-ink-200 px-1.5 text-[10px] text-ink-500 hover:border-element-enhancer hover:text-element-enhancer"
                  title="Remove this condition and re-run"
                >
                  drop
                </button>
              </div>
            ))}
          </div>
        </div>
      )}

      {result && (
        <div className="mt-3">
          <p className="mb-2 text-[12px] text-ink-700">
            <strong>{result.n_matched}</strong> gene{result.n_matched === 1 ? '' : 's'}{' '}
            match
          </p>
          {result.note && (
            <p className="mb-2 rounded border border-element-enhancer/40 bg-element-enhancer/5 px-2.5 py-1.5 text-[11px] text-ink-700">
              {result.note}
            </p>
          )}
          <div className="flex flex-wrap gap-1">
            {result.genes.map((g) => (
              <button
                key={g.gene_id}
                onClick={() => onPick?.(g.gene_symbol)}
                className="rounded border border-ink-200 px-2 py-1 text-[11px] hover:border-ink-400"
                title={
                  [g.group, g.in_sets?.length ? `in: ${g.in_sets.join(', ')}` : null]
                    .filter(Boolean)
                    .join(' · ') || undefined
                }
              >
                {g.gene_symbol}
                {g.in_sets && g.in_sets.length > 0 && (
                  <span className="ml-1 font-mono text-[9px] text-ink-400">
                    {g.in_sets.length}
                  </span>
                )}
              </button>
            ))}
          </div>

          {/* Overlap with external sets, counted AFTER selection. Reported as a
              tally rather than a filter, because these sets choosing the genes
              is exactly what the tool must not do. */}
          {(() => {
            const tally = new Map<string, number>()
            for (const g of result.genes)
              for (const s of g.in_sets ?? []) tally.set(s, (tally.get(s) ?? 0) + 1)
            const top = [...tally.entries()].sort((a, b) => b[1] - a[1]).slice(0, 6)
            if (!top.length) return null
            return (
              <div className="mt-2 rounded border border-ink-100 bg-ink-50/60 px-2.5 py-1.5">
                <p className="text-[10px] font-medium uppercase tracking-wide text-ink-500">
                  What these genes turned out to be
                  {result.n_matched > result.genes.length && (
                    <span className="ml-1 font-normal normal-case tracking-normal text-ink-400">
                      (over the {result.genes.length} shown, not all{' '}
                      {result.n_matched} matched)
                    </span>
                  )}
                </p>
                <p className="mt-1 flex flex-wrap gap-x-3 gap-y-0.5 text-[11px] text-ink-600">
                  {top.map(([s, n]) => (
                    <span key={s}>
                      <span className="font-mono">{n}</span>/{result.genes.length} {s}
                    </span>
                  ))}
                </p>
                <p className="mt-1 text-[10px] leading-relaxed text-ink-400">
                  Counted after the fact. These sets played no part in choosing the
                  genes, which is what makes the overlap worth reading.
                </p>
              </div>
            )
          })()}

          <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            Retrieval is a deterministic filter over the panel, the same query always
            returns the same genes. The model only produced the query. Selection uses
            contact architecture only.
          </p>
        </div>
      )}
    </div>
  )
}
