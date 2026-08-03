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
  const [err, setErr] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)

  useEffect(() => {
    api.vocabulary().then(setVocab).catch(() => setVocab(null))
  }, [])

  async function ask() {
    if (!question.trim()) return
    setBusy(true)
    setErr(null)
    try {
      const r = await api.ask(question.trim())
      setQuery(r.query ?? null)
      setInterpretation(r.interpretation ?? null)
      setResult(r.error ? null : r)
      if (r.error) setErr(r.error)
    } catch (e) {
      setErr(String((e as Error).message))
    } finally {
      setBusy(false)
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
                title={g.group ?? undefined}
              >
                {g.gene_symbol}
              </button>
            ))}
          </div>
          <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            Retrieval is a deterministic filter over the panel, the same query always
            returns the same genes. The model only produced the query.
          </p>
        </div>
      )}
    </div>
  )
}
