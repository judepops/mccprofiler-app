/**
 * Cohort view — where does a set of genes sit on each axis?
 *
 * The strongest non-circular evidence in the project. The archetypes are KMeans
 * clusters computed ON these features, so their separation needs a
 * within/pooled-ratio argument to be interpretable. These groups were defined
 * entirely outside the MCC data, so no such argument is needed.
 *
 * Two things this view refuses to do quietly:
 *  - plot a pasted list without saying how much of it is in the panel;
 *  - offer a set too small to say anything about.
 */

import { useEffect, useState } from 'react'
import { api, type CohortCompare, type CohortList } from './api'

export function CohortView({ onPick }: { onPick?: (symbol: string) => void }) {
  const [list, setList] = useState<CohortList | null>(null)
  const [selected, setSelected] = useState<string>('Eisenberg_HK')
  const [pasted, setPasted] = useState('')
  const [result, setResult] = useState<CohortCompare | null>(null)
  const [err, setErr] = useState<string | null>(null)

  useEffect(() => {
    api.cohorts().then(setList).catch((e) => setErr(String(e.message)))
  }, [])

  useEffect(() => {
    if (!selected) return
    setErr(null)
    api
      .cohortCompare({ group: selected })
      .then(setResult)
      .catch((e) => setErr(String(e.message)))
  }, [selected])

  function runPasted() {
    const syms = pasted.trim()
    if (!syms) return
    setErr(null)
    setSelected('')
    api
      .cohortCompare({ symbols: syms })
      .then(setResult)
      .catch((e) => setErr(String(e.message)))
  }

  const maxD = result ? Math.max(...result.axes.map((a) => Math.abs(a.cohen_d)), 0.2) : 1

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <h2 className="mb-1 text-xs font-semibold uppercase tracking-wide text-ink-500">
        Cohort view — externally-defined gene sets
      </h2>
      <p className="mb-4 text-[11px] leading-relaxed text-ink-500">
        {list?.notes.why_external}
      </p>

      <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
        {/* built-in sets ---------------------------------------------------- */}
        <div>
          <label className="mb-1 block text-[11px] font-medium text-ink-600">
            Reference set
          </label>
          <select
            value={selected}
            onChange={(e) => setSelected(e.target.value)}
            className="w-full rounded border border-ink-200 px-2 py-1.5 text-[13px]"
          >
            <option value="">— pick a set —</option>
            {list?.rows
              .filter((r) => r.usable)
              .map((r) => (
                <option key={r.group} value={r.group}>
                  {r.group} ({r.n_in_panel})
                  {r.is_positive_control ? ' — positive control' : ''}
                  {r.is_super_enhancer ? ' — SE, comparison only' : ''}
                </option>
              ))}
          </select>
          {list && list.n_filtered_out > 0 && (
            <p className="mt-1 text-[10px] text-ink-400">
              {list.n_filtered_out} set(s) hidden: fewer than {list.min_group_n} panel
              genes, so any enrichment would be noise.
            </p>
          )}
        </div>

        {/* pasted list ------------------------------------------------------ */}
        <div>
          <label className="mb-1 block text-[11px] font-medium text-ink-600">
            …or paste your own gene list
          </label>
          <div className="flex gap-2">
            <textarea
              value={pasted}
              onChange={(e) => setPasted(e.target.value)}
              placeholder="IL7R, CTCF, EEF1A1…"
              rows={2}
              className="flex-1 rounded border border-ink-200 px-2 py-1.5 font-mono text-[12px]"
            />
            <button
              onClick={runPasted}
              className="self-start rounded bg-ink-600 px-3 py-1.5 text-xs text-white hover:bg-ink-700"
            >
              Place
            </button>
          </div>
        </div>
      </div>

      {err && <p className="mt-3 text-[12px] text-element-enhancer">{err}</p>}

      {result && (
        <div className="mt-5">
          <div className="mb-3 flex flex-wrap items-baseline gap-x-3 gap-y-1">
            <span className="text-sm font-semibold">{result.label}</span>
            <span className="font-mono text-[11px] text-ink-500">
              {result.n_in_panel} panel genes
            </span>
          </div>

          {/* coverage — leads, because it decides whether the rest means anything */}
          {result.coverage && (
            <div className="mb-3 rounded border border-ink-200 bg-ink-50 px-3 py-2 text-[11px] leading-relaxed text-ink-700">
              <strong>
                {result.coverage.in_panel} of {result.coverage.requested} genes are in the
                panel.
              </strong>{' '}
              {result.coverage.note}
              {result.coverage.matched.length > 0 && (
                <span className="mt-1.5 flex flex-wrap gap-1">
                  {result.coverage.matched.map((sym) => (
                    <button
                      key={sym}
                      onClick={() => onPick?.(sym)}
                      className="rounded border border-ink-200 bg-white px-1.5 py-0.5 font-mono text-[10px] hover:border-ink-400"
                    >
                      {sym}
                    </button>
                  ))}
                </span>
              )}
              {result.coverage.missing.length > 0 && (
                <span className="mt-1 block font-mono text-ink-500">
                  not captured: {result.coverage.missing.join(', ')}
                </span>
              )}
            </div>
          )}

          {result.small_set_warning && (
            <div className="mb-3 rounded border border-element-enhancer/40 bg-element-enhancer/5 px-3 py-2 text-[11px] text-ink-700">
              {result.small_set_warning}
            </div>
          )}

          {/* diverging effect-size bars */}
          <div className="space-y-2">
            {result.axes.map((a) => {
              const frac = Math.abs(a.cohen_d) / maxD
              const pos = a.cohen_d >= 0
              return (
                <div key={a.axis} className="flex items-center gap-2">
                  <span className="w-56 shrink-0 truncate text-[11px] text-ink-600" title={a.label}>
                    {a.label}
                  </span>
                  <div className="relative h-4 flex-1">
                    <div className="absolute inset-y-0 left-1/2 w-px bg-ink-200" />
                    <div
                      className="absolute inset-y-0 rounded-sm"
                      style={{
                        left: pos ? '50%' : `${50 - frac * 50}%`,
                        width: `${frac * 50}%`,
                        background: pos ? '#c2703d' : '#2b5070',
                        opacity: 0.75,
                      }}
                    />
                  </div>
                  <span className="w-14 shrink-0 text-right font-mono text-[11px] text-ink-600">
                    {a.cohen_d >= 0 ? '+' : ''}
                    {a.cohen_d.toFixed(3)}
                  </span>
                </div>
              )
            })}
          </div>

          <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            Bars are standardised mean differences (Cohen's d) against the rest of the
            panel, sorted by magnitude. Effect sizes rather than p-values: with n in the
            hundreds almost any difference reaches significance, so significance carries
            no information here.
          </p>
        </div>
      )}

      {list?.notes && (
        <details className="mt-4 border-t border-ink-100 pt-3">
          <summary className="cursor-pointer text-[11px] font-medium text-ink-600">
            How to read this
          </summary>
          <div className="mt-2 space-y-2 text-[11px] leading-relaxed text-ink-500">
            <p>
              <strong>Positive control.</strong> {list.notes.positive_control}
            </p>
            <p>
              <strong>Display.</strong> {list.notes.display}
            </p>
            <p>
              <strong>Super-enhancers.</strong> {list.notes.super_enhancer}
            </p>
          </div>
        </details>
      )}
    </div>
  )
}
