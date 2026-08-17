/**
 * What this panel is and is not. Rendered above everything, on every page.
 *
 * Why it is a banner and not a footnote. `gw_cd4_1` is routinely called the
 * "genome-wide panel", and that phrase changes how every number below it should
 * be read: the panel is genome-DISTRIBUTED but not genome-REPRESENTATIVE. It was
 * built by TSS extraction then ATAC-accessibility selection, and accessible
 * human promoters are overwhelmingly CpG-island and TATA-less, so the gating
 * that chose the genes also chose the promoter class.
 *
 * Collapsed by default to one line, because a wall of caveat gets scrolled past
 * and stops working. Expanded it carries the per-claim status table, which is
 * the part that matters: this bias threatens some conclusions and not others,
 * and a blanket warning invites a reader to discount everything or nothing.
 */

import { useEffect, useState } from 'react'
import { api, type PanelBias } from './api'

const STATUS_STYLE: Record<string, string> = {
  robust: 'bg-ink-100 text-ink-600',
  QUALIFIED: 'bg-element-enhancer/15 text-element-enhancer',
  'CANNOT BE ASKED': 'bg-element-enhancer/25 text-element-enhancer',
  'PARTLY testable': 'bg-ink-100 text-ink-600',
}

export function PanelBiasBanner() {
  const [d, setD] = useState<PanelBias | null>(null)
  const [open, setOpen] = useState(false)

  useEffect(() => {
    api.panelBias().then(setD).catch(() => setD(null))
  }, [])

  if (!d) return null

  return (
    <div className="border-b border-element-enhancer/30 bg-element-enhancer/5">
      <div className="mx-auto max-w-[1400px] px-4 py-2">
        <button
          onClick={() => setOpen((v) => !v)}
          className="flex w-full items-baseline gap-2 text-left"
        >
          <span className="shrink-0 text-[11px] font-semibold uppercase tracking-wide text-element-enhancer">
            Panel scope
          </span>
          <span className="min-w-0 flex-1 truncate text-[11px] text-ink-700">
            {d.headline} {d.one_line}
          </span>
          <span className="shrink-0 text-[11px] text-ink-500">
            {open ? 'hide' : 'detail'}
          </span>
        </button>

        {open && (
          <div className="mt-2 space-y-3 border-t border-element-enhancer/20 pt-2">
            <table className="text-[11px]">
              <thead>
                <tr className="text-left text-ink-500">
                  <th className="pr-4 pb-1 font-medium">promoter class</th>
                  <th className="pr-4 pb-1 text-right font-medium">in panel</th>
                  <th className="pr-4 pb-1 text-right font-medium">not in panel</th>
                  <th className="pb-1 font-medium" />
                </tr>
              </thead>
              <tbody>
                {d.measured.map((m) => (
                  <tr key={m.what}>
                    <td className="pr-4 py-0.5 text-ink-700">{m.what}</td>
                    <td className="pr-4 py-0.5 text-right font-mono text-ink-800">
                      {m.panel}
                    </td>
                    <td className="pr-4 py-0.5 text-right font-mono text-ink-600">
                      {m.not_in_panel}
                    </td>
                    <td className="py-0.5 text-ink-500">{m.note}</td>
                  </tr>
                ))}
              </tbody>
            </table>

            <p className="max-w-4xl text-[11px] leading-relaxed text-ink-600">
              <strong className="text-ink-800">Cause.</strong> {d.cause}
            </p>
            <p className="max-w-4xl text-[11px] leading-relaxed text-ink-600">
              <strong className="text-ink-800">Framing.</strong> {d.framing}
            </p>
            <p className="max-w-4xl text-[11px] leading-relaxed text-ink-600">
              <strong className="text-ink-800">Not about.</strong> {d.not_about}
            </p>

            {/* The important part. A bias that threatens everything equally is a
                bias nobody can act on. */}
            <div>
              <p className="mb-1 text-[11px] font-medium text-ink-800">
                What this does and does not threaten
              </p>
              <ul className="space-y-1">
                {d.affects.map((a) => (
                  <li key={a.claim} className="flex gap-2 text-[11px] leading-relaxed">
                    <span
                      className={`mt-px shrink-0 rounded px-1.5 py-0.5 text-[9px] font-medium uppercase ${
                        STATUS_STYLE[a.status] ?? 'bg-ink-100 text-ink-600'
                      }`}
                    >
                      {a.status}
                    </span>
                    <span className="text-ink-600">
                      <span className="text-ink-800">{a.claim}</span> — {a.why}
                    </span>
                  </li>
                ))}
              </ul>
            </div>

            <div className="rounded border border-element-enhancer/30 bg-white/60 p-2">
              <p className="text-[11px] font-medium text-element-enhancer">
                {d.consequence}
              </p>
              <ul className="mt-1 list-inside list-disc space-y-0.5 text-[11px] leading-relaxed text-ink-600">
                {d.why_more_data.map((w, i) => (
                  <li key={i}>{w}</li>
                ))}
              </ul>
            </div>

            <p className="max-w-4xl text-[10px] leading-relaxed text-ink-500">
              {d.calibration}
            </p>
            <p className="max-w-4xl text-[10px] leading-relaxed text-ink-500">
              <strong>Scope.</strong> {d.scope}
            </p>
          </div>
        )}
      </div>
    </div>
  )
}
