/**
 * One peak, in detail, the J1 "tallest peak zoom" panel, made selectable.
 *
 * Everything shown is measured per peak in annotated.tsv. Position is the
 * signed offset computed in the store (peak_midpoint - viewpoint), not the
 * source table's unsigned distance.
 */

import { formatBp, type Peak } from './api'

const ELEMENT_COLOR: Record<string, string> = {
  enhancer: '#c2703d',
  ctcf: '#4a7c59',
  promoter: '#2b5070',
}

const FIELDS: [keyof Peak, string, (v: number) => string][] = [
  ['peak_max', 'height', (v) => v.toFixed(1)],
  ['peak_size', 'width', (v) => `${v.toFixed(0)} bp`],
  ['sharpness', 'sharpness', (v) => v.toFixed(2)],
  ['log2_enrichment', 'log2 enrichment', (v) => (v >= 0 ? '+' : '') + v.toFixed(2)],
  ['consensus_fraction', 'consensus', (v) => v.toFixed(2)],
]

export function PeakDetail({ peak, onClose }: { peak: Peak; onClose: () => void }) {
  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex items-baseline justify-between">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Selected peak
        </h2>
        <button onClick={onClose} className="text-[11px] text-ink-500 hover:text-ink-800">
          clear
        </button>
      </div>

      <div className="mb-3 flex items-center gap-2">
        <span
          className="inline-block h-3 w-3 rounded-full"
          style={{ background: ELEMENT_COLOR[peak.re] ?? '#888' }}
        />
        <span className="text-sm font-semibold">{peak.re}</span>
        <span className="font-mono text-[11px] text-ink-500">
          {formatBp(peak.offset_bp)} from viewpoint
        </span>
      </div>

      <dl className="grid grid-cols-2 gap-x-4 gap-y-1.5">
        {FIELDS.map(([key, label, fmt]) => {
          const v = peak[key]
          if (typeof v !== 'number' || !Number.isFinite(v)) return null
          return (
            <div key={String(key)} className="flex items-baseline justify-between">
              <dt className="text-[11px] text-ink-600">{label}</dt>
              <dd className="font-mono text-[12px] text-ink-900">{fmt(v)}</dd>
            </div>
          )
        })}
      </dl>

      <p className="mt-3 border-t border-ink-100 pt-2 font-mono text-[10px] leading-relaxed text-ink-400">
        {peak.chromosome}:{peak.start?.toLocaleString()}–{peak.end?.toLocaleString()}
        <br />
        midpoint {peak.peak_midpoint?.toLocaleString()} · {peak.viewpoint_id}
      </p>

      <p className="mt-2 text-[10px] leading-relaxed text-ink-500">
        Consensus is the fraction of replicates supporting the call. Element class comes
        from the pre-called peak table and is ATAC-intersected by construction, so
        ATAC-based validation of these classes is not independent.
      </p>
    </div>
  )
}
