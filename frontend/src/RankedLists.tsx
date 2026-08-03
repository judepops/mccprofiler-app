/**
 * Given an axis, which genes sit at its extremes?
 *
 * The inverse of the cohort view, and the answer to "I don't have a gene in
 * mind, show me something interesting". Both poles render side by side, named,
 * so the axis reads as a contrast rather than a ranking with a good end and a
 * bad end.
 */

import { useEffect, useState } from 'react'
import { api, type Ranked } from './api'

const GROUP_DISPLAY: Record<string, string> = {
  'arch-HK': 'dispersed',
  'arch-ME-constitutive': 'promoter-local',
  'arch-ME-effector': 'enhancer-focal',
  'arch-sparse': 'sparse',
  'arch-off': 'empty (QC)',
}

const AXES = ['pc1', 'pc2', 'pc3', 'pc4', 'pc5'] as const

function Column({
  genes,
  pole,
  onPick,
}: {
  genes: Ranked['top']
  pole: string | null | undefined
  onPick?: (s: string) => void
}) {
  return (
    <div>
      <h3 className="mb-1.5 text-[11px] font-medium text-ink-700">{pole ?? '-'}</h3>
      <div className="space-y-0.5">
        {genes?.map((g) => (
          <button
            key={g.gene_id}
            onClick={() => onPick?.(g.gene_symbol)}
            className="flex w-full items-baseline gap-2 rounded px-1.5 py-1 text-left hover:bg-ink-50"
          >
            <span className="w-24 shrink-0 truncate text-[12px] font-medium">
              {g.gene_symbol}
            </span>
            <span className="w-14 shrink-0 text-right font-mono text-[10px] text-ink-600">
              {g.value >= 0 ? '+' : ''}
              {g.value.toFixed(2)}
            </span>
            <span className="w-12 shrink-0 text-right font-mono text-[10px] text-ink-400">
              {g.percentile.toFixed(0)}%
            </span>
            <span className="truncate text-[10px] text-ink-500">
              {GROUP_DISPLAY[g.group ?? ''] ?? g.group}
            </span>
          </button>
        ))}
      </div>
    </div>
  )
}

export function RankedLists({ onPick }: { onPick?: (symbol: string) => void }) {
  const [axis, setAxis] = useState<string>('pc3')
  const [data, setData] = useState<Ranked | null>(null)

  useEffect(() => {
    api.ranked(axis, 12).then(setData).catch(() => setData(null))
  }, [axis])

  return (
    <div className="rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-3 flex flex-wrap items-baseline justify-between gap-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Extremes of an axis
        </h2>
        <select
          value={axis}
          onChange={(e) => setAxis(e.target.value)}
          className="rounded border border-ink-200 px-2 py-1 text-[11px]"
        >
          {AXES.map((a) => (
            <option key={a} value={a}>
              {a.toUpperCase()}
            </option>
          ))}
        </select>
      </div>

      {data && (
        <>
          <p className="mb-3 text-[11px] text-ink-600">{data.label}</p>
          <div className="grid grid-cols-1 gap-6 md:grid-cols-2">
            <Column genes={data.bottom} pole={data.bottom_pole} onPick={onPick} />
            <Column genes={data.top} pole={data.top_pole} onPick={onPick} />
          </div>
          <p className="mt-3 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-500">
            {data.caveat}
          </p>
        </>
      )}
    </div>
  )
}
