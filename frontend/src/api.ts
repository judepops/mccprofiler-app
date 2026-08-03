/**
 * Client for the local FastAPI backend.
 *
 * The API returns numbers, never images — every figure in this app is drawn
 * client-side from these payloads. See PLAN.md §1b.
 */

const BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8000'

export interface Health {
  ok: boolean
  store_version: number
  built: string
  panel: string
  n_genes: number
  channels: string[]
  scripts_cleaned_commit: string | null
}

export interface GeneSummary {
  gene_id: string
  gene_symbol: string
  group: string | null
  max_posterior: number | null
  confidence_class: string | null
  viewpoint_chrom: string | null
  viewpoint_pos: number | null
}

export interface Archetype {
  group: string | null
  /** Architecture-derived display name. Makes no functional claim. */
  display: string | null
  architecture: string | null
  /** Continuum caveat — render it wherever the label is rendered. */
  caveat: string | null
  max_posterior: number | null
  confidence_class: string | null
  is_core: boolean
  /** 44% of active genes are blends. Never render a bare label without this. */
  is_mixture: boolean
  entropy: number | null
  mixture: Record<string, number>
  /** arch-off is a QC class (near-empty signal), not a biological archetype. */
  is_qc_class: boolean
}

export interface PsFit {
  alpha_both?: number
  alpha_near?: number
  alpha_far?: number
  alpha_asymmetry?: number
  fit_r2_both?: number
  [k: string]: number | string | null | undefined
}

export interface Gene {
  gene_id: string
  gene_symbol: string
  viewpoint: { chrom: string | null; pos: number | null }
  archetype: Archetype
  ps?: PsFit
}

export interface Band {
  name: string
  start_bp: number
  end_bp: number
}

export interface Profile {
  gene_id: string
  channel: string
  mode: 'raw' | 'oe'
  level: number
  bp_per_bin: number
  start_bp: number
  end_bp: number
  n: number
  values: number[]
  bands: Band[]
  channel_role: string | null
}

export interface Peak {
  symbol_key: string
  viewpoint_id: string
  chromosome: string
  start: number
  end: number
  peak_midpoint: number
  /** SIGNED bp from the viewpoint. The source table's distance is unsigned. */
  offset_bp: number
  peak_max: number
  peak_size: number
  sharpness: number
  log2_enrichment: number
  consensus_fraction: number
  re: 'enhancer' | 'ctcf' | 'promoter'
}

export interface PeakSet {
  gene_id: string
  n: number
  by_class: Record<string, number>
  colors: Record<string, string>
  peaks: Peak[]
}

export interface FeatureBlock {
  block: string
  description: string
  features: { name: string; z: number; percentile: number }[]
}

export interface FeatureSet {
  gene_id: string
  n_features: number
  blocks: FeatureBlock[]
}

async function get<T>(path: string, params?: Record<string, unknown>): Promise<T> {
  const url = new URL(BASE + path)
  for (const [k, v] of Object.entries(params ?? {})) {
    if (v !== undefined && v !== null) url.searchParams.set(k, String(v))
  }
  const res = await fetch(url)
  if (!res.ok) {
    let detail = res.statusText
    try {
      detail = (await res.json()).detail ?? detail
    } catch {
      /* response was not JSON; keep the status text */
    }
    throw new Error(detail)
  }
  return res.json() as Promise<T>
}

export const api = {
  health: () => get<Health>('/api/health'),

  searchGenes: (q: string, limit = 25) =>
    get<{ n: number; genes: GeneSummary[] }>('/api/genes', { q, limit }),

  gene: (gene: string) => get<Gene>(`/api/genes/${encodeURIComponent(gene)}`),

  peaks: (gene: string) => get<PeakSet>(`/api/genes/${encodeURIComponent(gene)}/peaks`),

  features: (gene: string) =>
    get<FeatureSet>(`/api/genes/${encodeURIComponent(gene)}/features`),

  profile: (
    gene: string,
    opts: {
      channel?: string
      mode?: 'raw' | 'oe'
      level?: number
      start_bp?: number
      end_bp?: number
    } = {},
  ) => get<Profile>(`/api/genes/${encodeURIComponent(gene)}/profile`, opts),
}

/** Pick the coarsest pyramid level that still gives ~2 bins per output pixel. */
export function levelForSpan(spanBp: number, pixelWidth: number): number {
  const bpPerPixel = spanBp / Math.max(pixelWidth, 1)
  if (bpPerPixel < 125) return 0 // 50 bp bins
  if (bpPerPixel < 500) return 1 // 250 bp bins
  return 2 // 1 kb bins
}

export function formatBp(bp: number): string {
  const a = Math.abs(bp)
  const sign = bp < 0 ? '−' : bp > 0 ? '+' : ''
  if (a >= 1_000_000) return `${sign}${(a / 1_000_000).toFixed(a % 1_000_000 ? 2 : 0)} Mb`
  if (a >= 1_000) return `${sign}${(a / 1_000).toFixed(a % 1_000 ? 1 : 0)} kb`
  return `${sign}${a} bp`
}
