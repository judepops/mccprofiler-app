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

export interface EmbeddingPoint {
  gene_id: string
  gene_symbol: string
  x: number
  y: number
  group: string | null
  max_posterior: number | null
  is_core: boolean | null
}

export interface Axis {
  key: string
  label: string
  /** Which END is which. PCA sign is arbitrary, so the contrast in the label
   *  ("local vs long-range") does not say which side of the plot is local. */
  poles: { neg: string; pos: string } | null
}

export interface Embedding {
  x_axis: Axis
  y_axis: Axis
  axes_available: { key: string; label: string }[]
  highlight: string | null
  is_umap: boolean
  is_null: boolean
  umap_caveat: string
  continuum_caveat: string
  n: number
  points: EmbeddingPoint[]
}

export interface CohortRow {
  group: string
  n_in_panel: number
  usable: boolean
  is_positive_control: boolean
  is_super_enhancer: boolean
  stratification: Record<string, { pct_retained: number; signal_over_random: number }> | null
}

export interface CohortList {
  min_group_n: number
  n_offered: number
  n_filtered_out: number
  positive_control: string
  notes: {
    why_external: string
    positive_control: string
    display: string
    super_enhancer: string
  }
  rows: CohortRow[]
}

export interface CohortCompare {
  label: string
  n_in_panel: number
  coverage: {
    requested: number
    in_panel: number
    matched: string[]
    missing: string[]
    note: string
  } | null
  small_set_warning: string | null
  axes: {
    axis: string
    label: string
    n_in: number
    cohen_d: number
    in_quartiles: number[]
    out_quartiles: number[]
  }[]
}

export interface ScreeRow {
  pc: number
  variance_pct: number
  cumulative_pct: number
  noise_pct: number
  above_noise: boolean
}

export interface Scree {
  n_above_noise: number
  method: string
  note: string
  rows: ScreeRow[]
}

export interface Loadings {
  pc: number
  label: string
  variance_pct: number | null
  above_noise: boolean | null
  n_features: number
  loadings: { pc: number; feature: string; loading: number }[]
}

export interface Repro {
  caveat: string
  n_genes: number | null
  n_features: number
  median_rho: number
  n_above_0_7: number
  n_below_0_3: number
  note: string
  per_feature: { feature: string; pearson_r: number; spearman_rho: number; n: number }[]
  pairs?: {
    feature: string
    points: { symbol_key: string; gw: number; immune: number }[]
  }
}

export interface Ranked {
  axis: string
  label: string
  poles: { neg: string; pos: string } | null
  n: number
  caveat: string
  top?: { gene_id: string; gene_symbol: string; value: number; percentile: number; group: string | null }[]
  bottom?: { gene_id: string; gene_symbol: string; value: number; percentile: number; group: string | null }[]
  top_pole?: string | null
  bottom_pole?: string | null
}

export interface ExplainSection {
  claim: string
  detail: string
  citations?: string[]
  table?: Record<string, unknown>[]
}

export interface Explain {
  why_these_features: ExplainSection & { table?: Record<string, unknown>[] }
  why_not_clusters: ExplainSection
  why_name_regions_at_all: ExplainSection
  what_it_is_not: ExplainSection
  noise_floor: ExplainSection
  atac_is_not_a_feature: ExplainSection
  provenance: {
    panel: string
    n_genes: number
    built: string
    scripts_cleaned_commit: string | null
  }
  [k: string]: unknown
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

  labReproducibility: (feature?: string) =>
    get<Repro>('/api/lab/reproducibility', { feature }),

  ranked: (axis: string, limit = 25) =>
    get<Ranked>('/api/ranked', { axis, limit }),

  explain: () => get<Explain>('/api/explain'),

  scree: () => get<Scree>('/api/dimensions/scree'),

  loadings: (pc: number, top = 15) =>
    get<Loadings>(`/api/dimensions/${pc}/loadings`, { top }),

  cohorts: () => get<CohortList>('/api/cohorts'),

  cohortCompare: (opts: { group?: string; symbols?: string }) =>
    get<CohortCompare>('/api/cohorts/compare', opts),

  embedding: (x: string, y: string, highlight?: string) =>
    get<Embedding>('/api/embedding', { x, y, highlight }),

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
