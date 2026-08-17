/**
 * Client for the local FastAPI backend.
 *
 * The API returns numbers, never images, every figure in this app is drawn
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
  /** Continuum caveat, render it wherever the label is rendered. */
  caveat: string | null
  /** True since 2026-08-16: labels were fit on the earlier 91-feature substrate
   *  and have not been re-derived. Render the label demoted, never as a badge. */
  provisional?: boolean
  provisional_note?: string | null
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

/** Panel median / p10 / p90 for each P(s) exponent. An exponent without its
 *  reference distribution is unreadable: alpha_far 0.204 looks like a finding
 *  until you know the panel median is 1.291. */
export type PsReference = Record<
  string,
  { median: number; p10: number; p90: number }
>

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
  /** Panel median / p10 / p90 for each exponent, so the numbers are readable. */
  ps_reference?: PsReference
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
  /** Superseded 2026-07-21 archetype. Kept for provenance joins only. */
  group: string | null
  /** Current taxonomy region. Colour by this. */
  region?: string | null
  top_weight?: number | null
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

/** Loading vectors for the two displayed axes, for the biplot overlay.
 *  `available` is false over UMAP, which has no linear map back to features. */
export interface PlaneLoadings {
  available: boolean
  reason?: string
  x_axis?: string
  y_axis?: string
  n_features?: number
  vectors?: { feature: string; x: number; y: number; length: number }[]
  note?: string
}

/** One external set's footprint on the map. Post-hoc: these sets never entered
 *  the coordinates. */
export interface EnrichmentPanel {
  group: string
  n: number
  /** Row-major nx by ny, log2(observed rate / base rate). null below `min_n`. */
  cells: (number | null)[]
  max_abs: number
  /** Cells beating the label-permutation null. Zero is a real answer: the set is
   *  spread across the map rather than pooled in one place. */
  n_above_null: number | null
  null_threshold: number | null
  /** Standardised shift of members along each axis. Catches smooth gradients,
   *  which cell-wise testing cannot see. */
  axis_shift: Record<string, number>
  /** Displacement of the set's centroid in the FULL above-noise space, not just
   *  the two axes on screen, plus the overlap that displacement corresponds to.
   *  Every set clears significance against 1,846 genes, so `overlap_pct` is the
   *  number that says whether you could ever see the difference. */
  displacement: {
    z: number | null
    p: number
    strongest_dim: string
    cohens_d: number
    overlap_pct: number
  } | null
  /** Member vs rest densities along each displayed axis, shared bin edges.
   *  This is the "displaced but not separated" claim made visible: the curves
   *  shift but overlap almost completely. */
  distributions?: Record<string, {
    edges: number[]
    members: number[]
    rest: number[]
    member_mean: number
    rest_mean: number
    cohens_d: number
    overlap_pct: number
  }>
  /** Indices into `points` for this set's members. */
  members: number[]
  is_super_enhancer: boolean
  is_positive_control: boolean
  stratification: {
    stratifier: string
    pct_retained: number
    signal_over_random: number
  } | null
}

export interface EnrichmentGrid {
  x_axis: Axis
  y_axis: Axis
  nx: number
  ny: number
  min_n: number
  n_perm: number
  n_genes: number
  cell_totals: number[]
  is_umap: boolean
  is_null: boolean
  /** Every gene's [x, y] in real axis units, so each tile can draw the actual
   *  scatter rather than an abstraction of it. */
  points: [number, number][]
  extent: { x0: number; x1: number; y0: number; y1: number }
  panels: EnrichmentPanel[]
  caveats: Record<string, string>
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

export interface PanelBias {
  headline: string
  one_line: string
  not_about: string
  measured: { what: string; panel: string; genome: string
              not_in_panel: string; note: string }[]
  cause: string
  framing: string
  scope: string
  affects: { claim: string; status: string; why: string }[]
  calibration: string
  consequence: string
  why_more_data: string[]
}

export interface Taxonomy {
  regions: { region: string; reach: string; composition: string; n: number
             mean_top_weight: number | null }[]
  radar: { axis: string; feature: string; values: Record<string, number>
           q25: Record<string, number>; q75: Record<string, number> }[]
  radar_note: string
  mixture: { median_top_weight: number; pct_below_half: number; uniform: number }
  caveat: string
  derivation: { step: number; name: string; detail: string
                result?: string; evidence?: string }[]
  k_evidence: { half: string; n: number; k: number; smallest: number
                stability: number; chosen: boolean }[]
}

/** One gene's blend across the regions. Never a single label: every gene has a
 *  top weight below 0.5, so the nearest region is a coordinate not a category. */
export interface GeneTaxonomy {
  gene_id: string
  nearest_region: string | null
  top_weight: number | null
  mixedness: number | null
  /** Fewer than 3 called peaks: excluded from the taxonomy fit, no region. */
  under_evidenced?: boolean
  under_evidenced_note?: string | null
  /** Weights scaled by how much of the gene the scheme describes, plus one
   *  "outside the scheme" / "no peaks called" component. Sums to 1. */
  mixture: {
    region: string
    weight: number
    /** The conditional weight: share of the DESCRIBED portion. */
    weight_within_described: number | null
    /** Share of the gene's own peak signal in that band+class. */
    own_signal_pct: number | null
    /** False = the gene has no peaks of this kind, so the weight is proximity
     *  in feature space rather than architecture. */
    supported: boolean | null
    note?: string
  }[]
  /** Where the gene's peaks actually are: the raw evidence the weights
   *  summarise. */
  evidence: {
    n_peaks: number
    total_signal: number
    rows: { band: string; n_peaks: number
            promoter: number | null; enhancer: number | null; ctcf: number | null }[]
    classes: string[]
    proximal_pct: number | null
    described_pct: number
    /** Panel reference so a cell can be read as high or typical, not just
     *  large: 40% prox-enhancer is the 85th percentile, but the prox-promoter
     *  median is nearly the same number and far-enhancer's median is 0. */
    panel_median: {
      cells: Record<string, number>
      bands: Record<string, number>
      /** Share of genes with exactly zero in that cell. In the far row this is
       *  41-63%, so a percentile there is inflated by ties. */
      zero_pct: Record<string, number>
      /** This gene's percentile within the cell, strictly-below. */
      percentile: Record<string, number>
      note: string
    }
    coverage_note: string
    weakly_described: boolean
    no_peaks_called?: boolean
  }
  /** How decisive the group assignment is. Panel median margin over the
   *  runner-up is 8.1%, and 59% of genes have a runner-up within 10%, so a
   *  single label shown alone overstates for most genes. */
  assignment: {
    group: string
    runner_up: string | null
    margin_pct: number | null
    close_call: boolean
    panel_median_margin_pct: number
    /** Composition of contacts beyond 50 kb: the slice the clustering keys on,
     *  and therefore what explains an assignment that differs from the largest
     *  (usually proximal) cell in the table. */
    beyond_50kb: {
      total_pct: number
      shares: Record<string, number>
      dominant: string | null
    }
  }
  /** Whether the membership profile says anything at all. A softmax over 5
   *  regions cannot fall below 20% for a gene that resembles nothing in
   *  particular, so 20% is the FLOOR, not a fifth of the architecture. 63% of
   *  genes have their entire profile within 10 pp of that floor. */
  flatness: {
    uniform_pct: number | null
    n_regions: number
    /** Total variation distance from uniform. 0 = no information, 1 = a hard
     *  label. Panel median 0.119. */
    tvd: number | null
    panel_median_tvd: number
    /** Scale references for the one-bar membership summary. */
    panel_p90_tvd: number | null
    panel_max_tvd: number | null
    max_deviation_pp: number | null
    near_uniform: boolean
    panel_pct_near_uniform: number
    panel_pct_above_40: number
    note: string
    /** The opposite failure. The softmax uses exp(-d^2/tau), so a gene far from
     *  EVERY centroid gets a confident-looking profile out of a small relative
     *  gap. corr(nearest distance, top weight) is +0.436, and 29 of the 54 genes
     *  reaching any weight of 40% are in the top distance decile. */
    nearest_distance: number | null
    nearest_distance_pct: number | null
    panel_median_nearest_distance: number | null
    distance_inflated: boolean
    distance_note: string
  }
  /** Distance to each centroid and the softmax that turns it into a weight. */
  derivation: {
    formula: string
    tau: number | null
    tau_note: string
    nearest_distance: number | null
    reading: string
    steps: { region: string; distance: number; excess_over_nearest: number
             exp_term: number; weight: number
             /** Signed distance from the uniform floor, in percentage points.
              *  Negative is evidence AGAINST, which a share cannot express. */
             deviation_pp: number | null }[]
  }
  note: string
}

/** Which genes are worth comparing against, and why. The comparisons that teach
 *  something hold one taxonomy axis fixed and vary the other, which is not what
 *  a reader types into a free-text box. */
export interface CompareSuggestions {
  gene_id: string
  n_components: number
  distance_note: string
  panel_median_nearest: number
  groups: {
    key: string
    label: string
    note: string
    genes: {
      gene_symbol: string
      gene_id: string
      region: string | null
      distance: number | null
      why: string
    }[]
  }[]
}

export interface TaxonomyMap {
  x_axis: Axis
  y_axis: Axis
  regions: string[]
  n: number
  xs: number[]
  ys: number[]
  r: number[]
  w: number[]
  symbols: string[]
  note: string
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

export interface Vocabulary {
  axes: string[]
  cohorts: string[]
  groups: string[]
  features: string[]
  axis_labels: Record<string, string>
  axis_poles: Record<string, { neg: string; pos: string }>
  archetype_labels: Record<string, string>
}

export interface QueryStep {
  filter: Record<string, unknown>
  reads_as: string
  before: number
  after: number
}

export interface QueryResult {
  n_matched: number
  /** True when the answer is the top of a ranking rather than everything past a
   *  threshold. Changes how n_matched should be read: it is the pool that was
   *  ranked, not the size of the answer. */
  ranked?: boolean
  rank_by?: { feature: string; direction: string; weight: number }[] | null
  steps: QueryStep[]
  genes: {
    gene_id: string
    gene_symbol: string
    group: string | null
    /** External sets this gene belongs to. Annotation on an
     *  architecture-selected list; these sets never drive selection. */
    in_sets?: string[]
    /** Composite z across the ranked features. Present only when ranked. */
    score?: number
  }[]
  note: string | null
  query?: Record<string, unknown>
  interpretation?: string | null
  /** Part of the question the store cannot answer, for example TADs or
   *  insulation. Set by the translator, never inferred here. */
  unsupported?: string | null
  /** The model's stated working. A rationale it wrote alongside the query, not
   *  a trace of how it actually computed one, and labelled as such in the UI. */
  reasoning?: string[] | null
  annotation_note?: string
  error?: string
}

export interface FeatureExplain {
  feature: string
  drawable: boolean
  kind: string | null
  measures: string | null
  reads?: string
  why_not_drawable?: string
  quantile?: number | null
  element_class?: string | null
  regions?: {
    start_bp: number
    end_bp: number
    label: string
    mirrored?: boolean
    role?: string
  }[]
  distribution: { min: number; max: number; median: number }
  examples: {
    role: 'high' | 'low'
    gene_id: string
    gene_symbol: string
    group: string | null
    z: number
    percentile: number
    profile: Profile
  }[]
  caveat: string
}

async function post<T>(path: string, body: unknown): Promise<T> {
  const res = await fetch(BASE + path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
  if (!res.ok) {
    let detail = res.statusText
    try {
      detail = (await res.json()).detail ?? detail
    } catch {
      /* not JSON */
    }
    throw new Error(detail)
  }
  return res.json() as Promise<T>
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

  feature: (name: string) =>
    get<FeatureExplain>(`/api/feature/${encodeURIComponent(name)}`),

  vocabulary: () => get<Vocabulary>('/api/vocabulary'),

  ask: (question: string) => post<QueryResult>('/api/ask', { question }),

  query: (q: Record<string, unknown>) => post<QueryResult>('/api/query', q),

  ranked: (axis: string, limit = 25) =>
    get<Ranked>('/api/ranked', { axis, limit }),

  explain: () => get<Explain>('/api/explain'),

  panelBias: () => get<PanelBias>('/api/panel-bias'),
  taxonomy: () => get<Taxonomy>('/api/taxonomy'),
  taxonomyMap: (x: string, y: string) =>
    get<TaxonomyMap>('/api/taxonomy/map', { x, y }),
  geneTaxonomy: (g: string) => get<GeneTaxonomy>(`/api/genes/${g}/taxonomy`),
  psPercentile: (g: string) =>
    get<Record<string, number>>(`/api/ps-percentile/${g}`),
  compareSuggestions: (g: string) =>
    get<CompareSuggestions>(`/api/genes/${g}/compare-suggestions`),
  scree: (space: 'corrected' | 'raw' = 'corrected') =>
    get<Scree>('/api/dimensions/scree', { space }),

  loadings: (pc: number, top = 15, space: 'corrected' | 'raw' = 'corrected') =>
    get<Loadings>(`/api/dimensions/${pc}/loadings`, { top, space }),

  cohorts: () => get<CohortList>('/api/cohorts'),

  cohortCompare: (opts: { group?: string; symbols?: string }) =>
    get<CohortCompare>('/api/cohorts/compare', opts),

  embedding: (x: string, y: string, highlight?: string) =>
    get<Embedding>('/api/embedding', { x, y, highlight }),

  planeLoadings: (x: string, y: string, top = 8) =>
    get<PlaneLoadings>('/api/embedding/loadings', { x, y, top }),

  enrichmentGrid: (x: string, y: string, bins = 12, n_perm = 200) =>
    get<EnrichmentGrid>('/api/enrichment/grid', { x, y, bins, n_perm }),

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
