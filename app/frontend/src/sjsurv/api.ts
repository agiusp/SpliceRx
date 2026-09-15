export class ApiError extends Error {
  status: number;
  detail: unknown;
  constructor(status: number, detail: unknown) {
    super(typeof detail === "string" ? detail : JSON.stringify(detail));
    this.status = status;
    this.detail = detail;
  }
}

async function j<T>(resP: Response | Promise<Response>): Promise<T> {
  const res = await resP;
  const body = await res.json().catch(() => null);
  if (!res.ok) {
    const detail =
      body && typeof body === "object" && "detail" in (body as any)
        ? (body as any).detail
        : `HTTP ${res.status}${res.statusText ? " " + res.statusText : ""}`;
    throw new ApiError(res.status, detail);
  }
  return body as T;
}

const postJson = (url: string, body: unknown) =>
  fetch(url, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(body) });

export type SjdatKind = "junction_counts" | "rrs_scores" | "gene_matrix" | "pathway_matrix";

export interface SjdatOption {
  kind: SjdatKind;
  label: string;
  description: string;
  loaded: boolean;
  n_features: number;
  n_samples: number;
  sparse: boolean;
  // rows with at least one non-zero entry — only populated for pathway_matrix
  n_nonzero_rows: number | null;
}

export interface AgeBands {
  /** ascending; `null` = open-ended (only meaningful as the last edge) */
  edges: (number | null)[];
  include_lowest: boolean;
  labels: string[];
}

export interface HistologyCount {
  value: string;
  n: number;
}

export interface MetadataLoaded {
  n_matched: number;
  n_unmatched: number;
  n_with_age: number;
  age_min: number | null;
  age_max: number | null;
  age_bands: AgeBands | null;
  min_group_n: number | null;
  histology_counts: HistologyCount[];
  use_histology: boolean;
  histology_map: Record<string, string>;
}

export interface SessionState {
  has_raw_metadata: boolean;
  has_metadata: boolean;
  metadata: MetadataLoaded | null;
  sjdat_options: SjdatOption[];
  active_sjdat: SjdatKind | null;
  gencode_label: string | null;
  pathways_enabled: boolean;
  has_junction_metadata: boolean;
  has_geneset: boolean;
  has_features: boolean;
  selected_group: string | null;
  has_selection: boolean;
  has_model: boolean;
}

// gene-set feature resolution (Type genes / Upload list / Pathway tabs) —
// the exact shapes sjvc/api.ts's equivalents return, since the backend
// re-exports sjvc's own request/response models verbatim (see
// sjsurv/models.py)
export interface GeneSetResponse {
  matched: string[];
  unmatched: string[];
  n_genes: number;
  source: string;
  warnings: string[];
}

export interface FeaturesResponse {
  feature_kind: "junction" | "gene";
  n_features: number;
  n_samples: number;
  n_junctions: number;
  suggest_condense: boolean;
  feature_preview: string[];
  warnings: string[];
}

export interface JunctionMetadataLoaded {
  n_rows: number;
  n_genes: number;
}

export interface GroupCount {
  group: string;
  label: string;
  n_total: number;
  n_good: number;
  n_poor: number;
  n_labelled: number;
}

export interface SelectResponse {
  group: string;
  sjdat_kind: SjdatKind;
  n_group_samples: number;
  n_good: number;
  n_poor: number;
  n_candidates: number;
  n_after_coverage: number;
  n_selected: number;
  feature_preview: string[];
  warnings: string[];
}

export interface CVResponse {
  n_samples: number;
  n_good: number;
  n_poor: number;
  n_features: number;
  n_splits: number;
  auc: number;
  fold_aucs: number[];
  fold_auc_mean: number;
  fold_auc_sd: number;
  accuracy: number;
  sensitivity: number;
  specificity: number;
  confusion: number[][];
  baseline_accuracy: number;
  messages: string[];
}

export interface FeatureWeight {
  feature: string;
  weight: number;
  abs_weight: number;
  direction: string;
  mean_good: number;
  mean_poor: number;
}

export interface ModelResponse {
  sjdat_kind: SjdatKind;
  group: string;
  n_samples: number;
  n_good: number;
  n_poor: number;
  n_features: number;
  auc_resub: number;
  accuracy_resub: number;
  confusion_resub: number[][];
  features: FeatureWeight[];
  saved: boolean;
  messages: string[];
  warnings: string[];
}

export interface SelectParams {
  group: string;
  // omitted (null/undefined) skips the coverage prefilter entirely — what
  // the frontend sends whenever N/X are hidden (any sjdat but junction
  // counts / RRS scores)
  n_min: number | null;
  x_min: number | null;
  top_n: number;
}

export interface StratifyParams {
  edges: (number | null)[];
  include_lowest: boolean;
  min_group_n: number;
  use_histology: boolean;
  histology_map: Record<string, string>;
}

export const api = {
  async createSession(): Promise<string> {
    return (await j<{ session_id: string }>(await fetch("/api/sjsurv/session", { method: "POST" }))).session_id;
  },
  sessionState: (sid: string) => j<SessionState>(fetch(`/api/sjsurv/session/${sid}/state`)),
  groups: (sid: string) =>
    j<{ groups: GroupCount[]; warnings: string[] }>(fetch(`/api/sjsurv/session/${sid}/groups`)),
  suggestAgeBands: (sid: string, n_bands: number) =>
    j<AgeBands>(postJson(`/api/sjsurv/session/${sid}/age-bands/suggest`, { n_bands })),
  stratify: (sid: string, p: StratifyParams) =>
    j<MetadataLoaded>(postJson(`/api/sjsurv/session/${sid}/stratify`, p)),
  activateSjdat: (sid: string, kind: SjdatKind) =>
    j<SessionState>(postJson(`/api/sjsurv/session/${sid}/sjdat`, { kind })),
  select: (sid: string, p: SelectParams) =>
    j<SelectResponse>(postJson(`/api/sjsurv/session/${sid}/select`, {
      group: p.group, n_min: p.n_min, x_min: p.x_min, top_n: p.top_n,
    })),
  // gene-set feature resolution — mirrors sjvc's setGeneset/buildFeatures
  // exactly, just against SJSurv's own session/matrix
  setGeneset: (
    sid: string,
    body: { mode: string; text?: string; library?: string; term?: string; prefix?: boolean },
  ) => j<GeneSetResponse>(postJson(`/api/sjsurv/session/${sid}/geneset`, body)),
  buildFeatures: (sid: string, madTopN: number | null = null) =>
    j<FeaturesResponse>(postJson(`/api/sjsurv/session/${sid}/features`, { condense: false, mad_top_n: madTopN })),
  selectGeneset: (sid: string, group: string) =>
    j<SelectResponse>(postJson(`/api/sjsurv/session/${sid}/select-geneset`, { group })),
  // pathway search is session-independent — call 2D View's own endpoints
  // rather than duplicating them for SJSurv
  pathwayLibraries: (): Promise<string[]> =>
    j<{ libraries: string[] }>(fetch("/api/sjvc/pathways/libraries")).then((r) => r.libraries),
  pathwaySearch: (library: string, q: string) =>
    j<{ terms: { term: string; n_genes: number }[] }>(
      fetch(`/api/sjvc/pathways/search?library=${encodeURIComponent(library)}&q=${encodeURIComponent(q)}`),
    ).then((r) => r.terms),
  crossValidate: (sid: string, n_cv?: number) =>
    j<CVResponse>(postJson(`/api/sjsurv/session/${sid}/cross-validate`, n_cv ? { n_cv } : {})),
  trainModel: (sid: string) =>
    j<ModelResponse>(postJson(`/api/sjsurv/session/${sid}/model`, {})),
};
