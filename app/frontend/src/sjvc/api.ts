export class ApiError extends Error {
  status: number;
  detail: unknown;
  constructor(status: number, detail: unknown) {
    // `.message` (used everywhere errors are shown) is always a readable
    // string: pass the detail through as-is if it already is one, otherwise
    // fall back to its JSON form rather than the default `[object Object]`.
    super(typeof detail === "string" ? detail : JSON.stringify(detail));
    this.status = status;
    this.detail = detail;
  }
}

async function j<T>(res: Response): Promise<T> {
  const body = await res.json().catch(() => null);
  if (!res.ok) {
    const detail =
      body && typeof body === "object" && "detail" in (body as any)
        ? (body as any).detail
        : `HTTP ${res.status}${res.statusText ? " " + res.statusText : ""} — the server returned no details`;
    throw new ApiError(res.status, detail);
  }
  return body as T;
}

const upload = (url: string, file: File) => {
  const fd = new FormData();
  fd.append("file", file);
  return fetch(url, { method: "POST", body: fd });
};
const postJson = (url: string, body: unknown) =>
  fetch(url, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(body) });

export interface ClinicalColumn {
  name: string;
  type: "numeric" | "categorical";
  n_unique: number;
  n_missing: number;
  eligible: boolean;
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

export interface LegendColor {
  feature: string;
  kind: "categorical" | "sequential";
  items?: { value: string; color: string }[];
  min?: number;
  max?: number;
  stops?: string[];
}
export interface LegendShape {
  feature: string;
  items: { value: string; shape: string }[];
}
export interface BivariateLegend {
  features: [string, string];
  corners: { lo_lo: string; hi_lo: string; lo_hi: string; hi_hi: string };
  ranges: [[number, number], [number, number]];
}

export interface ProjectionPoint {
  sample: string;
  x: number;
  y: number;
  color: string;
  shape: string;
  // true when this point is drawn grey — missing the selected clinical
  // feature(s) — so it can be hidden without recomputing the projection.
  missing: boolean;
  clinical: Record<string, string | number | null>;
}

export interface ProjectionResponse {
  method: "pca" | "umap";
  axis_labels: string[];
  explained_variance: number[] | null;
  n_components: number;
  pc_x: number;
  pc_y: number;
  points: ProjectionPoint[];
  encoding: { color: { kind: string; feature?: string; features?: string[] }; shape: { kind: string; feature?: string } };
  legend: { color?: LegendColor; shape?: LegendShape; bivariate?: BivariateLegend };
  warnings: string[];
}

export interface HeatmapAnnotation {
  feature: string;
  type: "numeric" | "categorical";
  values: (string | number | null)[];
  colors: string[];
  legend?: { value: string; color: string }[];
  min?: number;
  max?: number;
  stops?: string[];
}

export interface HeatmapResponse {
  feature_kind: "junction" | "gene";
  row_labels: string[];
  // "<gene name>:<novel-splicing-event type>" per row, from the cohort's
  // junction metadata table — only present for a junction-level heatmap when
  // that table was loaded; `null` per-row where the table has no gene
  // annotation for that junction (fall back to `row_labels`/`row_ids` there).
  row_labels_gene_type: (string | null)[] | null;
  row_ids: string[];
  samples: string[];
  values: number[][];
  row_zscore: boolean;
  row_dendro: number[][][] | null;
  col_dendro: number[][][] | null;
  annotations: HeatmapAnnotation[];
  warnings: string[];
}

export interface SjdatOption {
  kind: "junction_counts" | "rrs_scores" | "gene_matrix" | "pathway_matrix";
  label: string;
  description: string;
  loaded: boolean;
  n_features: number;
  n_samples: number;
  sparse: boolean;
  // rows with at least one non-zero entry — only populated for pathway_matrix
  n_nonzero_rows: number | null;
}

export interface SessionState {
  has_junctions: boolean;
  feature_kind: "junction" | "gene" | null;
  samples: string[];
  n_junctions: number;
  clinical: { columns: ClinicalColumn[]; n_matched: number; warnings: string[] } | null;
  gencode_label: string | null;
  pathways_enabled: boolean;
  sjdat_options: SjdatOption[];
  active_sjdat: SjdatOption["kind"] | null;
  has_junction_metadata: boolean;
}

export const api = {
  async createSession(): Promise<string> {
    return (await j<{ session_id: string }>(await fetch("/api/sjvc/session", { method: "POST" }))).session_id;
  },
  async sessionState(sid: string) {
    return j<SessionState>(await fetch(`/api/sjvc/session/${sid}/state`));
  },
  async uploadJunctions(sid: string, file: File) {
    return j<{ samples: string[]; n_junctions: number; feature_kind: "junction" | "gene"; warnings: string[] }>(
      await upload(`/api/sjvc/session/${sid}/junctions`, file),
    );
  },
  async uploadClinical(sid: string, file: File) {
    return j<{ columns: ClinicalColumn[]; n_matched: number; warnings: string[] }>(
      await upload(`/api/sjvc/session/${sid}/clinical`, file),
    );
  },
  async activateSjdat(sid: string, kind: SjdatOption["kind"]) {
    return j<SessionState>(await postJson(`/api/sjvc/session/${sid}/sjdat`, { kind }));
  },
  async releases() {
    return (await j<{ releases: Record<string, string[]> }>(await fetch("/api/sjvc/gencode/releases"))).releases;
  },
  async selectGencode(sid: string, species: string, release: string) {
    return j<{ label: string; pathways_enabled: boolean; warnings: string[] }>(
      await postJson(`/api/sjvc/session/${sid}/gencode`, { species, release }),
    );
  },
  async uploadGtf(sid: string, file: File) {
    return j<{ label: string; pathways_enabled: boolean; warnings: string[] }>(
      await upload(`/api/sjvc/session/${sid}/gtf`, file),
    );
  },
  async suggestGenes(sid: string, q: string): Promise<string[]> {
    if (!q.trim()) return [];
    return (await j<{ names: string[] }>(await fetch(`/api/sjvc/session/${sid}/genes/suggest?q=${encodeURIComponent(q)}`))).names;
  },
  async pathwayLibraries(): Promise<string[]> {
    return (await j<{ libraries: string[] }>(await fetch("/api/sjvc/pathways/libraries"))).libraries;
  },
  async pathwaySearch(library: string, q: string) {
    return (
      await j<{ terms: { term: string; n_genes: number }[] }>(
        await fetch(`/api/sjvc/pathways/search?library=${encodeURIComponent(library)}&q=${encodeURIComponent(q)}`),
      )
    ).terms;
  },
  async setGeneset(
    sid: string,
    body: { mode: string; text?: string; library?: string; term?: string; prefix?: boolean },
  ) {
    return j<{ matched: string[]; unmatched: string[]; n_genes: number; source: string; warnings: string[] }>(
      await postJson(`/api/sjvc/session/${sid}/geneset`, body),
    );
  },
  async buildFeatures(sid: string, madTopN: number | null = null) {
    return j<FeaturesResponse>(
      await postJson(`/api/sjvc/session/${sid}/features`, { condense: false, mad_top_n: madTopN }),
    );
  },
  async buildFeaturesMad(
    sid: string, topN: number, proteinCodingOnly = false,
    nMin: number | null = null, xMin = 0,
  ) {
    return j<FeaturesResponse>(
      await postJson(`/api/sjvc/session/${sid}/features/mad`, {
        condense: false,
        top_n: topN,
        protein_coding_only: proteinCodingOnly,
        n_min: nMin,
        x_min: xMin,
      }),
    );
  },
  async featuresCsv(sid: string): Promise<Blob> {
    const res = await fetch(`/api/sjvc/session/${sid}/features.csv`);
    if (!res.ok) {
      const body = await res.json().catch(() => null);
      const detail =
        body && typeof body === "object" && "detail" in body
          ? (body as { detail: unknown }).detail
          : `HTTP ${res.status}`;
      throw new ApiError(res.status, detail);
    }
    return res.blob();
  },
  async projection(sid: string, body: Record<string, unknown>) {
    return j<ProjectionResponse>(await postJson(`/api/sjvc/session/${sid}/projection`, body));
  },
  async heatmap(sid: string, body: Record<string, unknown>) {
    return j<HeatmapResponse>(await postJson(`/api/sjvc/session/${sid}/heatmap`, body));
  },
};
