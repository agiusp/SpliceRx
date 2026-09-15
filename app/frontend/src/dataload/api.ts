export class ApiError extends Error {
  status: number;
  constructor(status: number, detail: unknown) {
    super(typeof detail === "string" ? detail : JSON.stringify(detail));
    this.status = status;
  }
}

async function j<T>(resP: Response | Promise<Response>): Promise<T> {
  const res = await resP;
  const body = await res.json().catch(() => null);
  if (!res.ok) {
    const detail =
      body && typeof body === "object" && "detail" in (body as any)
        ? (body as any).detail
        : `HTTP ${res.status} — the server returned no details`;
    throw new ApiError(res.status, detail);
  }
  return body as T;
}

const post = (url: string, body: unknown) =>
  fetch(url, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(body) });

export type FileStatus = "ok" | "caution" | "skip" | "unsupported";

export interface ScanFile {
  name: string;
  size: number;
  role: string;
  target: string | null;
  status: FileStatus;
  note: string;
}

export interface ScanResponse {
  cohort: string;
  dir: string;
  files: ScanFile[];
}

export interface MatrixLoaded {
  samples: string[];
  n_junctions: number;
  feature_kind?: "junction" | "gene";
  sparse?: boolean;
  warnings: string[];
}

export interface ColumnsLoaded {
  columns: { name: string }[];
  n_matched: number;
  warnings: string[];
}

export interface GencodeApplied {
  label: string;
  pathways_enabled: boolean;
  applied_to: string[];
  warnings: string[];
}

interface Sids {
  sjv_sid: string | null;
  sjvc_sid: string | null;
  sjsurv_sid: string | null;
}

export const dataloadApi = {
  scan: (path: string) => j<ScanResponse>(post("/api/dataload/scan", { path })),

  // GENCODE reference — one choice for both apps
  gencodeReleases: async () =>
    (await j<{ releases: Record<string, string[]> }>(await fetch("/api/dataload/gencode/releases")))
      .releases,
  applyGencode: (sids: Sids, species: string, release: string) =>
    j<GencodeApplied>(post("/api/dataload/gencode", { ...sids, species, release })),
  applyGtf: (sids: Sids, path: string) =>
    j<GencodeApplied>(post("/api/dataload/gtf", { ...sids, path })),

  // SJVC (2D View) — any of the 3 sjdat matrices, like SJSurv, plus the
  // optional junction-metadata table for fast gene lookup on a junction-level one
  loadGeneMatrix: (sid: string, path: string) =>
    j<MatrixLoaded>(post(`/api/dataload/sjvc/${sid}/gene-matrix`, { path })),
  loadClinical: (sid: string, path: string) =>
    j<ColumnsLoaded>(post(`/api/dataload/sjvc/${sid}/clinical`, { path })),
  loadSjvcSjdat: (sid: string, kind: string, path: string) =>
    j<SjdatLoaded>(post(`/api/dataload/sjvc/${sid}/sjdat/${kind}`, { path })),
  loadSjvcJunctionMetadata: (sid: string, path: string) =>
    j<JunctionMetadataLoaded>(post(`/api/dataload/sjvc/${sid}/junction-metadata`, { path })),
  // SJV (sashimi plot)
  loadJunctions: (sid: string, path: string) =>
    j<MatrixLoaded>(post(`/api/dataload/sjv/${sid}/junctions`, { path })),
  loadSampleMetadata: (sid: string, path: string) =>
    j<ColumnsLoaded>(post(`/api/dataload/sjv/${sid}/sample-metadata`, { path })),

  // SJSurv (survivor groups)
  loadSjsurvSjdat: (sid: string, kind: string, path: string) =>
    j<SjdatLoaded>(post(`/api/dataload/sjsurv/${sid}/sjdat/${kind}`, { path })),
  loadSjsurvMetadata: (sid: string, path: string) =>
    j<MetadataLoaded>(post(`/api/dataload/sjsurv/${sid}/metadata`, { path })),
  loadSjsurvJunctionMetadata: (sid: string, path: string) =>
    j<JunctionMetadataLoaded>(post(`/api/dataload/sjsurv/${sid}/junction-metadata`, { path })),

  // SJ Lookup (per-junction lookup) — only ever needs the junction metadata
  // table, never a matrix
  loadSjlookupJunctionMetadata: (sid: string, path: string) =>
    j<SjlookupJunctionMetadataLoaded>(post(`/api/dataload/sjlookup/${sid}/junction-metadata`, { path })),

  // push SJSurv's Group/SurviverGroup into the Sashimi-plot / 2D View sessions
  pushGroups: (sjsurv_sid: string, sjv_sid: string | null, sjvc_sid: string | null) =>
    j<PushGroupsResult>(post("/api/dataload/push-groups", { sjsurv_sid, sjv_sid, sjvc_sid })),
};

export interface SjdatLoaded {
  kind: string;
  n_features: number;
  n_samples: number;
  sparse: boolean;
  warnings: string[];
}

export interface MetadataLoaded {
  n_matched: number;
  n_unmatched: number;
  warnings: string[];
}

export interface JunctionMetadataLoaded {
  n_rows: number;
  n_genes: number;
}

export interface SjlookupJunctionMetadataLoaded {
  n_rows: number;
  n_duplicate_rownames: number;
  warnings: string[];
}

export interface PushGroupsResult {
  applied_to: string[];
  n_sjsurv_samples: number;
  n_group: number;
  n_labelled: number;
  matched: Record<string, number>;
  warnings: string[];
}
