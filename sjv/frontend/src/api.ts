export interface RdsUploaded {
  samples: string[];
  n_junctions: number;
  warnings: string[];
}

export interface GeneRecord {
  name: string;
  gene_id: string;
  chrom: string;
  start: number;
  end: number;
  strand: string;
}

export interface ExonModel {
  start: number;
  end: number;
  kind: "CDS" | "UTR";
}

export interface TranscriptModel {
  transcript_id: string;
  strand: string;
  exons: ExonModel[];
  gene_name: string;
}

export interface ArcModel {
  id: string;
  start: number;
  end: number;
  strand: string;
  count: number;
  height: number;
  category: string;
}

export interface LegendEntry {
  category: string;
  label: string;
  color: string;
}

export interface PlotResponse {
  genes: GeneRecord[];
  x_domain: [number, number];
  layout: { baseline_fraction: number; track_row_px: number };
  transcripts: TranscriptModel[];
  arcs: ArcModel[];
  legend: LegendEntry[];
  series_label: string;
  count_kind: "count" | "median";
  warnings: string[];
}

export interface StratColumn {
  name: string;
  n_values: number;
}

export interface SampleMetadataUploaded {
  columns: StratColumn[];
  n_matched: number;
  warnings: string[];
}

export interface StratValue {
  value: string;
  n_samples: number;
}

export class ApiError extends Error {
  status: number;
  detail: unknown;
  constructor(status: number, detail: unknown) {
    super(typeof detail === "string" ? detail : JSON.stringify(detail));
    this.status = status;
    this.detail = detail;
  }
}

async function j<T>(res: Response): Promise<T> {
  const body = await res.json().catch(() => ({}));
  if (!res.ok) throw new ApiError(res.status, (body as any).detail ?? body);
  return body as T;
}

export const api = {
  async createSession(): Promise<string> {
    const d = await j<{ session_id: string }>(await fetch("/api/session", { method: "POST" }));
    return d.session_id;
  },

  async uploadRds(sid: string, file: File): Promise<RdsUploaded> {
    const fd = new FormData();
    fd.append("file", file);
    return j(await fetch(`/api/session/${sid}/rds`, { method: "POST", body: fd }));
  },

  async releases(): Promise<Record<string, string[]>> {
    const d = await j<{ releases: Record<string, string[]> }>(await fetch("/api/gencode/releases"));
    return d.releases;
  },

  async selectGencode(sid: string, species: string, release: string): Promise<{ label: string; warnings: string[] }> {
    return j(
      await fetch(`/api/session/${sid}/gencode`, {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ species, release }),
      }),
    );
  },

  async uploadGtf(sid: string, file: File): Promise<{ label: string; warnings: string[] }> {
    const fd = new FormData();
    fd.append("file", file);
    return j(await fetch(`/api/session/${sid}/gtf`, { method: "POST", body: fd }));
  },

  async suggestGenes(sid: string, q: string): Promise<string[]> {
    if (!q.trim()) return [];
    const d = await j<{ names: string[] }>(
      await fetch(`/api/session/${sid}/genes/suggest?q=${encodeURIComponent(q)}`),
    );
    return d.names;
  },

  async uploadSampleMetadata(sid: string, file: File): Promise<SampleMetadataUploaded> {
    const fd = new FormData();
    fd.append("file", file);
    return j(await fetch(`/api/session/${sid}/sample-metadata`, { method: "POST", body: fd }));
  },

  async stratValues(sid: string, column: string): Promise<StratValue[]> {
    const d = await j<{ values: StratValue[] }>(
      await fetch(`/api/session/${sid}/sample-metadata/values?column=${encodeURIComponent(column)}`),
    );
    return d.values;
  },

  async plot(
    sid: string,
    series: { sample?: string; strat_column?: string; strat_value?: string },
    genes: string[],
    min_reads = 0,
  ): Promise<PlotResponse> {
    return j(
      await fetch(`/api/session/${sid}/plot`, {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ ...series, genes, min_reads }),
      }),
    );
  },
};
