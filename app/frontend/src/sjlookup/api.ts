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

export interface SessionState {
  has_index: boolean;
  n_rows: number;
}

export interface JunctionInfo {
  junction: string;
  found: boolean;
  error: string | null;
  gene_id: string | null;
  gene_name: string | null;
  width: number | null;
  annotated: boolean | null;
  left_motif: string | null;
  right_motif: string | null;
  left_annotated: string | null;
  right_annotated: string | null;
}

export interface LookupResponse {
  results: JunctionInfo[];
  n_found: number;
  n_not_found: number;
  n_invalid: number;
}

export const api = {
  async createSession(): Promise<string> {
    return (await j<{ session_id: string }>(await fetch("/api/sjlookup/session", { method: "POST" }))).session_id;
  },
  sessionState: (sid: string) => j<SessionState>(fetch(`/api/sjlookup/session/${sid}/state`)),
  lookup: (sid: string, junctions: string) =>
    j<LookupResponse>(postJson(`/api/sjlookup/session/${sid}/lookup`, { junctions })),
};
