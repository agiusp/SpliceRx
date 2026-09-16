import { useEffect, useMemo, useState } from "react";
import { api, ApiError, CATEGORY_LABELS, type LookupResponse, type SampleValue, type SessionState } from "./api";

interface Props {
  /** Session id, owned by the Shell so the Data tab can load files into it. */
  sessionId: string | null;
  /** Bump to re-hydrate from the server (after a Data-tab load). */
  reloadNonce: number;
}

function Panel({ n, title, done, disabled, children }: any) {
  return (
    <div className={"panel" + (disabled ? " disabled" : "")}>
      <h2>
        <span className="step-num">{n}</span>
        {title} {done && <span style={{ color: "var(--accent)" }}>✓</span>}
      </h2>
      {children}
    </div>
  );
}

const EXAMPLE = "chr3:148857542-149059801:+";

type SampleSortKey = "sample" | "count" | "rrs_score";
interface SampleSort {
  key: SampleSortKey;
  dir: "asc" | "desc";
}

/** Server sends `sample_values` sorted descending by RRS score; this
 *  re-sorts client-side by whichever column the user last clicked (nulls —
 *  a matrix that isn't loaded, or a sample missing from it — always sort to
 *  the bottom, in either direction) and filters by a sample-id search. */
function sortSamples(rows: SampleValue[], sort: SampleSort): SampleValue[] {
  const sign = sort.dir === "asc" ? 1 : -1;
  return [...rows].sort((a, b) => {
    const av = a[sort.key];
    const bv = b[sort.key];
    if (av == null && bv == null) return 0;
    if (av == null) return 1;
    if (bv == null) return -1;
    if (av < bv) return -sign;
    if (av > bv) return sign;
    return 0;
  });
}

export default function App({ sessionId, reloadNonce }: Props) {
  const [state, setState] = useState<SessionState | null>(null);
  const [text, setText] = useState("");
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [result, setResult] = useState<LookupResponse | null>(null);
  const [sampleQuery, setSampleQuery] = useState("");
  const [sampleSort, setSampleSort] = useState<SampleSort>({ key: "rrs_score", dir: "desc" });

  useEffect(() => {
    if (!sessionId) return;
    api
      .sessionState(sessionId)
      .then(setState)
      .catch((e) => setErr(e instanceof ApiError ? e.message : String(e)));
  }, [sessionId, reloadNonce]);

  async function doLookup() {
    if (!sessionId || !text.trim()) return;
    setBusy(true);
    setErr(null);
    setSampleQuery("");
    setSampleSort({ key: "rrs_score", dir: "desc" });
    try {
      setResult(await api.lookup(sessionId, text));
    } catch (e) {
      setResult(null);
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  }

  const dataReady = !!state?.has_index;

  const shownSampleValues = useMemo(() => {
    const rows = result?.sample_values ?? [];
    const q = sampleQuery.trim().toLowerCase();
    const filtered = q ? rows.filter((sv) => sv.sample.toLowerCase().includes(q)) : rows;
    return sortSamples(filtered, sampleSort);
  }, [result, sampleQuery, sampleSort]);

  function sortableHeader(key: SampleSortKey, label: string) {
    const active = sampleSort.key === key;
    return (
      <th
        style={{ cursor: "pointer", userSelect: "none" }}
        onClick={() =>
          setSampleSort((cur) =>
            cur.key === key ? { key, dir: cur.dir === "desc" ? "asc" : "desc" } : { key, dir: "desc" },
          )
        }
      >
        {label}
        {active ? (sampleSort.dir === "desc" ? " ▼" : " ▲") : ""}
      </th>
    );
  }

  return (
    <div className="app">
      <h1>SJ Lookup — per-junction metadata</h1>
      <p className="sub">
        Paste one or more splice junctions (<code>chr:start-end:strand</code>, e.g.{" "}
        <code>{EXAMPLE}</code>) and see whatever the cohort's junction metadata table
        (<code>prepTCGAdata::annotate_sj()</code> output) records for each — gene overlap,
        annotation status, and splice-site motifs. Query a single junction and, if junction
        counts and/or RRS scores are loaded too, you also get every sample's read count and RRS
        score for it, sorted by RRS score.
      </p>

      {err && <div className="err">{err}</div>}

      {/* 1 · Data --------------------------------------------------------- */}
      <Panel n={1} title="Data" done={dataReady}>
        <div className="row">
          <label style={{ display: "block" }}>Junction metadata table</label>
          {dataReady ? (
            <span className="chip">✓ {state!.n_rows.toLocaleString()} junction(s) loaded</span>
          ) : (
            <span className="muted">
              Load the junction metadata file (<code>TCGA_&lt;cohort&gt;_junction_metadata.rds</code>)
              on the Data tab.
            </span>
          )}
        </div>
        <div className="row" style={{ marginTop: 8 }}>
          <label style={{ display: "block" }}>Per-sample matrices (optional)</label>
          {state && state.sjdat_loaded.length > 0 ? (
            <span className="chip">
              ✓ {state.sjdat_loaded.map((k) => (k === "junction_counts" ? "junction counts" : "RRS scores")).join(", ")}
            </span>
          ) : (
            <span className="muted">
              Load junction counts and/or RRS scores on the Data tab to also get the per-sample
              table for a single-junction lookup.
            </span>
          )}
        </div>
      </Panel>

      {/* 2 · Look up junctions --------------------------------------------- */}
      <Panel n={2} title="Look up junctions" done={!!result} disabled={!dataReady}>
        <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
          One or more, separated by commas, spaces, or newlines.
        </p>
        <div style={{ display: "flex", gap: 10, alignItems: "flex-start" }}>
          <textarea
            rows={3}
            style={{
              flex: 1, minWidth: 320, font: "inherit", padding: 8, borderRadius: 6,
              border: "1px solid var(--border)", background: "var(--bg)", color: "var(--text)",
            }}
            placeholder={EXAMPLE}
            value={text}
            onChange={(e) => setText(e.target.value)}
          />
          <button disabled={!dataReady || busy || !text.trim()} onClick={doLookup}>
            {busy ? "Looking up…" : "Look up"}
          </button>
        </div>

        {result && (
          <div style={{ marginTop: 12 }}>
            <div className="muted" style={{ fontSize: 13, marginBottom: 8 }}>
              <b>{result.n_found}</b> found, <b>{result.n_not_found}</b> not in the table,{" "}
              <b>{result.n_invalid}</b> not a valid <code>chr:start-end:strand</code>.
            </div>
            <div className="scroll-x">
              <table className="grid">
                <thead>
                  <tr>
                    <th>Junction</th>
                    <th>Gene</th>
                    <th>Width</th>
                    <th>Annotated</th>
                    <th>Left motif</th>
                    <th>Right motif</th>
                    <th>Left end</th>
                    <th>Right end</th>
                  </tr>
                </thead>
                <tbody>
                  {result.results.map((r, i) => (
                    <tr key={i} className={!r.found ? (r.error ? "invalid" : "not-found") : ""}>
                      <td className="mono">{r.junction}</td>
                      {r.error ? (
                        <td colSpan={7}>{r.error}</td>
                      ) : !r.found ? (
                        <td colSpan={7} className="muted">not in the loaded junction metadata</td>
                      ) : (
                        <>
                          <td>
                            {r.gene_name ?? <span className="muted">—</span>}
                            {r.gene_id && (
                              <div className="muted mono" style={{ fontSize: 11 }}>
                                {r.gene_id}
                              </div>
                            )}
                          </td>
                          <td className="mono">{r.width ?? <span className="muted">—</span>}</td>
                          <td>
                            {r.annotated === null ? (
                              <span className="muted">—</span>
                            ) : r.annotated ? (
                              "yes"
                            ) : (
                              <>
                                no
                                {r.category && (
                                  <span className="muted"> — {CATEGORY_LABELS[r.category] ?? r.category}</span>
                                )}
                              </>
                            )}
                          </td>
                          <td className="mono">{r.left_motif ?? <span className="muted">—</span>}</td>
                          <td className="mono">{r.right_motif ?? <span className="muted">—</span>}</td>
                          <td>{r.left_annotated ?? <span className="muted">—</span>}</td>
                          <td>{r.right_annotated ?? <span className="muted">—</span>}</td>
                        </>
                      )}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            {result.sample_values !== null && (
              <div style={{ marginTop: 16 }}>
                <div className="muted" style={{ fontSize: 13, marginBottom: 8 }}>
                  Per sample — click a column to sort by it (click again to reverse).
                </div>
                {result.sample_values_warnings.map((w, i) => (
                  <div key={i} className="warn">{w}</div>
                ))}
                {result.sample_values.length > 0 && (
                  <>
                    <div className="row" style={{ marginBottom: 8 }}>
                      <input
                        type="text"
                        placeholder="find a sample…"
                        value={sampleQuery}
                        onChange={(e) => setSampleQuery(e.target.value)}
                        style={{ minWidth: 220 }}
                      />
                      {sampleQuery && (
                        <span className="muted" style={{ fontSize: 13 }}>
                          {shownSampleValues.length} of {result.sample_values.length} sample(s)
                        </span>
                      )}
                    </div>
                    <div className="scroll-y">
                      <table className="grid">
                        <thead>
                          <tr>
                            {sortableHeader("sample", "Sample")}
                            {sortableHeader("count", "Read count")}
                            {sortableHeader("rrs_score", "RRS score")}
                          </tr>
                        </thead>
                        <tbody>
                          {shownSampleValues.map((sv) => (
                            <tr key={sv.sample}>
                              <td className="mono">{sv.sample}</td>
                              <td className="mono">
                                {sv.count === null ? <span className="muted">—</span> : Math.round(sv.count).toLocaleString()}
                              </td>
                              <td className="mono">
                                {sv.rrs_score === null ? <span className="muted">—</span> : sv.rrs_score.toFixed(4)}
                              </td>
                            </tr>
                          ))}
                          {shownSampleValues.length === 0 && (
                            <tr>
                              <td colSpan={3} className="muted">no sample matches "{sampleQuery}"</td>
                            </tr>
                          )}
                        </tbody>
                      </table>
                    </div>
                  </>
                )}
              </div>
            )}
          </div>
        )}
      </Panel>
    </div>
  );
}
