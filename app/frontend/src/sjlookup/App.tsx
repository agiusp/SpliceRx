import { useEffect, useState } from "react";
import { api, ApiError, type LookupResponse, type SessionState } from "./api";

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

export default function App({ sessionId, reloadNonce }: Props) {
  const [state, setState] = useState<SessionState | null>(null);
  const [text, setText] = useState("");
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [result, setResult] = useState<LookupResponse | null>(null);

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

  return (
    <div className="app">
      <h1>SJ Lookup — per-junction metadata</h1>
      <p className="sub">
        Paste one or more splice junctions (<code>chr:start-end:strand</code>, e.g.{" "}
        <code>{EXAMPLE}</code>) and see whatever the cohort's junction metadata table
        (<code>prepTCGAdata::annotate_sj()</code> output) records for each — gene overlap,
        annotation status, and splice-site motifs.
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
                          <td>{r.annotated === null ? <span className="muted">—</span> : r.annotated ? "yes" : "no"}</td>
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
          </div>
        )}
      </Panel>
    </div>
  );
}
