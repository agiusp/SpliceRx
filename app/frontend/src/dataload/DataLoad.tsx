import { useEffect, useState } from "react";
import { ApiError, dataloadApi, type ScanFile, type ScanResponse } from "./api";
import { api as sjvApi } from "../sjv/api";

type Target = "sjv" | "sjvc" | "sjsurv" | "sjlookup";
const TARGET_LABEL: Record<Target, string> = {
  sjv: "Sashimi plot",
  sjvc: "2D View",
  sjsurv: "SJSurv",
  sjlookup: "SJ Lookup",
};

interface Props {
  sjvSessionId: string | null;
  sjvcSessionId: string | null;
  sjsurvSessionId: string | null;
  sjlookupSessionId: string | null;
  /** `full` is true only when every checked item for this target loaded
   *  without error. Always fires — even on a partial failure — so the target
   *  tab re-hydrates and reflects whatever *did* load into its session. */
  onLoaded: (target: Target, full: boolean) => void;
  /** GENCODE reference was set — both app tabs should re-hydrate. */
  onGencodeApplied: () => void;
}

const DEFAULT_PATH = "~/Work/SJ.Sep2026/Data/PAAD";

// Default the release dropdown to GENCODE v29 when the species offers it.
const pickRelease = (list: string[] = []) => (list.includes("v29") ? "v29" : list[0] ?? "");

function fmtSize(n: number): string {
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(0)} KB`;
  if (n < 1024 * 1024 * 1024) return `${(n / 1024 / 1024).toFixed(0)} MB`;
  return `${(n / 1024 / 1024 / 1024).toFixed(1)} GB`;
}

const STATUS_COLOR: Record<string, string> = {
  ok: "var(--accent)",
  caution: "#b8860b",
  skip: "var(--muted)",
  unsupported: "#c0392b",
};

const CLINICAL_ROLES = new Set(["clinical", "clinical_extra"]);
// sjdat matrix roles that 2D View / SJSurv can use, mapped to the sjdat kind
const SJDAT_ROLE_KIND: Record<string, string> = {
  junction_counts: "junction_counts",
  rrs_scores: "rrs_scores",
  gene_matrix: "gene_matrix",
  pathway_matrix: "pathway_matrix",
};
const loadable = (f: ScanFile) =>
  !!f.target && f.status !== "unsupported" && f.role !== "gtf";

export default function DataLoad({
  sjvSessionId,
  sjvcSessionId,
  sjsurvSessionId,
  sjlookupSessionId,
  onLoaded,
  onGencodeApplied,
}: Props) {
  const [path, setPath] = useState(DEFAULT_PATH);
  const [scan, setScan] = useState<ScanResponse | null>(null);
  const [checked, setChecked] = useState<Record<string, boolean>>({});
  const [scanning, setScanning] = useState(false);
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState<string[]>([]);
  const [warnings, setWarnings] = useState<string[]>([]);
  const [err, setErr] = useState<string | null>(null);
  const [loadedInto, setLoadedInto] = useState<Target[] | null>(null);

  // --- GENCODE reference (one choice for every app) ---
  const [releases, setReleases] = useState<Record<string, string[]>>({});
  const [species, setSpecies] = useState("human");
  const [release, setRelease] = useState("");
  const [gencodeLabel, setGencodeLabel] = useState<string | null>(null);
  const [gencodeBusy, setGencodeBusy] = useState(false);
  const [gencodeErr, setGencodeErr] = useState<string | null>(null);

  useEffect(() => {
    dataloadApi.gencodeReleases().then((r) => {
      setReleases(r);
      setRelease(pickRelease(r[species]));
    });
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  // show what's already set (e.g. after a refresh)
  useEffect(() => {
    if (sjvSessionId) {
      sjvApi
        .sessionState(sjvSessionId)
        .then((st) => st.gencode_label && setGencodeLabel(st.gencode_label))
        .catch(() => {});
    }
  }, [sjvSessionId]);

  async function applyReference(gtfPath?: string) {
    const sids = { sjv_sid: sjvSessionId, sjvc_sid: sjvcSessionId, sjsurv_sid: sjsurvSessionId };
    setGencodeBusy(true);
    setGencodeErr(null);
    try {
      const r = gtfPath
        ? await dataloadApi.applyGtf(sids, gtfPath)
        : await dataloadApi.applyGencode(sids, species, release);
      setGencodeLabel(r.label);
      onGencodeApplied();
    } catch (e) {
      setGencodeErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setGencodeBusy(false);
    }
  }

  async function doScan() {
    setScanning(true);
    setErr(null);
    setScan(null);
    setWarnings([]);
    setProgress([]);
    setLoadedInto(null);
    try {
      const r = await dataloadApi.scan(path);
      setScan(r);
      setChecked(Object.fromEntries(r.files.filter(loadable).map((f) => [f.name, f.status === "ok"])));
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setScanning(false);
    }
  }

  function pickedFor(target: Target): ScanFile[] {
    if (!scan) return [];
    return scan.files.filter((f) => {
      if (!checked[f.name] || !loadable(f)) return false;
      if (target === "sjlookup") return f.role === "junction_metadata"; // no matrix, no clinical table
      if (CLINICAL_ROLES.has(f.role)) return true; // sample/clinical table -> any app
      if (target === "sjv") return f.role === "junction_counts";
      // sjvc / sjsurv: any sjdat matrix, plus the optional junction-metadata
      // table (fast gene lookup on a junction-level one) for both now
      return f.role in SJDAT_ROLE_KIND || f.role === "junction_metadata";
    });
  }

  // Which apps a checked matrix makes relevant, purely from what's ticked —
  // no separate "load into X" choice needed. SJ Lookup is the exception: it
  // never needs a matrix, only the junction-metadata table itself.
  const wantSjv = pickedFor("sjv").some((f) => f.role === "junction_counts");
  const wantSjvc = pickedFor("sjvc").some((f) => f.role in SJDAT_ROLE_KIND);
  const wantSjsurv = pickedFor("sjsurv").some((f) => f.role in SJDAT_ROLE_KIND);
  const wantSjlookup = pickedFor("sjlookup").length > 0;
  const anyWanted = wantSjv || wantSjvc || wantSjsurv || wantSjlookup;

  async function loadSjsurv(sid: string, dir: string, p: (l: string) => void, w: (ws: string[]) => void) {
    const picked = pickedFor("sjsurv");
    const matrices = picked.filter((f) => f.role in SJDAT_ROLE_KIND);
    const jmeta = picked.find((f) => f.role === "junction_metadata");
    const table = picked.find((f) => CLINICAL_ROLES.has(f.role));
    for (const m of matrices) {
      p(`Loading ${m.name}…`);
      const r = await dataloadApi.loadSjsurvSjdat(sid, SJDAT_ROLE_KIND[m.role], `${dir}/${m.name}`);
      p(`✓ ${m.name} — ${r.n_features.toLocaleString()} features × ${r.n_samples} samples${r.sparse ? " (sparse)" : ""}`);
      w(r.warnings);
    }
    if (jmeta) {
      p(`Loading ${jmeta.name}…`);
      const r = await dataloadApi.loadSjsurvJunctionMetadata(sid, `${dir}/${jmeta.name}`);
      p(`✓ ${jmeta.name} — ${r.n_rows.toLocaleString()} junctions, ${r.n_genes.toLocaleString()} genes annotated`);
    }
    if (table) {
      p(`Loading ${table.name}…`);
      const r = await dataloadApi.loadSjsurvMetadata(sid, `${dir}/${table.name}`);
      p(`✓ ${table.name} — ${r.n_matched} sample(s) matched`);
      w(r.warnings);
    } else {
      w(["SJSurv: no sample-metadata table selected — it needs one with Group / SurviverGroup columns"]);
    }
  }

  async function loadSjlookup(sid: string, dir: string, p: (l: string) => void, w: (ws: string[]) => void) {
    const jmeta = pickedFor("sjlookup").find((f) => f.role === "junction_metadata");
    if (!jmeta) return; // shouldn't happen — caller only invokes when `wantSjlookup` is true
    p(`Loading ${jmeta.name}…`);
    const r = await dataloadApi.loadSjlookupJunctionMetadata(sid, `${dir}/${jmeta.name}`);
    p(`✓ ${jmeta.name} — ${r.n_rows.toLocaleString()} junction(s) loaded`);
    w(r.warnings);
  }

  async function loadSjv(sid: string, dir: string, p: (l: string) => void, w: (ws: string[]) => void) {
    const picked = pickedFor("sjv");
    const matrix = picked.find((f) => f.role === "junction_counts");
    const tables = picked.filter((f) => CLINICAL_ROLES.has(f.role));
    if (!matrix) return; // shouldn't happen — caller only invokes when `wantSjv` is true
    p(`Loading ${matrix.name}…`);
    const mres = await dataloadApi.loadJunctions(sid, `${dir}/${matrix.name}`);
    p(`✓ ${matrix.name} — ${mres.n_junctions.toLocaleString()} junctions × ${mres.samples.length} samples`);
    w(mres.warnings);
    for (const t of tables) {
      p(`Loading ${t.name}…`);
      const tres = await dataloadApi.loadSampleMetadata(sid, `${dir}/${t.name}`);
      p(`✓ ${t.name} — ${tres.n_matched} sample(s) matched, ${tres.columns.length} column(s)`);
      w(tres.warnings);
    }
  }

  // 2D View: any of the 3 sjdat matrices (like SJSurv), plus the optional
  // junction-metadata table (fast gene lookup on a junction-level matrix)
  // and the clinical table.
  async function loadSjvc(sid: string, dir: string, p: (l: string) => void, w: (ws: string[]) => void) {
    const picked = pickedFor("sjvc");
    const matrices = picked.filter((f) => f.role in SJDAT_ROLE_KIND);
    const jmeta = picked.find((f) => f.role === "junction_metadata");
    const tables = picked.filter((f) => CLINICAL_ROLES.has(f.role));
    for (const m of matrices) {
      p(`Loading ${m.name}…`);
      const r = await dataloadApi.loadSjvcSjdat(sid, SJDAT_ROLE_KIND[m.role], `${dir}/${m.name}`);
      p(`✓ ${m.name} — ${r.n_features.toLocaleString()} features × ${r.n_samples} samples${r.sparse ? " (sparse)" : ""}`);
      w(r.warnings);
    }
    if (jmeta) {
      p(`Loading ${jmeta.name}…`);
      const r = await dataloadApi.loadSjvcJunctionMetadata(sid, `${dir}/${jmeta.name}`);
      p(`✓ ${jmeta.name} — ${r.n_rows.toLocaleString()} junctions, ${r.n_genes.toLocaleString()} genes annotated`);
    }
    for (const t of tables) {
      p(`Loading ${t.name}…`);
      const tres = await dataloadApi.loadClinical(sid, `${dir}/${t.name}`);
      p(`✓ ${t.name} — ${tres.n_matched} sample(s) matched, ${tres.columns.length} column(s)`);
      w(tres.warnings);
    }
  }

  // One click: load every checked file into every tab it's relevant to. A big
  // matrix like the junction-count file is only parsed once per tab that
  // needs it, but a second parse of the same bytes hits the on-disk cache
  // (a few seconds instead of the ~1 min first pass).
  async function loadCohort() {
    if (!scan) return;
    const targets: { key: Target; sid: string | null; want: boolean }[] = [
      { key: "sjv", sid: sjvSessionId, want: wantSjv },
      { key: "sjvc", sid: sjvcSessionId, want: wantSjvc },
      { key: "sjsurv", sid: sjsurvSessionId, want: wantSjsurv },
      { key: "sjlookup", sid: sjlookupSessionId, want: wantSjlookup },
    ];
    const attempted = targets.filter((t) => t.want);
    if (attempted.length === 0) {
      setErr(
        "select at least one matrix (junction counts / gene matrix / RRS scores / pathway matrix), "
        + "or the junction metadata table for SJ Lookup, to load",
      );
      return;
    }

    setLoading(true);
    setErr(null);
    setProgress([]);
    setWarnings([]);
    setLoadedInto(null);
    const p = (line: string) => setProgress((cur) => [...cur, line]);
    const w = (ws: string[]) => setWarnings((cur) => [...cur, ...ws]);
    const succeeded: Target[] = [];

    for (const t of attempted) {
      if (!t.sid) continue;
      p(`— ${TARGET_LABEL[t.key]} —`);
      let full = false;
      try {
        if (t.key === "sjsurv") await loadSjsurv(t.sid, scan.dir, p, w);
        else if (t.key === "sjv") await loadSjv(t.sid, scan.dir, p, w);
        else if (t.key === "sjlookup") await loadSjlookup(t.sid, scan.dir, p, w);
        else await loadSjvc(t.sid, scan.dir, p, w);
        full = true;
        succeeded.push(t.key);
      } catch (e) {
        setErr(e instanceof ApiError ? e.message : String(e));
      } finally {
        // always re-hydrate the tab — even a partial failure may have landed
        // some files server-side (see the SJSurv metadata-column case).
        onLoaded(t.key, full);
      }
    }

    setLoading(false);
    setLoadedInto(succeeded);
  }

  return (
    <div className="app">
      <h1>Data — load a cohort</h1>
      <p className="sub">
        Point at a directory of prepared TCGA cohort files (from the <code>prepTCGAdata</code>{" "}
        package), check the ones you want, and load once — each file is sent to every tab that
        can use it: the junction matrix to the <b>Sashimi plot</b>; any of the feature matrices
        (junction counts, RRS scores, the gene matrix, or the pathway matrix) to <b>2D View</b>{" "}
        and <b>SJSurv</b>; the junction gene annotation to those two as well (a faster,
        no-GENCODE-needed gene lookup) and to <b>SJ Lookup</b>, which uses it on its own to look
        up individual junctions; a clinical / sample-metadata table to whichever of those are
        loaded.
      </p>

      <div className="panel">
        <h2 style={{ marginTop: 0 }}>GENCODE reference</h2>
        <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
          Chosen once here and used by every tab. First use of a release takes ~1 min to index.
        </p>
        <div className="row" style={{ alignItems: "flex-end" }}>
          <label>
            Species
            <select
              value={species}
              onChange={(e) => {
                setSpecies(e.target.value);
                setRelease(pickRelease(releases[e.target.value]));
              }}
            >
              {Object.keys(releases).map((s) => (
                <option key={s}>{s}</option>
              ))}
            </select>
          </label>
          <label>
            GENCODE release
            <select value={release} onChange={(e) => setRelease(e.target.value)}>
              {(releases[species] ?? []).map((r) => (
                <option key={r}>{r}</option>
              ))}
            </select>
          </label>
          <button
            disabled={
              !release || gencodeBusy || (!sjvSessionId && !sjvcSessionId && !sjsurvSessionId)
            }
            onClick={() => applyReference()}
          >
            {gencodeBusy ? "Indexing…" : "Use this reference"}
          </button>
          {gencodeLabel && <span className="chip">✓ {gencodeLabel}</span>}
        </div>
        {gencodeBusy && <div className="warn">Fetching / indexing annotation for every tab…</div>}
        {gencodeErr && <div className="err">{gencodeErr}</div>}
      </div>

      <div className="panel">
        <div className="row" style={{ alignItems: "flex-end" }}>
          <label style={{ flex: 1 }}>
            Cohort directory
            <input
              type="text"
              value={path}
              onChange={(e) => setPath(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && doScan()}
              style={{ width: "100%" }}
            />
          </label>
          <button disabled={scanning || !path.trim()} onClick={doScan}>
            {scanning ? "Scanning…" : "Scan"}
          </button>
        </div>
      </div>

      {err && <div className="err">{err}</div>}

      {scan && (
        <div className="panel">
          <h2 style={{ marginTop: 0 }}>
            {scan.cohort}{" "}
            <span className="muted" style={{ fontWeight: 400 }}>· {scan.files.length} files</span>
          </h2>

          <div style={{ overflowX: "auto" }}>
            <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 13 }}>
              <thead>
                <tr style={{ textAlign: "left", color: "var(--muted)" }}>
                  <th style={{ padding: "4px 8px", width: 28 }} />
                  <th style={{ padding: "4px 8px" }}>File</th>
                  <th style={{ padding: "4px 8px" }}>Size</th>
                  <th style={{ padding: "4px 8px" }}>Loads into</th>
                  <th style={{ padding: "4px 8px" }}>Notes</th>
                </tr>
              </thead>
              <tbody>
                {scan.files.map((f) => {
                  const can = loadable(f);
                  const dest =
                    f.role === "junction_counts"
                      ? "Sashimi plot, 2D View, SJSurv"
                      : f.role === "gene_matrix" || f.role === "rrs_scores" || f.role === "pathway_matrix"
                        ? "2D View, SJSurv"
                        : f.role === "junction_metadata"
                          ? "2D View, SJSurv, SJ Lookup"
                          : CLINICAL_ROLES.has(f.role)
                            ? "whichever tabs load"
                            : "—";
                  return (
                    <tr
                      key={f.name}
                      style={{ borderTop: "1px solid var(--border)", opacity: can ? 1 : 0.55 }}
                    >
                      <td style={{ padding: "4px 8px" }}>
                        {can && (
                          <input
                            type="checkbox"
                            checked={!!checked[f.name]}
                            onChange={(e) =>
                              setChecked((c) => ({ ...c, [f.name]: e.target.checked }))
                            }
                          />
                        )}
                      </td>
                      <td style={{ padding: "4px 8px", fontFamily: "ui-monospace, monospace" }}>
                        {f.name}
                      </td>
                      <td style={{ padding: "4px 8px", whiteSpace: "nowrap" }}>{fmtSize(f.size)}</td>
                      <td style={{ padding: "4px 8px", whiteSpace: "nowrap" }}>
                        <span style={{ color: STATUS_COLOR[f.status] ?? "var(--text)" }}>
                          {can ? dest : f.role.replace(/_/g, " ")}
                        </span>
                      </td>
                      <td style={{ padding: "4px 8px" }}>
                        {f.role === "gtf" ? (
                          <button
                            className="linklike"
                            disabled={
                              gencodeBusy || (!sjvSessionId && !sjvcSessionId && !sjsurvSessionId)
                            }
                            onClick={() => applyReference(`${scan.dir}/${f.name}`)}
                          >
                            use as GENCODE reference
                          </button>
                        ) : (
                          f.note
                        )}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          <div className="row" style={{ marginTop: 12, alignItems: "center" }}>
            <button disabled={loading || !anyWanted} onClick={loadCohort}>
              {loading ? "Loading…" : "Load cohort"}
            </button>
            {anyWanted && !loading && (
              <span className="muted">
                will load into:{" "}
                {[
                  wantSjv && "Sashimi plot",
                  wantSjvc && "2D View",
                  wantSjsurv && "SJSurv",
                  wantSjlookup && "SJ Lookup",
                ]
                  .filter(Boolean)
                  .join(", ")}
              </span>
            )}
            {(!sjvSessionId || !sjvcSessionId || !sjsurvSessionId) && (
              <span className="muted">connecting…</span>
            )}
          </div>

          {progress.map((line, i) => (
            <div key={i} className="muted" style={{ fontSize: 13, marginTop: 4 }}>
              {line}
            </div>
          ))}
          {[...new Set(warnings)].map((line, i) => (
            <div key={i} className="warn">
              {line}
            </div>
          ))}
          {loadedInto && loadedInto.length > 0 && (
            <div className="muted" style={{ marginTop: 8 }}>
              ✓ Loaded into {loadedInto.map((t) => TARGET_LABEL[t]).join(", ")} — switch tabs above
              to explore.
            </div>
          )}
        </div>
      )}
    </div>
  );
}
