import { useEffect, useState } from "react";
import { api, ApiError, type FeaturesResponse } from "../api";
import Typeahead from "./Typeahead";

export type GeneSetMode = "typed" | "list" | "pathway" | "mad";

interface Props {
  sessionId: string;
  pathwaysEnabled: boolean;
  matrixKind: "junction" | "gene";
  onResolved: (matched: string[], unmatched: string[], source: string) => void;
  onFeatures: (f: FeaturesResponse | null) => void;
  onModeChange: (mode: GeneSetMode) => void;
  onTopNChange: (n: number) => void;
  onCodingOnlyChange: (v: boolean) => void;
}

export default function GeneSetPanel({
  sessionId,
  pathwaysEnabled,
  matrixKind,
  onResolved,
  onFeatures,
  onModeChange,
  onTopNChange,
  onCodingOnlyChange,
}: Props) {
  // A gene-level matrix (row names are gene symbols): typed / list / pathway
  // select rows by matching those names directly, no GENCODE resolution.
  const geneLevel = matrixKind === "gene";
  const [tab, setTab] = useState<GeneSetMode>("typed");
  const noun = geneLevel ? "genes" : "splice junctions";
  const where = geneLevel ? "the matrix" : "the reference";
  const [text, setText] = useState("");
  const [prefix, setPrefix] = useState(false);
  const [libraries, setLibraries] = useState<string[]>([]);
  const [library, setLibrary] = useState("");
  const [term, setTerm] = useState<{ term: string; n_genes: number } | null>(null);
  const [topN, setTopN] = useState(50);
  const [codingOnly, setCodingOnly] = useState(false);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [result, setResult] = useState<{ matched: string[]; unmatched: string[]; warnings: string[] } | null>(null);
  const [listFileName, setListFileName] = useState<string | null>(null);

  const [features, setFeatures] = useState<FeaturesResponse | null>(null);
  const [featBusy, setFeatBusy] = useState(false);
  const [featErr, setFeatErr] = useState<string | null>(null);

  const pathwayOk = pathwaysEnabled || geneLevel;
  useEffect(() => {
    if (!pathwayOk) return;
    api.pathwayLibraries().then((l) => {
      setLibraries(l);
      setLibrary(l[0] ?? "");
    });
  }, [pathwayOk]);

  // Switching mode invalidates whatever was previously resolved/built — the
  // "top by MAD" mode in particular doesn't execute here at all; it only
  // reports its parameters upward and runs later, from the View step.
  useEffect(() => {
    onModeChange(tab);
    setResult(null);
    setFeatures(null);
    setErr(null);
    setFeatErr(null);
    onFeatures(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [tab]);

  useEffect(() => {
    onTopNChange(topN);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [topN]);

  useEffect(() => {
    onCodingOnlyChange(codingOnly);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [codingOnly]);

  // Genes are resolved first; building the feature matrix (junctions -> genes)
  // always happens right after.
  async function buildFeatures() {
    setFeatBusy(true);
    setFeatErr(null);
    try {
      const f = await api.buildFeatures(sessionId);
      setFeatures(f);
      onFeatures(f);
    } catch (e) {
      setFeatures(null);
      onFeatures(null);
      setFeatErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setFeatBusy(false);
    }
  }

  async function submit(body: { mode: string; text?: string; library?: string; term?: string; prefix?: boolean }) {
    setBusy(true);
    setErr(null);
    try {
      const r = await api.setGeneset(sessionId, body);
      setResult({ matched: r.matched, unmatched: r.unmatched, warnings: r.warnings ?? [] });
      onResolved(r.matched, r.unmatched, r.source);
      await buildFeatures();
    } catch (e) {
      setResult(null);
      setFeatures(null);
      onFeatures(null);
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  }

  return (
    <div>
      <div className="tabs">
        <button className={tab === "typed" ? "active" : ""} onClick={() => setTab("typed")}>
          Type genes
        </button>
        <button className={tab === "list" ? "active" : ""} onClick={() => setTab("list")}>
          Upload list
        </button>
        <button
          className={tab === "pathway" ? "active" : ""}
          onClick={() => setTab("pathway")}
          disabled={!pathwayOk}
          title={pathwayOk ? "" : "pathway libraries are human-only"}
        >
          Pathway
        </button>
        <button className={tab === "mad" ? "active" : ""} onClick={() => setTab("mad")}>
          Top by MAD
        </button>
      </div>

      {geneLevel && tab !== "mad" && (
        <div className="muted" style={{ marginBottom: 8, fontSize: 13 }}>
          This matrix is gene-level — the entries below are matched directly against its gene row
          names (no GENCODE needed).
        </div>
      )}

      {(tab === "typed" || tab === "list") && (
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6, marginBottom: 8 }}>
          <input type="checkbox" checked={prefix} onChange={(e) => setPrefix(e.target.checked)} />
          Match entries as name prefixes (MUC → MUC1, MUC2, MUC4, …)
        </label>
      )}

      {tab === "typed" && (
        <div style={{ display: "flex", gap: 10, alignItems: "flex-start" }}>
          <textarea
            rows={2}
            style={{ flex: 1, minWidth: 320, font: "inherit", padding: 8, borderRadius: 6, border: "1px solid var(--border)", background: "var(--bg)", color: "var(--text)" }}
            placeholder={prefix ? "MUC, KRT, HLA-" : "TP53, KRAS, MUC1"}
            value={text}
            onChange={(e) => setText(e.target.value)}
          />
          <button disabled={busy || !text.trim()} onClick={() => submit({ mode: "typed", text, prefix })}>
            Resolve
          </button>
        </div>
      )}

      {tab === "list" && (
        <div className="row">
          <label>
            Gene list (one per line / comma / space)
            <span className="file-trigger">
              <span className="file-btn">Choose File</span>
              <input
                type="file"
                accept=".txt,.csv,.tsv,.list"
                onChange={async (e) => {
                  const f = e.target.files?.[0];
                  if (!f) return;
                  setListFileName(f.name);
                  const t = await f.text();
                  submit({ mode: "list", text: t, prefix });
                }}
              />
            </span>
          </label>
          {listFileName && <span className="chip">✓ {listFileName}</span>}
        </div>
      )}

      {tab === "pathway" && (
        <div className="row">
          <label>
            Collection
            <select value={library} onChange={(e) => { setLibrary(e.target.value); setTerm(null); }}>
              {libraries.map((l) => (
                <option key={l}>{l}</option>
              ))}
            </select>
          </label>
          <label>
            Pathway
            <Typeahead
              placeholder="search e.g. apoptosis"
              fetchOptions={async (q) =>
                (await api.pathwaySearch(library, q)).map((t) => ({ label: t.term, sub: `${t.n_genes} genes` }))
              }
              onPick={(label) => setTerm({ term: label, n_genes: 0 })}
            />
          </label>
          {term && (
            <button disabled={busy} onClick={() => submit({ mode: "pathway", library, term: term.term })}>
              Use “{term.term}”
            </button>
          )}
        </div>
      )}

      {tab === "mad" && (
        <>
          <div className="row">
            <label>
              Number of top {noun} (by MAD)
              <input
                type="number"
                min={1}
                value={topN}
                style={{ width: 100 }}
                onChange={(e) => setTopN(Math.max(1, Number(e.target.value) || 1))}
              />
            </label>
            <span className="muted" style={{ maxWidth: 420 }}>
              Ranks every {geneLevel ? "gene" : "splice junction"} in the matrix by descending median
              absolute deviation. Runs when you hit <b>Run</b> in the View step below.
            </span>
          </div>
          {geneLevel && (
            <div className="row" style={{ marginTop: 8 }}>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input type="radio" checked={!codingOnly} onChange={() => setCodingOnly(false)} />
                All genes
              </label>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input type="radio" checked={codingOnly} onChange={() => setCodingOnly(true)} />
                Protein-coding genes only
              </label>
              <span className="muted" style={{ maxWidth: 360 }}>
                “Protein-coding only” uses the selected GENCODE reference’s gene biotypes.
              </span>
            </div>
          )}
        </>
      )}

      {tab !== "mad" && busy && <div className="warn">Resolving genes…</div>}
      {tab !== "mad" && err && <div className="err">{err}</div>}
      {tab !== "mad" && result && (
        <div style={{ marginTop: 8, fontSize: 13 }}>
          <b>{result.matched.length}</b> gene(s) matched
          {result.matched.length > 0 && (
            <span className="muted"> · {result.matched.slice(0, 20).join(", ")}{result.matched.length > 20 ? "…" : ""}</span>
          )}
          {result.unmatched.length > 0 && (
            <div className="muted">
              {result.unmatched.length} not in {where}: {result.unmatched.slice(0, 12).join(", ")}
              {result.unmatched.length > 12 ? "…" : ""}
            </div>
          )}
          {result.warnings.map((w, i) => (
            <div className="warn" key={i}>
              {w}
            </div>
          ))}
        </div>
      )}

      {tab !== "mad" && featBusy && <div className="warn">Building feature matrix…</div>}
      {tab !== "mad" && featErr && <div className="err">{featErr}</div>}
      {tab !== "mad" && features && (
        <div style={{ marginTop: 8, fontSize: 13 }}>
          {geneLevel ? (
            <b>{features.n_features} gene feature(s)</b>
          ) : (
            <>
              {features.n_junctions} junction(s) in the selected genes →{" "}
              <b>{features.n_features} junction feature(s)</b>
            </>
          )}{" "}
          × {features.n_samples} samples
          <div className="muted" style={{ marginTop: 2 }}>
            preprocessing: NA→0, drop constant features, log1p, then z-score
          </div>
          {features.warnings.map((w, i) => (
            <div className="warn" key={i}>
              {w}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
