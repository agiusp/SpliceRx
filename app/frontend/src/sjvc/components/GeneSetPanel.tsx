import { useEffect, useState } from "react";
import { api, ApiError, type FeaturesResponse } from "../api";
import Typeahead from "./Typeahead";

export type GeneSetMode = "typed" | "list" | "pathway" | "mad";

interface Props {
  sessionId: string;
  pathwaysEnabled: boolean;
  matrixKind: "junction" | "gene";
  /** Which sjdat matrix is active — switching it (even between two
   *  junction-level ones, e.g. junction counts -> RRS scores) invalidates
   *  whatever was previously resolved/built against the old one. */
  activeSjdat: string | null;
  /** Junction metadata is loaded — typed/list/pathway lookups against a
   *  junction-level matrix work without a GENCODE reference (see
   *  services/junction_metadata.py). */
  fastLookupEnabled: boolean;
  onResolved: (matched: string[], unmatched: string[], source: string) => void;
  onFeatures: (f: FeaturesResponse | null) => void;
  onModeChange: (mode: GeneSetMode) => void;
  onTopNChange: (n: number) => void;
  onCodingOnlyChange: (v: boolean) => void;
  /** Exclude genes in the curated low-mappability paralog-family list (HLA,
   *  immunoglobulin/TCR, MT-*, olfactory receptors, and other named
   *  segmental-duplication clusters) — independent of, and combinable with,
   *  the protein-coding-only radio above it. */
  onExcludeParalogsChange: (v: boolean) => void;
  /** N / X — the coverage prefilter, only meaningful (and only shown here)
   *  for a junction-level sjdat (junction counts / RRS scores). Owned by the
   *  parent so it survives a tab switch. */
  nMin: number;
  xMin: number;
  onNMinChange: (n: number) => void;
  onXMinChange: (n: number) => void;
  /** Rows of the active sjdat with at least one non-zero entry — only known
   *  (non-null) for pathway_matrix. Seeds "top n by MAD" the first time the
   *  matrix becomes active. */
  nNonzeroRows: number | null;
  /** Whether the junction_counts sjdat is also loaded in this session — the
   *  RRS min-supporting-reads filter needs it to look up each RRS
   *  junction's corresponding raw read count. */
  junctionCountsLoaded: boolean;
  minSupportingReads: number | null;
  onMinSupportingReadsChange: (n: number | null) => void;
}

export default function GeneSetPanel({
  sessionId,
  pathwaysEnabled,
  matrixKind,
  activeSjdat,
  fastLookupEnabled,
  onResolved,
  onFeatures,
  onModeChange,
  onTopNChange,
  onCodingOnlyChange,
  onExcludeParalogsChange,
  nMin,
  xMin,
  onNMinChange,
  onXMinChange,
  nNonzeroRows,
  junctionCountsLoaded,
  minSupportingReads,
  onMinSupportingReadsChange,
}: Props) {
  // A gene-level matrix (row names are gene symbols): typed / list / pathway
  // select rows by matching those names directly, no GENCODE resolution —
  // same for a junction-level one once the junction metadata table is loaded.
  const geneLevel = matrixKind === "gene";
  const noGencodeNeeded = geneLevel || fastLookupEnabled;
  // N/X (the coverage prefilter) only make sense for a raw per-junction
  // matrix — meaningless (and hidden) for an already-condensed gene/pathway
  // score matrix.
  const showCoverage = activeSjdat === "junction_counts" || activeSjdat === "rrs_scores";
  // A pathway row is a signature (e.g. "xCell:aDC%HPCA%1.txt"), not a gene —
  // matching typed/uploaded/pathway-searched gene names against its rows
  // doesn't mean anything, so only ranking by MAD is offered.
  const pathwayMatrixOnly = activeSjdat === "pathway_matrix";
  const [tab, setTab] = useState<GeneSetMode>("typed");
  const noun = pathwayMatrixOnly ? "pathways" : geneLevel ? "genes" : "splice junctions";
  const where = geneLevel ? "the matrix" : fastLookupEnabled ? "the junction metadata" : "the reference";
  // RRS scores are bounded [0, 1] and mostly zero, so a row's median is
  // almost always exactly 0 and MAD collapses to (near-)0 for most rows —
  // variance ranks them better; MAD still suits count-like gene/pathway
  // data with a few strong outliers.
  const rankAbbr = activeSjdat === "rrs_scores" ? "Var" : "MAD";
  const rankStat = activeSjdat === "rrs_scores" ? "variance" : "median absolute deviation";
  const [text, setText] = useState("");
  const [prefix, setPrefix] = useState(false);
  const [libraries, setLibraries] = useState<string[]>([]);
  const [library, setLibrary] = useState("");
  const [term, setTerm] = useState<{ term: string; n_genes: number } | null>(null);
  const [topN, setTopN] = useState(50);
  const [codingOnly, setCodingOnly] = useState(false);
  const [excludeParalogs, setExcludeParalogs] = useState(false);
  // RRS scores only — filters junctions to those whose corresponding row in
  // the (separately loaded) junction_counts matrix has a max supporting read
  // count above this value.
  const isRrs = activeSjdat === "rrs_scores";
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [result, setResult] = useState<{ matched: string[]; unmatched: string[]; warnings: string[] } | null>(null);
  const [listFileName, setListFileName] = useState<string | null>(null);

  const [features, setFeatures] = useState<FeaturesResponse | null>(null);
  const [featBusy, setFeatBusy] = useState(false);
  const [featErr, setFeatErr] = useState<string | null>(null);

  // optional refinement, applied on top of a resolved gene set (Type genes /
  // Upload list / Pathway) — narrows it to its top-N most variable members.
  // Meaningless for pathway_matrix (a resolved pathway row isn't a gene to
  // begin with), same restriction as the MAD-only tab there.
  const [madRefine, setMadRefine] = useState(false);
  const [madRefineTopN, setMadRefineTopN] = useState(50);

  const pathwayOk = pathwaysEnabled || noGencodeNeeded;
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

  // Switching the active matrix invalidates whatever was resolved/built
  // against the old one, same as switching mode — even between two
  // junction-level matrices (junction counts -> RRS scores), since they're
  // different data even though `matrixKind` alone doesn't change.
  useEffect(() => {
    setResult(null);
    setFeatures(null);
    setErr(null);
    setFeatErr(null);
    onFeatures(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeSjdat]);

  // pathway_matrix only offers "Top by MAD" — force off any other tab, and
  // seed the "top n" default to every row that isn't entirely zero (rather
  // than the generic default), since that's the natural "rank everything
  // that has signal" starting point for a small, few-hundred-row matrix.
  useEffect(() => {
    if (!pathwayMatrixOnly) return;
    setTab("mad");
    if (nNonzeroRows != null) setTopN(nNonzeroRows);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeSjdat]);

  useEffect(() => {
    onTopNChange(topN);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [topN]);

  useEffect(() => {
    onCodingOnlyChange(codingOnly);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [codingOnly]);

  useEffect(() => {
    onExcludeParalogsChange(excludeParalogs);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [excludeParalogs]);

  // the min-supporting-reads filter is meaningless outside RRS scores —
  // clear it when the active matrix changes away from RRS so a stale value
  // doesn't silently apply if the user switches back to it later on genes
  useEffect(() => {
    if (!isRrs) onMinSupportingReadsChange(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeSjdat]);

  // Genes are resolved first; building the feature matrix (junctions -> genes)
  // always happens right after.
  async function buildFeatures() {
    setFeatBusy(true);
    setFeatErr(null);
    try {
      const f = await api.buildFeatures(sessionId, madRefine ? madRefineTopN : null);
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
        <button
          className={tab === "typed" ? "active" : ""}
          onClick={() => setTab("typed")}
          disabled={pathwayMatrixOnly}
          title={pathwayMatrixOnly ? "a pathway matrix's rows are signatures, not genes — rank them by MAD instead" : ""}
        >
          Type genes
        </button>
        <button
          className={tab === "list" ? "active" : ""}
          onClick={() => setTab("list")}
          disabled={pathwayMatrixOnly}
          title={pathwayMatrixOnly ? "a pathway matrix's rows are signatures, not genes — rank them by MAD instead" : ""}
        >
          Upload list
        </button>
        <button
          className={tab === "pathway" ? "active" : ""}
          onClick={() => setTab("pathway")}
          disabled={pathwayMatrixOnly || !pathwayOk}
          title={
            pathwayMatrixOnly
              ? "a pathway matrix's rows are signatures, not genes — rank them by MAD instead"
              : pathwayOk ? "" : "pathway libraries are human-only"
          }
        >
          Pathway
        </button>
        <button className={tab === "mad" ? "active" : ""} onClick={() => setTab("mad")}>
          Top by {rankAbbr}
        </button>
      </div>

      {pathwayMatrixOnly && tab !== "mad" && (
        <div className="muted" style={{ marginBottom: 8, fontSize: 13 }}>
          This is the pathway-level matrix — each row is already a pathway signature, not a gene,
          so only "Top by MAD" applies.
        </div>
      )}
      {!pathwayMatrixOnly && geneLevel && tab !== "mad" && (
        <div className="muted" style={{ marginBottom: 8, fontSize: 13 }}>
          This matrix is gene-level — the entries below are matched directly against its gene row
          names (no GENCODE needed).
        </div>
      )}
      {!pathwayMatrixOnly && !geneLevel && fastLookupEnabled && tab !== "mad" && (
        <div className="muted" style={{ marginBottom: 8, fontSize: 13 }}>
          The junction metadata table is loaded — the entries below are looked up against its
          precomputed junction/gene annotation (no GENCODE needed).
        </div>
      )}

      {tab !== "mad" && !pathwayMatrixOnly && (
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6, marginBottom: 8 }}>
          <input type="checkbox" checked={madRefine} onChange={(e) => setMadRefine(e.target.checked)} />
          Also rank by MAD and keep only the top
          <input
            type="number"
            min={1}
            disabled={!madRefine}
            value={madRefineTopN}
            style={{ width: 80 }}
            onChange={(e) => setMadRefineTopN(Math.max(1, Math.round(Number(e.target.value) || 1)))}
          />
          most variable {noun}
        </label>
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
          {showCoverage && (
            <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
              A feature is a candidate only when at least <b>N</b> of the samples have a non-zero
              entry that is also <b>≥ X</b>; of those, the top <b>n</b> by {rankStat} are used.
            </p>
          )}
          <div className="row">
            {showCoverage && (
              <>
                <label>
                  N — min samples (a value &lt; 1 is a fraction)
                  <input
                    type="number"
                    min={0}
                    step={0.05}
                    value={nMin}
                    style={{ width: 130 }}
                    onChange={(e) => onNMinChange(Math.max(0, Number(e.target.value) || 0))}
                  />
                </label>
                <label>
                  X — min non-zero magnitude
                  <input
                    type="number"
                    min={0}
                    step={1}
                    value={xMin}
                    style={{ width: 120 }}
                    onChange={(e) => onXMinChange(Math.max(0, Number(e.target.value) || 0))}
                  />
                </label>
              </>
            )}
            <label>
              Number of top {noun} (by {rankAbbr})
              <input
                type="number"
                min={1}
                value={topN}
                style={{ width: 100 }}
                onChange={(e) => setTopN(Math.max(1, Number(e.target.value) || 1))}
              />
            </label>
            <span className="muted" style={{ maxWidth: 420 }}>
              Ranks every {pathwayMatrixOnly ? "pathway" : geneLevel ? "gene" : "splice junction"} in
              the matrix by descending {rankStat}
              {pathwayMatrixOnly && nNonzeroRows != null
                ? ` — defaults to all ${nNonzeroRows} with a non-zero value`
                : ""}
              . Runs when you hit <b>Run</b> in the View step below.
            </span>
          </div>
          {!pathwayMatrixOnly && (
            <div className="row" style={{ marginTop: 8 }}>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input type="radio" checked={!codingOnly} onChange={() => setCodingOnly(false)} />
                All {geneLevel ? "genes" : "junctions"}
              </label>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input type="radio" checked={codingOnly} onChange={() => setCodingOnly(true)} />
                {geneLevel
                  ? "Protein-coding, non-MT genes only"
                  : "Junctions of protein-coding, non-MT genes only"}
              </label>
              <span className="muted" style={{ maxWidth: 360 }}>
                Uses the selected GENCODE reference’s gene biotypes, excluding MT-* genes
                {geneLevel
                  ? ""
                  : " — a live reference is needed even when the junction metadata table is loaded, since that table doesn't carry biotypes"}
                .
              </span>
            </div>
          )}
          {!pathwayMatrixOnly && (
            <div className="row" style={{ marginTop: 8 }}>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input
                  type="checkbox"
                  checked={excludeParalogs}
                  onChange={(e) => setExcludeParalogs(e.target.checked)}
                />
                Exclude low-mappability paralog gene families
              </label>
              <span className="muted" style={{ maxWidth: 400 }}>
                HLA, immunoglobulin/TCR loci, MT-*, olfactory receptors, and other named
                segmental-duplication clusters known to multi-map — independent of, and combinable
                with, the protein-coding option above. Also uses the GENCODE reference.
              </span>
            </div>
          )}
          {isRrs && (
            <div className="row" style={{ marginTop: 8 }}>
              <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
                <input
                  type="checkbox"
                  checked={minSupportingReads != null}
                  disabled={!junctionCountsLoaded}
                  onChange={(e) => onMinSupportingReadsChange(e.target.checked ? 0 : null)}
                />
                Max supporting read count &gt;
                <input
                  type="number"
                  min={0}
                  step={1}
                  disabled={!junctionCountsLoaded || minSupportingReads == null}
                  value={minSupportingReads ?? 0}
                  style={{ width: 90 }}
                  onChange={(e) => onMinSupportingReadsChange(Math.max(0, Number(e.target.value) || 0))}
                />
              </label>
              <span className="muted" style={{ maxWidth: 400 }}>
                {junctionCountsLoaded
                  ? "Keeps an RRS junction only when its corresponding row in the loaded junction counts matrix has a max read count above this, across the matrix's samples."
                  : "Load the junction counts matrix (Data tab) alongside RRS scores to enable this filter."}
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
