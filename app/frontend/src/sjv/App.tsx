import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { api, ApiError, type AnnotationSource, type PlotResponse } from "./api";
import GeneInput from "./components/GeneInput";
import Legend from "./components/Legend";
import SeriesPicker, { type Series } from "./components/SeriesPicker";
import SashimiPlot from "./components/SashimiPlot";
import UploadPanel from "./components/UploadPanel";

interface Props {
  /** Session id, owned by the Shell so the Data tab can load files into it. */
  sessionId: string | null;
  /** Bump to re-hydrate step-1 state from the server (after a Data-tab load). */
  reloadNonce: number;
}

export default function App({ sessionId, reloadNonce }: Props) {
  const [samples, setSamples] = useState<string[]>([]);
  const [metaColumns, setMetaColumns] = useState<string[]>([]);
  const [loadedInfo, setLoadedInfo] = useState<{ rds?: string; ann?: string; jmeta?: string }>({});
  const [series, setSeries] = useState<Series>({
    mode: "sample",
    sample: "",
    stratColumn: "",
    stratValue: "",
  });
  const [annotationReady, setAnnotationReady] = useState(false);
  const [warnings, setWarnings] = useState<string[]>([]);
  const [plot, setPlot] = useState<PlotResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [nearMatches, setNearMatches] = useState<string[]>([]);
  const [drawing, setDrawing] = useState(false);
  const [minReads, setMinReads] = useState(0);
  const [lastGenes, setLastGenes] = useState<string[] | null>(null);
  // How each arc's category is decided: live against the loaded GENCODE
  // transcript models ("gencode", the long-standing behaviour), or from the
  // cohort's own junction-metadata table's recount3/STAR-aligner `annotated`
  // column, when loaded ("metadata") — falling back to "gencode" per-arc for
  // any junction the table has no row for.
  const [annotationSource, setAnnotationSource] = useState<AnnotationSource>("gencode");
  const [hasJunctionMetadata, setHasJunctionMetadata] = useState(false);

  const wrapRef = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(1000);

  // Hydrate step 1 from the server: on mount and whenever the Data tab reports it
  // loaded files into this session. Also means the tab survives a page refresh.
  useEffect(() => {
    if (!sessionId) return;
    api
      .sessionState(sessionId)
      .then((st) => {
        // GENCODE reference is set on the Data tab, independently of the matrix
        if (st.gencode_label) setAnnotationReady(true);
        setLoadedInfo((li) => ({ ...li, ann: st.gencode_label ?? li.ann }));
        setHasJunctionMetadata(st.has_junction_metadata);
        if (!st.has_junction_metadata) setAnnotationSource("gencode");
        setLoadedInfo((li) => ({
          ...li,
          jmeta: st.has_junction_metadata
            ? `loaded${st.junction_metadata_has_detail ? "" : " (no left/right site detail)"}`
            : undefined,
        }));
        if (!st.has_rds) return;
        setSamples(st.samples);
        setSeries((prev) => ({ ...prev, sample: prev.sample || st.samples[0] || "" }));
        setMetaColumns(st.sample_metadata_columns);
        setLoadedInfo((li) => ({
          ...li,
          rds: `${st.sparse ? "sparse " : ""}junction matrix · ${st.n_junctions.toLocaleString()} junctions · ${st.samples.length} samples`,
        }));
      })
      .catch((e) => setError(String(e)));
  }, [sessionId, reloadNonce]); // eslint-disable-line react-hooks/exhaustive-deps

  useLayoutEffect(() => {
    if (!wrapRef.current) return;
    const ro = new ResizeObserver((entries) => setWidth(entries[0].contentRect.width));
    ro.observe(wrapRef.current);
    return () => ro.disconnect();
  }, []);

  const seriesSel =
    series.mode === "group"
      ? series.stratColumn && series.stratValue
        ? { strat_column: series.stratColumn, strat_value: series.stratValue }
        : null
      : series.sample
        ? { sample: series.sample }
        : null;

  async function draw(genes: string[], cutoff = minReads, source = annotationSource) {
    if (!sessionId || !seriesSel || genes.length === 0) return;
    setDrawing(true);
    setError(null);
    setNearMatches([]);
    setLastGenes(genes);
    try {
      const p = await api.plot(sessionId, seriesSel, genes, cutoff, source);
      setPlot(p);
      setWarnings(p.warnings);
    } catch (e) {
      setPlot(null);
      if (e instanceof ApiError && e.detail && typeof e.detail === "object") {
        const d = e.detail as { message?: string; near_matches?: string[] };
        setError(d.message ?? "lookup failed");
        setNearMatches(d.near_matches ?? []);
      } else {
        setError(e instanceof ApiError ? String(e.detail) : String(e));
      }
    } finally {
      setDrawing(false);
    }
  }

  const ready = Boolean(sessionId && seriesSel && annotationReady);

  // re-draw the current gene(s) when the min-reads cutoff changes
  useEffect(() => {
    if (!plot || !lastGenes) return;
    const t = setTimeout(() => draw(lastGenes, minReads), 300);
    return () => clearTimeout(t);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [minReads]);

  // re-draw immediately (no debounce — it's a discrete choice, not typing)
  // when the annotation source changes
  useEffect(() => {
    if (!plot || !lastGenes) return;
    draw(lastGenes, minReads, annotationSource);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [annotationSource]);

  return (
    <div className="app" ref={wrapRef}>
      <h1>SJV — Splice Junction Visualizer</h1>
      <p className="sub">
        Load a cohort on the Data tab, pick a sample (or a group) and a gene, and draw a
        SCANVIS-style sashimi plot.
      </p>

      <UploadPanel loaded={loadedInfo} />

      <div className="panel">
        <h2>2 · Series &amp; gene</h2>
        {sessionId && (
          <SeriesPicker
            sessionId={sessionId}
            samples={samples}
            initialColumns={metaColumns}
            value={series}
            onChange={setSeries}
          />
        )}
        <div className="row" style={{ marginTop: 12 }}>
          {sessionId && (
            <GeneInput
              sessionId={sessionId}
              disabled={!ready || drawing}
              onSubmit={(gs) => draw(gs)}
            />
          )}
          <label>
            Min reads to display
            <input
              type="number"
              min={0}
              step={1}
              value={minReads}
              style={{ width: 90 }}
              onChange={(e) => setMinReads(Math.max(0, Number(e.target.value) || 0))}
            />
          </label>
        </div>
        <div className="row" style={{ marginTop: 8 }}>
          <fieldset style={{ border: "none", padding: 0, margin: 0 }}>
            <legend className="muted" style={{ fontSize: 13, padding: 0 }}>
              Junction annotation
            </legend>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6, display: "flex" }}>
              <input
                type="radio"
                name="annotation-source"
                checked={annotationSource === "gencode"}
                onChange={() => setAnnotationSource("gencode")}
              />
              Computed live from GENCODE
            </label>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6, display: "flex" }}>
              <input
                type="radio"
                name="annotation-source"
                checked={annotationSource === "metadata"}
                disabled={!hasJunctionMetadata}
                onChange={() => setAnnotationSource("metadata")}
              />
              "Annotated" column from the junction metadata table (STAR/recount3)
            </label>
          </fieldset>
          <span className="muted" style={{ maxWidth: 460 }}>
            {hasJunctionMetadata
              ? "The loaded table's own annotated/left_annotated/right_annotated columns — from " +
                "the aligner that originally called these junctions — take precedence over the " +
                "live per-gene GENCODE classification below; a junction missing from the table " +
                "falls back to it."
              : "Load the cohort's junction metadata table (TCGA_<cohort>_junction_metadata.rds) " +
                "on the Data tab to enable this."}
          </span>
        </div>
        {!annotationReady && (
          <div className="muted" style={{ marginTop: 8 }}>
            Choose a GENCODE reference on the Data tab.
          </div>
        )}
        {error && <div className="err">{error}</div>}
        {nearMatches.length > 0 && (
          <div className="warn">
            Did you mean:{" "}
            {nearMatches.map((n, i) => (
              <span key={n}>
                {i > 0 && ", "}
                <button className="ghost" style={{ padding: "1px 6px" }} onClick={() => draw([n])}>
                  {n}
                </button>
              </span>
            ))}
          </div>
        )}
      </div>

      {[...new Set(warnings)].map((w, i) => (
        <div className="warn" key={i}>
          {w}
        </div>
      ))}

      {plot && (
        <div className="panel">
          <h2>
            3 · {plot.genes.map((g) => g.name).join(" + ")} — {plot.series_label}
            <span className="muted" style={{ fontWeight: 400 }}>
              {plot.count_kind === "median" && " · per-junction medians, excluding NA / 0"}
              {minReads > 0 && ` · junctions with ≥ ${minReads} reads`}
            </span>
          </h2>
          <Legend entries={plot.legend} />
          <SashimiPlot plot={plot} width={width - 34} />
        </div>
      )}
    </div>
  );
}
