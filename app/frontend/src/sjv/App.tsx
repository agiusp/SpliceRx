import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { api, ApiError, type PlotResponse } from "./api";
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
  const [loadedInfo, setLoadedInfo] = useState<{ rds?: string; ann?: string }>({});
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

  async function draw(genes: string[], cutoff = minReads) {
    if (!sessionId || !seriesSel || genes.length === 0) return;
    setDrawing(true);
    setError(null);
    setNearMatches([]);
    setLastGenes(genes);
    try {
      const p = await api.plot(sessionId, seriesSel, genes, cutoff);
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
