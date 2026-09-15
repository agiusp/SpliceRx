import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { api, ApiError, type PlotResponse } from "./api";
import GeneInput from "./components/GeneInput";
import Legend from "./components/Legend";
import SeriesPicker, { type Series } from "./components/SeriesPicker";
import SashimiPlot from "./components/SashimiPlot";
import UploadPanel from "./components/UploadPanel";

export default function App() {
  const [sessionId, setSessionId] = useState<string | null>(null);
  const [samples, setSamples] = useState<string[]>([]);
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

  useEffect(() => {
    api.createSession().then(setSessionId).catch((e) => setError(String(e)));
  }, []);

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
        Upload a junction matrix and a GENCODE reference, pick a sample (or a group) and a gene,
        and draw a SCANVIS-style sashimi plot.
      </p>

      {sessionId && (
        <UploadPanel
          sessionId={sessionId}
          onRds={(s, w) => {
            setSamples(s);
            setSeries((prev) => ({ ...prev, sample: s[0] ?? "" }));
            setWarnings((prev) => [...prev, ...w]);
          }}
          onAnnotation={(_label, w) => {
            setAnnotationReady(true);
            setWarnings((prev) => [...prev, ...w]);
          }}
        />
      )}

      <div className="panel">
        <h2>2 · Series &amp; gene</h2>
        {sessionId && (
          <SeriesPicker
            sessionId={sessionId}
            samples={samples}
            value={series}
            onChange={setSeries}
            onWarnings={(w) => setWarnings((prev) => [...prev, ...w])}
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
            Select a GENCODE release or upload a GTF above.
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
