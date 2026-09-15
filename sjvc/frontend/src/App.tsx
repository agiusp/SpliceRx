import { useEffect, useRef, useState } from "react";
import { api, ApiError, type ClinicalColumn, type FeaturesResponse, type HeatmapResponse, type ProjectionResponse } from "./api";
import ClinicalPanel from "./components/ClinicalPanel";
import GeneSetPanel, { type GeneSetMode } from "./components/GeneSetPanel";
import Heatmap from "./components/Heatmap";
import Projection from "./components/Projection";
import UploadPanel from "./components/UploadPanel";

type View = "pca" | "umap" | "heatmap";

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

export default function App() {
  const [sid, setSid] = useState<string | null>(null);
  const [samples, setSamples] = useState<string[]>([]);
  const [clinicalCols, setClinicalCols] = useState<ClinicalColumn[]>([]);
  const [gencodeLabel, setGencodeLabel] = useState<string | null>(null);
  const [pathwaysEnabled, setPathwaysEnabled] = useState(false);
  const [matrixKind, setMatrixKind] = useState<"junction" | "gene">("junction");
  const [features, setFeatures] = useState<FeaturesResponse | null>(null);
  const [warnings, setWarnings] = useState<string[]>([]);

  const [geneSetMode, setGeneSetMode] = useState<GeneSetMode>("typed");
  const [madTopN, setMadTopN] = useState(50);
  const [madCodingOnly, setMadCodingOnly] = useState(false);
  const [madBusy, setMadBusy] = useState(false);
  const [madErr, setMadErr] = useState<string | null>(null);

  const [view, setView] = useState<View>("pca");
  const [selected, setSelected] = useState<string[]>([]);
  const [overrides, setOverrides] = useState<Record<string, "numeric" | "categorical">>({});
  const [order, setOrder] = useState<"cluster" | "group">("cluster");
  const [groupBy, setGroupBy] = useState("");
  const [rowZ, setRowZ] = useState(true);
  const [nNeighbors, setNNeighbors] = useState(15);
  const [minDist, setMinDist] = useState(0.1);

  const [proj, setProj] = useState<ProjectionResponse | null>(null);
  const [heat, setHeat] = useState<HeatmapResponse | null>(null);
  const [plotErr, setPlotErr] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    api.createSession().then(setSid).catch((e) => setPlotErr(String(e)));
  }, []);

  // an already-condensed gene-level matrix needs no GENCODE — MAD just ranks
  // its rows; a raw junction matrix still needs the reference to resolve genes
  const dataReady = Boolean(sid && samples.length && (matrixKind === "gene" || gencodeLabel));
  // For "top by MAD" the View step's Run button is the thing that actually
  // builds features — so the panel only needs a valid top-N to be usable, not
  // a finished build. Every other mode still gates on a finished build.
  const viewPanelReady = Boolean(dataReady && (geneSetMode === "mad" ? madTopN > 0 : true));
  const ready = Boolean(dataReady && features);
  const catCols = clinicalCols.filter((c) => c.eligible && (overrides[c.name] ?? c.type) === "categorical");

  // changing the MAD parameters after a run invalidates it — require an
  // explicit Run again rather than silently rebuilding
  useEffect(() => {
    if (geneSetMode === "mad") setFeatures(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [madTopN, madCodingOnly]);

  // clamp selection for projection
  useEffect(() => {
    if (view !== "heatmap" && selected.length > 2) setSelected(selected.slice(-2));
  }, [view]); // eslint-disable-line

  const timer = useRef<number>();
  useEffect(() => {
    if (!ready) return;
    window.clearTimeout(timer.current);
    timer.current = window.setTimeout(runPlot, 250);
    return () => window.clearTimeout(timer.current);
    // eslint-disable-next-line
  }, [ready, view, selected, overrides, order, groupBy, rowZ, nNeighbors, minDist, features]);

  async function runMad() {
    if (!sid) return;
    setMadBusy(true);
    setMadErr(null);
    try {
      const f = await api.buildFeaturesMad(sid, madTopN, matrixKind === "gene" && madCodingOnly);
      setFeatures(f);
    } catch (e) {
      setFeatures(null);
      setMadErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setMadBusy(false);
    }
  }

  async function downloadMadData() {
    if (!sid) return;
    try {
      const blob = await api.featuresCsv(sid);
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "heatmap_data.csv";
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      setMadErr(e instanceof ApiError ? e.message : String(e));
    }
  }

  async function runPlot() {
    if (!sid || !ready) return;
    setLoading(true);
    setPlotErr(null);
    try {
      if (view === "heatmap") {
        const h = await api.heatmap(sid, {
          clinical: selected,
          order,
          group_by: order === "group" ? groupBy || catCols[0]?.name : undefined,
          row_zscore: rowZ,
          overrides,
        });
        setHeat(h);
        setProj(null);
      } else {
        const p = await api.projection(sid, {
          method: view,
          clinical: selected.slice(0, 2),
          overrides,
          n_neighbors: nNeighbors,
          min_dist: minDist,
        });
        setProj(p);
        setHeat(null);
      }
    } catch (e) {
      setPlotErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setLoading(false);
    }
  }

  return (
    <div className="app">
      <h1>SJVC — Splice Junction Visualizer with Clinical</h1>
      <p className="sub">
        Upload a junction matrix and a clinical table, choose a gene set, then explore the cohort as a
        heatmap or a 2D projection with clinical variables layered on.
      </p>

      {sid && (
        <Panel n={1} title="Data" done={ready}>
          <UploadPanel
            sessionId={sid}
            onJunctions={(s, w, kind) => {
              setSamples(s);
              setMatrixKind(kind);
              if (kind !== "gene") setMadCodingOnly(false);
              setWarnings((p) => [...p, ...w]);
            }}
            onClinical={(cols, n, w) => {
              setClinicalCols(cols);
              setWarnings((p) => [...p, ...w, ...(n < samples.length ? [`${samples.length - n} matrix sample(s) have no clinical row`] : [])]);
            }}
            onGencode={(label, pe, w) => {
              setGencodeLabel(label);
              setPathwaysEnabled(pe);
              setWarnings((p) => [...p, ...w]);
            }}
          />
        </Panel>
      )}

      {sid && (
        <Panel
          n={2}
          title="Gene set"
          done={geneSetMode === "mad" ? madTopN > 0 : !!features}
          disabled={!samples.length || (matrixKind !== "gene" && !gencodeLabel)}
        >
          <GeneSetPanel
            sessionId={sid}
            pathwaysEnabled={pathwaysEnabled}
            matrixKind={matrixKind}
            onResolved={() => setFeatures(null)}
            onFeatures={setFeatures}
            onModeChange={setGeneSetMode}
            onTopNChange={setMadTopN}
            onCodingOnlyChange={setMadCodingOnly}
          />
        </Panel>
      )}

      <Panel n={3} title="View" disabled={!viewPanelReady}>
        <div className="tabs">
          {(["pca", "umap", "heatmap"] as View[]).map((v) => (
            <button key={v} className={view === v ? "active" : ""} onClick={() => setView(v)}>
              {v === "pca" ? "PCA" : v === "umap" ? "UMAP" : "Heatmap"}
            </button>
          ))}
        </div>

        {geneSetMode === "mad" && (() => {
          const needsRef = matrixKind === "gene" && madCodingOnly && !gencodeLabel;
          return (
            <div className="row" style={{ marginTop: 10 }}>
              <button disabled={madBusy || madTopN < 1 || needsRef} onClick={runMad}>
                {features ? "Re-run" : "Run"} — top {madTopN}{" "}
                {matrixKind === "gene" ? (madCodingOnly ? "protein-coding genes" : "genes") : "junctions"} by MAD
              </button>
              {madBusy && <span className="muted">computing…</span>}
              {needsRef && (
                <span className="muted">select a GENCODE release in step 1 to filter to protein-coding genes</span>
              )}
            </div>
          );
        })()}
        {geneSetMode === "mad" && madErr && <div className="err">{madErr}</div>}
        {geneSetMode === "mad" && features && (
          <div className="muted" style={{ marginTop: 6, fontSize: 13 }}>
            <b>
              {features.n_features} {features.feature_kind === "gene" ? "gene" : "junction"} feature(s)
            </b>{" "}
            × {features.n_samples} samples
            {" · "}
            <button className="linklike" onClick={downloadMadData}>
              download data (CSV)
            </button>
            {features.feature_preview.length > 0 && (
              <div>
                preview: {features.feature_preview.slice(0, 10).join(", ")}
                {features.feature_preview.length > 10 ? "…" : ""}
              </div>
            )}
            {features.warnings.map((w, i) => (
              <div className="warn" key={i}>
                {w}
              </div>
            ))}
          </div>
        )}

        <ClinicalPanel
          columns={clinicalCols}
          selected={selected}
          overrides={overrides}
          maxFeatures={view === "heatmap" ? Infinity : 2}
          onSelected={setSelected}
          onOverride={(name, t) => setOverrides((o) => ({ ...o, [name]: t }))}
        />

        {view === "heatmap" && (
          <div className="row" style={{ marginTop: 10 }}>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
              <input type="checkbox" checked={rowZ} onChange={(e) => setRowZ(e.target.checked)} />
              row z-score
            </label>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
              <input
                type="radio"
                checked={order === "cluster"}
                onChange={() => setOrder("cluster")}
              />
              cluster columns
            </label>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
              <input
                type="radio"
                checked={order === "group"}
                onChange={() => setOrder("group")}
                disabled={catCols.length === 0}
              />
              group columns by
            </label>
            {order === "group" && (
              <select value={groupBy || catCols[0]?.name || ""} onChange={(e) => setGroupBy(e.target.value)}>
                {catCols.map((c) => (
                  <option key={c.name}>{c.name}</option>
                ))}
              </select>
            )}
          </div>
        )}

        {view === "umap" && (
          <div className="row" style={{ marginTop: 10 }}>
            <label>
              n_neighbors
              <input
                type="number"
                min={2}
                max={Math.max(2, samples.length - 1)}
                value={nNeighbors}
                style={{ width: 80 }}
                onChange={(e) => setNNeighbors(Math.max(2, Number(e.target.value) || 15))}
              />
            </label>
            <label>
              min_dist
              <input
                type="number"
                min={0}
                max={1}
                step={0.05}
                value={minDist}
                style={{ width: 80 }}
                onChange={(e) => setMinDist(Math.max(0, Number(e.target.value) || 0.1))}
              />
            </label>
          </div>
        )}
      </Panel>

      {[...new Set(warnings)].map((w, i) => (
        <div className="warn" key={i}>
          {w}
        </div>
      ))}
      {plotErr && <div className="err">{plotErr}</div>}

      {ready && (
        <div className="panel">
          <h2>
            <span className="step-num">4</span>
            {view === "heatmap" ? "Heatmap" : view.toUpperCase()} {loading && <span className="muted">· updating…</span>}
          </h2>
          {proj && view !== "heatmap" && <Projection data={proj} />}
          {heat && view === "heatmap" && <Heatmap data={heat} />}
        </div>
      )}
    </div>
  );
}
