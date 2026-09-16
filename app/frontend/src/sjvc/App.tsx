import { useEffect, useRef, useState } from "react";
import {
  api,
  ApiError,
  type ClinicalColumn,
  type FeaturesResponse,
  type HeatmapResponse,
  type ProjectionResponse,
  type SessionState,
  type SjdatOption,
} from "./api";
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

interface Props {
  /** Session id, owned by the Shell so the Data tab can load files into it. */
  sessionId: string | null;
  /** Bump to re-hydrate step-1 state from the server (after a Data-tab load). */
  reloadNonce: number;
}

export default function App({ sessionId, reloadNonce }: Props) {
  const sid = sessionId;
  const [samples, setSamples] = useState<string[]>([]);
  const [clinicalCols, setClinicalCols] = useState<ClinicalColumn[]>([]);
  const [gencodeLabel, setGencodeLabel] = useState<string | null>(null);
  const [pathwaysEnabled, setPathwaysEnabled] = useState(false);
  // Any of the 4 sjdat matrices can be active (junction counts, RRS scores,
  // the gene-level matrix, or the pathway-level matrix) — same choice as SJSurv.
  const [sjdatOptions, setSjdatOptions] = useState<SjdatOption[]>([]);
  const [activeSjdat, setActiveSjdat] = useState<SjdatOption["kind"] | null>(null);
  const [matrixKind, setMatrixKind] = useState<"junction" | "gene">("gene");
  const [hasJunctionMetadata, setHasJunctionMetadata] = useState(false);
  const [features, setFeatures] = useState<FeaturesResponse | null>(null);
  const [warnings, setWarnings] = useState<string[]>([]);
  const [pickErr, setPickErr] = useState<string | null>(null);

  const [geneSetMode, setGeneSetMode] = useState<GeneSetMode>("typed");
  const [madTopN, setMadTopN] = useState(50);
  const [madCodingOnly, setMadCodingOnly] = useState(false);
  const [madNMin, setMadNMin] = useState(0.2);
  const [madXMin, setMadXMin] = useState(1);
  const [madBusy, setMadBusy] = useState(false);
  const [madErr, setMadErr] = useState<string | null>(null);

  const [view, setView] = useState<View>("pca");
  const [selected, setSelected] = useState<string[]>([]);
  const [overrides, setOverrides] = useState<Record<string, "numeric" | "categorical">>({});
  const [order, setOrder] = useState<"cluster" | "group">("cluster");
  const [groupBy, setGroupBy] = useState("");
  const [rowZ, setRowZ] = useState(true);
  // Heatmap-only: drop samples missing a selected clinical feature instead of
  // showing them with a grey annotation cell. Unlike the projection's
  // equivalent (a pure display filter — see Projection.tsx), this actually
  // changes which columns are clustered, so it has to go back to the server.
  const [dropMissingClinical, setDropMissingClinical] = useState(false);
  const [nNeighbors, setNNeighbors] = useState(15);
  const [minDist, setMinDist] = useState(0.1);

  const [proj, setProj] = useState<ProjectionResponse | null>(null);
  const [heat, setHeat] = useState<HeatmapResponse | null>(null);
  const [plotErr, setPlotErr] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // Hydrate step-1 state from the server: on mount, and whenever the Data tab
  // reports it has loaded files into this session (reloadNonce bump). Also means
  // the tab survives a browser refresh.
  function applyState(st: SessionState) {
    // GENCODE reference is set on the Data tab, independently of the matrix
    if (st.gencode_label) {
      setGencodeLabel(st.gencode_label);
      setPathwaysEnabled(st.pathways_enabled);
    }
    setSjdatOptions(st.sjdat_options);
    setActiveSjdat(st.active_sjdat);
    setHasJunctionMetadata(st.has_junction_metadata);
    if (!st.has_junctions) return;
    setSamples(st.samples);
    setMatrixKind(st.feature_kind === "gene" ? "gene" : "junction");
    if (st.feature_kind !== "gene") setMadCodingOnly(false);
    if (st.clinical) {
      setClinicalCols(st.clinical.columns);
      if (st.clinical.warnings.length) setWarnings((p) => [...p, ...st.clinical!.warnings]);
    }
  }

  useEffect(() => {
    if (!sid) return;
    api.sessionState(sid).then(applyState).catch((e) => setPlotErr(String(e)));
  }, [sid, reloadNonce]); // eslint-disable-line react-hooks/exhaustive-deps

  async function pickSjdat(kind: SjdatOption["kind"]) {
    if (!sid) return;
    setPickErr(null);
    setFeatures(null);
    try {
      applyState(await api.activateSjdat(sid, kind));
    } catch (e) {
      setPickErr(e instanceof ApiError ? e.message : String(e));
    }
  }

  // an already-condensed gene-level matrix needs no GENCODE — MAD just ranks
  // its rows; a raw junction matrix still needs the reference to resolve genes
  const dataReady = Boolean(
    sid && samples.length && (matrixKind === "gene" || gencodeLabel || hasJunctionMetadata),
  );
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
  }, [madTopN, madCodingOnly, madNMin, madXMin]);

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
  }, [ready, view, selected, overrides, order, groupBy, rowZ, dropMissingClinical, nNeighbors, minDist, features]);

  async function runMad() {
    if (!sid) return;
    setMadBusy(true);
    setMadErr(null);
    try {
      const showCoverage = activeSjdat === "junction_counts" || activeSjdat === "rrs_scores";
      const f = await api.buildFeaturesMad(
        sid, madTopN, matrixKind === "gene" && madCodingOnly,
        showCoverage ? madNMin : null, showCoverage ? madXMin : 0,
      );
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
          drop_missing_clinical: dropMissingClinical,
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
      <h1>2D View</h1>
      <p className="sub">
        Explore a cohort's samples in a 2D projection (PCA / UMAP) or a heatmap, built from a
        selection of splice-junction data. Pick any one of the cohort's matrices — junction
        counts, RRS scores, the gene-level matrix, or the pathway-level matrix — loaded on the
        Data tab along with its clinical table, choose a gene set, then explore with clinical
        variables layered on.
      </p>

      {pickErr && <div className="err">{pickErr}</div>}
      {sid && (
        <Panel n={1} title="Data" done={ready}>
          <UploadPanel
            options={sjdatOptions}
            activeSjdat={activeSjdat}
            onPick={pickSjdat}
            hasJunctionMetadata={hasJunctionMetadata}
            loaded={{
              clinical: clinicalCols.length ? `clinical table · ${clinicalCols.length} columns` : undefined,
              gencode: gencodeLabel ?? undefined,
            }}
          />
        </Panel>
      )}

      {sid && (
        <Panel
          n={2}
          title="Gene set"
          done={geneSetMode === "mad" ? madTopN > 0 : !!features}
          disabled={!samples.length || (matrixKind !== "gene" && !gencodeLabel && !hasJunctionMetadata)}
        >
          <GeneSetPanel
            sessionId={sid}
            pathwaysEnabled={pathwaysEnabled}
            matrixKind={matrixKind}
            activeSjdat={activeSjdat}
            fastLookupEnabled={hasJunctionMetadata}
            onResolved={() => setFeatures(null)}
            onFeatures={setFeatures}
            onModeChange={setGeneSetMode}
            onTopNChange={setMadTopN}
            onCodingOnlyChange={setMadCodingOnly}
            nMin={madNMin}
            xMin={madXMin}
            onNMinChange={setMadNMin}
            onXMinChange={setMadXMin}
            nNonzeroRows={sjdatOptions.find((o) => o.kind === "pathway_matrix")?.n_nonzero_rows ?? null}
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
                {activeSjdat === "pathway_matrix"
                  ? "pathways"
                  : matrixKind === "gene"
                    ? madCodingOnly ? "protein-coding genes" : "genes"
                    : "junctions"}{" "}
                by MAD
              </button>
              {madBusy && <span className="muted">computing…</span>}
              {needsRef && (
                <span className="muted">choose a GENCODE reference on the Data tab to filter to protein-coding genes</span>
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
            <label
              style={{ flexDirection: "row", alignItems: "center", gap: 6 }}
              title={
                selected.length === 0
                  ? "select a clinical feature above first"
                  : "rebuilds the heatmap without the dropped samples"
              }
            >
              <input
                type="checkbox"
                checked={dropMissingClinical}
                disabled={selected.length === 0}
                onChange={(e) => setDropMissingClinical(e.target.checked)}
              />
              drop samples missing selected clinical data
            </label>
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
