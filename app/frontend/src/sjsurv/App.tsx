import { useEffect, useState } from "react";
import { dataloadApi, type PushGroupsResult } from "../dataload/api";
import {
  ALL_GROUPS,
  api,
  ApiError,
  type AgeBands,
  type CovariateColumn,
  type CVResponse,
  type FeaturesResponse,
  type GroupCount,
  type HistologyCount,
  type MetadataLoaded,
  type ModelResponse,
  type SelectResponse,
  type SessionState,
  type SjdatKind,
} from "./api";
import GeneSetPanel, { type GeneSetMode } from "./components/GeneSetPanel";

// The default, exactly as the old R-side surv_cohort() had it: ages of 30 or
// under are left unbanded (excluded from Group) rather than forming a band.
const CLASSIC_EDGES = [30, 50, 70];
const CLASSIC_MIN_GROUP_N = 20;

/** Resize an edges array to `n` entries, extrapolating new ones by the last step. */
function resizeEdges(edges: number[], n: number): number[] {
  if (edges.length === n) return edges;
  if (edges.length > n) return edges.slice(0, n);
  const out = edges.slice();
  while (out.length < n) {
    const step = out.length >= 2 ? out[out.length - 1] - out[out.length - 2] : 10;
    out.push(Math.round((out[out.length - 1] + (step || 10)) * 10) / 10);
  }
  return out;
}

interface Props {
  /** Session id, owned by the Shell so the Data tab can load files into it. */
  sessionId: string | null;
  /** Bump to re-hydrate from the server (after a Data-tab load). */
  reloadNonce: number;
  /** The other two apps' session ids, so this tab's Group/SurviverGroup can be
   * pushed into their sample-metadata / clinical tables. */
  sjvSessionId: string | null;
  sjvcSessionId: string | null;
  /** Fired once per target ("sjv" | "sjvc") a push actually updated, so the
   * Shell can re-hydrate that tab from its now-changed session. */
  onGroupsPushed: (target: "sjv" | "sjvc") => void;
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

const pct = (x: number) => (Number.isFinite(x) ? `${(x * 100).toFixed(1)}%` : "–");
const f3 = (x: number) => (Number.isFinite(x) ? x.toFixed(3) : "–");

export default function App({ sessionId, reloadNonce, sjvSessionId, sjvcSessionId, onGroupsPushed }: Props) {
  const [state, setState] = useState<SessionState | null>(null);
  const [groups, setGroups] = useState<GroupCount[]>([]);
  // most Groups a cohort fragments into have nobody labelled (see the
  // heterogeneous-histology case below) — hide those by default
  const [showAllGroups, setShowAllGroups] = useState(false);
  const [warnings, setWarnings] = useState<string[]>([]);
  const [err, setErr] = useState<string | null>(null);

  // age-band stratification editor. `edges` always holds nBands+1 numbers;
  // when `openEnded` the last one is sent as "no upper bound" rather than its
  // literal value (so it survives being toggled back off).
  const [nBands, setNBands] = useState(3);
  const [edges, setEdges] = useState<number[]>([...CLASSIC_EDGES, CLASSIC_EDGES[CLASSIC_EDGES.length - 1]]);
  const [openEnded, setOpenEnded] = useState(true);
  const [includeLowest, setIncludeLowest] = useState(false);
  const [minGroupN, setMinGroupN] = useState(CLASSIC_MIN_GROUP_N);
  const [appliedLabels, setAppliedLabels] = useState<string[]>([]);
  const [appliedMinGroupN, setAppliedMinGroupN] = useState(CLASSIC_MIN_GROUP_N);
  const [suggestedLabels, setSuggestedLabels] = useState<string[] | null>(null);
  const [suggestBusy, setSuggestBusy] = useState(false);
  const [stratifyBusy, setStratifyBusy] = useState(false);
  const [pushBusy, setPushBusy] = useState(false);
  const [pushResult, setPushResult] = useState<PushGroupsResult | null>(null);

  // histology: drop it from Group entirely, or merge selected raw values
  // together (some cancers are heterogeneous enough that splitting on every
  // distinct histology leaves cohorts too small to label)
  const [useHistology, setUseHistology] = useState(true);
  const [histologyCounts, setHistologyCounts] = useState<HistologyCount[]>([]);
  const [histologyMergeAs, setHistologyMergeAs] = useState<Record<string, string>>({});

  // clinical covariates: sample aspects from the loaded metadata file added
  // alongside the selected molecular features when cross-validating/training
  // — defaults to Stage, Age At Diagnosis, Histology, and any MSI column
  // present (see services/metadata.py's available_covariates())
  const [covariateColumns, setCovariateColumns] = useState<CovariateColumn[]>([]);
  const [selectedCovariates, setSelectedCovariates] = useState<string[]>([]);
  const [covBusy, setCovBusy] = useState(false);
  // the full list can run long (real cohort files carry dozens of columns) —
  // collapsed by default, with a search box once expanded
  const [covariatesExpanded, setCovariatesExpanded] = useState(false);
  const [covariateSearch, setCovariateSearch] = useState("");

  // Feature selection: "Top by MAD" (nMin/xMin/topN, unchanged from before)
  // or one of the gene-set tabs (Type genes / Upload list / Pathway), whose
  // resolve + build happens inside <GeneSetPanel> — `features` tracks
  // whether it has produced a built feature matrix ready to select from.
  const [geneSetMode, setGeneSetMode] = useState<GeneSetMode>("mad");
  const [nMin, setNMin] = useState(0.2);
  const [xMin, setXMin] = useState(1);
  const [topN, setTopN] = useState(100);
  const [codingOnly, setCodingOnly] = useState(false);
  const [features, setFeatures] = useState<FeaturesResponse | null>(null);
  const [nCv, setNCv] = useState(5);

  const [sel, setSel] = useState<SelectResponse | null>(null);
  const [cv, setCv] = useState<CVResponse | null>(null);
  const [model, setModel] = useState<ModelResponse | null>(null);
  const [featureSearch, setFeatureSearch] = useState("");
  const [busy, setBusy] = useState<"select" | "cv" | "model" | null>(null);

  // hydrate from the server on mount / after a Data-tab load
  useEffect(() => {
    if (!sessionId) return;
    const sid = sessionId;
    api
      .sessionState(sid)
      .then(async (st) => {
        setState(st);
        const bands = st.metadata?.age_bands;
        if (bands) seedBandEditor(bands, st.metadata?.min_group_n ?? CLASSIC_MIN_GROUP_N);
        if (st.metadata) seedHistologyEditor(st.metadata);
        setCovariateColumns(st.covariate_columns);
        setSelectedCovariates(st.selected_covariates);
        if (st.has_metadata) {
          const g = await api.groups(sid);
          setGroups(g.groups);
          setWarnings(g.warnings);
        }
      })
      .catch((e) => setErr(String(e)));
  }, [sessionId, reloadNonce]);

  function seedBandEditor(bands: AgeBands, minGroup: number) {
    const last = bands.edges[bands.edges.length - 1];
    const open = last === null;
    const finite = bands.edges.map((e, i) => (e === null ? bands.edges[i - 1] ?? 100 : e)) as number[];
    setNBands(bands.edges.length - 1);
    setEdges(finite);
    setOpenEnded(open);
    setIncludeLowest(bands.include_lowest);
    setMinGroupN(minGroup);
    setAppliedLabels(bands.labels);
    setAppliedMinGroupN(minGroup);
    setSuggestedLabels(null);
  }

  function seedHistologyEditor(md: MetadataLoaded) {
    setHistologyCounts(md.histology_counts);
    setUseHistology(md.use_histology);
    setHistologyMergeAs(
      Object.fromEntries(md.histology_counts.map((c) => [c.value, md.histology_map[c.value] ?? c.value])),
    );
  }

  function clearDownstream() {
    setSel(null);
    setCv(null);
    setModel(null);
  }

  async function run<T>(tag: "select" | "cv" | "model", fn: () => Promise<T>) {
    setBusy(tag);
    setErr(null);
    try {
      return await fn();
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
      return null;
    } finally {
      setBusy(null);
    }
  }

  async function pickSjdat(kind: SjdatKind) {
    if (!sessionId) return;
    setErr(null);
    clearDownstream();
    setFeatures(null); // the backend also clears any resolved geneset/features on this switch
    try {
      setState(await api.activateSjdat(sessionId, kind));
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    }
  }

  function changeNBands(n: number) {
    const clamped = Math.max(1, Math.round(n) || 1);
    setNBands(clamped);
    setEdges((cur) => resizeEdges(cur, clamped + 1));
  }

  function resetClassic() {
    setNBands(3);
    setEdges([...CLASSIC_EDGES, CLASSIC_EDGES[CLASSIC_EDGES.length - 1]]);
    setOpenEnded(true);
    setIncludeLowest(false);
    setMinGroupN(CLASSIC_MIN_GROUP_N);
    setSuggestedLabels(null);
  }

  async function doSuggestBands() {
    if (!sessionId) return;
    setSuggestBusy(true);
    setErr(null);
    try {
      const bands = await api.suggestAgeBands(sessionId, nBands);
      const last = bands.edges[bands.edges.length - 1];
      const finite = bands.edges.map((e, i) => (e === null ? bands.edges[i - 1] ?? 100 : e)) as number[];
      setEdges(finite);
      setOpenEnded(last === null);
      setIncludeLowest(bands.include_lowest);
      setSuggestedLabels(bands.labels);
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setSuggestBusy(false);
    }
  }

  async function doApplyStratify() {
    if (!sessionId) return;
    setStratifyBusy(true);
    setErr(null);
    try {
      const wireEdges: (number | null)[] = edges.slice(0, nBands + 1);
      if (openEnded) wireEdges[wireEdges.length - 1] = null;
      // only send merges that actually change something
      const histologyMap = Object.fromEntries(
        Object.entries(histologyMergeAs).filter(([raw, as]) => as.trim() && as.trim() !== raw),
      );
      const md = await api.stratify(sessionId, {
        edges: wireEdges, include_lowest: includeLowest, min_group_n: minGroupN,
        use_histology: useHistology, histology_map: histologyMap,
      });
      setAppliedLabels(md.age_bands?.labels ?? []);
      setAppliedMinGroupN(md.min_group_n ?? minGroupN);
      seedHistologyEditor(md);
      clearDownstream();
      const g = await api.groups(sessionId);
      setGroups(g.groups);
      setWarnings(g.warnings);
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setStratifyBusy(false);
    }
  }

  async function doPushGroups() {
    if (!sessionId) return;
    setPushBusy(true);
    setErr(null);
    setPushResult(null);
    try {
      const r = await dataloadApi.pushGroups(sessionId, sjvSessionId, sjvcSessionId);
      setPushResult(r);
      for (const target of r.applied_to) onGroupsPushed(target as "sjv" | "sjvc");
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setPushBusy(false);
    }
  }

  async function doSelect() {
    if (!sessionId) return;
    clearDownstream();
    const showCoverage = state?.active_sjdat === "junction_counts" || state?.active_sjdat === "rrs_scores";
    const r = await run("select", () =>
      api.select(sessionId, {
        group: ALL_GROUPS, top_n: topN,
        n_min: showCoverage ? nMin : null, x_min: showCoverage ? xMin : null,
        protein_coding_only: codingOnly,
      }),
    );
    if (r) {
      setSel(r);
      setWarnings(r.warnings);
    }
  }

  async function doSelectGeneset() {
    if (!sessionId) return;
    clearDownstream();
    const r = await run("select", () => api.selectGeneset(sessionId, ALL_GROUPS));
    if (r) {
      setSel(r);
      setWarnings(r.warnings);
    }
  }

  async function doCv() {
    if (!sessionId) return;
    const r = await run("cv", () => api.crossValidate(sessionId, nCv));
    if (r) setCv(r);
  }

  async function doModel() {
    if (!sessionId) return;
    const r = await run("model", () => api.trainModel(sessionId));
    if (r) {
      setModel(r);
      setFeatureSearch("");
    }
  }

  async function doToggleCovariate(key: string) {
    if (!sessionId) return;
    const next = selectedCovariates.includes(key)
      ? selectedCovariates.filter((k) => k !== key)
      : [...selectedCovariates, key];
    setCovBusy(true);
    setErr(null);
    try {
      const st = await api.setCovariates(sessionId, next);
      setSelectedCovariates(st.selected_covariates);
      setCv(null);
      setModel(null);
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setCovBusy(false);
    }
  }

  const opts = state?.sjdat_options ?? [];
  const activeOpt = opts.find((o) => o.kind === state?.active_sjdat && o.loaded);
  const anyLoaded = opts.some((o) => o.loaded);
  const metaReady = !!state?.has_metadata;
  const allGroup = groups.find((g) => g.group === ALL_GROUPS);
  const groupChosen = !!allGroup && allGroup.n_labelled > 0;
  const selectReady = metaReady && groupChosen && !!activeOpt;

  return (
    <div className="app">
      <h1>SJSurv — Splice-junction survivor groups</h1>
      <p className="sub">
        Does a cohort's splice-junction signal separate <b>Good</b> from <b>Poor</b> survivors?
        SJSurv stratifies the cohort right here — by histology, pathologic stage and an age band
        you control — and labels each stratum's samples Good/Poor by their own survival against
        their stratum's median. The classifier itself always trains on every Good/Poor-labelled
        sample: pick a feature matrix, select molecular features and clinical covariates, then
        cross-validate a classifier of the survivor label.
      </p>

      {err && <div className="err">{err}</div>}
      {[...new Set(warnings)].map((w, i) => (
        <div className="warn" key={i}>
          {w}
        </div>
      ))}

      {/* 1 · Data ------------------------------------------------------- */}
      <Panel n={1} title="Data" done={anyLoaded && metaReady}>
        <div className="row">
          <label style={{ display: "block" }}>
            Sample metadata (histology, stage, age at diagnosis, survival)
          </label>
          {metaReady ? (
            <span className="chip">
              ✓ {state!.metadata!.n_matched} sample(s)
              {state!.metadata!.n_unmatched ? `, ${state!.metadata!.n_unmatched} unmatched` : ""}
              {state!.metadata!.n_with_age < state!.metadata!.n_matched
                ? `, ${state!.metadata!.n_matched - state!.metadata!.n_with_age} without an age`
                : ""}
            </span>
          ) : (
            <span className="muted">Load a cohort on the Data tab.</span>
          )}
        </div>
        <div className="row" style={{ marginTop: 10 }}>
          <label style={{ display: "block" }}>sjdat feature matrices</label>
          {anyLoaded ? (
            <span className="chip">
              ✓ {opts.filter((o) => o.loaded).map((o) => o.label).join(", ")}
            </span>
          ) : (
            <span className="muted">
              Load one or more of junction counts / RRS scores / novel junction counts per gene /
              novel junction counts per pathway on the Data tab.
            </span>
          )}
        </div>
        <div className="row" style={{ marginTop: 10 }}>
          <label style={{ display: "block" }}>Junction gene annotation</label>
          {state?.has_junction_metadata ? (
            <span className="chip">✓ loaded</span>
          ) : (
            <span className="muted">
              Optional — load the junction metadata file on the Data tab for fast typed-gene /
              pathway lookups against a junction-level matrix (junction counts or RRS scores)
              without needing a GENCODE reference. Not needed for a gene- or pathway-level matrix,
              which is already condensed.
            </span>
          )}
        </div>
        <div className="row" style={{ marginTop: 10 }}>
          <label style={{ display: "block" }}>GENCODE reference</label>
          {state?.gencode_label ? (
            <span className="chip">✓ {state.gencode_label}</span>
          ) : (
            <span className="muted">
              Optional — choose a GENCODE release on the Data tab. Needed for a junction-level
              matrix's typed-gene / pathway lookups only when no junction metadata is loaded.
            </span>
          )}
        </div>
      </Panel>

      {/* 2 · Group & age-band stratification --------------------------- */}
      <Panel
        n={2}
        title="Group and Age-Band Stratification for Good/Poor survivor labels"
        done={groupChosen}
        disabled={!metaReady}
      >
        <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
          Each sample's <code>Group</code> is {useHistology ? "Histology | " : ""}Stage | age band;
          within a <code>Group</code> with at least <b>{appliedMinGroupN}</b> samples, a sample is
          labelled <b>Good</b> if its survival is at or above that group's median, <b>Poor</b>{" "}
          otherwise. Currently stratifying by age into: <b>{appliedLabels.join(", ") || "…"}</b>.
          Adjust the age bands and histology grouping below to try to grow the number of Good/Poor
          labels — but every classifier downstream always trains on <b>every</b> labelled sample
          (the "{groups[0]?.label ?? "All labelled samples"}" row below), regardless of which
          stratification Group a sample falls into.
        </p>
        {(() => {
          // the aggregate "All labelled samples" row is always first (see
          // services/metadata.py's group_counts()) and always shown; every
          // other Group is hidden by default once it has nobody labelled
          const alwaysVisible = groups[0]?.group;
          const visibleGroups = showAllGroups
            ? groups
            : groups.filter((g) => g.group === alwaysVisible || g.n_good > 0 || g.n_poor > 0);
          const hiddenCount = groups.length - visibleGroups.length;
          return (
            <>
              <div className="scroll-x">
                <table className="grid">
                  <thead>
                    <tr>
                      <th>Group</th>
                      <th className="num">Samples</th>
                      <th className="num">Good</th>
                      <th className="num">Poor</th>
                      <th className="num">Labelled</th>
                    </tr>
                  </thead>
                  <tbody>
                    {visibleGroups.map((g) => (
                      <tr key={g.group} className={g.group === ALL_GROUPS ? "sel" : ""}>
                        <td>{g.label}</td>
                        <td className="num">{g.n_total}</td>
                        <td className="num">{g.n_good}</td>
                        <td className="num">{g.n_poor}</td>
                        <td className="num">{g.n_labelled}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              {hiddenCount > 0 || showAllGroups ? (
                <button
                  type="button"
                  className="ghost"
                  style={{ marginTop: 8 }}
                  onClick={() => setShowAllGroups((v) => !v)}
                >
                  {showAllGroups
                    ? "Show only groups with labelled samples"
                    : `Show all groups (${hiddenCount} hidden with no labelled samples)`}
                </button>
              ) : null}
            </>
          );
        })()}
        {allGroup && allGroup.n_labelled > 0 && allGroup.n_labelled < 12 && (
          <div className="warn">
            Only {allGroup.n_labelled} labelled sample(s) total — the classifier will be unstable.
            Try coarser age bands, a lower per-group minimum, or dropping/merging histology below.
          </div>
        )}
        {groups.length > 0 && groups[0].n_labelled === 0 && (
          <div className="warn">
            No sample cleared the {appliedMinGroupN}-sample minimum under the current age bands —
            nobody has a survivor label yet. Try fewer, coarser age bands or a lower minimum below.
          </div>
        )}

        <div style={{ marginTop: 14, paddingTop: 12, borderTop: "1px solid var(--border)" }}>
          <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
            This <code>Group</code> / <code>SurviverGroup</code> is computed here, from the age
            bands above — not read from whatever the metadata file itself happens to carry. Push
            it to the other tabs to use as a sample-group series in the Sashimi plot, or as a
            clinical feature (coloring, PCA/UMAP) in 2D View. Pushing again after changing the age
            bands replaces what's there.
          </p>
          <div className="row" style={{ alignItems: "center" }}>
            <button
              className="ghost"
              disabled={!metaReady || pushBusy || (!sjvSessionId && !sjvcSessionId)}
              onClick={doPushGroups}
            >
              {pushBusy ? "Pushing…" : "Push groups to Sashimi plot & 2D View"}
            </button>
            {!sjvSessionId && !sjvcSessionId && (
              <span className="muted">connecting to the other tabs…</span>
            )}
          </div>
          {pushResult && (
            <div className="muted" style={{ fontSize: 13, marginTop: 6 }}>
              {pushResult.applied_to.length > 0 ? (
                <>
                  ✓ pushed to {pushResult.applied_to.map((t) => (t === "sjv" ? "Sashimi plot" : "2D View")).join(", ")}
                  {" — "}
                  {Object.entries(pushResult.matched)
                    .map(([t, n]) => `${n} sample(s) in ${t === "sjv" ? "Sashimi plot" : "2D View"}`)
                    .join(", ")}
                  {" received a Group/SurviverGroup value"}
                  {pushResult.n_labelled < pushResult.n_sjsurv_samples
                    ? ` (only ${pushResult.n_labelled} of this tab's ${pushResult.n_sjsurv_samples} samples are Good/Poor-labelled)`
                    : ""}
                  .
                </>
              ) : (
                "nothing was pushed — see the warning above."
              )}
            </div>
          )}
          {pushResult?.warnings.map((w, i) => (
            <div className="warn" key={i}>
              {w}
            </div>
          ))}
        </div>

        <div style={{ marginTop: 14, paddingTop: 12, borderTop: "1px solid var(--border)" }}>
          <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
            Some cancers are heterogeneous enough that splitting on every distinct histology value
            leaves each <code>Group</code> too small to label. Drop histology from <code>Group</code>{" "}
            entirely below, or give two or more values the same "merge as" label to combine them —
            either way, hit "Apply stratification" to recompute.
          </p>
          <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
            <input
              type="checkbox"
              checked={useHistology}
              onChange={(e) => setUseHistology(e.target.checked)}
            />
            Include histology in Group
          </label>
          {useHistology && histologyCounts.length > 0 && (
            <div className="scroll-x" style={{ marginTop: 8 }}>
              <table className="grid">
                <thead>
                  <tr>
                    <th>Histology value</th>
                    <th className="num">Samples</th>
                    <th>Merge as</th>
                  </tr>
                </thead>
                <tbody>
                  {histologyCounts.map((c) => (
                    <tr key={c.value}>
                      <td className="mono">{c.value}</td>
                      <td className="num">{c.n}</td>
                      <td>
                        <input
                          type="text"
                          value={histologyMergeAs[c.value] ?? c.value}
                          style={{ width: "100%" }}
                          onChange={(e) =>
                            setHistologyMergeAs((cur) => ({ ...cur, [c.value]: e.target.value }))
                          }
                        />
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <button
                className="ghost"
                style={{ marginTop: 6 }}
                onClick={() =>
                  setHistologyMergeAs(Object.fromEntries(histologyCounts.map((c) => [c.value, c.value])))
                }
              >
                Reset merges
              </button>
            </div>
          )}
        </div>

        <div style={{ marginTop: 14, paddingTop: 12, borderTop: "1px solid var(--border)" }}>
          <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
            The fixed 30 / 50 / 70 cut points above are a reasonable default, but they can fragment
            a cohort — especially a small one — into groups too small to ever clear the label
            minimum (as just happened here for PAAD). There may be a better way to stratify by age:
            a split based on the age distribution <i>actually observed in this cohort</i>{" "}
            {state?.metadata?.age_min != null && state?.metadata?.age_max != null && (
              <>({state.metadata.age_min.toFixed(0)}–{state.metadata.age_max.toFixed(0)} years,{" "}
                {state.metadata.n_with_age} sample(s) with an age)</>
            )}
            .
          </p>
          <div className="row" style={{ alignItems: "flex-end" }}>
            <label>
              Number of age bands
              <input
                type="number"
                min={1}
                step={1}
                value={nBands}
                style={{ width: 90 }}
                onChange={(e) => changeNBands(Number(e.target.value))}
              />
            </label>
            <button className="ghost" disabled={suggestBusy} onClick={doSuggestBands}>
              {suggestBusy ? "Computing…" : `Suggest a ${nBands}-band split from this cohort`}
            </button>
            <button className="ghost" onClick={resetClassic}>
              Reset to classic (30 / 50 / 70)
            </button>
          </div>
          {suggestedLabels && (
            <div className="muted" style={{ fontSize: 12, marginTop: 4 }}>
              suggested: {suggestedLabels.join(", ")} — edit below before applying if you like
            </div>
          )}

          <div className="row" style={{ marginTop: 10 }}>
            {Array.from({ length: nBands + 1 }, (_, i) => {
              const isOpenEnd = openEnded && i === nBands;
              return (
                <label key={i}>
                  {i === 0 ? "from" : i === nBands ? "to" : `edge ${i}`}
                  <input
                    type="number"
                    step={1}
                    value={isOpenEnd ? "" : edges[i] ?? ""}
                    placeholder={isOpenEnd ? "∞" : undefined}
                    disabled={isOpenEnd}
                    style={{ width: 90 }}
                    onChange={(e) => {
                      const v = Number(e.target.value);
                      setEdges((cur) => cur.map((x, j) => (j === i ? v : x)));
                    }}
                  />
                </label>
              );
            })}
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
              <input type="checkbox" checked={openEnded} onChange={(e) => setOpenEnded(e.target.checked)} />
              no upper bound on the last band
            </label>
            <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
              <input
                type="checkbox"
                checked={includeLowest}
                onChange={(e) => setIncludeLowest(e.target.checked)}
              />
              include the very lowest age (otherwise ages at or below "from" are excluded)
            </label>
          </div>

          <div className="row" style={{ marginTop: 10 }}>
            <label>
              Min samples per Group for a survivor label
              <input
                type="number"
                min={1}
                step={1}
                value={minGroupN}
                style={{ width: 90 }}
                onChange={(e) => setMinGroupN(Math.max(1, Math.round(Number(e.target.value) || 1)))}
              />
            </label>
            <button disabled={stratifyBusy} onClick={doApplyStratify}>
              {stratifyBusy ? "Applying…" : "Apply stratification"}
            </button>
          </div>
          <p className="muted" style={{ fontSize: 12, marginTop: 6 }}>
            e.g. for just two age bands (under 50 / 50 and over), set "Number of age bands" to 2 and
            edit the single edge to 50.
          </p>
        </div>
      </Panel>

      {/* 3 · sjdat ------------------------------------------------------ */}
      <Panel n={3} title="Splice-junction data (sjdat)" done={!!activeOpt} disabled={!metaReady}>
        <div className="radio-row">
          {opts.map((o) => (
            <label key={o.kind} className={o.loaded ? "" : "disabled"}>
              <input
                type="radio"
                disabled={!o.loaded}
                checked={state?.active_sjdat === o.kind && o.loaded}
                onChange={() => pickSjdat(o.kind)}
              />
              <span>
                <b>{o.label}</b>
                {o.loaded ? (
                  <span className="mono">
                    {" "}
                    — {o.n_features.toLocaleString()} {o.kind === "gene_matrix" ? "genes" : "features"} ×{" "}
                    {o.n_samples} samples{o.sparse ? " (sparse)" : ""}
                  </span>
                ) : (
                  <span className="muted"> — not loaded on the Data tab</span>
                )}
                <span className="desc">{o.description}</span>
              </span>
            </label>
          ))}
        </div>
      </Panel>

      {/* 4 · Feature selection --------------------------------------- */}
      <Panel n={4} title="Feature selection" done={!!sel} disabled={!selectReady}>
        <div style={{ marginBottom: 14, paddingBottom: 12, borderBottom: "1px solid var(--border)" }}>
          <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
            Clinical covariates — sample aspects from the loaded metadata file, added to the
            classifier alongside the molecular features selected below. Stage, Age At Diagnosis,
            Histology, and any MSI column are on by default; toggle any others the file has.
          </p>
          <div className="row" style={{ alignItems: "center" }}>
            <button
              type="button"
              className="ghost"
              onClick={() => setCovariatesExpanded((v) => !v)}
              disabled={covariateColumns.length === 0}
            >
              {covariatesExpanded
                ? "Hide covariate list"
                : `Edit covariates (${covariateColumns.length} available)`}
            </button>
            <span className="muted" style={{ fontSize: 13 }}>
              {covariateColumns.length === 0
                ? "no covariate columns available"
                : selectedCovariates.length > 0
                  ? `${selectedCovariates.length} selected: ${covariateColumns
                      .filter((c) => selectedCovariates.includes(c.key))
                      .map((c) => c.label)
                      .join(", ")}`
                  : "none selected"}
            </span>
          </div>
          {covariatesExpanded && covariateColumns.length > 0 && (
            <>
              <label style={{ display: "block", maxWidth: 320, marginTop: 10 }}>
                Search
                <input
                  type="text"
                  placeholder="filter covariates by name…"
                  value={covariateSearch}
                  onChange={(e) => setCovariateSearch(e.target.value)}
                />
              </label>
              {(() => {
                const q = covariateSearch.trim().toLowerCase();
                const filtered = q
                  ? covariateColumns.filter((c) => c.label.toLowerCase().includes(q))
                  : covariateColumns;
                return (
                  <div className="scroll-x" style={{ marginTop: 8 }}>
                    <table className="grid">
                      <thead>
                        <tr>
                          <th style={{ width: 28 }} />
                          <th>Covariate</th>
                          <th>Type</th>
                          <th className="num">Available</th>
                        </tr>
                      </thead>
                      <tbody>
                        {filtered.map((c) => (
                          <tr key={c.key}>
                            <td>
                              <input
                                type="checkbox"
                                checked={selectedCovariates.includes(c.key)}
                                disabled={covBusy}
                                onChange={() => doToggleCovariate(c.key)}
                              />
                            </td>
                            <td>{c.label}</td>
                            <td className="muted">{c.kind}</td>
                            <td className="num">{c.n_available}</td>
                          </tr>
                        ))}
                        {filtered.length === 0 && (
                          <tr>
                            <td colSpan={4} className="muted">
                              no covariate matches "{covariateSearch}"
                            </td>
                          </tr>
                        )}
                      </tbody>
                    </table>
                  </div>
                );
              })()}
            </>
          )}
        </div>
        {sessionId && (
          <GeneSetPanel
            sessionId={sessionId}
            pathwaysEnabled={!!state?.pathways_enabled}
            matrixKind={
              state?.active_sjdat === "junction_counts" || state?.active_sjdat === "rrs_scores"
                ? "junction"
                : "gene"
            }
            activeSjdat={(state?.active_sjdat as SjdatKind | undefined) ?? null}
            fastLookupEnabled={!!state?.has_junction_metadata}
            onResolved={() => setFeatures(null)}
            onFeatures={setFeatures}
            onModeChange={setGeneSetMode}
            onTopNChange={setTopN}
            onCodingOnlyChange={setCodingOnly}
            nMin={nMin}
            xMin={xMin}
            onNMinChange={setNMin}
            onXMinChange={setXMin}
            nNonzeroRows={opts.find((o) => o.kind === "pathway_matrix")?.n_nonzero_rows ?? null}
          />
        )}

        <div className="row" style={{ marginTop: 10, alignItems: "center" }}>
          <button
            disabled={
              !selectReady || busy === "select" || (geneSetMode !== "mad" && !features)
              || (geneSetMode === "mad" && codingOnly && !state?.gencode_label)
            }
            onClick={geneSetMode === "mad" ? doSelect : doSelectGeneset}
          >
            {busy === "select" ? "Selecting…" : "Select features"}
          </button>
          {geneSetMode === "mad" && codingOnly && !state?.gencode_label && (
            <span className="muted">choose a GENCODE reference on the Data tab to filter to protein-coding genes</span>
          )}
        </div>

        {sel && (
          <div className="muted" style={{ marginTop: 10, fontSize: 13 }}>
            <div>
              <b>{sel.n_selected}</b> feature(s) selected
              {geneSetMode === "mad" ? (
                <>
                  {" "}from {sel.n_candidates.toLocaleString()} candidates (
                  {sel.n_after_coverage.toLocaleString()} passed the N / X filter).
                </>
              ) : (
                "."
              )}
              {selectedCovariates.length > 0 && (
                <> Plus {selectedCovariates.length} clinical covariate(s) below.</>
              )}
            </div>
            <div>
              {sel.n_group_samples} samples — Good {sel.n_good}, Poor {sel.n_poor}.
            </div>
            {sel.feature_preview.length > 0 && (
              <div>
                preview: {sel.feature_preview.join(", ")}
                {sel.n_selected > sel.feature_preview.length ? " …" : ""}
              </div>
            )}
          </div>
        )}
      </Panel>

      {/* 5 · Cross-validation -------------------------------------------- */}
      <Panel n={5} title="Cross-validated classification" done={!!cv} disabled={!sel}>
        <p className="muted" style={{ marginTop: 0, fontSize: 13 }}>
          The fold split uses a fixed random seed, so re-running with the same selected features
          and ncv always gives the same AUC — it only changes if the selection, group, or ncv
          itself changes.
        </p>
        <div className="row" style={{ alignItems: "flex-end" }}>
          <label>
            ncv — cross-validation folds
            <input
              type="number"
              min={2}
              step={1}
              value={nCv}
              style={{ width: 120 }}
              onChange={(e) => setNCv(Math.max(2, Math.round(Number(e.target.value) || 2)))}
            />
          </label>
          <button disabled={!sel || busy === "cv"} onClick={doCv}>
            {busy === "cv" ? "Building models…" : `Run ${nCv}-fold cross-validation`}
          </button>
          {sel && (
            <span className="muted">
              predicting SurviverGroup on {sel.n_group_samples} samples, {sel.n_selected} molecular
              feature(s)
              {selectedCovariates.length > 0
                ? ` plus ${selectedCovariates.length} clinical covariate(s) (a categorical one can `
                  + `contribute more than one model feature)`
                : ""}
            </span>
          )}
        </div>

        {cv && (
          <>
            <div className="log">{cv.messages.join("\n")}</div>
            <div className="statgrid">
              <div className="stat">
                <div className="k">Out-of-fold AUC</div>
                <div className="v">{f3(cv.auc)}</div>
              </div>
              <div className="stat">
                <div className="k">Per-fold AUC</div>
                <div className="v">
                  {Number.isFinite(cv.fold_auc_mean) ? `${f3(cv.fold_auc_mean)} ± ${f3(cv.fold_auc_sd)}` : "–"}
                </div>
              </div>
              <div className="stat">
                <div className="k">Accuracy</div>
                <div className="v">{pct(cv.accuracy)}</div>
              </div>
              <div className="stat">
                <div className="k">Baseline (majority)</div>
                <div className="v">{pct(cv.baseline_accuracy)}</div>
              </div>
              <div className="stat">
                <div className="k">Sensitivity (Good)</div>
                <div className="v">{pct(cv.sensitivity)}</div>
              </div>
              <div className="stat">
                <div className="k">Specificity (Poor)</div>
                <div className="v">{pct(cv.specificity)}</div>
              </div>
            </div>
            <table className="confusion">
              <thead>
                <tr>
                  <th />
                  <th>pred Poor</th>
                  <th>pred Good</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <th>actual Poor</th>
                  <td>{cv.confusion[0][0]}</td>
                  <td>{cv.confusion[0][1]}</td>
                </tr>
                <tr>
                  <th>actual Good</th>
                  <td>{cv.confusion[1][0]}</td>
                  <td>{cv.confusion[1][1]}</td>
                </tr>
              </tbody>
            </table>
            <div className="muted" style={{ marginTop: 8, fontSize: 12 }}>
              {cv.auc > 0.6
                ? "The features carry survivor-group signal. Train a model on all samples below to see which ones."
                : "AUC is close to chance — try another sjdat, a bigger group, or looser feature filters."}
            </div>
          </>
        )}
      </Panel>

      {/* 6 · Full model + feature importance ---------------------------- */}
      <Panel n={6} title="Model on all samples — feature importance" done={!!model} disabled={!sel}>
        <div className="row">
          <button disabled={!sel || busy === "model"} onClick={doModel}>
            {busy === "model" ? "Training…" : "Train on all samples"}
          </button>
          <span className="muted">
            trained and tested on every selected sample (resubstitution — optimistic); the fitted
            model is saved to this session as MODEL
          </span>
        </div>

        {model && (
          <>
            <div className="log">{model.messages.join("\n")}</div>
            <div className="statgrid">
              <div className="stat">
                <div className="k">Resub. AUC</div>
                <div className="v">{f3(model.auc_resub)}</div>
              </div>
              <div className="stat">
                <div className="k">Resub. accuracy</div>
                <div className="v">{pct(model.accuracy_resub)}</div>
              </div>
              <div className="stat">
                <div className="k">Features</div>
                <div className="v">{model.n_features}</div>
              </div>
              <div className="stat">
                <div className="k">Samples</div>
                <div className="v">
                  {model.n_samples} <span className="muted" style={{ fontSize: 12 }}>({model.n_good}G / {model.n_poor}P)</span>
                </div>
              </div>
            </div>

            {(() => {
              const q = featureSearch.trim().toLowerCase();
              const ranked = model.features.map((f, i) => ({ ...f, rank: i + 1 }));
              const filtered = q
                ? ranked.filter(
                    (f) => f.feature.toLowerCase().includes(q) || f.direction.toLowerCase().includes(q),
                  )
                : ranked;
              return (
                <>
                  <div className="row" style={{ marginTop: 10, alignItems: "center" }}>
                    <label style={{ flex: 1, maxWidth: 320 }}>
                      Search
                      <input
                        type="text"
                        placeholder="filter by feature or direction…"
                        value={featureSearch}
                        onChange={(e) => setFeatureSearch(e.target.value)}
                      />
                    </label>
                    {q && (
                      <span className="muted">
                        {filtered.length} of {ranked.length} feature(s)
                      </span>
                    )}
                  </div>
                  <div className="tbl-scroll scroll-x">
                    <table className="grid">
                      <thead>
                        <tr>
                          <th className="num">#</th>
                          <th>Feature</th>
                          <th className="num">Weight</th>
                          <th>Direction</th>
                          <th className="num">mean Good</th>
                          <th className="num">mean Poor</th>
                        </tr>
                      </thead>
                      <tbody>
                        {filtered.map((f) => (
                          <tr key={f.feature}>
                            <td className="num">{f.rank}</td>
                            <td className="mono">{f.feature}</td>
                            <td className="num">{f.weight.toFixed(3)}</td>
                            <td>{f.direction}</td>
                            <td className="num">{f.mean_good.toFixed(2)}</td>
                            <td className="num">{f.mean_poor.toFixed(2)}</td>
                          </tr>
                        ))}
                        {filtered.length === 0 && (
                          <tr>
                            <td colSpan={6} className="muted">
                              no feature matches "{featureSearch}"
                            </td>
                          </tr>
                        )}
                      </tbody>
                    </table>
                  </div>
                </>
              );
            })()}
          </>
        )}
      </Panel>
    </div>
  );
}
