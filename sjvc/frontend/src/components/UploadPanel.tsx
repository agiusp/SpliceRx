import { useEffect, useState } from "react";
import { api, ApiError, type ClinicalColumn } from "../api";

interface Props {
  sessionId: string;
  onJunctions: (samples: string[], warnings: string[], featureKind: "junction" | "gene") => void;
  onClinical: (cols: ClinicalColumn[], nMatched: number, warnings: string[]) => void;
  onGencode: (label: string, pathwaysEnabled: boolean, warnings: string[]) => void;
}

export default function UploadPanel({ sessionId, onJunctions, onClinical, onGencode }: Props) {
  const [releases, setReleases] = useState<Record<string, string[]>>({});
  const [species, setSpecies] = useState("human");
  const [release, setRelease] = useState("");
  const [customGtf, setCustomGtf] = useState(false);
  const [names, setNames] = useState<{ j?: string; c?: string; g?: string }>({});
  const [busy, setBusy] = useState<string | null>(null);
  const [err, setErr] = useState<string | null>(null);

  useEffect(() => {
    api.releases().then((r) => {
      setReleases(r);
      setRelease(r[species]?.[0] ?? "");
    });
  }, []); // eslint-disable-line

  const run = async (label: string, fn: () => Promise<void>) => {
    setBusy(label);
    setErr(null);
    try {
      await fn();
    } catch (e) {
      setErr(e instanceof ApiError ? e.message : String(e));
    } finally {
      setBusy(null);
    }
  };

  // Show the picked filename the instant it's chosen — before the (possibly
  // slow) upload/parse even starts, and regardless of whether it succeeds —
  // so it's never ambiguous which file is in flight or which one errored.
  function pick(slot: keyof typeof names, file: File, label: string, fn: (f: File) => Promise<void>) {
    setNames((n) => ({ ...n, [slot]: file.name }));
    run(`${label} “${file.name}”…`, () => fn(file));
  }

  return (
    <>
      <div className="row">
        <label>
          Junction matrix (.rds)
          <span className="file-trigger">
            <span className="file-btn">Choose File</span>
            <input
              type="file"
              accept=".rds"
              onChange={(e) =>
                e.target.files?.[0] &&
                pick("j", e.target.files[0], "Reading", async (f) => {
                  const r = await api.uploadJunctions(sessionId, f);
                  onJunctions(r.samples, r.warnings, r.feature_kind);
                })
              }
            />
          </span>
        </label>
        {names.j && <span className="chip">✓ {names.j}</span>}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label>
          Clinical table — needs a sample_id column (.csv / .tsv / .rds)
          <span className="file-trigger">
            <span className="file-btn">Choose File</span>
            <input
              type="file"
              accept=".csv,.tsv,.txt,.rds"
              onChange={(e) =>
                e.target.files?.[0] &&
                pick("c", e.target.files[0], "Reading", async (f) => {
                  const r = await api.uploadClinical(sessionId, f);
                  onClinical(r.columns, r.n_matched, r.warnings);
                })
              }
            />
          </span>
        </label>
        {names.c && <span className="chip">✓ {names.c}</span>}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
          <input type="checkbox" checked={customGtf} onChange={(e) => setCustomGtf(e.target.checked)} />
          Upload a custom GTF instead of a GENCODE release
        </label>
      </div>

      {!customGtf ? (
        <div className="row" style={{ marginTop: 8 }}>
          <label>
            Species
            <select
              value={species}
              onChange={(e) => {
                setSpecies(e.target.value);
                setRelease(releases[e.target.value]?.[0] ?? "");
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
            disabled={!release || !!busy}
            onClick={() =>
              run("Fetching / indexing annotation… (first use of a release can take a minute)", async () => {
                const r = await api.selectGencode(sessionId, species, release);
                setNames((n) => ({ ...n, g: r.label }));
                onGencode(r.label, r.pathways_enabled, r.warnings);
              })
            }
          >
            Use this release
          </button>
          {names.g && <span className="chip">✓ {names.g}</span>}
        </div>
      ) : (
        <div className="row" style={{ marginTop: 8 }}>
          <label>
            Annotation (.gtf / .gtf.gz)
            <span className="file-trigger">
              <span className="file-btn">Choose File</span>
              <input
                type="file"
                accept=".gtf,.gz"
                onChange={(e) =>
                  e.target.files?.[0] &&
                  pick("g", e.target.files[0], "Indexing", async (f) => {
                    const r = await api.uploadGtf(sessionId, f);
                    onGencode(r.label, r.pathways_enabled, r.warnings);
                  })
                }
              />
            </span>
          </label>
          {names.g && <span className="chip">✓ {names.g}</span>}
        </div>
      )}

      {busy && <div className="warn">{busy}</div>}
      {err && <div className="err">{err}</div>}
    </>
  );
}
