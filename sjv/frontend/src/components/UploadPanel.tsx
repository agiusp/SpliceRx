import { useEffect, useState } from "react";
import { api, ApiError } from "../api";

interface Props {
  sessionId: string;
  onRds: (samples: string[], warnings: string[]) => void;
  onAnnotation: (label: string, warnings: string[]) => void;
}

export default function UploadPanel({ sessionId, onRds, onAnnotation }: Props) {
  const [releases, setReleases] = useState<Record<string, string[]>>({});
  const [species, setSpecies] = useState("human");
  const [release, setRelease] = useState("");
  const [rdsName, setRdsName] = useState<string | null>(null);
  const [annLabel, setAnnLabel] = useState<string | null>(null);
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [useCustomGtf, setUseCustomGtf] = useState(false);

  useEffect(() => {
    api.releases().then((r) => {
      setReleases(r);
      const first = r[species]?.[0] ?? "";
      setRelease(first);
    });
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  async function handleRds(file: File) {
    setBusy("Reading RDS…");
    setError(null);
    try {
      const r = await api.uploadRds(sessionId, file);
      setRdsName(file.name);
      onRds(r.samples, r.warnings);
    } catch (e) {
      setError(e instanceof ApiError ? String(e.detail) : String(e));
    } finally {
      setBusy(null);
    }
  }

  async function handleGencode() {
    setBusy("Fetching / indexing annotation… (first use of a release can take a minute)");
    setError(null);
    try {
      const r = await api.selectGencode(sessionId, species, release);
      setAnnLabel(r.label);
      onAnnotation(r.label, r.warnings);
    } catch (e) {
      setError(e instanceof ApiError ? String(e.detail) : String(e));
    } finally {
      setBusy(null);
    }
  }

  async function handleGtf(file: File) {
    setBusy("Indexing GTF…");
    setError(null);
    try {
      const r = await api.uploadGtf(sessionId, file);
      setAnnLabel(r.label);
      onAnnotation(r.label, r.warnings);
    } catch (e) {
      setError(e instanceof ApiError ? String(e.detail) : String(e));
    } finally {
      setBusy(null);
    }
  }

  return (
    <div className="panel">
      <h2>1 · Inputs</h2>
      <div className="row">
        <label>
          Junction matrix (.rds)
          <input
            type="file"
            accept=".rds"
            onChange={(e) => e.target.files?.[0] && handleRds(e.target.files[0])}
          />
        </label>
        {rdsName && <span className="chip">✓ {rdsName}</span>}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
          <input
            type="checkbox"
            checked={useCustomGtf}
            onChange={(e) => setUseCustomGtf(e.target.checked)}
          />
          Upload a custom GTF instead of a GENCODE release
        </label>
      </div>

      {!useCustomGtf ? (
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
                <option key={s} value={s}>
                  {s}
                </option>
              ))}
            </select>
          </label>
          <label>
            GENCODE release
            <select value={release} onChange={(e) => setRelease(e.target.value)}>
              {(releases[species] ?? []).map((r) => (
                <option key={r} value={r}>
                  {r}
                </option>
              ))}
            </select>
          </label>
          <button onClick={handleGencode} disabled={!release || !!busy}>
            Use this release
          </button>
          {annLabel && <span className="chip">✓ {annLabel}</span>}
        </div>
      ) : (
        <div className="row" style={{ marginTop: 8 }}>
          <label>
            Annotation (.gtf / .gtf.gz)
            <input
              type="file"
              accept=".gtf,.gz"
              onChange={(e) => e.target.files?.[0] && handleGtf(e.target.files[0])}
            />
          </label>
          {annLabel && <span className="chip">✓ {annLabel}</span>}
        </div>
      )}

      {busy && <div className="warn">{busy}</div>}
      {error && <div className="err">{error}</div>}
    </div>
  );
}
