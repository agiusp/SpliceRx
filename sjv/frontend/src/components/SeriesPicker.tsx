import { useState } from "react";
import { api, ApiError, type StratColumn, type StratValue } from "../api";
import SampleSelect from "./SampleSelect";

export interface Series {
  mode: "sample" | "group";
  sample: string;
  stratColumn: string;
  stratValue: string;
}

interface Props {
  sessionId: string;
  samples: string[];
  value: Series;
  onChange: (s: Series) => void;
  onWarnings: (w: string[]) => void;
}

export default function SeriesPicker({ sessionId, samples, value, onChange, onWarnings }: Props) {
  const [columns, setColumns] = useState<StratColumn[]>([]);
  const [values, setValues] = useState<StratValue[]>([]);
  const [fileName, setFileName] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function handleFile(file: File) {
    setBusy(true);
    setError(null);
    setValues([]);
    try {
      const r = await api.uploadSampleMetadata(sessionId, file);
      setColumns(r.columns);
      setFileName(file.name);
      onWarnings(r.warnings);
      onChange({ ...value, mode: "group", stratColumn: "", stratValue: "" });
    } catch (e) {
      setColumns([]);
      setError(e instanceof ApiError ? String(e.detail) : String(e));
    } finally {
      setBusy(false);
    }
  }

  async function pickColumn(column: string) {
    onChange({ ...value, stratColumn: column, stratValue: "" });
    setValues([]);
    if (!column) return;
    try {
      const vs = await api.stratValues(sessionId, column);
      setValues(vs);
      if (vs.length) onChange({ ...value, stratColumn: column, stratValue: vs[0].value });
    } catch (e) {
      setError(e instanceof ApiError ? String(e.detail) : String(e));
    }
  }

  return (
    <div>
      <div className="radio-row">
        <label>
          <input
            type="radio"
            checked={value.mode === "sample"}
            onChange={() => onChange({ ...value, mode: "sample" })}
          />
          Single sample
        </label>
        <label>
          <input
            type="radio"
            checked={value.mode === "group"}
            onChange={() => onChange({ ...value, mode: "group" })}
          />
          Sample group (per-junction median, excluding NA / 0)
        </label>
      </div>

      {value.mode === "sample" ? (
        samples.length > 0 ? (
          <SampleSelect
            samples={samples}
            value={value.sample}
            onChange={(s) => onChange({ ...value, sample: s })}
          />
        ) : (
          <span className="muted">Upload an RDS to choose a sample.</span>
        )
      ) : (
        <div className="row">
          <label>
            Sample metadata — .csv / .tsv / .rds with a sample_id column
            <input
              type="file"
              accept=".csv,.tsv,.txt,.rds"
              onChange={(e) => e.target.files?.[0] && handleFile(e.target.files[0])}
            />
          </label>
          {fileName && <span className="chip">✓ {fileName}</span>}

          {columns.length > 0 && (
            <label>
              Stratify by
              <select value={value.stratColumn} onChange={(e) => pickColumn(e.target.value)}>
                <option value="">— choose a column —</option>
                {columns.map((c) => (
                  <option key={c.name} value={c.name}>
                    {c.name} ({c.n_values})
                  </option>
                ))}
              </select>
            </label>
          )}

          {values.length > 0 && (
            <label>
              Group
              <select
                value={value.stratValue}
                onChange={(e) => onChange({ ...value, stratValue: e.target.value })}
              >
                {values.map((v) => (
                  <option key={v.value} value={v.value}>
                    {v.value} ({v.n_samples})
                  </option>
                ))}
              </select>
            </label>
          )}
        </div>
      )}

      {busy && <div className="warn">Reading metadata table…</div>}
      {error && <div className="err">{error}</div>}
    </div>
  );
}
