import { useEffect, useState } from "react";
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
  /** Stratification columns in the session (the sample metadata is loaded on
   *  the Data tab). */
  initialColumns?: string[];
  value: Series;
  onChange: (s: Series) => void;
}

export default function SeriesPicker({
  sessionId,
  samples,
  initialColumns,
  value,
  onChange,
}: Props) {
  const [columns, setColumns] = useState<StratColumn[]>([]);
  const [values, setValues] = useState<StratValue[]>([]);
  const [error, setError] = useState<string | null>(null);

  // mirror the stratification columns of the session populated by the Data tab
  useEffect(() => {
    setColumns((initialColumns ?? []).map((name) => ({ name, n_values: 0 })));
  }, [initialColumns]);

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
          <span className="muted">Load a cohort on the Data tab to choose a sample.</span>
        )
      ) : (
        <div className="row">
          {columns.length > 0 ? (
            <label>
              Stratify by
              <select value={value.stratColumn} onChange={(e) => pickColumn(e.target.value)}>
                <option value="">— choose a column —</option>
                {columns.map((c) => (
                  <option key={c.name} value={c.name}>
                    {c.name}
                    {c.n_values ? ` (${c.n_values})` : ""}
                  </option>
                ))}
              </select>
            </label>
          ) : (
            <span className="muted">
              Load a cohort with sample metadata on the Data tab to stratify into groups.
            </span>
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

      {error && <div className="err">{error}</div>}
    </div>
  );
}
