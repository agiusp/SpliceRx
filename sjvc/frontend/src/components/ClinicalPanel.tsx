import type { ClinicalColumn } from "../api";

interface Props {
  columns: ClinicalColumn[];
  selected: string[];
  overrides: Record<string, "numeric" | "categorical">;
  maxFeatures: number;
  onSelected: (s: string[]) => void;
  onOverride: (name: string, type: "numeric" | "categorical") => void;
}

export default function ClinicalPanel({
  columns,
  selected,
  overrides,
  maxFeatures,
  onSelected,
  onOverride,
}: Props) {
  const eligible = columns.filter((c) => c.eligible);
  const available = eligible.filter((c) => !selected.includes(c.name));
  const typeOf = (c: ClinicalColumn) => overrides[c.name] ?? c.type;
  const full = selected.length >= maxFeatures;

  function add(name: string) {
    if (!name) return;
    if (selected.length < maxFeatures) onSelected([...selected, name]);
  }

  return (
    <div>
      <div className="row" style={{ alignItems: "center" }}>
        <label style={{ maxWidth: 320 }}>
          Clinical feature{maxFeatures === 2 ? "s (up to 2)" : "s"}
          <select
            value=""
            disabled={eligible.length === 0 || (full && maxFeatures !== Infinity)}
            onChange={(e) => {
              add(e.target.value);
              e.currentTarget.value = "";
            }}
          >
            <option value="">
              {eligible.length === 0
                ? "no eligible columns"
                : full && maxFeatures !== Infinity
                  ? "remove one to add another"
                  : "+ add feature…"}
            </option>
            {available.map((c) => (
              <option key={c.name} value={c.name}>
                {c.name} ({typeOf(c)})
              </option>
            ))}
          </select>
        </label>
      </div>

      {selected.length > 0 && (
        <div className="pill-list" style={{ marginTop: 8 }}>
          {selected.map((name) => {
            const col = columns.find((c) => c.name === name);
            if (!col) return null;
            const t = typeOf(col);
            return (
              <span className="feat-chip" key={name}>
                {name}
                <span
                  className="type"
                  title="click to switch numeric ⇄ categorical"
                  onClick={() => onOverride(name, t === "numeric" ? "categorical" : "numeric")}
                >
                  {t}
                </span>
                <span className="x" onClick={() => onSelected(selected.filter((s) => s !== name))}>
                  ×
                </span>
              </span>
            );
          })}
        </div>
      )}
    </div>
  );
}
