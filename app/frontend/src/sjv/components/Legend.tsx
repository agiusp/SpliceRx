import type { LegendEntry } from "../api";

export default function Legend({ entries }: { entries: LegendEntry[] }) {
  if (!entries.length) return null;
  return (
    <div className="legend">
      {entries.map((e) => (
        <span className="item" key={e.category}>
          <span className="swatch" style={{ background: e.color }} />
          {e.label}
        </span>
      ))}
    </div>
  );
}
