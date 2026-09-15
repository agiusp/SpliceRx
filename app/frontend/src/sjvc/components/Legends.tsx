import type { BivariateLegend, LegendColor, LegendShape } from "../api";
import { shapePath } from "../plot/shapes";

function Glyph({ shape, color = "var(--text)" }: { shape: string; color?: string }) {
  return (
    <svg width={16} height={16} viewBox="-8 -8 16 16">
      <path d={shapePath(shape, 5)} fill={color} />
    </svg>
  );
}

export function ColorLegend({ legend }: { legend: LegendColor }) {
  if (legend.kind === "categorical") {
    return (
      <div className="grp">
        <span className="ttl">{legend.feature}</span>
        {(legend.items ?? []).map((it) => (
          <span className="item" key={it.value}>
            <span className="sw" style={{ background: it.color }} />
            {it.value}
          </span>
        ))}
      </div>
    );
  }
  const stops = legend.stops ?? [];
  const grad = `linear-gradient(to right, ${stops.join(",")})`;
  return (
    <div className="grp">
      <span className="ttl">{legend.feature}</span>
      <div style={{ width: 140, height: 10, borderRadius: 3, background: grad }} />
      <div className="mono" style={{ display: "flex", justifyContent: "space-between", width: 140, fontSize: 10 }}>
        <span>{legend.min?.toFixed(1)}</span>
        <span>{legend.max?.toFixed(1)}</span>
      </div>
    </div>
  );
}

export function ShapeLegend({ legend }: { legend: LegendShape }) {
  return (
    <div className="grp">
      <span className="ttl">{legend.feature}</span>
      {legend.items.map((it) => (
        <span className="item" key={it.value}>
          <Glyph shape={it.shape} />
          {it.value}
        </span>
      ))}
    </div>
  );
}

function hex2rgb(h: string) {
  const n = h.replace("#", "");
  return [parseInt(n.slice(0, 2), 16), parseInt(n.slice(2, 4), 16), parseInt(n.slice(4, 6), 16)];
}
function blend(corners: BivariateLegend["corners"], f1: number, f2: number) {
  const ll = hex2rgb(corners.lo_lo);
  const hl = hex2rgb(corners.hi_lo);
  const lh = hex2rgb(corners.lo_hi);
  const hh = hex2rgb(corners.hi_hi);
  const top = ll.map((c, i) => c + (hl[i] - c) * f1);
  const bot = lh.map((c, i) => c + (hh[i] - c) * f1);
  const rgb = top.map((c, i) => Math.round(c + (bot[i] - c) * f2));
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

export function BivariateKey({ legend }: { legend: BivariateLegend }) {
  const { corners, features, ranges } = legend;
  const N = 6;
  const cell = 14;
  return (
    <div className="grp">
      <span className="ttl">
        {features[0]} × {features[1]}
      </span>
      <div style={{ display: "flex", gap: 4 }}>
        <span style={{ writingMode: "vertical-rl", transform: "rotate(180deg)", fontSize: 10 }}>
          {features[1]} →
        </span>
        <div>
          <svg width={N * cell} height={N * cell} style={{ border: "1px solid var(--border)", borderRadius: 3 }}>
            {Array.from({ length: N }).map((_, r) =>
              Array.from({ length: N }).map((_, c) => (
                <rect
                  key={`${r}-${c}`}
                  x={c * cell}
                  y={(N - 1 - r) * cell}
                  width={cell}
                  height={cell}
                  fill={blend(corners, c / (N - 1), r / (N - 1))}
                />
              )),
            )}
          </svg>
          <div style={{ fontSize: 10, textAlign: "center" }}>{features[0]} →</div>
        </div>
      </div>
      <div className="mono" style={{ fontSize: 10 }}>
        {features[0]}: {ranges[0][0].toFixed(1)}–{ranges[0][1].toFixed(1)} · {features[1]}: {ranges[1][0].toFixed(1)}–
        {ranges[1][1].toFixed(1)}
      </div>
    </div>
  );
}
