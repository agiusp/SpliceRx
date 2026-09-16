import { useMemo, useRef, useState } from "react";
import type { ProjectionResponse } from "../api";
import { extent, linear } from "../plot/scales";
import { shapePath } from "../plot/shapes";
import { BivariateKey, ColorLegend, ShapeLegend } from "./Legends";

const W = 640;
const H = 520;
const M = { top: 16, right: 16, bottom: 44, left: 52 };

export default function Projection({ data }: { data: ProjectionResponse }) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [hover, setHover] = useState<{ x: number; y: number; lines: string[] } | null>(null);
  // Points the backend drew grey — missing the selected clinical feature(s).
  // Hiding them is a pure display filter: the projection itself (computed
  // from every sample) is unchanged, so no new request is made.
  const [hideMissing, setHideMissing] = useState(false);
  const nMissing = data.points.reduce((n, p) => n + (p.missing ? 1 : 0), 0);
  const points = hideMissing ? data.points.filter((p) => !p.missing) : data.points;

  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);
  const xScale = useMemo(() => linear(extent(xs), [M.left, W - M.right]), [points]);
  const yScale = useMemo(() => linear(extent(ys), [H - M.bottom, M.top]), [points]);

  function download() {
    if (!svgRef.current) return;
    const blob = new Blob([new XMLSerializer().serializeToString(svgRef.current)], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${data.method}.svg`;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  return (
    <div>
      <div className="legend-block">
        {data.legend.color && <ColorLegend legend={data.legend.color} />}
        {data.legend.shape && <ShapeLegend legend={data.legend.shape} />}
        {data.legend.bivariate && <BivariateKey legend={data.legend.bivariate} />}
        {!data.legend.color && !data.legend.bivariate && (
          <span className="muted">no clinical feature selected</span>
        )}
      </div>
      <div className="toolbar">
        <button className="ghost" onClick={download}>
          Download SVG
        </button>
        {nMissing > 0 && (
          <label style={{ display: "flex", flexDirection: "row", alignItems: "center", gap: 6 }}>
            <input type="checkbox" checked={hideMissing} onChange={(e) => setHideMissing(e.target.checked)} />
            hide {nMissing} sample{nMissing === 1 ? "" : "s"} missing this data
          </label>
        )}
        <span className="chip">
          {data.method.toUpperCase()} · {points.length} samples
          {data.explained_variance &&
            ` · PC${data.pc_x} + PC${data.pc_y} capture ${(
              (data.explained_variance[data.pc_x - 1] + data.explained_variance[data.pc_y - 1]) *
              100
            ).toFixed(0)}%`}
        </span>
      </div>

      <div className="plot-card">
        <svg
          ref={svgRef}
          width={W}
          height={H}
          viewBox={`0 0 ${W} ${H}`}
          xmlns="http://www.w3.org/2000/svg"
          style={{ fontFamily: "system-ui, sans-serif", display: "block" }}
          onMouseLeave={() => setHover(null)}
        >
          <line x1={M.left} x2={W - M.right} y1={H - M.bottom} y2={H - M.bottom} stroke="currentColor" strokeOpacity={0.35} />
          <line x1={M.left} x2={M.left} y1={M.top} y2={H - M.bottom} stroke="currentColor" strokeOpacity={0.35} />
          <text x={(W) / 2} y={H - 10} textAnchor="middle" fontSize={11} fill="currentColor" fillOpacity={0.7}>
            {data.axis_labels[0]}
          </text>
          <text
            x={14}
            y={H / 2}
            textAnchor="middle"
            fontSize={11}
            fill="currentColor"
            fillOpacity={0.7}
            transform={`rotate(-90 14 ${H / 2})`}
          >
            {data.axis_labels[1]}
          </text>

          {points.map((p) => {
            const cx = xScale(p.x);
            const cy = yScale(p.y);
            const lines = [
              p.sample,
              ...Object.entries(p.clinical).map(([k, v]) => `${k}: ${v ?? "—"}`),
            ];
            return (
              <path
                key={p.sample}
                d={shapePath(p.shape, 5)}
                transform={`translate(${cx},${cy})`}
                fill={p.color}
                stroke="var(--surface)"
                strokeWidth={1}
                onMouseMove={(e) => setHover({ x: e.clientX, y: e.clientY, lines })}
                onMouseLeave={() => setHover(null)}
              />
            );
          })}
        </svg>
      </div>

      {data.warnings.map((w, i) => (
        <div className="warn" key={i}>
          {w}
        </div>
      ))}

      {hover && (
        <div className="tooltip" style={{ left: hover.x + 12, top: hover.y + 12 }}>
          {hover.lines.map((l, i) => (
            <div key={i}>{l}</div>
          ))}
        </div>
      )}
    </div>
  );
}
