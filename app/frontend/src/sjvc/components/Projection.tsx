import { useMemo, useRef, useState } from "react";
import type { ProjectionResponse } from "../api";
import { extent, linear } from "../plot/scales";
import { shapePath } from "../plot/shapes";
import { BivariateKey, ColorLegend, ShapeLegend, blend } from "./Legends";

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

  // The colour/shape/bivariate key renders as plain HTML above the plot (see
  // <ColorLegend>/<ShapeLegend>/<BivariateKey> below) so it looks right on
  // screen, but that means it isn't part of `svgRef`'s DOM tree and a plain
  // serialize-and-save would silently drop it from the downloaded file.
  // Redraw the same legend as SVG onto a clone, stacked below the plot, from
  // `data.legend` directly — mirrors Heatmap.tsx's `exportSvg`'s approach of
  // drawing its own legend from the response data rather than cloning DOM.
  function download() {
    const src = svgRef.current;
    if (!src) return;
    const NS = "http://www.w3.org/2000/svg";
    const ink = getComputedStyle(src).color || "#111";
    const bg =
      getComputedStyle(document.documentElement).getPropertyValue("--bg").trim() ||
      getComputedStyle(document.body).backgroundColor ||
      "#ffffff";

    const svg = src.cloneNode(true) as SVGSVGElement;
    const defs = document.createElementNS(NS, "defs");
    svg.insertBefore(defs, svg.firstChild);

    const add = (tag: string, attrs: Record<string, string | number>, text?: string) => {
      const el = document.createElementNS(NS, tag);
      for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
      if (text != null) el.textContent = text;
      svg.appendChild(el);
      return el;
    };

    const x0 = 12;
    let y = H + 14;

    function title(t: string) {
      add("text", { x: x0, y: y + 11, "font-size": 12, "font-weight": 600, fill: "currentColor" }, t);
      y += 18;
    }

    const lc = data.legend.color;
    if (lc) {
      title(lc.feature);
      if (lc.kind === "categorical") {
        for (const it of lc.items ?? []) {
          add("rect", { x: x0, y, width: 12, height: 12, fill: it.color, stroke: ink, "stroke-opacity": 0.2 });
          add("text", { x: x0 + 17, y: y + 10, "font-size": 11, fill: "currentColor" }, it.value);
          y += 16;
        }
      } else {
        const stops = lc.stops ?? [];
        const gid = "projlegendgrad";
        const grad = document.createElementNS(NS, "linearGradient");
        grad.setAttribute("id", gid);
        stops.forEach((c, k) => {
          const stop = document.createElementNS(NS, "stop");
          stop.setAttribute("offset", `${(k / Math.max(stops.length - 1, 1)) * 100}%`);
          stop.setAttribute("stop-color", c);
          grad.appendChild(stop);
        });
        defs.appendChild(grad);
        add("rect", { x: x0, y, width: 150, height: 10, fill: `url(#${gid})`, stroke: ink, "stroke-opacity": 0.2 });
        add("text", { x: x0, y: y + 22, "font-size": 10, fill: "currentColor" }, (lc.min ?? 0).toFixed(1));
        add(
          "text", { x: x0 + 150, y: y + 22, "font-size": 10, "text-anchor": "end", fill: "currentColor" },
          (lc.max ?? 0).toFixed(1),
        );
        y += 26;
      }
      y += 10;
    }

    const ls = data.legend.shape;
    if (ls) {
      title(ls.feature);
      for (const it of ls.items) {
        const g = document.createElementNS(NS, "g");
        g.setAttribute("transform", `translate(${x0 + 6},${y + 6})`);
        const path = document.createElementNS(NS, "path");
        path.setAttribute("d", shapePath(it.shape, 5));
        path.setAttribute("fill", "currentColor");
        g.appendChild(path);
        svg.appendChild(g);
        add("text", { x: x0 + 17, y: y + 10, "font-size": 11, fill: "currentColor" }, it.value);
        y += 16;
      }
      y += 10;
    }

    const bv = data.legend.bivariate;
    if (bv) {
      title(`${bv.features[0]} × ${bv.features[1]}`);
      const N = 6;
      const cell = 12;
      for (let r = 0; r < N; r++) {
        for (let c = 0; c < N; c++) {
          add("rect", {
            x: x0 + c * cell, y: y + (N - 1 - r) * cell, width: cell, height: cell,
            fill: blend(bv.corners, c / (N - 1), r / (N - 1)),
          });
        }
      }
      y += N * cell + 14;
      add(
        "text", { x: x0, y, "font-size": 10, fill: "currentColor" },
        `${bv.features[0]}: ${bv.ranges[0][0].toFixed(1)}–${bv.ranges[0][1].toFixed(1)} · ` +
          `${bv.features[1]}: ${bv.ranges[1][0].toFixed(1)}–${bv.ranges[1][1].toFixed(1)}`,
      );
      y += 16;
    }

    const w = Math.max(W, x0 + 300);
    const h = lc || ls || bv ? y + 4 : H;
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
    svg.setAttribute("width", String(w));
    svg.setAttribute("height", String(h));
    svg.setAttribute("style", `font-family: system-ui, sans-serif; color: ${ink}`);

    const bgRect = document.createElementNS(NS, "rect");
    bgRect.setAttribute("width", String(w));
    bgRect.setAttribute("height", String(h));
    bgRect.setAttribute("fill", bg);
    svg.insertBefore(bgRect, svg.firstChild);

    const blob = new Blob([new XMLSerializer().serializeToString(svg)], { type: "image/svg+xml" });
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
