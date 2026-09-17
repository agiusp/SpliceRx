import { useEffect, useMemo, useRef, useState } from "react";
import type { HeatmapResponse } from "../api";
import { diverging, sequential } from "../plot/scales";

const LABEL_W = 130;
const ANNO_H = 14;
const DENDRO = 46;
const MAX_LABEL_ROWS = 60;

export default function Heatmap({ data }: { data: HeatmapResponse }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const overlayRef = useRef<SVGSVGElement>(null);
  const [hover, setHover] = useState<{ x: number; y: number; lines: string[] } | null>(null);
  const [compact, setCompact] = useState(false);
  const [msg, setMsg] = useState<string | null>(null);
  // "junction" = the raw chr:start-end:strand row id (always available);
  // "gene_type" = "<gene name>:<novel-splicing-event type>" from the
  // cohort's junction metadata table — only offered when the backend sent
  // it. Toggling is purely a label swap: row order/values are unaffected,
  // so no new request is made.
  const [labelMode, setLabelMode] = useState<"junction" | "gene_type">("junction");
  const hasGeneTypeLabels = !!data.row_labels_gene_type;
  const rowLabels =
    labelMode === "gene_type" && data.row_labels_gene_type
      ? data.row_labels_gene_type.map((lab, i) => lab ?? data.row_labels[i])
      : data.row_labels;

  const nRows = data.values.length;
  const nCols = data.samples.length;

  // Fit the cell grid inside a slide-friendly box: a big cohort must not blow
  // the width out to thousands of px (which pastes into slides distorted).
  // Cells may go sub-pixel for a very large matrix; `compact` shrinks further
  // and drops labels. Non-destructive — just a view mode.
  const box = compact ? { w: 560, h: 300 } : { w: 760, h: 520 };
  const cellW = Math.min(compact ? 12 : 30, box.w / Math.max(nCols, 1));
  const cellH = Math.min(compact ? 7 : 20, box.h / Math.max(nRows, 1));
  // integer cell boundaries so cells tile seamlessly even at a fractional pitch
  const colX = (i: number) => Math.round(i * cellW);
  const rowY = (i: number) => Math.round(i * cellH);

  const rowDendW = data.row_dendro ? DENDRO : 0;
  const colDendH = data.col_dendro ? DENDRO : 0;
  const annoH = data.annotations.length * (ANNO_H + 2);
  const showRowLabels = !compact && nRows <= MAX_LABEL_ROWS && cellH >= 7;
  const labelW = showRowLabels ? LABEL_W : data.annotations.length ? 84 : 10;
  // Feature labels sit immediately left of the grid; the feature dendrogram
  // (if any) is pushed to the grid's right instead of the more conventional
  // left, so the labels stay flush against the axis they annotate and read
  // top-to-bottom without a dendrogram gap in between. Sample ids are never
  // drawn under the grid — a wide cohort made them illegible anyway, and the
  // per-cell hover tooltip already names the sample.
  const gridLeft = labelW;
  const gridTop = colDendH + annoH + 4;
  const gridW = colX(nCols);
  const gridH = rowY(nRows);
  const rowDendLeft = gridLeft + gridW;

  const [vmin, vmax] = useMemo(() => {
    let lo = Infinity;
    let hi = -Infinity;
    for (const row of data.values) for (const v of row) {
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    return [lo, hi];
  }, [data]);

  const cellColor = (v: number) =>
    data.row_zscore ? diverging(v) : sequential((v - vmin) / (vmax - vmin || 1), ["#f2f5fb", "#9ec5f4", "#3987e5", "#184f95"]);

  useEffect(() => {
    const cv = canvasRef.current;
    if (!cv) return;
    const dpr = window.devicePixelRatio || 1;
    cv.width = gridW * dpr;
    cv.height = gridH * dpr;
    cv.style.width = `${gridW}px`;
    cv.style.height = `${gridH}px`;
    const ctx = cv.getContext("2d")!;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    for (let r = 0; r < nRows; r++) {
      const y0 = rowY(r);
      const h = rowY(r + 1) - y0;
      for (let c = 0; c < nCols; c++) {
        const x0 = colX(c);
        ctx.fillStyle = cellColor(data.values[r][c]);
        ctx.fillRect(x0, y0, colX(c + 1) - x0, h);
      }
    }
  }, [data, cellW, cellH]); // eslint-disable-line

  function dendroLines(
    segments: number[][][],
    kind: "col" | "row",
  ): string[] {
    let maxD = 0;
    for (const [, ys] of segments) for (const y of ys) if (y > maxD) maxD = y;
    maxD = maxD || 1;
    return segments.map(([xs, ys]) => {
      const pts = xs.map((ic, i) => {
        const along = (ic / 10) * (kind === "col" ? cellW : cellH);
        const depth = (ys[i] / maxD) * (DENDRO - 4);
        return kind === "col"
          ? `${gridLeft + along},${colDendH - depth}`
          : `${rowDendLeft + depth},${gridTop + along}`;
      });
      return `M${pts.join("L")}`;
    });
  }

  const totalW = gridLeft + gridW + rowDendW + 8;
  const totalH = gridTop + gridH + 8;

  // --- export ------------------------------------------------------------- //
  function themeColors() {
    const ink = getComputedStyle(overlayRef.current ?? document.body).color || "#111";
    const bg =
      getComputedStyle(document.documentElement).getPropertyValue("--bg").trim() ||
      getComputedStyle(document.body).backgroundColor ||
      "#ffffff";
    return { ink, bg };
  }

  /** One self-contained SVG (string + dimensions): the cell grid as an
   *  embedded raster under the dendrogram / label / annotation overlay,
   *  with the annotation legend drawn in below. */
  function exportSvg(): { markup: string; w: number; h: number } | null {
    const overlay = overlayRef.current;
    const canvas = canvasRef.current;
    if (!overlay || !canvas) return null;
    const { ink, bg } = themeColors();
    const NS = "http://www.w3.org/2000/svg";

    const svg = overlay.cloneNode(true) as SVGSVGElement;

    const img = document.createElementNS(NS, "image");
    img.setAttribute("x", String(gridLeft));
    img.setAttribute("y", String(gridTop));
    img.setAttribute("width", String(gridW));
    img.setAttribute("height", String(gridH));
    img.setAttribute("href", canvas.toDataURL("image/png"));
    svg.insertBefore(img, svg.firstChild);

    const defs = document.createElementNS(NS, "defs");
    svg.insertBefore(defs, svg.firstChild);

    const add = (tag: string, attrs: Record<string, string | number>, text?: string) => {
      const el = document.createElementNS(NS, tag);
      for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
      if (text != null) el.textContent = text;
      svg.appendChild(el);
      return el;
    };

    // legend, below the plot
    const x0 = 8;
    let y = totalH + 4;
    data.annotations.forEach((a, ai) => {
      add("text", { x: x0, y: y + 11, "font-size": 12, "font-weight": 600, fill: "currentColor" }, a.feature);
      y += 18;
      if (a.type === "categorical" && a.legend) {
        for (const it of a.legend) {
          add("rect", { x: x0, y, width: 12, height: 12, fill: it.color, stroke: ink, "stroke-opacity": 0.2 });
          add("text", { x: x0 + 17, y: y + 10, "font-size": 11, fill: "currentColor" }, String(it.value));
          y += 16;
        }
      } else {
        const stops = a.stops ?? [];
        const gid = `hmgrad${ai}`;
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
        add("text", { x: x0, y: y + 22, "font-size": 10, fill: "currentColor" }, (a.min ?? 0).toFixed(1));
        add("text", { x: x0 + 150, y: y + 22, "font-size": 10, "text-anchor": "end", fill: "currentColor" }, (a.max ?? 0).toFixed(1));
        y += 26;
      }
      y += 10;
    });

    const w = Math.max(totalW, x0 + 270);
    const h = data.annotations.length ? y + 4 : totalH;

    svg.setAttribute("xmlns", NS);
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
    svg.setAttribute("width", String(w));
    svg.setAttribute("height", String(h));
    svg.setAttribute("style", `font-family: system-ui, sans-serif; color: ${ink}`);

    const bgRect = document.createElementNS(NS, "rect");
    bgRect.setAttribute("width", String(w));
    bgRect.setAttribute("height", String(h));
    bgRect.setAttribute("fill", bg);
    svg.insertBefore(bgRect, svg.firstChild);

    return { markup: new XMLSerializer().serializeToString(svg), w, h };
  }

  function flash(m: string) {
    setMsg(m);
    window.setTimeout(() => setMsg((cur) => (cur === m ? null : cur)), 2500);
  }

  function saveBlob(blob: Blob, name: string) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = name;
    a.click();
    URL.revokeObjectURL(url);
  }

  function downloadSvg() {
    const e = exportSvg();
    if (e) saveBlob(new Blob([e.markup], { type: "image/svg+xml" }), "heatmap.svg");
  }

  async function rasterize(scale = 2): Promise<Blob | null> {
    const e = exportSvg();
    if (!e) return null;
    const img = new Image();
    await new Promise<void>((resolve, reject) => {
      img.onload = () => resolve();
      img.onerror = () => reject(new Error("render failed"));
      img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(e.markup);
    });
    const cv = document.createElement("canvas");
    cv.width = Math.ceil(e.w * scale);
    cv.height = Math.ceil(e.h * scale);
    const ctx = cv.getContext("2d")!;
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0);
    return await new Promise<Blob | null>((resolve) => cv.toBlob(resolve, "image/png"));
  }

  async function downloadPng() {
    try {
      const blob = await rasterize(2);
      if (blob) saveBlob(blob, "heatmap.png");
    } catch {
      flash("PNG export failed");
    }
  }

  async function copyImage() {
    try {
      const blob = await rasterize(2);
      if (!blob) return;
      await navigator.clipboard.write([new ClipboardItem({ "image/png": blob })]);
      flash("copied to clipboard");
    } catch {
      flash("copy failed — use Download PNG");
    }
  }

  return (
    <div>
      <div className="toolbar">
        <button className="ghost" onClick={copyImage}>
          Copy image
        </button>
        <button className="ghost" onClick={downloadPng}>
          Download PNG
        </button>
        <button className="ghost" onClick={downloadSvg}>
          Download SVG
        </button>
        <label style={{ display: "flex", flexDirection: "row", alignItems: "center", gap: 6 }}>
          <input type="checkbox" checked={compact} onChange={(e) => setCompact(e.target.checked)} />
          compact
        </label>
        {hasGeneTypeLabels && (
          <label style={{ display: "flex", flexDirection: "row", alignItems: "center", gap: 6 }}>
            row labels
            <select value={labelMode} onChange={(e) => setLabelMode(e.target.value as typeof labelMode)}>
              <option value="junction">junction id</option>
              <option value="gene_type">gene : splice type</option>
            </select>
          </label>
        )}
        <span className="chip">
          {data.feature_kind === "gene" ? "gene signatures" : "junctions"} · {nRows} × {nCols}
          {data.row_zscore ? " · row z-score" : ""}
          {compact ? " · compact" : ""}
        </span>
        {msg && <span className="muted">{msg}</span>}
      </div>
      {data.warnings.map((w, i) => (
        <div className="warn" key={i}>
          {w}
        </div>
      ))}

      <div className="plot-card">
        <div style={{ position: "relative", width: totalW, height: totalH }}>
          {/* dendrograms + annotation bars as SVG overlay */}
          <svg
            ref={overlayRef}
            width={totalW}
            height={totalH}
            style={{ position: "absolute", inset: 0, fontFamily: "system-ui, sans-serif", pointerEvents: "none" }}
          >
            {data.col_dendro &&
              dendroLines(data.col_dendro, "col").map((d, i) => (
                <path key={i} d={d} fill="none" stroke="currentColor" strokeOpacity={0.4} />
              ))}
            {data.row_dendro &&
              dendroLines(data.row_dendro, "row").map((d, i) => (
                <path key={i} d={d} fill="none" stroke="currentColor" strokeOpacity={0.4} />
              ))}

            {/* clinical annotation bars */}
            {data.annotations.map((a, ai) => {
              const y = colDendH + ai * (ANNO_H + 2);
              return (
                <g key={a.feature}>
                  {a.colors.map((c, ci) => (
                    <rect key={ci} x={gridLeft + colX(ci)} y={y} width={colX(ci + 1) - colX(ci)} height={ANNO_H} fill={c} />
                  ))}
                  <text x={gridLeft - 6} y={y + ANNO_H - 3} textAnchor="end" fontSize={10} fill="currentColor" fillOpacity={0.75}>
                    {a.feature}
                  </text>
                </g>
              );
            })}

            {/* row labels */}
            {showRowLabels &&
              rowLabels.map((lab, ri) => (
                <text
                  key={ri}
                  x={gridLeft - 6}
                  y={gridTop + ri * cellH + cellH / 2 + 3}
                  textAnchor="end"
                  fontSize={Math.min(11, cellH)}
                  fill="currentColor"
                  fillOpacity={0.8}
                >
                  {lab}
                </text>
              ))}
          </svg>

          <canvas
            ref={canvasRef}
            style={{ position: "absolute", left: gridLeft, top: gridTop, imageRendering: "pixelated" }}
            onMouseMove={(e) => {
              const rect = (e.target as HTMLCanvasElement).getBoundingClientRect();
              const c = Math.min(nCols - 1, Math.floor(((e.clientX - rect.left) / rect.width) * nCols));
              const r = Math.min(nRows - 1, Math.floor(((e.clientY - rect.top) / rect.height) * nRows));
              if (r < 0 || r >= nRows || c < 0 || c >= nCols) return setHover(null);
              const annoLines = data.annotations.map((a) => `${a.feature}: ${a.values[c] ?? "—"}`);
              setHover({
                x: e.clientX,
                y: e.clientY,
                lines: [`${rowLabels[r]} · ${data.samples[c]}`, `value ${data.values[r][c].toFixed(2)}`, ...annoLines],
              });
            }}
            onMouseLeave={() => setHover(null)}
          />
        </div>
      </div>

      {/* annotation legends */}
      <div className="legend-block" style={{ marginTop: 8 }}>
        {data.annotations.map((a) =>
          a.type === "categorical" && a.legend ? (
            <div className="grp" key={a.feature}>
              <span className="ttl">{a.feature}</span>
              {a.legend.map((it) => (
                <span className="item" key={it.value}>
                  <span className="sw" style={{ background: it.color }} />
                  {it.value}
                </span>
              ))}
            </div>
          ) : (
            <div className="grp" key={a.feature}>
              <span className="ttl">{a.feature}</span>
              <div
                style={{
                  width: 120,
                  height: 10,
                  borderRadius: 3,
                  background: `linear-gradient(to right, ${(a.stops ?? []).join(",")})`,
                }}
              />
              <div className="mono" style={{ display: "flex", justifyContent: "space-between", width: 120, fontSize: 10 }}>
                <span>{a.min?.toFixed(1)}</span>
                <span>{a.max?.toFixed(1)}</span>
              </div>
            </div>
          ),
        )}
      </div>

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
