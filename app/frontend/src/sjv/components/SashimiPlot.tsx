import { useMemo, useRef, useState } from "react";
import type { ExonModel, PlotResponse, TranscriptModel } from "../api";
import { arcApex, arcPath, chevronXs } from "../plot/arcs";
import { categoryColor, CATEGORY_LABELS } from "../plot/palette";
import { bpTicks, formatBp, linearScale } from "../plot/scales";

const M = { left: 12, right: 12, top: 10 };
const ARC_AREA = 240;
const ARC_STROKE = 2; // uniform line width — read support is shown by arc height
const AXIS_PX = 26;
const ROW_PX = 22;
const CDS_H = 12;
const UTR_H = 6;

// expanded (scroll) mode: on-screen px per kb at zoom 1
const PX_PER_KB = 90;

function mergeUnion(transcripts: TranscriptModel[]): ExonModel[] {
  const segs = transcripts
    .flatMap((t) => t.exons)
    .slice()
    .sort((a, b) => a.start - b.start);
  const out: ExonModel[] = [];
  for (const s of segs) {
    const last = out[out.length - 1];
    if (last && s.start <= last.end + 1) {
      last.end = Math.max(last.end, s.end);
      if (s.kind === "CDS") last.kind = "CDS";
    } else {
      out.push({ ...s });
    }
  }
  return out;
}

export default function SashimiPlot({ plot, width }: { plot: PlotResponse; width: number }) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [collapsed, setCollapsed] = useState(true);
  const [fit, setFit] = useState(true);
  const [zoom, setZoom] = useState(1);
  const [hover, setHover] = useState<{ x: number; y: number; text: string[] } | null>(null);

  const multiGene = plot.genes.length > 1;

  const rows: (TranscriptModel & { gene_name: string })[] = (
    collapsed
      ? plot.genes.map((g) => ({
          transcript_id: `${g.name} (gene model)`,
          strand: g.strand,
          gene_name: g.name,
          exons: mergeUnion(plot.transcripts.filter((t) => t.gene_name === g.name)),
        }))
      : plot.transcripts
  ).filter((r) => r.exons.length > 0);

  const fitW = Math.max(width, 640);
  const bpSpan = Math.max(1, plot.x_domain[1] - plot.x_domain[0]);
  const innerW = fit
    ? fitW
    : Math.min(
        40000,
        Math.max(
          Math.round(fitW * 1.25 * zoom),
          Math.round((bpSpan / 1000) * PX_PER_KB * zoom),
        ),
      );

  const x = useMemo(
    () => linearScale(plot.x_domain, [M.left, innerW - M.right]),
    [plot.x_domain, innerW],
  );

  const baselineY = M.top + ARC_AREA;
  const maxArcPx = ARC_AREA - 26;
  const axisY = baselineY + 6;
  const trackTop = baselineY + AXIS_PX + 8;
  const height = trackTop + rows.length * ROW_PX + 12;

  // Draw order: all annotated (grey) arcs first / underneath, then the colored
  // arcs in decreasing order of read support (count, or group median) so the
  // lower-support arcs and their labels land on top and stay readable.
  const arcs = plot.arcs
    .map((a) => ({ ...a, archPx: (a.height / 1.1) * maxArcPx }))
    .sort((p, q) => {
      const pa = p.category === "annotated" ? 0 : 1;
      const qa = q.category === "annotated" ? 0 : 1;
      return pa - qa || q.count - p.count;
    });

  const topArc = arcs.length ? arcs.reduce((a, b) => (b.archPx > a.archPx ? b : a)) : null;
  const maxCount = topArc ? topArc.count : 0;

  const guideXs = Array.from(
    new Set(plot.arcs.flatMap((a) => [x(a.start), x(a.end)]).map((v) => Math.round(v))),
  );
  const ticks = bpTicks(plot.x_domain);

  const chrom = plot.genes[0]?.chrom ?? "";
  const lo = Math.min(...plot.genes.map((g) => g.start));
  const hi = Math.max(...plot.genes.map((g) => g.end));
  const strands = Array.from(new Set(plot.genes.map((g) => g.strand))).join("/");
  const fileBase = plot.genes.map((g) => g.name).join("-");

  function downloadSvg() {
    if (!svgRef.current) return;
    const blob = new Blob([new XMLSerializer().serializeToString(svgRef.current)], {
      type: "image/svg+xml",
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${fileBase}.svg`;
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <div>
      <div className="toolbar">
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
          <input type="checkbox" checked={collapsed} onChange={(e) => setCollapsed(e.target.checked)} />
          Collapse transcripts
        </label>
        <label style={{ flexDirection: "row", alignItems: "center", gap: 6 }}>
          <input type="checkbox" checked={fit} onChange={(e) => setFit(e.target.checked)} />
          Fit width
        </label>
        {!fit && (
          <span style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
            <button
              className="ghost"
              style={{ padding: "1px 8px" }}
              onClick={() => setZoom((z) => Math.max(0.25, +(z / 1.6).toFixed(3)))}
            >
              −
            </button>
            <span className="chip">{Math.round(innerW)} px</span>
            <button
              className="ghost"
              style={{ padding: "1px 8px" }}
              onClick={() => setZoom((z) => Math.min(24, +(z * 1.6).toFixed(3)))}
            >
              +
            </button>
          </span>
        )}
        <button className="ghost" onClick={downloadSvg}>
          Download SVG
        </button>
        <span className="chip">
          {plot.genes.map((g) => g.name).join(" + ")} · {chrom}:{formatBp(lo)}–{formatBp(hi)} ({strands}) ·{" "}
          {plot.arcs.length} junction{plot.arcs.length === 1 ? "" : "s"}
        </span>
      </div>

      <div className="plot-wrap">
        <svg
          ref={svgRef}
          width={innerW}
          height={height}
          viewBox={`0 0 ${innerW} ${height}`}
          xmlns="http://www.w3.org/2000/svg"
          style={{ display: "block", fontFamily: "system-ui, sans-serif" }}
          onMouseLeave={() => setHover(null)}
        >
          {/* alignment guides */}
          {guideXs.map((gx, i) => (
            <line
              key={`g${i}`}
              x1={gx}
              x2={gx}
              y1={baselineY}
              y2={height - 8}
              stroke="currentColor"
              strokeOpacity={0.12}
              strokeDasharray="2 3"
            />
          ))}

          {/* baseline */}
          <line x1={M.left} x2={innerW - M.right} y1={baselineY} y2={baselineY} stroke="currentColor" strokeOpacity={0.35} />

          {/* arcs */}
          {arcs.map((a) => {
            const color = categoryColor(a.category);
            const apex = arcApex(x, a.start, a.end, baselineY, a.archPx);
            const tip = [
              a.id,
              `${plot.count_kind === "median" ? "median" : "count"} ${a.count.toLocaleString("en-US")}`,
              CATEGORY_LABELS[a.category] ?? a.category,
            ];
            return (
              <g key={a.id}>
                <path
                  d={arcPath(x, a.start, a.end, baselineY, a.archPx)}
                  fill="none"
                  stroke={color}
                  strokeWidth={ARC_STROKE}
                  strokeLinecap="round"
                  opacity={0.9}
                  onMouseMove={(e) => setHover({ x: e.clientX, y: e.clientY, text: tip })}
                  onMouseLeave={() => setHover(null)}
                />
                {topArc && a.id === topArc.id ? (
                  <text
                    x={apex.x}
                    y={Math.max(apex.y - 4, 9)}
                    textAnchor="middle"
                    fontSize={10}
                    fontWeight={600}
                    fill="currentColor"
                  >
                    max {maxCount.toLocaleString("en-US")}
                    {plot.count_kind === "median" ? " median reads" : " reads"}
                  </text>
                ) : (
                  a.category !== "annotated" && (
                    <text x={apex.x} y={apex.y - 3} textAnchor="middle" fontSize={9} fill={color}>
                      {a.count.toLocaleString("en-US")}
                    </text>
                  )
                )}
              </g>
            );
          })}

          {/* axis */}
          <line x1={M.left} x2={innerW - M.right} y1={axisY} y2={axisY} stroke="currentColor" strokeOpacity={0.3} />
          {ticks.map((t) => (
            <g key={t} transform={`translate(${x(t)},${axisY})`}>
              <line y1={0} y2={4} stroke="currentColor" strokeOpacity={0.4} />
              <text y={15} textAnchor="middle" fontSize={10} fill="currentColor" fillOpacity={0.65}>
                {formatBp(t)}
              </text>
            </g>
          ))}

          {/* annotation track */}
          {rows.map((tx, ri) => {
            const midY = trackTop + ri * ROW_PX + ROW_PX / 2;
            const xs = tx.exons.map((e) => [x(e.start), x(e.end)]).flat();
            const left = Math.min(...xs);
            const right = Math.max(...xs);
            const chevGlyph = tx.strand === "-" ? "‹" : "›";
            const label =
              multiGene && !collapsed ? `${tx.gene_name} · ${tx.transcript_id}` : tx.transcript_id;
            return (
              <g key={`${tx.gene_name}:${tx.transcript_id}`}>
                <line x1={left} x2={right} y1={midY} y2={midY} stroke="var(--track-blue)" strokeWidth={1} />
                {chevronXs(left, right).map((cx, ci) => (
                  <text
                    key={ci}
                    x={cx}
                    y={midY + 3}
                    textAnchor="middle"
                    fontSize={9}
                    fill="var(--track-blue)"
                    fillOpacity={0.7}
                  >
                    {chevGlyph}
                  </text>
                ))}
                {tx.exons.map((e, ei) => {
                  const h = e.kind === "CDS" ? CDS_H : UTR_H;
                  const w = Math.max(x(e.end) - x(e.start), 1);
                  return (
                    <rect
                      key={ei}
                      x={x(e.start)}
                      y={midY - h / 2}
                      width={w}
                      height={h}
                      fill={e.kind === "CDS" ? "var(--track-blue)" : "var(--track-blue-utr)"}
                    />
                  );
                })}
                <text x={Math.max(M.left, left)} y={midY - CDS_H / 2 - 3} fontSize={9} fill="currentColor" fillOpacity={0.55}>
                  {label}
                </text>
              </g>
            );
          })}
        </svg>
      </div>

      {hover && (
        <div className="tooltip" style={{ left: hover.x + 12, top: hover.y + 12 }}>
          {hover.text.map((line, i) => (
            <div key={i}>{line}</div>
          ))}
        </div>
      )}
    </div>
  );
}
