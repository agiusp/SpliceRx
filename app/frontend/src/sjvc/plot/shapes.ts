// SVG path for a point glyph centred at (0,0), roughly area-matched. r ~ radius.
export function shapePath(shape: string, r: number): string {
  const s = r;
  switch (shape) {
    case "square":
      return `M${-s},${-s}h${2 * s}v${2 * s}h${-2 * s}z`;
    case "triangle-up":
      return `M0,${-s * 1.2}L${s * 1.1},${s * 0.9}L${-s * 1.1},${s * 0.9}z`;
    case "triangle-down":
      return `M0,${s * 1.2}L${s * 1.1},${-s * 0.9}L${-s * 1.1},${-s * 0.9}z`;
    case "diamond":
      return `M0,${-s * 1.3}L${s * 1.3},0L0,${s * 1.3}L${-s * 1.3},0z`;
    case "plus": {
      const t = s * 0.42;
      return `M${-t},${-s}h${2 * t}v${s - t}h${s - t}v${2 * t}h${-(s - t)}v${s - t}h${-2 * t}v${-(s - t)}h${-(s - t)}v${-2 * t}h${s - t}z`;
    }
    case "cross": {
      const t = s * 0.42;
      const a = s;
      return `M${-a},${-a + t}L${-a + t},${-a}L0,${-t}L${a - t},${-a}L${a},${-a + t}L${t},0L${a},${a - t}L${a - t},${a}L0,${t}L${-a + t},${a}L${-a},${a - t}L${-t},0z`;
    }
    case "star": {
      let d = "";
      for (let i = 0; i < 10; i++) {
        const rr = i % 2 === 0 ? s * 1.35 : s * 0.55;
        const ang = (Math.PI / 5) * i - Math.PI / 2;
        d += (i === 0 ? "M" : "L") + `${(rr * Math.cos(ang)).toFixed(2)},${(rr * Math.sin(ang)).toFixed(2)}`;
      }
      return d + "z";
    }
    case "dot-small":
      return `M0,0m${-s * 0.5},0a${s * 0.5},${s * 0.5} 0 1,0 ${s},0a${s * 0.5},${s * 0.5} 0 1,0 ${-s},0`;
    default: // circle
      return `M0,0m${-s},0a${s},${s} 0 1,0 ${2 * s},0a${s},${s} 0 1,0 ${-2 * s},0`;
  }
}
