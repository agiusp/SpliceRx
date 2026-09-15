export function linear(domain: [number, number], range: [number, number]) {
  const [d0, d1] = domain;
  const [r0, r1] = range;
  const m = (r1 - r0) / (d1 - d0 || 1);
  return (x: number) => r0 + (x - d0) * m;
}

export function extent(xs: number[]): [number, number] {
  let lo = Infinity;
  let hi = -Infinity;
  for (const x of xs) {
    if (Number.isFinite(x)) {
      if (x < lo) lo = x;
      if (x > hi) hi = x;
    }
  }
  if (!Number.isFinite(lo)) return [0, 1];
  if (lo === hi) return [lo - 1, hi + 1];
  const pad = (hi - lo) * 0.06;
  return [lo - pad, hi + pad];
}

/** diverging blue-white-red for z-scored heatmap values, t in [-1,1] after /clip */
export function diverging(v: number, clip = 3): string {
  const t = Math.max(-1, Math.min(1, v / clip));
  const lerp = (a: number[], b: number[], f: number) => a.map((x, i) => Math.round(x + (b[i] - x) * f));
  const rgb =
    t < 0
      ? lerp([32, 92, 171], [245, 245, 242], 1 + t) // blue -> white
      : lerp([245, 245, 242], [208, 59, 59], t); // white -> red
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

/** sequential (blue ramp) for non-z-scored heatmap values, t in [0,1] */
export function sequential(t: number, stops: string[]): string {
  const p = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  return stops[Math.round(p)];
}
