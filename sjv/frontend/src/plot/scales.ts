export interface LinearScale {
  (x: number): number;
  invert(px: number): number;
  domain: [number, number];
  range: [number, number];
}

export function linearScale(domain: [number, number], range: [number, number]): LinearScale {
  const [d0, d1] = domain;
  const [r0, r1] = range;
  const m = (r1 - r0) / (d1 - d0 || 1);
  const f = ((x: number) => r0 + (x - d0) * m) as LinearScale;
  f.invert = (px: number) => d0 + (px - r0) / m;
  f.domain = domain;
  f.range = range;
  return f;
}

/** "Nice" tick positions for a genomic bp axis. */
export function bpTicks(domain: [number, number], target = 6): number[] {
  const span = domain[1] - domain[0];
  const raw = span / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm >= 5 ? 5 : norm >= 2 ? 2 : 1) * mag;
  const first = Math.ceil(domain[0] / step) * step;
  const out: number[] = [];
  for (let t = first; t <= domain[1]; t += step) out.push(Math.round(t));
  return out;
}

export function formatBp(x: number): string {
  return x.toLocaleString("en-US");
}
