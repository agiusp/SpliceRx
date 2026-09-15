import type { LinearScale } from "./scales";

/**
 * SVG path for an upper half-ellipse arc sitting on a horizontal baseline.
 * Bulges upward (toward smaller y). Endpoint order does not matter.
 */
export function arcPath(
  xScale: LinearScale,
  gStart: number,
  gEnd: number,
  baselineY: number,
  archPx: number,
): string {
  const a = xScale(gStart);
  const b = xScale(gEnd);
  const xL = Math.min(a, b);
  const xR = Math.max(a, b);
  const rx = Math.max((xR - xL) / 2, 0.5);
  const ry = Math.max(archPx, 1);
  // sweep-flag 1 => clockwise on screen (y-down) => arch goes up and over
  return `M ${xL.toFixed(2)} ${baselineY.toFixed(2)} A ${rx.toFixed(2)} ${ry.toFixed(2)} 0 0 1 ${xR.toFixed(
    2,
  )} ${baselineY.toFixed(2)}`;
}

export function arcApex(
  xScale: LinearScale,
  gStart: number,
  gEnd: number,
  baselineY: number,
  archPx: number,
): { x: number; y: number } {
  return {
    x: (xScale(gStart) + xScale(gEnd)) / 2,
    y: baselineY - Math.max(archPx, 1),
  };
}

/** Evenly spaced positions along [x0, x1] for strand chevrons. */
export function chevronXs(x0: number, x1: number, spacing = 34): number[] {
  const lo = Math.min(x0, x1);
  const hi = Math.max(x0, x1);
  const out: number[] = [];
  for (let x = lo + spacing / 2; x < hi; x += spacing) out.push(x);
  return out;
}
