// Printability readouts (a subset of the desktop app's printcheck.py).

import type { TriMesh } from "./kernel";

export const FILAMENT_DENSITY = 1.24; // g/cm3, PLA
export const FILL_FACTOR = 0.3; // printed mass / solid volume
export const PRINT_MIN_PER_GRAM = 300 / 72; // "a 72 g print in 5 hours"
export const MAX_OVERHANG_FRACTION = 0.15;
export const ALUMINIUM_DENSITY = 2.7; // g/cm3, 6061-T6

/**
 * Fraction of the surface that needs support when printed back-down
 * (build direction +x): faces whose normal points within 40 deg of
 * straight down, excluding the mounting face on the bed.
 */
export function overhangFraction(m: TriMesh, maxAngleDeg = 50, bedTol = 0.75): number {
  const { pos: p, idx } = m;
  const cosLim = Math.cos(((90 - maxAngleDeg) * Math.PI) / 180);
  let bad = 0, total = 0;
  for (let i = 0; i < idx.length; i += 3) {
    const a = 3 * idx[i], b = 3 * idx[i + 1], c = 3 * idx[i + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2];
    const vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const A = Math.hypot(nx, ny, nz);
    if (!A) continue;
    total += A;
    const cx = (p[a] + p[b] + p[c]) / 3;
    if (cx >= bedTol && -nx / A > cosLim) bad += A;
  }
  return total ? bad / total : 0;
}

export function printEstimate(volumeMm3: number) {
  const grams = (volumeMm3 / 1000) * FILAMENT_DENSITY * FILL_FACTOR;
  return { grams, minutes: grams * PRINT_MIN_PER_GRAM };
}
