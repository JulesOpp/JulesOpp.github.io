// What a hold is, as far as the page is concerned: a set of sliders and a
// function that turns their values into a body plus the screw voids.
//
// Body and screws are returned separately so the texture can be applied
// to the body first and the screws cut afterwards -- a textured bore or
// a roughened countersink seat would be useless.

import type { Manifold } from "./kernel";

export interface Param {
  key: string;
  label: string;
  lo: number;
  hi: number;
  def: number;
  step?: number; // default 1
  unit?: string; // default " mm"
  tip?: string;
}

export type Values = Record<string, number>;

export interface Built {
  body: Manifold;
  screws: Manifold[];
}

export interface HoldDef {
  id: string;
  name: string;
  group: "Parametric" | "Sculpted" | "Machinable" | "Scanned (CC BY-SA)";
  blurb: string;
  params: Param[];
  build(v: Values, screws: boolean): Built;
  /** Holds with an analytic B-rep twin can be exported as STEP. */
  step?: boolean;
  /** Fixed-mesh holds: credit line shown under the viewer and in downloads. */
  credit?: string;
}

export const P = (key: string, label: string, lo: number, hi: number, def: number, step = 1, unit = " mm", tip?: string): Param =>
  ({ key, label, lo, hi, def, step, unit, tip });

export const deg = (key: string, label: string, lo: number, hi: number, def: number, tip?: string) => P(key, label, lo, hi, def, 1, "°", tip);
export const frac = (key: string, label: string, lo: number, hi: number, def: number, step = 0.01, tip?: string) => P(key, label, lo, hi, def, step, "", tip);

export function defaults(h: HoldDef): Values {
  return Object.fromEntries(h.params.map((p) => [p.key, p.def]));
}

/** Snap a raw slider value to the parameter's grid and range. */
export function quantize(p: Param, v: number): number {
  const step = p.step ?? 1;
  const q = Math.round((v - p.lo) / step) * step + p.lo;
  return +Math.min(p.hi, Math.max(p.lo, q)).toFixed(6);
}
