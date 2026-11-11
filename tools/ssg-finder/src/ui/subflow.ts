// Subgroup page logic without the DOM: transformation text, the "Send to finder" record
// and its sessionStorage hand-over. subgroups.ts (subgroup page) and main.ts (finder page)
// both use it; the Node tests call the same functions.
//
// Licence rule (D1): nothing here contains BYU data. The record handed to the finder is
// computed from the viewer's own paste and lives only in their browser's sessionStorage.

import { Q } from "../rational";
import { buildOps, type GroupRecord } from "../ssg";
import type { SubRow } from "../subfinder";
import type { Standardised } from "../standardise";

/** sessionStorage key of a subgroup record on its way from the subgroup page to the finder. */
export const SEND_KEY = "ssgsub:send";

/** Linear expression Σ c_i v_i + c0 in the usual xyz notation ("x-y", "-z+1/4", "0"). */
function linText(coef: Q[], vars: string[], c0: Q): string {
  let s = "";
  coef.forEach((c, i) => {
    if (c.isZero()) return;
    const neg = c.n < 0n;
    const a = neg ? c.neg() : c;
    const mag = a.n === 1n && a.d === 1n ? "" : a.toString();
    s += (neg ? "-" : s ? "+" : "") + mag + vars[i];
  });
  if (!c0.isZero()) s += (c0.n < 0n ? "-" : s ? "+" : "") + (c0.n < 0n ? c0.neg() : c0).toString();
  return s || "0";
}

/**
 * The map from the parent's axes to the type's setting, as new coordinates in terms of the
 * old: x_std = P x + p and t' = m·x + η t + p4, written "(-y, x-y, -z-1/4; -t)".
 */
export function affineText(t: Pick<Standardised, "P" | "p" | "eta" | "m" | "p4">): string {
  const xyz = ["x", "y", "z"];
  const ext = t.P.map((row, i) => linText(row.map((x) => Q.parse(x)), xyz, Q.parse(t.p[i])));
  const int = linText([...t.m.map((x) => Q.of(x)), Q.of(t.eta)], [...xyz, "t"], Q.parse(t.p4));
  return `(${ext.join(", ")}; ${int})`;
}

/** The subgroup as a finder record: the parent's lattice, q and centrings, the subgroup's operators (parent axes). */
export function subgroupRecord(parent: GroupRecord, row: SubRow, parentSgNumber: number): GroupRecord {
  const name = row.inParentAxes ?? (row.type ? row.type.symbol.full : "subgroup");
  const rec: GroupRecord = {
    sym: `${name} < ${parent.sym}`,
    cls: parent.cls,
    d: parent.d,
    q: [...parent.q],
    params: [...parent.params],
    cen: [...parent.cen],
    gen: [...row.generators],
  };
  // The operators are written in the PARENT's axes, so the finder must pick the crystal system
  // of those axes (the parent's space-group number), not that of the subgroup's own type.
  rec.sgNumber = parentSgNumber;
  rec.nops = buildOps(rec).length;
  return rec;
}

/** Parse and check a record from sessionStorage; null when it is missing or malformed. */
export function decodeSent(text: string | null): GroupRecord | null {
  if (!text) return null;
  try {
    const o = JSON.parse(text);
    const strs = (a: unknown) => Array.isArray(a) && a.every((x) => typeof x === "string");
    if (!o || typeof o.sym !== "string" || typeof o.cls !== "string" || !Number.isInteger(o.d)) return null;
    if (!strs(o.q) || !strs(o.params) || !strs(o.cen) || !strs(o.gen) || o.q.length !== o.d) return null;
    const rec: GroupRecord = { sym: o.sym, cls: o.cls, d: o.d, q: o.q, params: o.params, cen: o.cen, gen: o.gen };
    if (Number.isInteger(o.sgNumber) && o.sgNumber >= 1 && o.sgNumber <= 230) rec.sgNumber = o.sgNumber;
    rec.nops = buildOps(rec).length; // throws on malformed operators
    return rec;
  } catch {
    return null;
  }
}
