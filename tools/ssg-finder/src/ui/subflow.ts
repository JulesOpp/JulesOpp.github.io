// Subgroup page logic without the DOM: transformation text, the "Send to finder" record
// and its sessionStorage hand-over. subgroups.ts (subgroup page) and main.ts (finder page)
// both use it; the Node tests call the same functions.
//
// Licence rule (D1): nothing here contains BYU data. The record handed to the finder is
// computed from the viewer's own paste and lives only in their browser's sessionStorage.

import { Q } from "../rational";
import { buildOps, type GroupRecord } from "../ssg";
import type { AffineMap, SubRow } from "../subfinder";

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
export function affineText(t: AffineMap): string {
  const xyz = ["x", "y", "z"];
  const ext = t.P.map((row, i) => linText(row.map((x) => Q.parse(x)), xyz, Q.parse(t.p[i])));
  const int = t.M.map((row, i) => linText([...row.map(Q.parse), ...t.E[i].map(x=>typeof x === "string"?Q.parse(x):Q.of(x))], [...xyz, ...["t", "u", "v"].slice(0, t.E.length)], Q.parse(t.p4[i]))).join(", ");
  return `(${ext.join(", ")}; ${int})`;
}

/** Display an equivalent map with origin shifts in [0, 1): integer shifts lie in Σ_B = Z^n. */
export function displayMap(t: AffineMap): typeof t {
  return { ...t, p: t.p.map((x) => Q.parse(x).mod1().toString()), p4: t.p4.map((x) => Q.parse(x).mod1().toString()) };
}

/** Build the hand-over before attempting browser storage, so build errors remain distinct. */
export function prepareSend(parent: GroupRecord, row: SubRow, sg: number): { ok: true; json: string } | { ok: false; reason: string } {
  try {
    return { ok: true, json: JSON.stringify(subgroupRecord(parent, row, sg)) };
  } catch (e) {
    return { ok: false, reason: e instanceof Error ? e.message : String(e) };
  }
}

/** Finder record: parent axes at k=1, primitive subgroup basis at k>1. */
export function subgroupRecord(parent: GroupRecord, row: SubRow, parentSgNumber: number): GroupRecord {
  if (row.rebased) {
    const rec = {...row.rebased, sym:`${row.type?.symbol.full ?? "unidentified subgroup"} < ${parent.sym}`};
    // Primitive rebased axes may differ from the parent's crystal system.
    // When naming failed, let the finder infer the axes from the operations.
    if(row.type) rec.sgNumber=row.type.sgNumber;
    return rec;
  }
  if(row.kIndex>1) throw new Error("The subgroup lattice could not be rebased; finder handoff is unavailable.");
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
