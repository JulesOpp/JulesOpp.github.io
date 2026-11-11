// Superspace groups: operator parsing/formatting, canonical operator list (closure),
// label q-vectors and parameters. Mirrors reference/exact.py (CONTRACTS.md §1–2).

import { Q, QMat, QVec } from "./rational";

export const VARS = "xyztuv";

/** Group record, CONTRACTS.md §2. */
export interface GroupRecord {
  label?: string;
  sym: string;
  cls: string;
  d: number;
  q: string[];
  params: string[];
  cen: string[];
  gen: string[];
  nops?: number;
  /** Basic space-group number, for records without a label (e.g. a subgroup sent from the subgroup page). */
  sgNumber?: number;
}

/** A superspace operation x' = W x + w, W integer (3+d)×(3+d), w rational in [0,1). */
export interface SOp {
  W: number[][];
  w: Q[];
}

// ------------------------------------------------------------------ linear expressions

/**
 * Parse a linear expression over named variables, e.g. "-x+1/2", "2a1-b1+1/4", "z-u".
 * Returns coefficients (rational) per variable name and the constant.
 */
export function parseLinear(expr: string, names: string[]): { coef: Q[]; c: Q } {
  let e = expr.replace(/\s+/g, "").replace(/\*/g, "");
  if (!e) throw new SyntaxError("empty expression");
  if (e[0] !== "+" && e[0] !== "-") e = "+" + e;
  const coef = names.map(() => Q.ZERO);
  let c = Q.ZERO;
  const re = /([+-])([^+-]+)/g;
  let m: RegExpExecArray | null;
  let consumed = 0;
  while ((m = re.exec(e))) {
    consumed += m[0].length;
    const sgn = m[1] === "-" ? -1 : 1;
    const body = m[2];
    const t = /^(\d+(?:\/\d+)?)?([a-zA-Z][a-zA-Z]*\d*)$/.exec(body);
    if (t) {
      const idx = names.indexOf(t[2]);
      if (idx < 0) throw new SyntaxError(`unknown variable "${t[2]}" in "${expr}"`);
      const k = t[1] ? Q.parse(t[1]) : Q.ONE;
      coef[idx] = coef[idx].add(sgn < 0 ? k.neg() : k);
    } else {
      const k = Q.parse(body);
      c = c.add(sgn < 0 ? k.neg() : k);
    }
  }
  if (consumed !== e.length) throw new SyntaxError(`cannot parse "${expr}"`);
  return { coef, c };
}

// ------------------------------------------------------------------ operators

export function parseOp(s: string, nv: number): SOp {
  const parts = s.split(",");
  if (parts.length !== nv) throw new SyntaxError(`operator "${s}" has ${parts.length} components, expected ${nv}`);
  const names = VARS.slice(0, nv).split("");
  const W: number[][] = [];
  const w: Q[] = [];
  for (const p of parts) {
    const { coef, c } = parseLinear(p, names);
    W.push(coef.map((x) => {
      if (!x.isInt()) throw new SyntaxError(`non-integer rotation coefficient in "${s}"`);
      return Number(x.n);
    }));
    w.push(c);
  }
  return { W, w };
}

export function formatOp(op: SOp): string {
  const n = op.W.length;
  const out: string[] = [];
  for (let i = 0; i < n; i++) {
    let s = "";
    for (let j = 0; j < n; j++) {
      const c = op.W[i][j];
      if (!c) continue;
      const sign = c < 0 ? "-" : s ? "+" : "";
      s += sign + (Math.abs(c) === 1 ? "" : String(Math.abs(c))) + VARS[j];
    }
    if (!op.w[i].isZero()) s += (s ? "+" : "") + op.w[i].toString();
    out.push(s || "0");
  }
  return out.join(",");
}

const opKey = (op: SOp) => op.W.map((r) => r.join(",")).join(";") + "|" + op.w.map(String).join(",");

export function opMul(a: SOp, b: SOp): SOp {
  const n = a.W.length;
  const W = a.W.map((row) => Array.from({ length: n }, (_, j) => row.reduce((s, x, k) => s + x * b.W[k][j], 0)));
  const w = a.w.map((wi, i) => a.W[i].reduce((s, x, k) => s.add(b.w[k].mulInt(x)), wi).mod1());
  return { W, w };
}

function opCmp(a: SOp, b: SOp): number {
  const fa = a.W.flat(), fb = b.W.flat();
  for (let i = 0; i < fa.length; i++) if (fa[i] !== fb[i]) return fa[i] - fb[i];
  for (let i = 0; i < a.w.length; i++) {
    const c = a.w[i].cmp(b.w[i]);
    if (c) return c;
  }
  return 0;
}

export function identityOp(nv: number): SOp {
  return {
    W: Array.from({ length: nv }, (_, i) => Array.from({ length: nv }, (_, j) => (i === j ? 1 : 0))),
    w: Array.from({ length: nv }, () => Q.ZERO),
  };
}

/** Canonical operator list (CONTRACTS §1): closure mod Z^{3+d}, identity first, rest sorted. */
export function buildOps(rec: GroupRecord, maxOrder = 100000): SOp[] {
  const nv = 3 + rec.d;
  const gens: SOp[] = rec.gen.map((g) => {
    const op = parseOp(g, nv);
    return { W: op.W, w: op.w.map((x) => x.mod1()) };
  });
  for (const c of rec.cen) gens.push({ W: identityOp(nv).W, w: c.split(",").map((x) => Q.parse(x).mod1()) });
  const e = identityOp(nv);
  const seen = new Map<string, SOp>([[opKey(e), e]]);
  let frontier = [e];
  while (frontier.length) {
    const next: SOp[] = [];
    for (const a of frontier)
      for (const g of gens) {
        const b = opMul(a, g);
        const k = opKey(b);
        if (!seen.has(k)) {
          seen.set(k, b);
          next.push(b);
          if (seen.size > maxOrder) throw new Error("operator closure is not finite (bad generators?)");
        }
      }
    frontier = next;
  }
  const ek = opKey(e);
  const rest = [...seen.entries()].filter(([k]) => k !== ek).map(([, v]) => v).sort(opCmp);
  return [e, ...rest];
}

/** Split an op into R (3×3), M (d×3), E (d×d), tE (3), tI (d). */
export function splitOp(op: SOp, d: number) {
  const R = op.W.slice(0, 3).map((r) => r.slice(0, 3));
  const upperRight = op.W.slice(0, 3).flatMap((r) => r.slice(3));
  if (upperRight.some((x) => x !== 0)) throw new Error("operator is not in basic-setting form (upper-right block ≠ 0)");
  const M = op.W.slice(3).map((r) => r.slice(0, 3));
  const E = op.W.slice(3).map((r) => r.slice(3, 3 + d));
  return { R, M, E, tE: op.w.slice(0, 3), tI: op.w.slice(3, 3 + d) };
}

// ------------------------------------------------------------------ q-vectors and parameters

/** Free parameters of the q strings, in first-appearance order (e.g. a1, b1, g2). */
export function paramsOf(q: string[]): string[] {
  const out: string[] = [];
  for (const s of q) for (const tok of s.match(/[a-z]+\d*/g) ?? []) if (!out.includes(tok)) out.push(tok);
  return out;
}

/** Evaluate the q-vectors at rational parameter values: d×3 matrix Q (rows q_j). */
export function evalQ(rec: GroupRecord, values: Record<string, Q>): QMat {
  const names = rec.params;
  for (const p of names) if (!(p in values)) throw new Error(`missing value for parameter ${p}`);
  return rec.q.map((qs) =>
    qs.split(",").map((comp) => {
      const { coef, c } = parseLinear(comp, names);
      return coef.reduce((s, k, i) => s.add(k.mul(values[names[i]])), c);
    }),
  );
}

/** Check M + E·Q = Q·R exactly for every operator; returns the first violation or null. */
export function checkConsistency(ops: SOp[], Qm: QMat, d: number): string | null {
  for (const op of ops) {
    const { R, M, E } = splitOp(op, d);
    for (let i = 0; i < d; i++)
      for (let j = 0; j < 3; j++) {
        let lhs = Q.of(M[i][j]);
        for (let k = 0; k < d; k++) lhs = lhs.add(Qm[k][j].mulInt(E[i][k]));
        let rhs = Q.ZERO;
        for (let k = 0; k < 3; k++) rhs = rhs.add(Qm[i][k].mulInt(R[k][j]));
        if (!lhs.eq(rhs)) return `operator ${formatOp(op)}: M+EQ ≠ QR at (${i},${j})`;
      }
  }
  return null;
}

export type { QVec };
