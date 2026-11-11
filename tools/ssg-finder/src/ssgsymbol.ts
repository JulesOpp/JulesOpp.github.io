// (3+1)D superspace-group symbol in the GIVEN axes: basic HM symbol + q + phase-shift letters,
// written the way the BYU ISO(3+d)D tables write it, e.g. "P-3c1(0,0,g)000", "P321(0,0,g)t00",
// "P2_1/b(1/2,0,g)00", "Pbnm(0,0,g)s00". Validated against all 775 (3+1)D BYU symbols
// (test/ssgsymbol.test.ts, local-only).
//
// basic: the compact short HM of the basic space group in the given (possibly non-standard)
//   axes: settingSymbol, then the IT A setting-symbol conventions it does not apply itself
//   (itConventions), then the itShortForm decision with screw underscores kept.
// q: the q text of the family, verbatim, in parentheses.
// shifts: one letter per element written in the basic symbol, including "1" positions and
//   both halves of an "n/m" token ("P4/nbm" -> 4, n, b, m; "P321" -> 3, 2, 1; "P2/m" -> 2, m).
//
// Phase letter of an element (0 = 0, s = 1/2, t = 1/3, q = 1/4, h = 1/6; sign ignored).
// In the internal coordinate t~ = t - q_r.x (q_r = rational part of q) an operation
// (R, t3; t -> M.x + eps t + w4) acts as t~ -> eps t~ + tau~ with tau~ = w4 - q_r.t3, where t3 is
// its exact (unreduced) 3D translation. tau~ does not change under 3D origin shifts, but it does
// change by q_r.n between the parallel copies T(n)g of an element, so the copy matters:
//  1. operations: those of the group (incl. centrings) that realise the written element about
//     any direction of its symbol position (rotation of the written order, positive sense and
//     written screw; reflection with the written glide letter; rotoinversion), allowing a
//     lattice translation;
//  2. copies: their translates T(n)g, n in {-1,0,1}^3, that are themselves symbol elements
//     (screw lambda in [0,1); glide vector in [0,1)^3 with a glide letter, but not the horizontal
//     glide (a+b)/2 of a diagonal plane, which IT calls g);
//  3. keep the copies nearest the origin (squared fractional offset of the axis/plane, 0 when it
//     passes through the origin) and take the smallest letter among them;
//  4. eps = -1 always gives 0 (an internal origin shift removes the internal translation).
// Step 3 is what makes P4nc(1/2,1/2,g) read 0 (the n-glide through the origin) while P4bm(1/2,1/2,g)
// reads the m-plane at x-y = 1/2, and Pmn2_1(0,1/2,g) reads the 2_1 axis at (1/4,0,z). With q_r = 0
// every copy gives the same letter.

import { Q, type QVec } from "./rational";
import { parseLinear, paramsOf, type SOp } from "./ssg";
import type { Op3 } from "./identify";
import { _internal, settingSymbol, shortTokens, SYMMETRY_DIRECTIONS, type CrystalSystem } from "./symbol";

const { glideName, sense, shortestAlong, parallel } = _internal;

export interface SsgSymbol {
  basic: string;
  q: string;
  shifts: string;
  full: string;
}

const Z = Q.ZERO;
const qi = (n: number | bigint) => Q.of(BigInt(n));
const RANK = "0stqh";
const LETTERS: [string, Q][] = [
  ["0", Z],
  ["s", Q.of(1n, 2n)],
  ["t", Q.of(1n, 3n)],
  ["q", Q.of(1n, 4n)],
  ["h", Q.of(1n, 6n)],
];

function letterOf(x: Q): string | null {
  const r = x.mod1();
  for (const [l, v] of LETTERS) if (r.eq(v) || r.eq(Q.ONE.sub(v).mod1())) return l;
  return null;
}

// ------------------------------------------------------------------ small integer-matrix helpers

const R3 = (op: SOp) => op.W.slice(0, 3).map((r) => r.slice(0, 3));

function det3(R: number[][]): number {
  return (
    R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
    R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
    R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0])
  );
}

/** Order of the proper part (R or -R). */
function properOrder(R: number[][]): number {
  const tr = (R[0][0] + R[1][1] + R[2][2]) * det3(R);
  return { 3: 1, 2: 6, 1: 4, 0: 3, [-1]: 2 }[tr]!;
}

/** Axis of R (kernel of R - det(R) I), null for 1 and -1. */
function axisOf(R: number[][]): QVec | null {
  const s = det3(R);
  const M = R.map((r, i) => r.map((x, j) => qi(i === j ? x - s : x)));
  const cross = (a: QVec, b: QVec): QVec => [
    a[1].mul(b[2]).sub(a[2].mul(b[1])),
    a[2].mul(b[0]).sub(a[0].mul(b[2])),
    a[0].mul(b[1]).sub(a[1].mul(b[0])),
  ];
  for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) {
    const c = cross(M[i], M[j]);
    if (c.some((x) => !x.isZero())) return c;
  }
  return null;
}

function isIdentityW(W: number[][]): boolean {
  return W.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0)));
}

/** 3D intrinsic (screw/glide) translation (1/n) sum_{k<n} R^k t, n = order of R. */
function intrinsic(R: number[][], t: Q[]): Q[] {
  let acc = t.map(() => Z);
  let v = t;
  let P = R;
  let n = 1;
  for (;;) {
    acc = acc.map((x, i) => x.add(v[i]));
    v = R.map((row) => row.reduce((s, x, k) => s.add(v[k].mulInt(BigInt(x))), Z));
    if (isIdentityW(P)) break;
    P = P.map((row) => [0, 1, 2].map((j) => row.reduce((s, x, k) => s + x * R[k][j], 0)));
    n++;
  }
  return acc.map((x) => x.div(qi(n)));
}

// ------------------------------------------------------------------ symbol elements

type Elem =
  | { kind: "one" }
  | { kind: "rot"; n: number; screw: number }
  | { kind: "mirror"; glide: string }
  | { kind: "rotoinv"; n: number };

function parseToken(tok: string): Elem[] {
  return tok.split("/").map((part): Elem => {
    let m: RegExpExecArray | null;
    if (part === "1") return { kind: "one" };
    if ((m = /^-([1346])$/.exec(part))) return { kind: "rotoinv", n: Number(m[1]) };
    if ((m = /^([2346])(?:_(\d))?$/.exec(part))) return { kind: "rot", n: Number(m[1]), screw: m[2] ? Number(m[2]) : 0 };
    if (/^[mabcnd]$/.test(part)) return { kind: "mirror", glide: part };
    throw new Error(`cannot parse symbol token "${tok}"`);
  });
}

interface Ctx {
  ops: SOp[];
  cents: QVec[]; // centring vectors incl. 0 (3D, mod 1)
  axisLen: Q[];
}

/** A lattice translate T(n)·g of an operation that is a symbol element. */
interface Copy {
  t3: Q[]; // exact 3D translation
  w4: Q; // internal translation (unchanged by a 3D lattice translation)
  eps: number;
  dist: Q; // squared fractional offset of the element's axis/plane from the origin
  exact: boolean; // is the written element itself (screw index / glide letter)
}

/** Integer shifts {-1,0,1}^3. */
const SHIFTS3: number[][] = [-1, 0, 1].flatMap((a) => [-1, 0, 1].flatMap((b) => [-1, 0, 1].map((c) => [a, b, c])));

/**
 * Translates T(n)·op, n in {-1,0,1}^3, that are symbol elements of the same kind as `e`
 * about/perpendicular to d (see the header, step 2). Empty if op's rotation part does not fit.
 */
function copies(e: Elem, d: QVec, op: SOp, ctx: Ctx): Copy[] {
  if (e.kind === "one") return [];
  const R = R3(op);
  const dt = det3(R);
  const ord = properOrder(R);
  if (ord === 1) {
    if (!(e.kind === "rotoinv" && e.n === 1 && dt < 0)) return [];
  } else {
    const ax = axisOf(R);
    if (!ax || !parallel(ax, d)) return [];
    if (e.kind === "rotoinv" && !(dt < 0 && ord === e.n)) return []; // -n = -(n-fold), either sense
    if (e.kind === "mirror" && !(dt < 0 && ord === 2)) return [];
    if (e.kind === "rot" && !(dt > 0 && ord === e.n && (ord === 2 || sense(R.map((r) => r.map(qi)), d) > 0))) return [];
  }
  const axial = d.filter((x) => !x.isZero()).length === 1;
  const out: Copy[] = [];
  for (const s of SHIFTS3) {
    const t3 = op.w.slice(0, 3).map((x, i) => x.add(qi(s[i])));
    const v3 = intrinsic(R, t3);
    if (!v3.every((x) => x.cmp(Z) >= 0 && x.cmp(Q.ONE) < 0)) continue;
    let exact = true;
    if (e.kind === "mirror") {
      const name = glideName(v3.map((x, i) => x.div(ctx.axisLen[i])));
      if (name === "g?") continue;
      if (!axial && v3[2].isZero() && !(v3[0].isZero() && v3[1].isZero())) continue; // horizontal g
      exact = name === e.glide;
    } else if (e.kind === "rot") {
      const u = d.map((x) => x.mul(shortestAlong(d, ctx.cents)));
      const k = u.findIndex((x) => !x.isZero());
      const lam = v3[k].div(u[k]);
      if (!v3.every((x, i) => x.eq(u[i].mul(lam)))) continue;
      exact = lam.eq(Q.of(BigInt(e.screw), BigInt(e.n)));
    } else if (!v3.every((x) => x.isZero())) continue;
    const off = t3.map((x, i) => x.sub(v3[i]));
    out.push({ t3, w4: op.w[3], eps: op.W[3][3], dist: off.reduce((a, x) => a.add(x.mul(x)), Z), exact });
  }
  return out;
}

/** Does op (up to a lattice translation) realise exactly element e about/perpendicular to d? */
function realises(e: Elem, d: QVec, op: SOp, ctx: Ctx): boolean {
  return copies(e, d, op, ctx).some((c) => c.exact);
}

function planeAvailable(letter: string, d: QVec, ctx: Ctx): boolean {
  const e: Elem = { kind: "mirror", glide: letter };
  return ctx.ops.some((op) => realises(e, d, op, ctx));
}

const AXES: QVec[] = [0, 1, 2].map((i) => [0, 1, 2].map((j) => (i === j ? Q.ONE : Z)));

// ------------------------------------------------------------------ IT setting-symbol conventions

// Old-style (IT A before 2002, as BYU writes them) standard short symbols of the orthorhombic
// types whose glide letters are not fixed by a local rule: e-planes (two axial glides in one
// plane) and I-centred planes with two axial glides at different heights. Keyed by space-group
// number; tokens are [centring, ⟂a, ⟂b, ⟂c] ("2" = no plane at that position).
const ORTHO_STANDARD: Record<number, string[]> = {
  39: ["A", "b", "m", "2"],
  41: ["A", "b", "a", "2"],
  45: ["I", "b", "a", "2"],
  46: ["I", "m", "a", "2"],
  64: ["C", "m", "c", "a"],
  67: ["C", "m", "m", "a"],
  68: ["C", "c", "c", "a"],
  72: ["I", "b", "a", "m"],
  73: ["I", "b", "c", "a"],
  74: ["I", "m", "m", "a"],
};

// Axis permutations, new axis i = old (standard) axis PERMS[k][i], in order of preference:
// q (along the new c axis) along the standard c, then a, then b (abc, ba-c, -cba, bca, cab, a-cb).
const PERMS = [[0, 1, 2], [1, 0, 2], [2, 1, 0], [1, 2, 0], [2, 0, 1], [0, 2, 1]];

/** IT setting name of a group of one of the ORTHO_STANDARD types; null if none fits. */
function orthoSetting(cen: string, tokens: string[], ctx: Ctx): string[] | null {
  const plane = (t: string) => (t.includes("/") ? t.split("/")[1] : /^[mabcnde]$/.test(t) ? t : null);
  for (const std of Object.values(ORTHO_STANDARD)) {
    for (const perm of PERMS) {
      const newIdx = (old: number) => perm.indexOf(old);
      const c = /^[ABC]$/.test(std[0]) ? "ABC"[newIdx("ABC".indexOf(std[0]))] : std[0];
      if (c !== cen) continue;
      const cand = [0, 1, 2].map((i) => {
        const L = std[1 + perm[i]];
        return /^[abc]$/.test(L) ? "abc"[newIdx("abc".indexOf(L))] : L;
      });
      const ok = cand.every((L, i) => {
        const p = plane(tokens[i]);
        if (L === "2") return p === null;
        return p !== null && planeAvailable(L, AXES[i], ctx);
      });
      if (ok) return tokens.map((t, i) => (cand[i] === "2" ? t : t.includes("/") ? `${t.split("/")[0]}/${cand[i]}` : cand[i]));
    }
  }
  return null;
}

/** I222 vs I2_12_12_1: do pure 2-fold axes along a, b, c meet in a common point? */
function pureTwofoldsIntersect(ctx: Ctx): boolean {
  // locations mod 1/2 of pure 2-fold axes along axis i: (x_j, x_k) = (w_j/2, w_k/2)
  const half = Q.of(1n, 2n);
  const locs = [0, 1, 2].map((i) => {
    const out: Q[][] = [];
    for (const op of ctx.ops) {
      const R = R3(op);
      if (!(det3(R) > 0 && R[i][i] === 1 && R.every((r, a) => r.every((x, b) => x === (a === b ? (a === i ? 1 : -1) : 0))))) continue;
      if (!op.w[i].mod1().isZero()) continue; // screw
      out.push([0, 1, 2].map((j) => (j === i ? Z : op.w[j].mul(half).modQ(half))));
    }
    return out;
  });
  const vals = (j: number) => [...new Map(locs.flat().map((l) => [l[j].toString(), l[j]])).values()];
  for (const X of vals(0)) for (const Y of vals(1)) for (const Zv of vals(2)) {
    const P = [X, Y, Zv];
    if ([0, 1, 2].every((i) => locs[i].some((l) => [0, 1, 2].every((j) => j === i || l[j].eq(P[j]))))) return true;
  }
  return false;
}

/** The IT A setting conventions settingSymbol leaves out. Tokens without the centring. */
function itConventions(system: CrystalSystem, cen: string, tokens: string[], ctx: Ctx): string[] | null {
  // a position whose members carry glides relabelled by the main axis ("[b;a]") is written by
  // its first member, the [100] direction (P4bm, P-4b2, P4/nbm)
  tokens = tokens.map((t) => (t.startsWith("[") ? t.slice(1, -1).split(";")[0] : t));
  // "<cen> 1" / "<cen> -1": no symmetry direction carries an element
  if (tokens.length !== SYMMETRY_DIRECTIONS[system].length) return tokens.some((t) => /[e?\[]/.test(t)) ? null : tokens;
  if (system === "tetragonal" && cen === "I") {
    // I-centred tetragonal: the glide perpendicular to [100] is called c, not b (I4cm, I4/mcm)
    const t = tokens[1];
    const p = t.includes("/") ? t.split("/")[1] : t;
    if (p === "b" && planeAvailable("c", AXES[0], ctx)) tokens[1] = t.includes("/") ? `${t.split("/")[0]}/c` : "c";
  }
  if (system === "orthorhombic") {
    if (cen === "I" && tokens.every((t) => /^2(_1)?$/.test(t)) && !pureTwofoldsIntersect(ctx)) return ["2_1", "2_1", "2_1"];
    const it = orthoSetting(cen, tokens, ctx);
    if (it) return it;
  }
  if (tokens.some((t) => /[e?\[]/.test(t))) return null;
  return tokens;
}

// ------------------------------------------------------------------ main entry

/** Rational part of q (constant terms of the q text). */
function rationalPart(q: string): Q[] {
  const names = paramsOf([q]);
  return q.split(",").map((c) => parseLinear(c, names).c);
}

/**
 * BYU-style (3+1)D symbol of the group `ops` (canonical (3+1)D operation list incl. centrings)
 * in its own axes. `qFamily` is the q text (one string, d = 1). Null if the basic symbol cannot
 * be written in these axes.
 */
export function ssgSymbolInAxes(ops: SOp[], qFamily: string[], system: CrystalSystem): SsgSymbol | null {
  if (qFamily.length !== 1 || ops.some((o) => o.W.length !== 4)) throw new Error("only d = 1 is supported for now");
  const ops3: Op3[] = ops.map((o) => ({ R: R3(o), t: o.w.slice(0, 3) }));
  const setting = settingSymbol(ops3, system);
  if (!setting) return null;
  const [cen0, ...tok0] = setting.split(" ");
  if (!/^[PABCIFR]$/.test(cen0)) return null;

  const cents: QVec[] = [[Z, Z, Z]];
  for (const o of ops) if (isIdentityW(R3(o)) && o.w.slice(0, 3).some((x) => !x.isZero())) cents.push(o.w.slice(0, 3));
  const ctx: Ctx = { ops, cents, axisLen: AXES.map((a) => shortestAlong(a, cents)) };

  const fixed = system === "triclinic" ? tok0 : itConventions(system, cen0, tok0, ctx);
  if (!fixed) return null;
  const { cen, tokens } = shortTokens([cen0, ...fixed].join(" "));
  const positions = SYMMETRY_DIRECTIONS[system];
  // one letter per element (header, steps 1-4)
  const qr = rationalPart(qFamily[0]);
  let shifts = "";
  for (const { pos, token } of tokens) {
    const dirs: QVec[] = (positions[pos] ?? [[0, 0, 1]]).map((v) => v.map(qi));
    for (const e of parseToken(token)) {
      if (e.kind === "one") {
        shifts += "0";
        continue;
      }
      const cs = dirs.flatMap((d) => ops.filter((op) => realises(e, d, op, ctx)).flatMap((op) => copies(e, d, op, ctx)));
      if (!cs.length) return null;
      const near = cs.reduce((a, c) => (c.dist.cmp(a) < 0 ? c.dist : a), cs[0].dist);
      let best: string | null = null;
      for (const c of cs) {
        if (!c.dist.eq(near)) continue;
        const l = c.eps === -1 ? "0" : letterOf(qr.reduce((acc, x, i) => acc.sub(x.mul(c.t3[i])), c.w4));
        if (l !== null && (best === null || RANK.indexOf(l) < RANK.indexOf(best))) best = l;
      }
      if (best === null) return null;
      shifts += best;
    }
  }

  const basic = cen + tokens.map((t) => t.token).join("");
  const q = `(${qFamily[0]})`;
  return { basic, q, shifts, full: basic + q + shifts };
}
