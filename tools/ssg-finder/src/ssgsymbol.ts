// d>1 follows Stokes/Campbell/van Smaalen 2011 §5.2: choose the globally nicest
// whole phase assignment from the scg/diag copy set (nearest copies plus SCG-lattice
// equivalents with zero or equal diagonal q offsets; both half-half tetragonal mirror classes).
// Glide phases use unsigned representatives; rotation phases retain their orientation signs.
// About 0.75% of d=2 symbols differ from BYU because of its step-3 degeneracy order, in
// P222, C222, A222, F222, Pm2m, Pmm2, Pb2_1m, Pc2m, Pn2_1m, Am2m, Ac2m, Pmmb,
// I422, I4_122 and P4_2/mmc.
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
//     This is squared offset in lattice coordinates, not Cartesian distance. For
//     hexagonal 3/6-fold axes and diagonal planes it is a tie-break order; the lattice
//     coordinate norm can rank copies differently from the Cartesian metric.
//     Retain the smallest letter on disagreement: the 775-record d = 1 sweep finds
//     two such positions, in basic groups 126 and 133, both with q = (1/2,1/2,g).
//     This is empirical support for these tables, not a guarantee for arbitrary inputs.
//  4. eps = -1 always gives 0 (an internal origin shift removes the internal translation).
// Step 3 is what makes P4nc(1/2,1/2,g) read 0 (the n-glide through the origin) while P4bm(1/2,1/2,g)
// reads the m-plane at x-y = 1/2, and Pmn2_1(0,1/2,g) reads the 2_1 axis at (1/4,0,z). With q_r = 0
// every copy gives the same letter.

import { Q, denLcm, identity, matEq, matInv, matMul, matVec, toSafeInt, transpose, type QMat, type QVec } from "./rational";
import { hnf, imatDet } from "./intmat";
import { parseLinear, paramsOf, type SOp } from "./ssg";
import type { Op3 } from "./identify";
import { _internal, settingSymbol, shortTokens, SYMMETRY_DIRECTIONS, type CrystalSystem } from "./symbol";

const { glideName, sense, shortestAlong, parallel } = _internal;

export interface SsgSymbol {
  basic: string;
  q: string;
  shifts: string[];
  full: string;
}

/** Test-only number of symbol positions whose nearest copies give different letters. */
export const _stats = { disagree: 0 };

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
  wI: Q[]; // internal translation (unchanged by a 3D lattice translation)
  E: number[][];
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
    const v3raw = intrinsic(R, t3);
    // In a threefold star a c glide can use a lattice-equivalent z copy.
    const starGlide = op.W.length > 4 && e.kind === "mirror" && e.glide === "c"
      && op.W.slice(3).some(row => row.slice(3).filter(x => x !== 0).length > 1);
    const v3 = starGlide ? v3raw.map(x => x.mod1()) : v3raw;
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
    const off = t3.map((x, i) => x.sub(v3raw[i]));
    out.push({ t3, wI: op.w.slice(3), E: op.W.slice(3).map(row => row.slice(3)), dist: off.reduce((a, x) => a.add(x.mul(x)), Z), exact });
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
  const alternatives: string[][] = [];
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
      if (ok) {
        const choice = tokens.map((t, i) => (cand[i] === "2" ? t : t.includes("/") ? `${t.split("/")[0]}/${cand[i]}` : cand[i]));
        if (ctx.ops[0].w.length === 4) return choice;
        alternatives.push(choice);
      }
    }
  }
  // When several conventional axial-glide descriptions describe the same
  // higher-dimensional basic setting, prefer the alphabetically first one.
  // This leaves the established d1 permutation preference unchanged.
  // In the I-centred m-glide-m setting, the conventional middle glide
  // is c rather than a. Keep the scalar setting preference above intact.
  if (cen === "I") {
    const mcm = alternatives.find(a => a.map(t => t.split("/").at(-1)).join("") === "mcm");
    if (mcm) return mcm;
  }
  alternatives.sort((a, b) => a.join("") < b.join("") ? -1 : a.join("") > b.join("") ? 1 : 0);
  return alternatives[0] ?? null;
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
  const d = qFamily.length, n = 3 + d;
  if (d < 1 || d > 3 || ops.some(o => o.W.length !== n || o.w.length !== n || o.W.some(row => row.length !== n)))
    throw new RangeError("inconsistent superspace dimension");
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
  const shortened = shortTokens([cen0, ...fixed].join(" "));
  const cen = shortened.cen;
  let tokens = shortened.tokens;
  if (d > 1 && system === "cubic") {
    // Cubic23 andm-3 omit the absent third direction in the short HM symbol.
    if (tokens.length === 3 && tokens[2].token === "1") tokens = tokens.slice(0, 2);
    if (tokens[1]?.token === "-3") tokens = tokens.map(t => ({ ...t, token: t.token.includes("/") ? t.token.split("/")[1] : t.token }));
    // The same nonsymmorphic axis-intersection distinction used forI222.
    if (cen === "I" && tokens[0]?.token === "2" && tokens[1]?.token === "3" && !pureTwofoldsIntersect(ctx))
      tokens = tokens.map((t, i) => i === 0 ? { ...t, token: "2_1" } : t);
    // Conventional primitive cubic diagonal glides use the fulln vector.
    if (cen === "P" && tokens[2]?.token === "c" && planeAvailable("n", [qi(1), qi(-1), Z], ctx))
      tokens = tokens.map((t, i) => i === 2 ? { ...t, token: "n" } : t);
  }
  const positions = SYMMETRY_DIRECTIONS[system];
  // one letter per element (header, steps 1-4)
  const qr = rationalPart(qFamily[0]);
  // EXPERIMENT ONLY. Coefficient-lattice SCG basis and copy selection not validated.
  const { D, Q0 } = d === 1
    ? { D: identity(1), Q0: [qr] }
    : hasConventionalInternalMetric(qFamily, system) ? { D: identity(d), Q0: qFamily.map(rationalPart) } : proposedScgBasis(qFamily);
  const shifts = qFamily.map(() => "");
  if (d > 1) {
    const g = globalShifts(tokens, positions, system, ops, ctx, Q0, D, d);
    if (!g) return null;
    const basic = cen + tokens.map(t => t.token).join("");
    return { basic, q: qFamily.map(q => `(${q})`).join(""), shifts: g,
      full: basic + qFamily.map((q, i) => `(${q})${g[i]}`).join("") };
  }
  for (const { pos, token } of tokens) {
    const dirs: QVec[] = (positions[pos] ?? [[0, 0, 1]]).map((v) => v.map(qi));
    for (const e of parseToken(token)) {
      if (e.kind === "one") {
        for (let i = 0; i < d; i++) shifts[i] += "0";
        continue;
      }
      const point = d > 1 ? fixedGenerator(e, pos, system) : null;
      const element: Elem = e;
      const matching = point === null ? ops : ops.filter(op => R3(op).every((row, i) => row.every((x, j) => x === point[i][j])));
      const cs = dirs.flatMap((direction) => matching.filter((op) => realises(element, direction, op, ctx)).flatMap((op) => copies(element, direction, op, ctx)));
      if (!cs.length) return null;
      const near = cs.reduce((a, c) => (c.dist.cmp(a) < 0 ? c.dist : a), cs[0].dist);
      if (d === 1) {
        // Preserve the Phase A scalar rule and its two measured disagreements.
        let best: string | null = null;
        const nearestLetters = new Set<string>();
        for (const c of cs) {
          if (!c.dist.eq(near)) continue;
          const l = c.E[0][0] === -1 ? "0" : letterOf(qr.reduce((acc, x, i) => acc.sub(x.mul(c.t3[i])), c.wI[0]));
          if (l !== null) nearestLetters.add(l);
          if (l !== null && (best === null || RANK.indexOf(l) < RANK.indexOf(best))) best = l;
        }
        if (nearestLetters.size > 1) _stats.disagree++;
        if (best === null) return null;
        shifts[0] += best;
      }
    }
  }

  const basic = cen + tokens.map((t) => t.token).join("");
  const q = qFamily.map(q => `(${q})`).join("");
  return { basic, q, shifts, full: basic + qFamily.map((q, i) => `(${q})${shifts[i]}`).join("") };
}

/** Test-only opt-in candidate capture for local diagnostics; null disables collection. */
export const _testOnly = { candidates: null as string[][] | null };

/**
 * d > 1: Stokes, Campbell & van Smaalen (2011) §5.2. For each symbol generator, every lattice
 * translate whose external intrinsic translation matches the written element (step 2) gives a
 * candidate internal intrinsic translation (IIT, SCG coordinates); the symbol is the "nicest"
 * whole assignment (step 4, compareIntrinsicAssignments).
 */
function globalShifts(
  tokens: { pos: number; token: string }[], positions: number[][][], system: CrystalSystem,
  ops: SOp[], ctx: Ctx, Q0: QMat, D: QMat, d: number,
): string[] | null {
  const options: Q[][][] = []; // per generator: candidate IIT vectors (signed, in (-1/2, 1/2])
  for (const { pos, token } of tokens) {
    const dirs: QVec[] = (positions[pos] ?? [[0, 0, 1]]).map((v) => v.map(qi));
    for (const e of parseToken(token)) {
      if (e.kind === "one") { options.push([Array.from({ length: d }, () => Z)]); continue; }
      const point = fixedGenerator(e, pos, system);
      const matching = point === null ? ops : ops.filter(op => R3(op).every((row, i) => row.every((x, j) => x === point[i][j])));
      const seen = new Map<string, Q[]>();
      const all = dirs.flatMap(direction => matching.flatMap(op => copies(e, direction, op, ctx))).filter(c => c.exact);
      const near = all.reduce((a: Q | null, c) => (a === null || c.dist.cmp(a) < 0 ? c.dist : a), null);
      const bases = all.filter(c => c.dist.eq(near!));
      const scgStep = (diff: Q[]) => ctx.cents.some(cv => {
        const n = diff.map((x, i) => x.sub(cv[i]));
        if (!n.every(x => x.isInt())) return false;
        const qn = Q0.map(row => row.reduce((s, x, i) => s.add(x.mul(n[i])), Z).mod1());
        if (qn.every(x => x.isZero())) return true;
        return qn.every(x => x.eq(qn[0]));
      });
      // Half-half tetragonal axial mirrors need both parallel plane classes.
      const halfHalfMirror = system === "tetragonal" && pos === 1 && e.kind === "mirror"
        && Q0.some(row => row[0].mod1().eq(Q.of(1, 2)) && row[1].mod1().eq(Q.of(1, 2)));
      for (const c of all) {
        if (!halfHalfMirror && !bases.some(b => scgStep(c.t3.map((x, i) => x.sub(b.t3[i]))))) continue;
        for (const v of intrinsicVectors(c.E, c.wI, c.t3, Q0, D)) {
          if (!v.every(x => signedLetter(x) !== null)) continue;
          // Reversing a glide reverses its intrinsic translation; choose
          // the unsigned representative, retaining oriented rotation signs.
          const s = v.map(x => e.kind === "mirror" ? (signed(x).sign() < 0 ? signed(x).neg() : signed(x)) : signed(x));
          seen.set(s.map(String).join(","), s);
        }
      }
      if (!seen.size) return null;
      options.push([...seen.values()]);
    }
  }
  // order of the flat sequence used by the "nicest" comparison
  const flat = (choice: Q[][]): Q[] => Array.from({ length: d }, (_, i) => choice.map(v => v[i])).flat();
  let best: Q[][] | null = null, bestFlat: Q[] | null = null;
  const pick: Q[][] = [];
  const rec = (k: number) => {
    if (k === options.length) {
      const f = flat(pick);
      _testOnly.candidates?.push(Array.from({ length: d }, (_, i) => pick.map(v => signedLetter(v[i])!).join("")));
      if (!bestFlat || compareIntrinsicAssignments(f, bestFlat) < 0) { best = pick.slice(); bestFlat = f; }
      return;
    }
    for (const v of options[k]) { pick.push(v); rec(k + 1); pick.pop(); }
  };
  rec(0);
  if (!best) return null;
  const chosen: Q[][] = best;
  return Array.from({ length: d }, (_, i) => chosen.map(v => signedLetter(v[i])!).join(""));
}

// Exact projection and finite translation quotient helpers.

const add = (a: Q[], b: Q[]) => a.map((x, i) => x.add(b[i]));
const mod = (a: Q[]) => a.map(x => x.mod1());
const key = (a: Q[]) => a.map(String).join(",");

/** P_E = (I + E + ... + E^(n-1))/n, the origin-independent projection. */
export function internalProjection(E: number[][]): QMat {
  const d = E.length;
  if (d < 1 || d > 3 || E.some(row => row.length !== d || row.some(x => !Number.isSafeInteger(x))))
    throw new RangeError("invalid internal representation dimension or integer entry");
  const R = E.map(row => row.map(x => Q.of(x)));
  let power = identity(d), sum = identity(d).map(row => row.map(() => Z));
  for (let n = 1; n <= 12; n++) {
    sum = sum.map((row, i) => row.map((x, j) => x.add(power[i][j])));
    power = matMul(power, R);
    if (matEq(power, identity(d))) return sum.map(row => row.map(x => x.div(Q.of(n))));
  }
  throw new RangeError("internal representation does not have crystallographic finite order");
}

/** All A Z^d modulo Z^d, without a denominator^d Cartesian enumeration. */
export function projectedIntegerTranslations(A: QMat): Q[][] {
  const columns = transpose(A).map(mod);
  const zero = A.map(() => Z), queue = [zero], seen = new Set([key(zero)]);
  for (let at = 0; at < queue.length; at++) for (const v of columns) {
    const next = mod(add(queue[at], v)), k = key(next);
    if (!seen.has(k)) { seen.add(k); queue.push(next); }
  }
  return queue;
}

/** Coupled vectors attainable by changing the integer internal translate of one copy. */
export function intrinsicVectors(
  E: number[][], wI: Q[], t3: Q[], Q0: QMat, D: QMat,
): Q[][] {
  // SCG coordinates T = D^-1 (t - Q0 x). D is a direct basis matrix.
  const A = matMul(matInv(D), internalProjection(E));
  const tau = wI.map((x, i) => x.sub(matVec(Q0, t3)[i]));
  const base = matVec(A, tau);
  return projectedIntegerTranslations(A).map(v => mod(add(base, v)));
}

export function signedLetter(x: Q): string | null {
  const r = x.mod1();
  const values: [string, Q][] = [
    ["0", Z], ["s", Q.of(1, 2)], ["t", Q.of(1, 3)], ["q", Q.of(1, 4)], ["h", Q.of(1, 6)],
    ["-t", Q.of(2, 3)], ["-q", Q.of(3, 4)], ["-h", Q.of(5, 6)],
  ];
  return values.find(([, v]) => r.eq(v))?.[0] ?? null;
}

const signed = (x: Q) => { const r = x.mod1(); return r.cmp(Q.of(1, 2)) > 0 ? r.sub(Q.ONE) : r; };
const cmpInt = (a: bigint, b: bigint) => a < b ? -1 : a > b ? 1 : 0;

/** Paper section 5.2 step 4. Caller supplies a WHOLE ordered phase assignment. */
export function compareIntrinsicAssignments(a0: Q[], b0: Q[]): number {
  const a = a0.map(signed), b = b0.map(signed);
  if (a.length !== b.length) throw new RangeError("different phase assignment dimensions");
  const negatives = (v: Q[]) => v.filter(x => x.sign() < 0).length;
  const zeros = (v: Q[]) => v.filter(x => x.isZero()).length;
  const maxden = (v: Q[]) => v.reduce((m, x) => x.d > m ? x.d : m, 1n);
  let c = negatives(a) - negatives(b) || zeros(b) - zeros(a) || cmpInt(maxden(a), maxden(b));
  if (c) return c;
  for (let i = 0; i < a.length; i++) if ((c = cmpInt(a[i].d, b[i].d))) return c;
  for (let i = 0; i < a.length; i++) if ((c = cmpInt(a[i].n, b[i].n))) return c;
  return 0;
}

/**
 * EXPERIMENTAL ONLY: coefficient-lattice SCG basis proposal from source notes.
 * Preserve primitive coefficient lattices, so hexagonal 120-degree star axes
 * are not automatically changed into a different conventional setting.
 * This basis has not been checked against the local corpus or published SCG
 * transforms. It must not be enabled in production before those checks.
 */
export function proposedScgBasis(qFamily: string[]): { D: QMat; Q0: QMat; index: Q } {
  const d = qFamily.length, names = paramsOf(qFamily);
  const parsed = qFamily.map(q => q.split(",").map(c => parseLinear(c, names)));
  if (parsed.some(row => row.length !== 3)) throw new RangeError("q dimension");
  const Q0 = parsed.map(row => row.map(x => x.c));
  const C = parsed.map(row => names.flatMap((_, p) => row.map(x => x.coef[p])));
  const den = denLcm(C.flat());
  const columns = transpose(C).map(row => row.map(x => x.mulInt(den).n));
  const H = hnf(columns).H;
  if (H.length !== d) throw new RangeError("q-family coefficient lattice has deficient rank");
  const B = transpose(H).map(row => row.map(x => Q.of(x, den)));
  const determinant = imatDet(transpose(H));
  const index = Q.of(determinant < 0n ? -determinant : determinant, den ** BigInt(d));
  if (index.eq(Q.ONE) && B.flat().every(x => x.isInt())) return { D: identity(d), Q0, index };
  const Cprimitive = matMul(matInv(B), C).map(row => row.map(x => {
    // Keep bigint here; no accidental rational-to-number truncation.
    if (!x.isInt()) throw new RangeError("coefficient lattice basis is not exact");
    return x.n;
  }));
  const { U } = hnf(Cprimitive);
  const D = matMul(B, matInv(U.map(row => row.map(x => Q.of(x)))));
  return { D, Q0, index };
}

/** Diagnostic: require the candidate basis to keep E crystallographic/integral. */
export function internalInScg(E: number[][], D: QMat): number[][] {
  return matMul(matMul(matInv(D), E.map(row => row.map(x => Q.of(x)))), D)
    .map(row => row.map(x => toSafeInt(x, "SCG internal operator entry")));
}

/**
 * Fixed point parts for the symbol generators, from Table2 of Stokes,
 * Campbell & van Smaalen, Acta Cryst.A67(2011),45–55, p51.
 * https://iso.byu.edu/2011%20Stokes.pdf
 * These are published crystallographic generator conventions, not ISO rows.
 * Apply only to d>1, preserving the validated scalar convention.
 */
function fixedGenerator(e: Elem, pos: number, system: CrystalSystem): number[][] | null {
  let xyz: string | null = null;
  if (system === "cubic") {
    if (pos === 1) xyz = e.kind === "rotoinv" ? "-z,-x,-y" : "z,x,y";
    else if (pos === 2) xyz = e.kind === "mirror" ? "y,x,z" : "-y,-x,-z";
    else if (e.kind === "mirror") xyz = "x,y,-z";
    else if (e.kind === "rot") xyz = e.n === 4 ? "-y,x,z" : "-x,-y,z";
    else if (e.kind === "rotoinv" && e.n === 4) xyz = "y,-x,-z";
  } else if (system === "tetragonal") {
    if (pos === 0) {
      if (e.kind === "mirror") xyz = "x,y,-z";
      else if (e.kind === "rot") xyz = "-y,x,z";
      else if (e.kind === "rotoinv") xyz = "y,-x,-z";
    } else if (pos === 1) xyz = e.kind === "mirror" ? "-x,y,z" : "x,-y,-z";
    else xyz = e.kind === "mirror" ? "y,x,z" : "-y,-x,-z";
  } else if (system === "hexagonal" || system === "trigonal") {
    if (pos === 0) {
      if (e.kind === "mirror") xyz = "x,y,-z";
      else if (e.kind === "rot") xyz = e.n === 6 ? "x-y,x,z" : "-y,x-y,z";
      else if (e.kind === "rotoinv") xyz = e.n === 6 ? "-x+y,-x,-z" : "y,-x+y,-z";
    } else if (pos === 1) xyz = e.kind === "mirror" ? "x,x-y,z" : "-x,-x+y,-z";
    else xyz = e.kind === "mirror" ? "y,x,z" : "-y,-x,-z";
  }
  return xyz === null ? null : xyz.split(",").map(c => parseLinear(c, ["x", "y", "z"]).coef.map(x => toSafeInt(x, "published point-generator coefficient")));
}

/** Whether the formal q rows already have orthogonal or conventional120-degree
 * internal axes for every parameter value. Coefficients are compared exactly;
 * evaluating all independent parameters at1 would spuriously turn an ordinary
 * orthorhombic reflected pair into a square pair.
 * Diagnostic basis-selection hypothesis; validate against every local symbol.
 */
function hasConventionalInternalMetric(qFamily: string[], system: CrystalSystem): boolean {
  const names = paramsOf(qFamily), d = qFamily.length;
  const parsed = qFamily.map(q => q.split(",").map(c => parseLinear(c, names)));
  const coeffs = names.map((_, k) => parsed.map(row => row.map(x => x.coef[k])));
  const metric = identity(3);
  if (system === "hexagonal" || system === "trigonal") metric[0][1] = metric[1][0] = Q.of(1, 2);
  const terms: QMat[] = [];
  for (let p = 0; p < names.length; p++) for (let q = p; q < names.length; q++) {
    const a = matMul(matMul(coeffs[p], metric), transpose(coeffs[q]));
    const b = matMul(matMul(coeffs[q], metric), transpose(coeffs[p]));
    terms.push(p === q ? a : a.map((row, i) => row.map((x, j) => x.add(b[i][j]))));
  }
  const orthogonal = (i: number, j: number) => terms.every(t => t[i][j].isZero());
  if (Array.from({ length: d }, (_, i) => i).every(i => Array.from({ length: i }, (_, j) => j).every(j => orthogonal(i, j)))) return true;
  for (let i = 0; i < d; i++) for (let j = i + 1; j < d; j++) {
    if (!terms.every(t => t[i][i].eq(t[j][j]) && t[i][j].mulInt(2).add(t[i][i]).isZero())) continue;
    let ok = true;
    for (let a = 0; a < d; a++) for (let b = a + 1; b < d; b++) if (!(a === i && b === j) && !orthogonal(a, b)) ok = false;
    if (ok) return true;
  }
  return false;
}
