// Standard (BYU-convention) setting and type name of a (3+1)D superspace group (d = 1).
//
// Input: the group's operations (canonical (3+1)D list incl. centrings, any setting) and its q
// text in those axes. Output: the type symbol (ssgSymbolInAxes in the standard setting, q written
// as its family), the actual q in that setting, and one map x' = Px + p, x4' = m.x + eta x4 + p4
// that takes the given operations there.
//
// Algorithm.
//  1. moyo standardises the basic group: x_std = T x + t (IT A standard setting, unique axis b,
//     hexagonal axes for R).
//  2. The affine normaliser of the standard basic group: linear parts N (integer, entries in
//     {-1,0,1}, det +-1, normalising the point group and the centred lattice; one per coset N.K)
//     with every translation n solving the congruences (exactly, Smith normal form).
//  3. BYU's axes for the crystal system (B below): monoclinic unique axis c (b-unique C2/c becomes
//     B2/b); the 6 orthorhombic settings of IT A; identity otherwise. With x' = B((N,n)(T x + t)):
//     P = B N T, p = B(N t + n). The basic group is then the same for every (N, n); only q and the
//     4th rows change.
//  4. For eta = +-1: q' = eta q P^-1 must lie in one of the (3+1)D Bravais-class q forms of that
//     crystal system and centring (Q_FORMS: the IT C list of q forms, (0,0,g), (1/2,1/2,g), ...)
//     after adding an integer row m' (so m = m'P): fixed components must match exactly, free
//     components may hold anything. m' on a free component is not unique, and it shifts the
//     internal translations by m'.t3 (e.g. P6_2(0,0,g)h = P6_2(0,0,g)s via x4' = x4 + z), so
//     every residue modulo the translation denominators is tried.
//  5. The transformed group must follow BYU's lattice convention (no centring translation has an
//     internal component); its symbol is ssgSymbolInAxes. Among the candidates (all the same
//     group) the choice rule in `choose` picks one; it was fitted to the BYU tables (test (a):
//     every record in 3 random settings; test (b): all 8204 k-index-1 subgroup classes).
//     Known gap: for 14 orthorhombic setting pairs related by a <-> b (q along c), BYU's choice
//     follows no rule found here (e.g. it writes Pcn2 rather than Pnc2, but Pmc2_1 rather than
//     Pcm2_1), and the IT A order of settings is used.

import { gcd, matInv, matMul, Q, type QMat, type QVec } from "./rational";
import { identify, numToQ, type Op3 } from "./identify";
import { paramsOf, parseLinear, type SOp } from "./ssg";
import { crystalSystemOf, type CrystalSystem } from "./symbol";
import { ssgSymbolInAxes, type SsgSymbol } from "./ssgsymbol";
import { hnf, snf, type IMat } from "./intmat";

export interface Standardised {
  ok: true;
  sgNumber: number;
  symbol: SsgSymbol; // type symbol, q written as the family (a, b, g letters)
  qActual: string; // the group's own q in that setting, e.g. "0,0,g" or "-g,g,0"
  P: string[][];
  p: string[]; // 3D: x_std = P x + p
  eta: 1 | -1;
  m: number[];
  p4: string;
}

const Z = Q.ZERO;
const ONE = Q.ONE;
const qi = (n: number | bigint) => Q.of(BigInt(n));

// ------------------------------------------------------------------ data (IT conventions)

/**
 * (3+1)D Bravais-class q forms per crystal system and centring of the conventional cell, in
 * the axes used below (monoclinic unique axis c, q along c in orthorhombic and higher systems);
 * IT C Table 9.8.3.x. A number is a fixed component (modulo the reciprocal lattice), a letter
 * a free one.
 */
const Q_FORMS: Record<string, string[]> = {
  "triclinic P": ["a,b,g"],
  "monoclinic P": ["a,b,0", "a,b,1/2", "0,0,g", "1/2,0,g"],
  "monoclinic B": ["a,b,0", "0,0,g", "0,1/2,g"],
  "orthorhombic P": ["0,0,g", "0,1/2,g", "1/2,1/2,g"],
  "orthorhombic C": ["0,0,g", "1,0,g"],
  "orthorhombic A": ["0,0,g", "1/2,0,g"],
  "orthorhombic F": ["0,0,g", "1,0,g"],
  "orthorhombic I": ["0,0,g"],
  "tetragonal P": ["0,0,g", "1/2,1/2,g"],
  "tetragonal I": ["0,0,g"],
  "trigonal P": ["0,0,g", "1/3,1/3,g"],
  "trigonal R": ["0,0,g"],
  "hexagonal P": ["0,0,g"],
};

/**
 * From the IT A standard setting (moyo) to the axes BYU writes symbols in, in order of
 * preference. Coordinate maps x' = B x. Monoclinic: unique axis b -> c (a' = a, b' = c,
 * c' = -b; C 1 2/c 1 -> B 1 1 2/b). Orthorhombic: the six IT A settings abc, ba-c, -cba, bca,
 * cab, a-cb ("q along the standard c, then a, then b").
 */
const B_SETTINGS: Record<CrystalSystem, number[][][]> = {
  triclinic: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  monoclinic: [[[1, 0, 0], [0, 0, 1], [0, -1, 0]]],
  orthorhombic: [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]], // abc
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]], // ba-c
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]], // -cba
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]], // bca
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]], // cab
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]], // a-cb
  ],
  tetragonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  trigonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  hexagonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  cubic: [],
};

// ------------------------------------------------------------------ small exact helpers

const qmat = (A: number[][]): QMat => A.map((r) => r.map((x) => qi(x)));
const mv = (A: QMat, v: QVec): QVec => A.map((r) => r.reduce((s, a, j) => s.add(a.mul(v[j])), Z));
const vadd = (a: QVec, b: QVec): QVec => a.map((x, i) => x.add(b[i]));
const vsub = (a: QVec, b: QVec): QVec => a.map((x, i) => x.sub(b[i]));
const vkey = (v: QVec) => v.map(String).join(",");
const mkey = (A: QMat | number[][]) => A.map((r) => r.map(String).join(",")).join(";");
const isIntV = (v: QVec) => v.every((x) => x.isInt());
const imul = (A: number[][], B: number[][]) => A.map((r) => B[0].map((_, j) => r.reduce((s, x, k) => s + x * B[k][j], 0)));
const idet = (A: number[][]) =>
  A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
/** Inverse of an integer det +-1 matrix (adjugate). */
function iinv(A: number[][]): number[][] {
  const d = idet(A);
  const c = (i: number, j: number) => {
    const r = [0, 1, 2].filter((x) => x !== i), s = [0, 1, 2].filter((x) => x !== j);
    return A[r[0]][s[0]] * A[r[1]][s[1]] - A[r[0]][s[1]] * A[r[1]][s[0]];
  };
  return [0, 1, 2].map((i) => [0, 1, 2].map((j) => (((i + j) % 2 ? -1 : 1) * c(j, i)) / d + 0));
}

/** Linear form over the q parameters. */
interface Lin {
  coef: Q[];
  c: Q;
}

/** q text, e.g. "-g,g,0" or "1/2,0,g": components as linear forms; zero -> "0". */
export function formatQ(q: Lin[], names: string[]): string {
  return q
    .map(({ coef, c }) => {
      let s = "";
      coef.forEach((k, i) => {
        if (k.isZero()) return;
        const a = k.sign() < 0 ? k.neg() : k;
        s += (k.sign() < 0 ? "-" : s ? "+" : "") + (a.eq(ONE) ? "" : a.toString()) + names[i];
      });
      if (!c.isZero()) s += (c.sign() < 0 ? "-" : s ? "+" : "") + (c.sign() < 0 ? c.neg() : c).toString();
      return s || "0";
    })
    .join(",");
}

// ------------------------------------------------------------------ (3+1)D change of setting

const opKey = (o: SOp) => o.W.map((r) => r.join(",")).join(";") + "|" + vkey(o.w);
const isIdW = (W: number[][]) => W.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0)));

/**
 * The operations of a (3+1)D group (canonical list incl. centrings, translations mod 1) in the
 * setting x' = Px + p, x4' = m.x + eta x4 + p4. Null if the map does not keep the superspace
 * lattice integral (a new basis vector that is not a lattice vector, or a non-integral W').
 */
export function transformOps(ops: SOp[], P: QMat, p: QVec, eta: 1 | -1, m: QVec, p4: Q): SOp[] | null {
  const S: QMat = [...P.map((r) => [...r, Z]), [...m, qi(eta)]];
  const Pi = matInv(P);
  const mPi = [0, 1, 2].map((j) => m.reduce((s, x, i) => s.add(x.mul(Pi[i][j])), Z));
  const Si: QMat = [...Pi.map((r) => [...r, Z]), [...mPi.map((x) => x.mulInt(-eta)), qi(eta)]];
  const s = [...p, p4];
  // lattice: the old translations (Z^4 + centrings) in new coordinates must contain Z^4
  const oldCen = ops.filter((o) => isIdW(o.W)).map((o) => o.w);
  for (let j = 0; j < 4; j++) {
    const v = Si.map((r) => r[j]);
    if (!oldCen.some((c) => isIntV(vsub(v, c)))) return null;
  }
  // new centrings: closure of the images of the old lattice generators, mod 1
  const gens = [...[0, 1, 2, 3].map((j) => S.map((r) => r[j])), ...oldCen.map((c) => mv(S, c))].map((v) => v.map((x) => x.mod1()));
  const cen = new Map<string, QVec>([[vkey([Z, Z, Z, Z]), [Z, Z, Z, Z]]]);
  let frontier = [...cen.values()];
  while (frontier.length) {
    const next: QVec[] = [];
    for (const a of frontier)
      for (const g of gens) {
        const b = vadd(a, g).map((x) => x.mod1());
        const k = vkey(b);
        if (!cen.has(k)) {
          cen.set(k, b);
          next.push(b);
        }
      }
    frontier = next;
    if (cen.size > 1000) return null;
  }
  const out = new Map<string, SOp>();
  for (const o of ops) {
    const Wq = matMul(matMul(S, qmat(o.W)), Si);
    if (!Wq.every((r) => isIntV(r))) return null;
    const W = Wq.map((r) => r.map((x) => Number(x.n)));
    const w0 = vsub(vadd(mv(S, o.w), s), mv(Wq, s));
    for (const c of cen.values()) {
      const op = { W, w: vadd(w0, c).map((x) => x.mod1()) };
      out.set(opKey(op), op);
    }
  }
  const all = [...out.values()];
  const id = all.find((o) => isIdW(o.W) && o.w.every((x) => x.isZero()))!;
  return [id, ...all.filter((o) => o !== id).sort((a, b) => (opKey(a) < opKey(b) ? -1 : 1))];
}

// ------------------------------------------------------------------ lattices

/** Primitive basis (rows, conventional fractions) of Z^3 + centrings. */
function primitiveBasis(cents: QVec[]): QMat {
  let D = 1n;
  for (const v of cents) for (const x of v) D = (D * x.d) / gcd(D, x.d);
  const rows: IMat = [[D, 0n, 0n], [0n, D, 0n], [0n, 0n, D], ...cents.map((v) => v.map((x) => (x.n * D) / x.d))];
  const H = hnf(rows).H.slice(0, 3);
  return H.map((r) => r.map((x) => Q.of(x, D)));
}

/** Centring letter of a conventional cell from its centring vectors (mod 1, zero included). */
function centringLetter(cents: QVec[]): string | null {
  const set = new Set(cents.map((v) => vkey(v.map((x) => x.mod1()))));
  const h = "1/2";
  const sets: Record<string, string[]> = {
    P: ["0,0,0"],
    A: ["0,0,0", `0,${h},${h}`],
    B: ["0,0,0", `${h},0,${h}`],
    C: ["0,0,0", `${h},${h},0`],
    I: ["0,0,0", `${h},${h},${h}`],
    F: ["0,0,0", `0,${h},${h}`, `${h},0,${h}`, `${h},${h},0`],
    R: ["0,0,0", "2/3,1/3,1/3", "1/3,2/3,2/3"],
  };
  for (const [l, s] of Object.entries(sets)) if (s.length === set.size && s.every((k) => set.has(k))) return l;
  return null;
}

// ------------------------------------------------------------------ affine normaliser

/** Integer 3x3 matrices with entries in {-1,0,1} and det +-1. */
const SMALL: number[][][] = (() => {
  const out: number[][][] = [];
  const vals = [0, 1, -1];
  const rows: number[][] = [];
  for (const a of vals) for (const b of vals) for (const c of vals) if (a || b || c) rows.push([a, b, c]);
  for (const r0 of rows) for (const r1 of rows) for (const r2 of rows) {
    const M = [r0, r1, r2];
    if (Math.abs(idet(M)) === 1) out.push(M);
  }
  return out;
})();

interface Group3 {
  /** one translation per rotation (mod 1), conventional coordinates */
  reps: { R: number[][]; t: QVec }[];
  cents: QVec[];
}

interface NormEl {
  N: number[][];
  n: QVec;
}

const normCache = new Map<string, NormEl[]>();

/** Elements (N, n) of the affine normaliser of g: N from SMALL, one per coset N.K; all n mod L. */
function normaliser(g: Group3, onlyIdentity: boolean): NormEl[] {
  const key = (onlyIdentity ? "1|" : "") + g.reps.map((o) => mkey(o.R) + ":" + vkey(o.t)).sort().join(" ") + "|" + g.cents.map(vkey).sort().join(" ");
  const hit = normCache.get(key);
  if (hit) return hit;
  const Kset = new Map(g.reps.map((o) => [mkey(o.R), o]));
  const cenSet = new Set(g.cents.map((c) => vkey(c.map((x) => x.mod1()))));
  const Bp = primitiveBasis(g.cents); // rows
  const BpT = Bp[0].map((_, j) => Bp.map((r) => r[j])); // x = BpT y
  const BpTi = matInv(BpT);
  const out: NormEl[] = [];
  const seenCosets = new Set<string>();
  const cands = onlyIdentity ? [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] : SMALL;
  for (const N of cands) {
    const Ni = iinv(N);
    let ok = true;
    for (const o of g.reps) if (!Kset.has(mkey(imul(imul(N, o.R), Ni)))) { ok = false; break; }
    if (!ok) continue;
    const Nq = qmat(N);
    for (const c of g.cents) if (!cenSet.has(vkey(mv(Nq, c).map((x) => x.mod1())))) { ok = false; break; }
    if (!ok) continue;
    const coset = g.reps.map((o) => mkey(imul(N, o.R))).sort().join(" ");
    if (seenCosets.has(coset)) continue;
    // translations: (I - R') n = t_target(R') - N t  (mod L), in primitive coordinates y
    const A: IMat = [];
    const b: Q[] = [];
    for (const o of g.reps) {
      const Rp = imul(imul(N, o.R), Ni);
      const target = Kset.get(mkey(Rp))!.t;
      const Ry = matMul(matMul(BpTi, qmat(Rp)), BpT);
      const rhs = mv(BpTi, vsub(target, mv(Nq, o.t)));
      for (let i = 0; i < 3; i++) {
        A.push([0, 1, 2].map((j) => BigInt((i === j ? 1 : 0) - Number(Ry[i][j].n))));
        b.push(rhs[i]);
      }
    }
    const sols = solveMod1(A, b);
    if (!sols.length) continue;
    seenCosets.add(coset);
    const seenN = new Set<string>();
    for (const y of sols) {
      const n = mv(BpT, y).map((x) => x.mod1());
      // modulo centrings
      const canon = g.cents.map((c) => vkey(vadd(n, c).map((x) => x.mod1()))).sort()[0];
      if (seenN.has(canon)) continue;
      seenN.add(canon);
      out.push({ N, n });
    }
  }
  normCache.set(key, out);
  return out;
}

/** All y in [0,1)^n (mod Z^n) with A y = b (mod Z^m); free directions set to 0. */
function solveMod1(A: IMat, b: Q[]): QVec[] {
  const n = A[0].length;
  const { U, V, diag } = snf(A);
  const Ub = U.map((r) => r.reduce((s, x, j) => s.add(b[j].mulInt(x)), Z));
  for (let i = diag.length; i < Ub.length; i++) if (!Ub[i].isInt()) return [];
  // y = V z, d_i z_i = Ub_i + k
  let zs: QVec[] = [[]];
  for (let i = 0; i < n; i++) {
    const next: QVec[] = [];
    for (const z of zs) {
      if (i < diag.length) {
        const d = diag[i];
        for (let k = 0n; k < d; k++) next.push([...z, Ub[i].add(qi(k)).div(Q.of(d))]);
      } else next.push([...z, Z]);
    }
    zs = next;
  }
  return zs.map((z) => V.map((r) => r.reduce((s, x, j) => s.add(z[j].mulInt(x)), Z)).map((x) => x.mod1()));
}

// ------------------------------------------------------------------ q forms

/**
 * Integer row m' that brings q into the form `form`: fixed components equal exactly, free
 * components unrestricted (m' chosen to bring their constants nearest 0). Null if impossible.
 * Whether q + m' suits the superspace lattice is checked on the transformed centrings.
 */
function matchForm(q: Lin[], form: string): QVec | null {
  const f = form.split(",");
  const fixed = f.map((s) => (/^-?\d/.test(s) ? Q.parse(s) : null));
  const out: Q[] = [];
  for (let i = 0; i < 3; i++) {
    const F = fixed[i];
    if (F) {
      if (q[i].coef.some((x) => !x.isZero())) return null;
      const d = F.sub(q[i].c);
      if (!d.isInt()) return null;
      out.push(d);
    } else out.push(qi(q[i].c.add(Q.of(1n, 2n)).floor()).neg()); // constant nearest 0
  }
  return out;
}

// ------------------------------------------------------------------ main entry

const toM = (mp: QVec, P: QMat): QVec => [0, 1, 2].map((j) => [0, 1, 2].reduce((s, i) => s.add(mp[i].mul(P[i][j])), Z));
const lcmB = (a: bigint, b: bigint) => (a / gcd(a, b)) * b;
/** All integer vectors k with 0 <= k_j < den_j. */
function grid(den: number[]): number[][] {
  let out: number[][] = [[]];
  for (const d of den) out = out.flatMap((v) => Array.from({ length: d }, (_, i) => [...v, i]));
  return out;
}

interface Candidate {
  symbol: SsgSymbol;
  bIndex: number;
  P: QMat;
  p: QVec;
  eta: 1 | -1;
  m: QVec;
  qActual: string;
}

/** All candidate BYU-axes forms of the group (exported for tests/debugging). */
export function candidates(ops: SOp[], Qtext: string[]): { sgNumber: number; list: Candidate[] } | { ok: false; reason: string } {
  if (Qtext.length !== 1 || ops.some((o) => o.W.length !== 4)) return { ok: false, reason: "only d = 1 is supported for now" };
  const names = paramsOf(Qtext);
  const q: Lin[] = Qtext[0].split(",").map((s) => parseLinear(s, names));
  if (q.length !== 3) return { ok: false, reason: `bad q "${Qtext[0]}"` };

  const ops3: Op3[] = ops.map((o) => ({ R: o.W.slice(0, 3).map((r) => r.slice(0, 3)), t: o.w.slice(0, 3) }));
  const id = identify(ops3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
  if (!id.ok) return { ok: false, reason: id.reason };
  const system = crystalSystemOf(id.number);
  if (system === "cubic") return { ok: false, reason: "cubic basic groups have no incommensurate (3+1)D form" };
  // x_conv = L x_std + s  =>  x_std = T x + t,  T = L^-1, t = -T s
  const L = id.stdLinear.map((r) => r.map(numToQ));
  const T = matInv(L);
  const t = mv(T, id.stdOriginShift.map(numToQ)).map((x) => x.neg());

  // the standard basic group
  const reps = new Map<string, { R: number[][]; t: QVec }>();
  const cents: QVec[] = [];
  const Ti = L;
  for (const o of ops3) {
    const Rq = matMul(matMul(T, qmat(o.R)), Ti);
    const R = Rq.map((r) => r.map((x) => Number(x.n)));
    const tt = vsub(vadd(mv(T, o.t), t), mv(Rq, t)).map((x) => x.mod1());
    if (isIdW(R)) cents.push(tt);
    const k = mkey(R);
    if (!reps.has(k) || vkey(tt) < vkey(reps.get(k)!.t)) reps.set(k, { R, t: tt });
  }
  // the centring set may come from the ops (centrings) and from T (primitive -> conventional)
  const centSet = new Map<string, QVec>();
  for (const c of cents) centSet.set(vkey(c), c);
  const tgens = [0, 1, 2].map((j) => T.map((r) => r[j]).map((x) => x.mod1()));
  for (let grew = true; grew; ) {
    grew = false;
    for (const c of [...centSet.values()])
      for (const v of tgens) {
        const s = vadd(c, v).map((x) => x.mod1());
        if (!centSet.has(vkey(s))) {
          centSet.set(vkey(s), s);
          grew = true;
        }
      }
  }
  const g: Group3 = { reps: [...reps.values()], cents: [...centSet.values()] };
  const norm = normaliser(g, system === "triclinic");

  const list: Candidate[] = [];
  const seen = new Set<string>();
  const symbols = new Map<string, SsgSymbol | null>();
  B_SETTINGS[system].forEach((Bn, bIndex) => {
    const B = qmat(Bn);
    const bc = g.cents.map((c) => mv(B, c).map((x) => x.mod1()));
    const letter = centringLetter(bc);
    const forms = letter ? Q_FORMS[`${system} ${letter}`] : undefined;
    if (!forms) return;
    for (const { N, n } of norm) {
      const BN = matMul(B, qmat(N));
      const P = matMul(BN, T);
      const p = vadd(mv(BN, t), mv(B, n));
      const Pi = matInv(P);
      for (const eta of [1, -1] as const) {
        const q0: Lin[] = [0, 1, 2].map((j) => ({
          coef: names.map((_, k) => [0, 1, 2].reduce((s, i) => s.add(q[i].coef[k].mulInt(eta).mul(Pi[i][j])), Z)),
          c: [0, 1, 2].reduce((s, i) => s.add(q[i].c.mulInt(eta).mul(Pi[i][j])), Z),
        }));
        for (const form of forms) {
          const mp0 = matchForm(q0, form);
          if (!mp0) continue;
          const base = transformOps(ops, P, p, eta, toM(mp0, P), Z);
          if (!base) continue;
          // The free components of q' may also take any integer shift: it changes the internal
          // translations by m'.t3, so try every residue mod the translation denominators.
          const free = form.split(",").map((x) => !/^-?\d/.test(x));
          const den = [0, 1, 2].map((j) => (free[j] ? Number(base.reduce((l, o) => lcmB(l, o.w[j].d), 1n)) : 1));
          for (const k of grid(den)) {
            const mp = mp0.map((x, j) => x.add(qi(k[j])));
            const m = toM(mp, P);
            if (!isIntV(m)) continue;
            const newOps = k.every((x) => x === 0) ? base : transformOps(ops, P, p, eta, m, Z);
            // BYU convention: no centring translation involves the internal coordinate
            if (!newOps || newOps.some((o) => isIdW(o.W) && !o.w[3].isZero())) continue;
            const key = form + "#" + newOps.map(opKey).join(" ");
            if (!symbols.has(key)) symbols.set(key, ssgSymbolInAxes(newOps, [form], system));
            const symbol = symbols.get(key);
            if (!symbol) continue;
            const qActual = formatQ(q0.map((l, i) => ({ coef: l.coef, c: l.c.add(mp[i]) })), names);
            if (seen.has(key + "#" + qActual)) continue;
            seen.add(key + "#" + qActual);
            list.push({ symbol, bIndex, P, p, eta, m, qActual });
          }
        }
      }
    }
  });
  return { sgNumber: id.number, list };
}

// Phase letters by value: 0 < h (1/6) < q (1/4) < t (1/3) < s (1/2).
const VALUE = "0hqts";

/**
 * Order in which the letters are compared (indices into the shift string), fitted to the BYU
 * tables: orthorhombic compares the 2nd position first, then the 1st, then the 3rd; tetragonal
 * compares the last letter (the [110] position) first, then the rest in order; every other
 * system compares in symbol order.
 */
function letterOrder(system: CrystalSystem, n: number): number[] {
  const nat = Array.from({ length: n }, (_, i) => i);
  if (system === "orthorhombic" && n === 3) return [1, 0, 2];
  if (system === "tetragonal" && n >= 3) return [n - 1, ...nat.slice(0, n - 1)];
  return nat;
}

/**
 * Choice rule among the candidates (all are the same group in BYU axes):
 *  1. setting: lowest index in B_SETTINGS (IT A order of the orthorhombic settings);
 *  2. phase letters: smallest by value, compared in letterOrder;
 *  3. then the simplest q' and map: fewest minus signs and non-zero constants in q', eta = +1,
 *     smallest |m|.
 */
function choose(list: Candidate[], system: CrystalSystem): Candidate {
  const nz = (c: Candidate) => c.qActual.split(",").filter((x) => /\d/.test(x) && x !== "0").length;
  const cmp = (a: Candidate, b: Candidate) => {
    if (a.bIndex !== b.bIndex) return a.bIndex - b.bIndex;
    const sa = a.symbol.shifts, sb = b.symbol.shifts;
    if (sa.length === sb.length)
      for (const i of letterOrder(system, sa.length)) {
        const d = VALUE.indexOf(sa[i]) - VALUE.indexOf(sb[i]);
        if (d) return d;
      }
    if (a.symbol.full !== b.symbol.full) return a.symbol.full < b.symbol.full ? -1 : 1;
    const neg = (c: Candidate) => c.qActual.split("-").length;
    return neg(a) - neg(b) || nz(a) - nz(b) || b.eta - a.eta || a.m.reduce((s, x) => s + Math.abs(x.toNumber()), 0) - b.m.reduce((s, x) => s + Math.abs(x.toNumber()), 0);
  };
  return [...list].sort(cmp)[0];
}

/**
 * Standard (BYU-convention) setting and type of a (3+1)D group given by its operations (canonical
 * list incl. centrings) and q text in the same axes. Requires loadMoyo().
 */
export function standardise(ops: SOp[], Qtext: string[]): Standardised | { ok: false; reason: string } {
  const c = candidates(ops, Qtext);
  if ("ok" in c) return c;
  if (!c.list.length) return { ok: false, reason: "no standard setting with a (3+1)D Bravais-class q form" };
  const best = choose(c.list, crystalSystemOf(c.sgNumber));
  return {
    ok: true,
    sgNumber: c.sgNumber,
    symbol: best.symbol,
    qActual: best.qActual,
    P: best.P.map((r) => r.map(String)),
    p: best.p.map(String),
    eta: best.eta,
    m: best.m.map((x) => Number(x.n)),
    p4: "0",
  };
}
