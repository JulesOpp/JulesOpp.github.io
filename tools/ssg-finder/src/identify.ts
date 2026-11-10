// Identification of a 3D space group given as a list of operations (plan §3, milestone M4).
//
// moyo has no "type from operations" entry point, so we build a pseudo-structure that has
// exactly the given symmetry and let moyo's analyze_cell find its type:
//   1. reduce to a primitive cell of the translation lattice T = Z^3 + <centrings>,
//   2. take a metric invariant under the point group (average R^T G R if necessary),
//   3. put 4 species at one generic rational point each (1/97 grid, far from all symmetry
//      elements and from each other) and expand each under G/T,
//   4. analyze_cell, then check that moyo found exactly our operations (count AND set).
// Any inconsistency returns {ok: false, reason}; we never return an unchecked answer.

import { Q, gcd, matInv, matMul, matVec, transpose, type QMat, type QVec } from "./rational";
import { fromColMajor, getMoyo } from "./moyo";

/** A 3D operation x -> R x + t in fractional coordinates of the conventional supercell. */
export interface Op3 {
  R: number[][];
  t: Q[];
}

export type Identification =
  | {
      ok: true;
      number: number;
      hallNumber: number;
      hmShort: string;
      hmFull: string;
      /** x_conv = stdLinear · x_std + stdOriginShift (conventional-supercell fractions → moyo standard setting). */
      stdLinear: number[][];
      stdOriginShift: number[];
      /** |G/T|: number of operations modulo the full translation lattice (= moyo's count). */
      nOps: number;
      /** Number of lattice translations per conventional supercell (1 + centrings). */
      nCentring: number;
      /** Reduced primitive basis of T, rows in conventional-supercell fractions. */
      primitiveBasis: Q[][];
    }
  | { ok: false; reason: string; nOps?: number };

// ------------------------------------------------------------------ small exact helpers

/** Exact rational from a JS number (integers, or floats close to a small-denominator fraction). */
export function numToQ(x: number): Q {
  if (Number.isInteger(x)) return Q.of(BigInt(x));
  // continued fractions, max denominator 10^6
  let h0 = 0n, h1 = 1n, k0 = 1n, k1 = 0n;
  let v = x;
  for (let i = 0; i < 40; i++) {
    const a = Math.floor(v);
    const ab = BigInt(a);
    [h0, h1] = [h1, ab * h1 + h0];
    [k0, k1] = [k1, ab * k1 + k0];
    if (k1 > 1000000n) break;
    if (Math.abs(Number(h1) / Number(k1) - x) < 1e-12 * Math.max(1, Math.abs(x))) return Q.of(h1, k1);
    const f = v - a;
    if (f < 1e-15) break;
    v = 1 / f;
  }
  throw new RangeError(`cannot rationalise matrix entry ${x}`);
}

const qmatOfNum = (A: number[][]): QMat => A.map((r) => r.map(numToQ));
const qkey = (v: QVec): string => v.map((x) => x.toString()).join(",");
const mkey = (A: QMat): string => A.map(qkey).join(";");
const isIntMat = (A: QMat): boolean => A.every((r) => r.every((x) => x.isInt()));
const qdet3 = (A: QMat): Q =>
  A[0][0].mul(A[1][1].mul(A[2][2]).sub(A[1][2].mul(A[2][1])))
    .sub(A[0][1].mul(A[1][0].mul(A[2][2]).sub(A[1][2].mul(A[2][0]))))
    .add(A[0][2].mul(A[1][0].mul(A[2][1]).sub(A[1][1].mul(A[2][0]))));
const det3 = (A: number[][]): number =>
  A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
  A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
  A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
const nmul = (A: number[][], B: number[][]): number[][] =>
  A.map((r) => B[0].map((_, j) => r.reduce((s, a, k) => s + a * B[k][j], 0)));
const ntr = (A: number[][]): number[][] => A[0].map((_, j) => A.map((r) => r[j]));

/** Row-style Hermite normal form basis (3 rows) of the Z-span of integer rows of rank 3. */
function hnfRows3(rows: bigint[][]): bigint[][] {
  const M = rows.map((r) => [...r]);
  const out: bigint[][] = [];
  let r0 = 0;
  for (let c = 0; c < 3; c++) {
    // gcd-eliminate column c among rows r0..end
    for (;;) {
      let piv = -1;
      for (let i = r0; i < M.length; i++)
        if (M[i][c] !== 0n && (piv < 0 || abs(M[i][c]) < abs(M[piv][c]))) piv = i;
      if (piv < 0) throw new Error("translation lattice is not of full rank");
      [M[r0], M[piv]] = [M[piv], M[r0]];
      let done = true;
      for (let i = r0 + 1; i < M.length; i++) {
        if (M[i][c] === 0n) continue;
        const k = M[i][c] / M[r0][c];
        M[i] = M[i].map((x, j) => x - k * M[r0][j]);
        if (M[i][c] !== 0n) done = false;
      }
      if (done) break;
    }
    if (M[r0][c] < 0n) M[r0] = M[r0].map((x) => -x);
    out.push(M[r0]);
    r0++;
  }
  // reduce above-diagonal entries
  for (let c = 1; c < 3; c++)
    for (let i = 0; i < c; i++) {
      const k = fdiv(out[i][c], out[c][c]);
      if (k) out[i] = out[i].map((x, j) => x - k * out[c][j]);
    }
  return out;
}
const abs = (a: bigint): bigint => (a < 0n ? -a : a);
const fdiv = (a: bigint, b: bigint): bigint => {
  const q = a / b;
  return (a % b !== 0n && (a < 0n) !== (b < 0n)) ? q - 1n : q;
};

interface QOp {
  W: QMat;
  w: QVec;
}

/** Express ops in a new basis whose rows P are given in old coordinates: x_old = P^T x_new. */
function changeBasis(ops: QOp[], P: QMat): QOp[] {
  const Pt = transpose(P);
  const Pti = matInv(Pt);
  return ops.map(({ W, w }) => ({ W: matMul(matMul(Pti, W), Pt), w: matVec(Pti, w) }));
}

/** Pairwise (Gauss) reduction of a basis w.r.t. Gram matrix G; returns unimodular U (rows). */
function reduceBasis(G0: number[][]): number[][] {
  let U = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  let G = G0.map((r) => [...r]);
  for (let iter = 0; iter < 200; iter++) {
    let changed = false;
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++) {
        if (i === j) continue;
        const k = Math.round(G[i][j] / G[j][j]);
        if (k === 0) continue;
        // b_i -= k b_j only if it strictly shortens b_i
        const newLen = G[i][i] - 2 * k * G[i][j] + k * k * G[j][j];
        if (newLen >= G[i][i] * (1 - 1e-12)) continue;
        const E = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        E[i][j] = -k;
        U = nmul(E, U);
        G = nmul(nmul(E, G), ntr(E));
        changed = true;
      }
    if (!changed) break;
  }
  if (det3(U) < 0) U = [U[0].map((x) => -x), U[1], U[2]];
  return U;
}

function cholesky(G: number[][]): number[][] {
  const L = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j <= i; j++) {
      let s = G[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 0)) throw new Error("metric is not positive definite");
        L[i][i] = Math.sqrt(s);
      } else L[i][j] = s / L[j][j];
    }
  return L;
}

/** Deterministic pseudo-random integers (LCG), so identification is reproducible. */
function* lcg(seed: number): Generator<number> {
  let s = seed >>> 0;
  for (;;) {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    yield s;
  }
}

// ------------------------------------------------------------------ main entry

/**
 * Identify the space-group type of `ops` (conventional-supercell fractions, centrings
 * included as R = I operations) with Cartesian supercell basis `basis` (rows a, b, c).
 * Requires loadMoyo() to have completed.
 */
export function identify(ops: Op3[], basis: number[][]): Identification {
  try {
    return identifyInner(ops, basis);
  } catch (e) {
    return { ok: false, reason: `identification failed: ${e instanceof Error ? e.message : String(e)}` };
  }
}

function identifyInner(ops: Op3[], basis: number[][]): Identification {
  if (!ops.length) return { ok: false, reason: "no operations" };
  // 0. exact ops, de-duplicated modulo Z^3
  const seen = new Map<string, QOp>();
  for (const o of ops) {
    const W = qmatOfNum(o.R);
    const w = o.t.map((x) => x.mod1());
    const k = mkey(W) + "|" + qkey(w);
    if (!seen.has(k)) seen.set(k, { W, w });
  }
  const conv = [...seen.values()];

  // 1. translation lattice T = Z^3 + centrings, primitive basis by HNF
  const I3 = mkey(qmatOfNum([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
  const cents = conv.filter((o) => mkey(o.W) === I3).map((o) => o.w);
  if (!cents.length) return { ok: false, reason: "identity operation missing" };
  let D = 1n;
  for (const v of cents) for (const x of v) D = (D * x.d) / gcd(D, x.d);
  const gens: bigint[][] = [[D, 0n, 0n], [0n, D, 0n], [0n, 0n, D]];
  for (const v of cents) gens.push(v.map((x) => (x.n * D) / x.d));
  const H = hnfRows3(gens);
  const B: QMat = H.map((r) => r.map((x) => Q.of(x, D)));
  const detB = qdet3(B); // = 1 / (number of lattice points per conventional cell)
  const nCentring = Number(Q.ONE.div(detB).n);
  if (!Q.ONE.div(detB).isInt()) return { ok: false, reason: "centring set is not a lattice" };
  if (cents.length !== nCentring)
    return { ok: false, reason: `centring translations not closed: ${cents.length} given, lattice has ${nCentring}` };

  // 2. point group and coset representatives in the HNF primitive basis
  let prim = changeBasis(conv, B);
  const reps = new Map<string, QOp>();
  for (const o of prim) {
    if (!isIntMat(o.W)) return { ok: false, reason: "a rotation does not preserve the translation lattice" };
    const k = mkey(o.W);
    const r = reps.get(k);
    if (!r) reps.set(k, { W: o.W, w: o.w.map((x) => x.mod1()) });
    else if (!o.w.every((x, i) => x.sub(r.w[i]).isInt()))
      return { ok: false, reason: "two operations with the same rotation differ by a non-lattice translation (centring missing)" };
  }
  const nOps = reps.size;
  if (conv.length !== nOps * nCentring)
    return { ok: false, reason: `operation count ${conv.length} != ${nOps} cosets x ${nCentring} centrings`, nOps };
  if (nOps > 48) return { ok: false, reason: `point group of order ${nOps} > 48`, nOps };
  // closure modulo T
  const repList = [...reps.values()];
  for (const a of repList)
    for (const b of repList) {
      const W = matMul(a.W, b.W);
      const r = reps.get(mkey(W));
      if (!r) return { ok: false, reason: "operations are not closed under composition (rotation part)", nOps };
      const w = matVec(a.W, b.w).map((x, i) => x.add(a.w[i]));
      if (!w.every((x, i) => x.sub(r.w[i]).isInt()))
        return { ok: false, reason: "operations are not closed under composition (translation part)", nOps };
      if (!qdet3(a.W).isInt() || Math.abs(qdet3(a.W).toNumber()) !== 1)
        return { ok: false, reason: "rotation with determinant other than +-1", nOps };
    }

  // 3. metric of the primitive cell, invariant under the point group
  const Bn = B.map((r) => r.map((x) => x.toNumber()));
  const Ap = nmul(Bn, basis); // rows: primitive vectors, Cartesian
  let G = nmul(Ap, ntr(Ap));
  const Wn = repList.map((o) => o.W.map((r) => r.map((x) => x.toNumber())));
  const scale = Math.max(G[0][0], G[1][1], G[2][2]);
  const invariant = Wn.every((W) => {
    const G2 = nmul(nmul(ntr(W), G), W);
    return G2.every((r, i) => r.every((x, j) => Math.abs(x - G[i][j]) <= 1e-9 * scale));
  });
  if (!invariant) {
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const W of Wn) {
      const G2 = nmul(nmul(ntr(W), G), W);
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) S[i][j] += G2[i][j] / Wn.length;
    }
    G = S;
  }

  // 4. reduce the primitive basis, re-express everything in it
  const U = reduceBasis(G);
  const Uq = U.map((r) => r.map((x) => Q.of(BigInt(x))));
  const B2 = matMul(Uq, B); // rows in conventional fractions
  G = nmul(nmul(U, G), ntr(U));
  prim = changeBasis(repList, Uq).map((o) => ({ W: o.W, w: o.w.map((x) => x.mod1()) }));
  let L = cholesky(G); // rows: Cartesian primitive vectors (right-handed)
  if (det3(basis) < 0) L = L.map((r) => [-r[0], r[1], r[2]]); // keep the input handedness
  const vol = Math.abs(det3(L));
  const Lc = Math.cbrt(vol);
  const symprec = 1e-4 * Lc;

  // 5. pseudo-structure: 4 species at generic points on a 1/97 grid
  const nAtoms = 4 * nOps;
  const dmin = Math.max(0.08 * Math.cbrt(vol / nAtoms), 200 * symprec);
  const cart = (f: number[]): number[] => [0, 1, 2].map((k) => f[0] * L[0][k] + f[1] * L[1][k] + f[2] * L[2][k]);
  const minImageDist = (a: number[], b: number[]): number => {
    const d0 = a.map((x, i) => x - b[i]).map((x) => x - Math.round(x));
    let best = Infinity;
    for (let i = -1; i <= 1; i++)
      for (let j = -1; j <= 1; j++)
        for (let k = -1; k <= 1; k++) {
          const c = cart([d0[0] + i, d0[1] + j, d0[2] + k]);
          const r = Math.hypot(c[0], c[1], c[2]);
          if (r < best) best = r;
        }
    return best;
  };
  const positions: number[][] = [];
  const numbers: number[] = [];
  const centres: number[][] = [];
  const rnd = lcg(0x5eed + nOps);
  let tries = 0;
  while (centres.length < 4) {
    if (++tries > 20000) return { ok: false, reason: "could not place generic points", nOps };
    const p = [0, 1, 2].map(() => Q.of(BigInt((rnd.next().value as number) % 96) + 1n, 97n));
    const orbit = prim.map((o) => matVec(o.W, p).map((x, i) => x.add(o.w[i]).mod1().toNumber()));
    // generic + well separated: all images of p distinct and far from each other and from earlier atoms
    let good = true;
    for (let i = 0; i < orbit.length && good; i++) {
      for (let j = 0; j < i && good; j++) if (minImageDist(orbit[i], orbit[j]) < dmin) good = false;
      for (const q of positions) if (good && minImageDist(orbit[i], q) < dmin) good = false;
    }
    if (!good) continue;
    const pc = cart(p.map((x) => x.toNumber()));
    if (centres.length === 3) {
      const d = [0, 1, 2].map((i) => [0, 1, 2].map((k) => centres[i][k] - pc[k]));
      if (Math.abs(det3(d)) < 1e-3 * vol) continue; // non-coplanar
    }
    centres.push(pc);
    for (const x of orbit) {
      positions.push(x);
      numbers.push(centres.length);
    }
  }

  // 6. moyo
  const moyo = getMoyo();
  const cell = {
    lattice: { basis: [...L[0], ...L[1], ...L[2]] },
    positions,
    numbers,
  };
  let ds;
  try {
    ds = moyo.analyze_cell(JSON.stringify(cell), symprec, "Standard");
  } catch (e) {
    return { ok: false, reason: `moyo analyze_cell failed: ${String(e)}`, nOps };
  }
  if (ds.operations.length !== nOps)
    return { ok: false, reason: `moyo found ${ds.operations.length} operations, expected ${nOps}`, nOps };
  // each moyo operation must be one of ours (rotation exact, translation mod 1)
  const ours = new Map(prim.map((o) => [o.W.map((r) => r.map((x) => x.toNumber()).join(",")).join(";"), o.w]));
  for (const op of ds.operations) {
    const R = fromColMajor(op.rotation);
    const w = ours.get(R.map((r) => r.join(",")).join(";"));
    if (!w) return { ok: false, reason: "moyo found a rotation that is not in the group", nOps };
    const dw = op.translation.map((x, i) => x - w[i].toNumber());
    if (dw.some((x) => Math.abs(x - Math.round(x)) > 1e-6))
      return { ok: false, reason: "moyo found a translation that is not in the group", nOps };
  }
  const entry = moyo.hall_symbol_entry(ds.hall_number);
  if (entry.number !== ds.number) return { ok: false, reason: "inconsistent Hall entry", nOps };

  // 7. transformation to the standard setting, relative to the conventional supercell
  const Pm = fromColMajor(ds.std_linear); // x_prim = Pm x_std + p
  const B2t = ntr(B2.map((r) => r.map((x) => x.toNumber()))); // x_conv = B2^T x_prim
  const stdLinear = nmul(B2t, Pm).map((r) => r.map(clean));
  const stdOriginShift = B2t.map((r) => clean(r.reduce((s, a, k) => s + a * ds.std_origin_shift[k], 0)));

  return {
    ok: true,
    number: ds.number,
    hallNumber: ds.hall_number,
    hmShort: entry.hm_short,
    hmFull: entry.hm_full,
    stdLinear,
    stdOriginShift,
    nOps,
    nCentring,
    primitiveBasis: B2,
  };
}

/** Snap floating noise to the nearest multiple of 1/48 when within 1e-8 (keeps -0 out). */
function clean(x: number): number {
  const r = Math.round(x * 48) / 48;
  const y = Math.abs(x - r) < 1e-8 ? r : x;
  return y === 0 ? 0 : y;
}
