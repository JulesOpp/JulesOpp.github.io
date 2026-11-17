// Deterministic canonical setting and type name of a (3+d)D superspace group (d = 1, 2, 3).
// d = 1 follows BYU conventions; d = 2, 3 use this tool's own canonical setting.
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
//     group) the choice rule in `rank` picks one; it was fitted to the BYU tables (test (a):
//     every record in 3 random settings; test (b): all 8204 k-index-1 subgroup classes).
//     Ties: with q along c, an orthorhombic group may be written equally well in two settings
//     related by a <-> b (e.g. Pnc2 / Pcn2). The primary name uses the IT A Table 4.3.2.1 order of
//     settings; the other names are returned in `alternates`. BYU's own choice between such
//     names follows no rule found here, so callers should accept any of primary + alternates.

// The {-1,0,1} normaliser search is a finite conventional-setting search: signed axial
// permutations and the hexagonal-axis changes needed by the supported d = 1 tables fit
// this bound. Low-symmetry point groups can have unbounded integral normalisers; this is
// not an enumeration of their full normaliser and gives no completeness guarantee there.
//
// For k = 1 we search with p4 = 0: the integer row m is varied through its translation
// denominator residues, which supplies the internal-phase choices used by this d = 1
// convention. Internal-origin shifts of eps = -1 elements do not affect their printed
// zero letter. This suffices in the tested conventional settings, not a general proof
// of equivalence for arbitrary superspace lattices or higher d.
//
// The letter-position comparison order is empirical, fixed by the 775-record d = 1
// sweep. It is a naming convention fitted to those data, not a derived universal order.

import { standardise as legacyStandardise } from "./standardise-legacy";
import { gcd, toSafeInt, matInv, Q, type QMat, type QVec } from "./rational";
import { identify, numToQ, type Op3 } from "./identify";
import { paramsOf, parseLinear, type SOp } from "./ssg";
import { crystalSystemOf, type CrystalSystem } from "./symbol";
import { ssgSymbolInAxes, type SsgSymbol } from "./ssgsymbol";
import { hnf, snf, imatDet, type IMat } from "./intmat";

export interface Standardised {
  ok: true;
  /** Exact compatibility result from the frozen pre-Phase-D writer. */
  namingFallback?: "d244096";
  sgNumber: number;
  symbol: SsgSymbol; // type symbol, q written as the family (a, b, g letters)
  qActual: string[]; // the group's own q in that setting, e.g. "0,0,g" or "-g,g,0"
  P: string[][];
  p: string[]; // 3D: x_std = P x + p
  E: number[][];
  M: string[][];
  p4: string[];
  /** Full symbols of the same type in the other tied settings (sorted, without the primary). */
  alternates: string[];
}

export interface Fail { ok: false; reason: string }

export function isFail(x: Standardised | Fail): x is Fail {
  return x.ok === false;
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
 * c' = -b; C 1 2/c 1 -> B 1 1 2/b). Orthorhombic: the six IT A settings in the column order of
 * IT A Table 4.3.2.1: abc, ba-c, cab, -cba, bca, a-cb. This order decides the primary name when
 * two settings give the same type (see `alternates`).
 */
const B_SETTINGS: Record<CrystalSystem, number[][][]> = {
  triclinic: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  monoclinic: [[[1, 0, 0], [0, 0, 1], [0, -1, 0]]],
  orthorhombic: [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]], // abc
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]], // ba-c
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]], // cab
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]], // -cba
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]], // bca
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]], // a-cb
  ],
  tetragonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  trigonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  hexagonal: [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  cubic: [],
};

// ------------------------------------------------------------------ small exact helpers

const qmat = (A: number[][]): QMat => A.map((r) => r.map((x) => qi(x)));
// Exact sparse products avoid allocating rationals for zero summands.
const matMul = (A:QMat,B:QMat):QMat => A.map(row=>Array.from({length:B[0]?.length??0},(_,j)=>row.reduce((s,a,k)=>a.isZero()||B[k][j].isZero()?s:s.add(a.mul(B[k][j])),Z)));
const mv = (A: QMat, v: QVec): QVec => A.map((r) => r.reduce((s, a, j) => a.isZero()||v[j].isZero()?s:s.add(a.mul(v[j])), Z));
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
export function transformOps(ops: SOp[], P: QMat, p: QVec, E: number[][], M: QMat, p4: QVec): SOp[] | null {
  const d = E.length, n = 3 + d;
  const S = [...P.map(r => [...r, ...Array(d).fill(Z)]), ...M.map((r,i) => [...r, ...E[i].map(qi)])];
  const Si = matInv(S), shift = [...p, ...p4];
  const oldCen = ops.filter(o => isIdW(o.W)).map(o => o.w);
  for (let j=0;j<n;j++) if (!oldCen.some(c => isIntV(vsub(Si.map(r=>r[j]),c)))) return null;
  const gens = [...Array.from({length:n},(_,j)=>S.map(r=>r[j])), ...oldCen.map(c=>mv(S,c))].map(v=>v.map(x=>x.mod1()));
  const zero = Array(n).fill(Z), cen = new Map([[vkey(zero),zero]]);
  for(let frontier=[zero];frontier.length;) {
    const next: QVec[]=[];
    for(const a of frontier) for(const g of gens) {const b=vadd(a,g).map(x=>x.mod1()), k=vkey(b); if(!cen.has(k)){cen.set(k,b);next.push(b);}}
    frontier=next; if(cen.size>1000)return null;
  }
  const out = new Map<string,SOp>();
  for(const o of ops) {
    const Wq=matMul(matMul(S,qmat(o.W)),Si); if(!Wq.every(isIntV))return null;
    const W=Wq.map(r=>r.map(x=>toSafeInt(x,'operator matrix entry'))), w0=vsub(vadd(mv(S,o.w),shift),mv(Wq,shift));
    for(const c of cen.values()) {const op={W,w:vadd(w0,c).map(x=>x.mod1())};out.set(opKey(op),op);}
  }
  const all=[...out.values()],id=all.find(o=>isIdW(o.W)&&o.w.every(x=>x.isZero()))!;
  return [id,...[...out.entries()].filter(([,o])=>o!==id).sort(([a],[b])=>a<b?-1:a>b?1:0).map(([,o])=>o)];
}

/** Factor a second, integral internal change after a physical conjugation.
 * y=P x+p, s=U t; s'=B s+mp y-mp p, so the composed internal
 * origin is exactly zero. This is an automorphism of Z^(3+d): no new lattice
 * closure is needed. Integer matrix products use BigInt throughout. */
function internalTransformer(ops:SOp[],p:QVec,resultCache?:Map<string,Map<string,SOp[]>>):(B:number[][],mp:QMat)=>SOp[] {
  // Cache only exact operation sets and the full internal affine change.
  // Eight groups with 256 entries each bound memory within one naming call.
  const baseKey=resultCache?ops.map(opKey).join(' '):'';
  const results=resultCache?.get(baseKey)??new Map<string,SOp[]>();
  if(resultCache&&!resultCache.has(baseKey)){
    if(resultCache.size>=8)resultCache.delete(resultCache.keys().next().value!);
    resultCache.set(baseKey,results);
  }
  const d=ops[0].W.length-3;
  const mul=(A:IMat,B:IMat):IMat=>A.map(row=>B[0].map((_,j)=>row.reduce((s,x,k)=>s+x*B[k][j],0n)));
  const data=ops.map(o=>({o,R:o.W.slice(0,3).map(r=>r.slice(0,3).map(BigInt)),F:o.W.slice(3).map(r=>r.slice(0,3).map(BigInt)),D:o.W.slice(3).map(r=>r.slice(3).map(BigInt))}));
  const inverses=new Map<string,IMat>();
  return (basis,mp)=>{
    const shift=mv(mp,p).map(x=>x.neg()),key=mkey(basis),resultKey=key+'|'+mkey(mp)+'|'+vkey(shift);
    const hit=results.get(resultKey);if(hit)return hit;
    const B=basis.map(r=>r.map(BigInt));
    let Bi=inverses.get(key);if(!Bi){Bi=matInv(qmat(basis)).map(r=>r.map(x=>BigInt(toSafeInt(x,'inverse internal basis'))));inverses.set(key,Bi);}
    const m=mp.map(r=>r.map(x=>BigInt(toSafeInt(x,'internal residue'))));
    const out=data.map(({o,R,F,D})=>{
      const Dn=mul(mul(B,D),Bi!),BF=mul(B,F),mR=mul(m,R),Dm=mul(Dn,m);
      const Fn=BF.map((r,i)=>r.map((x,j)=>x+mR[i][j]-Dm[i][j]));
      const W=[...o.W.slice(0,3).map(r=>[...r.slice(0,3),...Array(d).fill(0)]),...Fn.map((r,i)=>[...r,...Dn[i]].map(x=>toSafeInt(x,'operator matrix entry')))];
      const wI=vsub(vadd(vadd(mv(qmat(basis),o.w.slice(3)),mv(mp,o.w.slice(0,3))),shift),mv(Dn.map(r=>r.map(qi)),shift));
      return {W,w:[...o.w.slice(0,3),...wI.map(x=>x.mod1())]};
    });
    const id=out.find(o=>isIdW(o.W)&&o.w.every(x=>x.isZero()))!;
    const sorted=[id,...out.filter(o=>o!==id).map(o=>({o,key:opKey(o)})).sort((a,b)=>a.key<b.key?-1:a.key>b.key?1:0).map(x=>x.o)];
    if(results.size>=256)results.clear();results.set(resultKey,sorted);
    return sorted;
  };
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
function normaliser(g: Group3, onlyIdentity: boolean, T: QMat, anchor?:number[][]): NormEl[] {
  const key = T.map((r) => r.map(String).join(",")).join(";") + "|" + (onlyIdentity ? "1|" : "") + g.reps.map((o) => mkey(o.R) + ":" + vkey(o.t)).sort().join(" ") + "|" + g.cents.map(vkey).sort().join(" ") + (anchor ? "|anchor:"+mkey(anchor) : "");
  const hit = normCache.get(key);
  if (hit) return hit;
  const Kset = new Map(g.reps.map((o) => [mkey(o.R), o]));
  const cenSet = new Set(g.cents.map((c) => vkey(c.map((x) => x.mod1()))));
  const Bp = primitiveBasis(g.cents); // rows
  const BpT = Bp[0].map((_, j) => Bp.map((r) => r[j])); // x = BpT y
  const BpTi = matInv(BpT);
  const out: NormEl[] = [];
  const bounded = onlyIdentity ? [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] : SMALL;
  const cands = anchor ? bounded.map(N=>imul(N,anchor)) : bounded;
  // Each coset N.K is represented by its simplest member (identity first, then fewest off-diagonal
  // entries, then fewest minus signs); cosets keep the order in which SMALL first meets them.
  // Simplicity is judged on the map N.T that is actually applied to the given axes (T: given -> standard).
  const cost = (N: number[][]) => {
    const M = matMul(qmat(N), T);
    let off = 0, neg = 0;
    M.forEach((r, i) => r.forEach((x, j) => { if (i !== j && !x.isZero()) off++; if (i === j && x.isZero()) off++; if (x.n < 0n) neg++; }));
    return off * 100 + neg;
  };
  const repOf = new Map<string, number[][]>();
  const order: string[] = [];
  for (const N of cands) {
    const Ni = iinv(N);
    let ok = true;
    for (const o of g.reps) if (!Kset.has(mkey(imul(imul(N, o.R), Ni)))) { ok = false; break; }
    if (!ok) continue;
    const Nq = qmat(N);
    for (const c of g.cents) if (!cenSet.has(vkey(mv(Nq, c).map((x) => x.mod1())))) { ok = false; break; }
    if (!ok) continue;
    const coset = g.reps.map((o) => mkey(imul(N, o.R))).sort().join(" ");
    const cur = repOf.get(coset);
    if (!cur) order.push(coset);
    if (!cur || cost(N) < cost(cur)) repOf.set(coset, N);
  }
  for (const coset of order) {
    const N = repOf.get(coset)!;
    const Ni = iinv(N);
    const Nq = qmat(N);
    // translations: (I - R') n = t_target(R') - N t  (mod L), in primitive coordinates y
    const A: IMat = [];
    const b: Q[] = [];
    for (const o of g.reps) {
      const Rp = imul(imul(N, o.R), Ni);
      const target = Kset.get(mkey(Rp))!.t;
      const Ry = matMul(matMul(BpTi, qmat(Rp)), BpT);
      const rhs = mv(BpTi, vsub(target, mv(Nq, o.t)));
      for (let i = 0; i < 3; i++) {
        A.push([0, 1, 2].map((j) => BigInt((i === j ? 1 : 0) - toSafeInt(Ry[i][j], "primitive rotation entry"))));
        b.push(rhs[i]);
      }
    }
    const sols = solveMod1(A, b);
    if (!sols.length) continue;
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

/** Complexity of the affine map (smaller is simpler): identity first, then fewest non-zero and negative entries. */
function mapCost(c: Candidate): number {
  const nzv = (v: Q[]) => v.filter((x) => !x.isZero()).length;
  const offP = c.P.reduce((s, r, i) => s + r.filter((x, j) => i !== j && !x.isZero()).length, 0);
  const negP = c.P.reduce((s, r) => s + r.filter((x) => x.n < 0n).length, 0);
  const diagP = c.P.reduce((s, r, i) => s + (r[i].isZero() ? 1 : 0), 0);
  return (offP + diagP) * 1e6 + negP * 1e4 + nzv(c.p) * 100 + nzv(c.m) * 10 + (c.eta === 1 ? 0 : 1);
}

/** All candidate BYU-axes forms of the group (exported for tests/debugging). */
export function candidates(ops: SOp[], Qtext: string[]): { sgNumber: number; list: Candidate[] } | Fail {
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
    const R = Rq.map((r) => r.map((x) => toSafeInt(x, "operator matrix entry")));
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
  const norm = normaliser(g, system === "triclinic", T);

  const list: Candidate[] = [];
  const seen = new Map<string, number>();
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
          const base = transformOps(ops, P, p, [[eta]], [toM(mp0, P)], [Z]);
          if (!base) continue;
          // The free components of q' may also take any integer shift: it changes the internal
          // translations by m'.t3, so try every residue mod the translation denominators.
          const free = form.split(",").map((x) => !/^-?\d/.test(x));
          const den = [0, 1, 2].map((j) => (free[j] ? toSafeInt(base.reduce((l, o) => lcmB(l, o.w[j].d), 1n), "translation denominator") : 1));
          for (const k of grid(den)) {
            const mp = mp0.map((x, j) => x.add(qi(k[j])));
            const m = toM(mp, P);
            if (!isIntV(m)) continue;
            const newOps = k.every((x) => x === 0) ? base : transformOps(ops, P, p, [[eta]], [m], [Z]);
            // BYU convention: no centring translation involves the internal coordinate
            if (!newOps || newOps.some((o) => isIdW(o.W) && !o.w[3].isZero())) continue;
            const key = form + "#" + newOps.map(opKey).sort().join(" "); // the operation set, order-independent
            if (!symbols.has(key)) symbols.set(key, ssgSymbolInAxes(newOps, [form], system));
            const symbol = symbols.get(key);
            if (!symbol) continue;
            const qActual = formatQ(q0.map((l, i) => ({ coef: l.coef, c: l.c.add(mp[i]) })), names);
            const cand: Candidate = { symbol, bIndex, P, p, eta, m, qActual };
            const sk = key + "#" + qActual;
            const prev = seen.get(sk);
            if (prev !== undefined) {
              // same operations in the same axes: keep the simplest affine map
              if (mapCost(cand) < mapCost(list[prev])) list[prev] = cand;
              continue;
            }
            seen.set(sk, list.length);
            list.push(cand);
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
 *  1. setting: lowest index in B_SETTINGS (IT A Table 4.3.2.1 order of the orthorhombic settings);
 *  2. phase letters: smallest by value, compared in letterOrder;
 *  3. then the simplest q' and map: fewest minus signs and non-zero constants in q', eta = +1,
 *     smallest |m|.
 */
function rank(list: Candidate[], system: CrystalSystem): Candidate[] {
  const nz = (c: Candidate) => c.qActual.split(",").filter((x) => /\d/.test(x) && x !== "0").length;
  const cmp = (a: Candidate, b: Candidate) => {
    if (a.bIndex !== b.bIndex) return a.bIndex - b.bIndex;
    const sa = a.symbol.shifts[0], sb = b.symbol.shifts[0];
    if (sa.length === sb.length)
      for (const i of letterOrder(system, sa.length)) {
        const d = VALUE.indexOf(sa[i]) - VALUE.indexOf(sb[i]);
        if (d) return d;
      }
    if (a.symbol.full !== b.symbol.full) return a.symbol.full < b.symbol.full ? -1 : 1;
    const neg = (c: Candidate) => c.qActual.split("-").length;
    const magnitude = (c: Candidate) => c.m.reduce((s, x) => s.add(x.sign() < 0 ? x.neg() : x), Z);
    return neg(a) - neg(b) || nz(a) - nz(b) || b.eta - a.eta || magnitude(a).cmp(magnitude(b));
  };
  return [...list].sort(cmp);
}

/**
 * Standard (BYU-convention) setting and type of a (3+1)D group given by its operations (canonical
 * list incl. centrings) and q text in the same axes. Requires loadMoyo().
 */
export function standardise(ops: SOp[], Qtext: string[]): Standardised | Fail {
  if (Qtext.length !== 1) return standardiseTypeFamily(ops, Qtext);
  const c = candidates(ops, Qtext);
  if ("ok" in c) return c;
  if (!c.list.length) return { ok: false, reason: "no standard setting with a (3+1)D Bravais-class q form" };
  const ranked = rank(c.list, crystalSystemOf(c.sgNumber));
  const best = ranked[0];
  // the best name in every other setting that writes the same type with the same q form
  const firstPerSetting = new Map<number, Candidate>();
  for (const x of ranked) if (!firstPerSetting.has(x.bIndex)) firstPerSetting.set(x.bIndex, x);
  const alternates = [...new Set([...firstPerSetting.values()]
    .filter((x) => x.symbol.q === best.symbol.q && x.symbol.full !== best.symbol.full)
    .map((x) => x.symbol.full))].sort();
  return {
    ok: true,
    sgNumber: c.sgNumber,
    symbol: best.symbol,
    qActual: [best.qActual],
    P: best.P.map((r) => r.map(String)),
    p: best.p.map(String),
    E: [[best.eta]],
    M: [best.m.map(String)],
    p4: ["0"],
    alternates,
  };
}

/** Anchor the monoclinic ac-plane shear to its coefficient lattice instead
 * of the supplied axes. Internal row HNF first removes input GL(d,Z). Physical
 * row HNF then reduces the two equal-character coordinate directions. Its
 * unimodular map may have arbitrarily large entries; the subsequent finite
 * normaliser representatives are composed with this exact anchor. */
function monoclinicAnchor(qs:Lin[][],names:string[],Ti:QMat):number[][] {
  const coeff=qs.map(row=>names.flatMap((_,k)=>[0,1,2].map(j=>row.reduce((v,x,h)=>v.add(x.coef[k].mul(Ti[h][j])),Z))));
  const den=coeff.flat().reduce((v,x)=>lcmB(v,x.d),1n);
  const H=hnf(coeff.map(row=>row.map(x=>x.n*(den/x.d)))).H;
  const axes=[0,2],columns=axes.map(j=>names.flatMap((_,k)=>H.map(row=>row[k*3+j])));
  const U=hnf(columns).U.map(row=>row.map(x=>toSafeInt(x,'physical coefficient anchor')));
  const A=[[U[0][0],0,U[0][1]],[0,1,0],[U[1][0],0,U[1][1]]],Ai=matInv(qmat(A));
  return Ai[0].map((_,j)=>Ai.map(row=>toSafeInt(row[j],'physical coefficient anchor')));
}

// Pass-3 higher-dimensional naming rule; BYU agreement is measured, not required.
/** Solve e H = target exactly in the full-rank row lattice. */
function coordinateSolver(H:QMat):(target:QVec)=>QVec|null {
  if(!H.length)return target=>target.every(x=>x.isZero())?[]:null;
  const d=H.length, pivots:number[]=[],A=H.map(r=>[...r]);let rank=0;
  for(let j=0;j<H[0].length&&rank<d;j++) {
    const k=A.findIndex((r,i)=>i>=rank&&!r[j].isZero());if(k<0)continue;
    [A[rank],A[k]]=[A[k],A[rank]];const v=A[rank][j];A[rank]=A[rank].map(x=>x.div(v));
    for(let i=rank+1;i<d;i++){const f=A[i][j];A[i]=A[i].map((x,h)=>x.sub(f.mul(A[rank][h])));}pivots.push(j);rank++;
  }
  if(rank!==d)return ()=>null;
  const inv=matInv(H.map(r=>pivots.map(j=>r[j])));
  return target=>{
    const e=H.map((_,j)=>pivots.reduce((v,k,i)=>v.add(target[k].mul(inv[i][j])),Z));
    if(target.some((x,j)=>!x.eq(H.reduce((v,r,i)=>v.add(r[j].mul(e[i])),Z))))return null;
    return e;
  };
}

/** Plain seed/orbit and independent-row solutions, without an entry bound on E.
 * Each seed coordinate is zero or one named parameter, with partial supports. Solve its
 * coordinates in the HNF lattice, then solve all point-group images in that same
 * lattice and retain every unimodular ordered tuple. */
function orbitBases(H:QMat,names:string[],rotations:number[][][]):number[][][] {
  const d=H.length,out=new Map<string,number[][]>();
  // H is fixed throughout this search: factor it once, then verify every
  // target in all columns. No tolerance or approximate lattice solve.
  const coordinates=coordinateSolver(H);
  const solve=(target:QVec)=>{
    const e=coordinates(target);
    return !e||e.some(x=>!x.isInt())?null:e.map(x=>toSafeInt(x,'internal lattice coordinate'));
  };
  const targets:QVec[]=[];
  for(let a=-1;a<names.length;a++)for(let b=-1;b<names.length;b++)for(let c=-1;c<names.length;c++) {
    const axes=[a,b,c];if(axes.every(k=>k<0))continue;
    targets.push(names.flatMap((_,k)=>axes.map(v=>v===k?ONE:Z)));
  }
  const independent=new Map<string,number[]>();
  for(const seed of targets) {
    const first=solve(seed);if(!first)continue;
    independent.set(first.join(','),first);
    const rows=new Map<string,number[]>();rows.set(first.join(','),first);
    for(const R of rotations) {
      const target=names.flatMap((_,k)=>[0,1,2].map(j=>[0,1,2].reduce((v,h)=>v.add(seed[k*3+h].mulInt(R[h][j])),Z)));
      const row=solve(target);if(row)rows.set(row.join(','),row);
    }
    const visit=(A:number[][])=>{if(A.length===d){if([1n,-1n].includes(imatDet(A.map(r=>r.map(BigInt)))))out.set(mkey(A),A);return;}for(const row of rows.values())visit([...A,row]);};
    visit([first]);
  }
  // Independent vectors may omit parameters belonging to an earlier star.
  // Solve in the exact HNF lattice; entries of the resulting GL map are unbounded.
  const visit=(A:number[][])=>{
    if(A.length===d){if([1n,-1n].includes(imatDet(A.map(r=>r.map(BigInt)))))out.set(mkey(A),A);return;}
    for(const row of independent.values())visit([...A,row]);
  };
  visit([]);
  return [...out.values()];
}

const reciprocalOffsetCache=new Map<string,QVec>();
/** Canonical representative modulo the reciprocal lattice of Z^3 + centrings.
 * Order: L1 absolute size, componentwise absolute sizes, [0,1) membership,
 * then nonnegative signs and numeric order. The axis periods below belong to
 * the reciprocal lattice, so the finite nearest-coordinate box is exhaustive. */
export function reciprocalOffsetRepresentative(offset:QVec,cents:QVec[]):QVec {
  const periods=[0,1,2].map(j=>qi(cents.reduce((v,c)=>lcmB(v,c[j].d),1n)));
  const residue=offset.map((x,j)=>x.modQ(periods[j]));
  const key=cents.map(c=>vkey(c.map(x=>x.mod1()))).sort().join(';')+'|'+vkey(residue);
  const cached=reciprocalOffsetCache.get(key);if(cached)return [...cached];
  const abs=(x:Q)=>x.sign()<0?x.neg():x;
  const compare=(a:QVec,b:QVec)=>{
    const aa=a.map(abs),bb=b.map(abs);
    const size=aa.reduce((v,x)=>v.add(x),Z).cmp(bb.reduce((v,x)=>v.add(x),Z));if(size)return size;
    for(let j=0;j<3;j++){const c=aa[j].cmp(bb[j]);if(c)return c;}
    const outside=(v:QVec)=>v.filter(x=>x.sign()<0||x.cmp(ONE)>=0).length;
    const range=outside(a)-outside(b);if(range)return range;
    const signs=a.filter(x=>x.sign()<0).length-b.filter(x=>x.sign()<0).length;if(signs)return signs;
    for(let j=0;j<3;j++){const c=b[j].cmp(a[j]);if(c)return c;}return 0;
  };
  let candidates:QVec[]=[[]];
  for(let j=0;j<3;j++){
    const m=periods[j],fraction=residue[j].mod1();
    const choices=Array.from({length:toSafeInt(m,'reciprocal axis period')},(_,k)=>{
      const x=fraction.add(qi(k));return x.mulInt(2).cmp(m)>0?x.sub(m):x;
    });
    candidates=candidates.flatMap(row=>choices.map(x=>[...row,x]));
  }
  const allowed=candidates.filter(row=>{
    const delta=row.map((x,j)=>x.sub(residue[j]));
    return delta.every(x=>x.isInt())&&cents.every(c=>c.reduce((v,x,j)=>v.add(x.mul(delta[j])),Z).isInt());
  });
  allowed.sort(compare);const best=allowed[0];
  if(!best)throw new Error('no reciprocal offset representative');
  if(reciprocalOffsetCache.size>=4096)reciprocalOffsetCache.clear();reciprocalOffsetCache.set(key,best);
  return [...best];
}

/** One shared row-parameter set uses a,b,g. Distinct sets are numbered by
 * first appearance, so a star followed by an independent row uses sets 1,2.
 * A parameter shared across different row sets retains its first owner's name. */
function parameterSetNames(sets:number[][],labels:string[],owners:number[]):string[] {
  const groups=new Map<string,number>();
  const blockGroup=sets.map(set=>{const key=set.join(',');if(!groups.has(key))groups.set(key,groups.size+1);return groups.get(key)!;});
  return labels.map((label,k)=>label[0]+(groups.size===1?'':String(blockGroup[owners[k]])));
}

/** Rename independent component forms, preserving every exact shared relation.
 * An invertible parameter change is separate from the unimodular internal basis.
 * Opposite components share a name; dependent sums require another q basis. */
export function simpleParameterForm(q:Lin[][],names:string[]):string[]|null {
  const vectors:QVec[]=[], labels:string[]=[], reduced:QVec[]=[];
  for(const [i,row] of q.entries())for(const [j,x] of row.entries()) {
    if(x.coef.every(v=>v.isZero()))continue;
    if(vectors.some(v=>v.every((a,k)=>a.eq(x.coef[k]))||v.every((a,k)=>a.eq(x.coef[k].neg()))))continue;
    const v=[...x.coef];
    let r=[...v];
    for(const old of reduced){const at=old.findIndex(a=>!a.isZero());const f=r[at];r=r.map((a,k)=>a.sub(f.mul(old[k])));}
    const at=r.findIndex(a=>!a.isZero());if(at<0)return null;
    const pivot=r[at];reduced.push(r.map(a=>a.div(pivot)));
    vectors.push(v);labels.push(['a','b','g'][j]+String(i+1));
  }
  const sets=q.map(row=>vectors.flatMap((v,k)=>row.some(x=>v.every((a,h)=>a.eq(x.coef[h]))||v.every((a,h)=>a.eq(x.coef[h].neg())))?[k]:[]));
  const displayNames=parameterSetNames(sets,labels,vectors.map((_,k)=>sets.findIndex(set=>set.includes(k))));
  return q.map(row=>row.map(x=>{
    if(x.coef.every(v=>v.isZero()))return x.c.toString();
    const k=vectors.findIndex(v=>v.every((a,h)=>a.eq(x.coef[h]))||v.every((a,h)=>a.eq(x.coef[h].neg())));
    const sign=vectors[k].every((a,h)=>a.eq(x.coef[h]))?'':'-';
    return sign+displayNames[k]+(x.c.isZero()?'':(x.c.sign()>0?'+':'')+x.c.toString());
  }).join(','));
}

/** Permit dependent component forms only in later exact star images. The
 * coefficient comparison is symbolic; the offset lies in the basic reciprocal
 * lattice precisely when it is integral and annihilates every centring. */
export function starParameterForm(q:Lin[][],names:string[],rotations:number[][][],cents:QVec[]):string[]|null {
  const vectors:QVec[]=[],labels:string[]=[],owners:number[]=[];
  let coordinates=coordinateSolver(vectors);
  for(const [i,row] of q.entries())for(const [j,x] of row.entries()) {
    if(coordinates(x.coef)!==null)continue;
    vectors.push(x.coef);labels.push(['a','b','g'][j]+String(i+1));owners.push(i);
    coordinates=coordinateSolver(vectors);
  }
  const renamed=q.map(row=>row.map(x=>({...x,coef:coordinates(x.coef)!})));
  for(const [i,row] of renamed.entries()) {
    // New parameters must have a plain introduction in their first block.
    for(const [k,owner] of owners.entries())if(owner===i&&!row.some(x=>x.coef.every((v,h)=>v.eq(h===k?ONE:Z))))return null;
    if(row.every(x=>plainComponent(x)))continue;
    if(owners.some((owner,k)=>owner===i&&row.some(x=>!x.coef[k].isZero())))return null;
    const star=q.slice(0,i).some(seed=>rotations.some(R=>{
      const delta=row.map((_,j)=>q[i][j].c.sub(seed.reduce((v,x,h)=>v.add(x.c.mulInt(R[h][j])),Z)));
      return delta.every(x=>x.isInt())&&cents.every(c=>c.reduce((v,x,j)=>v.add(x.mul(delta[j])),Z).isInt())&&
        q[i].every((x,j)=>names.every((_,k)=>x.coef[k].eq(seed.reduce((v,y,h)=>v.add(y.coef[k].mulInt(R[h][j])),Z))));
    }));
    if(!star)return null;
  }
  const sets=renamed.map(row=>vectors.flatMap((_,k)=>row.some(x=>!x.coef[k].isZero())?[k]:[]));
  const displayNames=parameterSetNames(sets,labels,owners);
  return renamed.map(row=>formatQ(row,displayNames));
}
const plainComponent=(x:Lin)=>{
  const nz=x.coef.filter(v=>!v.isZero());return !nz.length||nz.length===1&&(nz[0].eq(ONE)||nz[0].eq(ONE.neg()));
};

/** q-block comparison ignores parameter spelling: free < nonzero rational < 0,
 * then larger rational first. Later exact symbol/map keys make the order total. */
export function compareQBlocks(a:string,b:string):number {
  const aa=a.split(','),bb=b.split(',');
  for(let j=0;j<3;j++){
    const free=(s:string)=>/[a-z]/i.test(s),af=free(aa[j]),bf=free(bb[j]);
    const av=af?Z:Q.parse(aa[j]),bv=bf?Z:Q.parse(bb[j]);
    const at=af?0:av.isZero()?2:1,bt=bf?0:bv.isZero()?2:1;
    if(at!==bt)return at-bt;
    if(at===1){const c=bv.cmp(av);if(c)return c;}
  }
  return 0;
}
const compareQForms=(a:string[],b:string[])=>{
  for(let i=0;i<a.length;i++){const c=compareQBlocks(a[i],b[i]);if(c)return c;}return 0;
};
/** Prefer nontrivial phases first among blocks with equal semantic q keys.
 * Count only out-of-order pairs: all-zero assignments keep zero cost too.
 * This orders blocks without introducing a preference for new phase content. */
export function phaseOrderCost(symbol:SsgSymbol):number {
  const form=(symbol.q.match(/\([^)]*\)/g)??[]).map(row=>row.slice(1,-1));
  let cost=0;
  for(let i=0;i<form.length;i++)for(let j=i+1;j<form.length;j++)
    if(compareQBlocks(form[i],form[j])===0&&!/[^0]/.test(symbol.shifts[i])&&/[^0]/.test(symbol.shifts[j]))cost++;
  return cost;
}
const settingCost=(form:string[],bIndex:number,lowerForD?:number)=> {
  const components=form.flatMap(row=>row.split(',')).filter(Boolean);
  const compound=components.filter(t=>/[a-z]/i.test(t)&&!/^[-]?[abg][0-9]*(?:[+-][0-9]+(?:\/[0-9]+)?)?$/.test(t)).length;
  const outside=components.filter(t=>!/[a-z]/i.test(t)&&Q.parse(t).cmp(ONE)>=0).length;
  const negative=components.filter(t=>/^-[0-9]/.test(t)).length;
  // Retain the fitted rational setting preference, minimised over block
  // permutations so it cannot choose their order. Branches use the smallest
  // possible positional weight as an admissible lower bound, never a prune
  // based on the current block's eventual position.
  const fixed=form.map(row=>row.split(',').reduce((v,t,j)=>v+(/^[-0-9]/.test(t)&&t!=='0'?(j+1):0),0));
  const constants=lowerForD===undefined?fixed.sort((a,b)=>a-b).reduce((v,x,i)=>v+x*10**(3-i),0):fixed.reduce((v,x)=>v+x,0)*10**(4-lowerForD);
  return compound*1e12+negative*1e9+outside*1e8+bIndex*1e6+constants;
};

/** Additive descriptor cost, also an admissible branch bound. Phase content
 * is known only after exact operation transformation and is ranked later. */
const nontrivialPhaseCount=(symbol:SsgSymbol)=>symbol.shifts.reduce((n,s)=>n+(s.match(/[stqh]/g)??[]).length,0);
const multiCost=(form:string[],_bIndex:number,_lowerForD?:number)=> {
  const components=form.flatMap(row=>row.split(','));
  const occurrences=components.reduce((n,x)=>n+(x.match(/[abg][0-9]*/g)??[]).length,0);
  const names=paramsOf(form),offsets=components.filter(x=>!parseLinear(x,names).c.isZero()).length;
  return occurrences*100+offsets;
};

/** All candidate BYU-axes forms of the group (exported for tests/debugging). */
type Prepared={basis:number[][];q0:Lin[][];form:string[];conventional:boolean;mp:QMat;cost:number};
type PreparedFamily={basis:number[][];q0:Lin[][];variants:{mp:QVec;text:string|undefined;preferred:boolean;cost:number}[][];remaining:number[];minimum:number};
const multiPrepared=new Map<string,Prepared[]>();
const multiBases=new Map<number,number[][][]>();
interface MultiCandidate { symbol: SsgSymbol; bIndex: number; P: QMat; p: QVec; E: number[][]; M: QMat; qActual: string[]; basis: number[][]; conventional: boolean; operators: string }
export function candidatesMulti(ops: SOp[], Qtext: string[]): { sgNumber: number; list: MultiCandidate[] } | Fail {
  if (Qtext.length < 2 || Qtext.length > 3 || ops.some((o) => o.W.length !== 3+Qtext.length)) return { ok: false, reason: "expected d = 2 or 3 and matching operations" };
  const names = paramsOf(Qtext).sort();
  const qs: Lin[][] = Qtext.map(row => row.split(",").map(s => parseLinear(s,names)));
  const q = qs[0];
  const paramCoefficients=names.map((_,k)=>qs.flatMap(row=>row.map(x=>x.coef[k])));
  const paramDen=paramCoefficients.flat().reduce((v,x)=>lcmB(v,x.d),1n);
  const parameterRank=hnf(paramCoefficients.map(row=>row.map(x=>x.n*(paramDen/x.d)))).H.length;
  if (q.length !== 3) return { ok: false, reason: `bad q "${Qtext[0]}"` };

  const ops3: Op3[] = ops.map((o) => ({ R: o.W.slice(0, 3).map((r) => r.slice(0, 3)), t: o.w.slice(0, 3) }));
  const id = identify(ops3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
  if (!id.ok) return { ok: false, reason: id.reason };
  const system = crystalSystemOf(id.number);
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
    const R = Rq.map((r) => r.map((x) => toSafeInt(x, "operator matrix entry")));
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
  // The d=1 quotient by the point group suppresses sign choices needed
  // for fixed components of a second q. Retain each exact point-group lift.
  const anchor=system==="monoclinic"?monoclinicAnchor(qs,names,L):undefined;
  const norm = normaliser(g, system === "triclinic", T,anchor).flatMap(({N,n}) => g.reps.map(({R,t:rt})=>({N:imul(R,N),n:vadd(mv(qmat(R),n),rt)})));

  const d=Qtext.length, list:MultiCandidate[]=[];
  const bases:number[][][]=multiBases.get(d)??[];
  const enumerate=(entries:number[])=>{
    if(entries.length===d*d){const E=Array.from({length:d},(_,i)=>entries.slice(i*d,(i+1)*d));const det=imatDet(E.map(r=>r.map(BigInt)));if(det===1n||det===-1n)bases.push(E);return;}
    for(const v of [0,1,-1])enumerate([...entries,v]);
  };
  if(!multiBases.has(d))enumerate([]);
  const basisCost=(E:number[][])=>E.reduce((s,r,i)=>s+r.reduce((v,x,j)=>v+(i!==j&&x!==0?100:0)+(i===j&&x<=0?10:0),0),0);
  bases.sort((A,B)=>basisCost(A)-basisCost(B)||(mkey(A)<mkey(B)?-1:mkey(A)>mkey(B)?1:0));
  multiBases.set(d,bases);

  const settings = system==='cubic' ? [ [[1,0,0],[0,1,0],[0,0,1]] ] : system==='monoclinic' ? [...B_SETTINGS.monoclinic, [[0,0,1],[1,0,0],[0,1,0]]] : B_SETTINGS[system];
  const basisCache=multiPrepared,symbols=new Map<string,SsgSymbol|null>();
  const familyCache=new Map<string,PreparedFamily[]>();
  const transformedCache=new Map<string,Map<string,SOp[]>>();
  const operatorKeys=new WeakMap<SOp[],string>();
  let winning:number[]|undefined;
  const independentZeroBases=new Map<string,number>();
  // Plain components precede exact star compounds. Conventional descriptors
  // compete at the same sparsity level; they cannot hide a reduced basis.
  for(const simpleOnly of [true,false]) {
  for(const preferredOnly of [false]) {
  // Seed the fallback with the exact lowest eligible descriptor cost of each
  // affine lift until one is valid. Its cost is then a proven upper bound.
  // The full traversal restarts below with that bound and an empty list, so
  // even equal rank keys retain the original affine-lift/tie order.
  for(const minimumOnly of preferredOnly?[false]:[true,false]) {
  if(!minimumOnly)list.length=0;
  for(const [bIndex,Bn] of settings.entries()) for(const {N,n} of norm) {
    // Every independent parameter needs an occurrence. At that absolute lower
    // bound with no offsets/phases, a worse basic setting cannot improve.
    if(winning&&winning[0]===parameterRank*100&&winning[1]===0&&winning[2]===0&&winning[3]<bIndex*1e6)continue;

    const B=qmat(Bn),BN=matMul(B,qmat(N)),P=matMul(BN,T),p=vadd(mv(BN,t),mv(B,n)),Pi=matInv(P),pk=mkey(P);
    // Cache by the exact HNF q lattice, not the input-relative affine map.
    const coeff=qs.map(row=>names.flatMap((_,k)=>[0,1,2].map(j=>row.reduce((v,x,h)=>v.add(x.coef[k].mul(Pi[h][j])),Z))));
    const D=coeff.flat().reduce((v,x)=>lcmB(v,x.d),1n),H=hnf(coeff.map(row=>row.map(x=>x.n*(D/x.d))));
    if(H.H.length!==d)continue;
    // GL row changes preserve the rank of each physical-component matrix.
    // Rank r needs at least r nonzero components, hence at least r parameter
    // occurrences. This exact bound also rejects sheared dense families before
    // preparing thousands of descriptors that cannot beat a sparse winner.
    const occurrenceLower=[0,1,2].reduce((v,j)=>v+hnf(H.H.map(row=>names.map((_,k)=>row[k*3+j]))).H.length,0);
    if(winning&&occurrenceLower*100>winning[0])continue;
    const U=H.U.map(row=>row.map(x=>toSafeInt(x,'internal basis')));
    const rotations=ops.map(o=>matMul(matMul(P,qmat(o.W.slice(0,3).map(r=>r.slice(0,3)))),Pi).map(r=>r.map(x=>toSafeInt(x,"rotation"))));
    const canonical=H.H.map(row=>row.map(x=>Q.of(x,D)));
    const constants=U.map(row=>[0,1,2].map(j=>qs.reduce((v,r,l)=>v.add(r.reduce((w,x,h)=>w.add(x.c.mul(Pi[h][j]).mulInt(row[l])),Z)),Z)));
    const centred=g.cents.map(c=>mv(B,c));
    const active=[0,1,2].filter(j=>canonical.some(row=>names.some((_,k)=>!row[k*3+j].isZero())));
    // Full component rank means every GL basis has the same independent
    // displayed components. With zero offsets its formatted q is fixed too.
    const independent=parameterRank===d*active.length&&constants.flat().every(x=>x.isZero());
    const zeroKeyForSetting=pk+'|'+vkey(p),zeroLimit=independentZeroBases.get(zeroKeyForSetting);
    const impossibleConventional=active.length===1&&active[0]!==2||active.length===2&&active.join(',')!=='0,1';
    if(independent&&impossibleConventional&&winning&&winning[0]===parameterRank*100&&winning[1]===0&&winning[2]===0)continue;
    const costLimit=winning?.[0];
    const familyKey=[String(costLimit),String(minimumOnly),String(zeroLimit),String(simpleOnly),String(preferredOnly),names.join(','),system,mkey(canonical),mkey(constants),[...new Set(rotations.map(mkey))].sort().join(' '),centred.map(vkey).sort().join(' ')].join('|');
    const cacheKey=[String(minimumOnly),String(costLimit),familyKey].join('|');
    let prepared=basisCache.get(cacheKey);
    if(!prepared){
      if([...basisCache.values()].reduce((n,a)=>n+a.length,0)>1_000_000)basisCache.clear();
      prepared=[];
      let families=familyCache.get(familyKey);
      if(!families){
      families=[];
      // A small seed pass establishes valid upper bounds in every basic
      // setting; the full pass still considers every potentially winning basis.
      const seeds=minimumOnly?bases.filter(A=>A.every(r=>r.filter(x=>x!==0).length===1)):bases;
      const allBases=new Map(seeds.map(A=>[mkey(A),A]));
      if(!minimumOnly)for(const A of orbitBases(canonical,names,rotations))allBases.set(mkey(A),A);
      // All coefficient entries share D. Accumulate integral numerators once
      // instead of normalising every multiply/add in each trial basis.
      const constantDen=constants.flat().reduce((v,x)=>lcmB(v,x.d),1n);
      const constantNumerators=constants.map(row=>row.map(x=>x.n*(constantDen/x.d)));
      for(const basis of allBases.values()){
        if(independent&&zeroLimit!==undefined&&basisCost(basis)>zeroLimit)continue;
        const q0=basis.map(row=>{
          const integers=row.map(BigInt);
          return [0,1,2].map(j=>({
            coef:names.map((_,k)=>Q.of(integers.reduce((v,x,l)=>v+x*H.H[l][k*3+j],0n),D)),
            c:Q.of(integers.reduce((v,x,l)=>v+x*constantNumerators[l][j],0n),constantDen)
          }));
        });
        // A nonzero component needs at least one displayed occurrence under
        // any invertible parameter renaming. Keep equality for every later tie.
        if(costLimit!==undefined&&q0.reduce((n,row)=>n+row.filter(x=>x.coef.some(v=>!v.isZero())).length,0)*100>costLimit)continue;
        // Reject ineligible coefficient patterns before preparing offset and
        // conventional-setting variants; these predicates are independent.
        const simple=simpleParameterForm(q0,names);
        if(simpleOnly!==(simple!==null))continue;
        const parameterForm=simple??starParameterForm(q0,names,rotations,centred);
        if(!parameterForm)continue;
        // Constants cannot repair a failed coefficient pattern. This is only
        // pruning within the preferred score level; the fallback pass includes it.
        const coefficientQ=q0.map(row=>row.map(x=>({...x,c:Z}))),rotationOps=rotations.map(W=>({W,w:[]}));
        const coefficientForm=conventional(coefficientQ,names,system,rotationOps);
        if(preferredOnly&&!coefficientForm)continue;
        // Reduce the whole offset vector in the centred reciprocal lattice.
        // The resulting shift is applied to M and the operations as well as
        // the displayed q: phase letters must be recomputed in that gauge.
        // Every formerly enumerated shift differs by a reciprocal vector and
        // therefore has this same representative; none needs a separate branch.
        const choices=q0.map(row=>{
          const offsets=row.map(x=>x.c),canonical=reciprocalOffsetRepresentative(offsets,centred);
          return [canonical.map((x,j)=>x.sub(offsets[j]))];
        });
        // Coefficients are fixed here. All remaining conventional() tests on
        // constants are row-local (including the shared seed's zero constants).
        // Test each row against zero constants elsewhere, preserving that exact
        // predicate and caching its formatting instead of repeating the product.
        const simpleNames=paramsOf(parameterForm);
        const simpleQ=parameterForm.map(row=>row.split(',').map(x=>parseLinear(x,simpleNames)));
        const variants=choices.map((row,i)=>row.map(mp=>{
          const shifted=q0[i].map((x,j)=>({...x,c:mp[j].isZero()?x.c:x.c.add(mp[j])}));
          const form=coefficientForm?conventional(coefficientQ.map((r,k)=>k===i?shifted:r),names,system,rotationOps):null;
          const display=simpleQ[i].map((x,j)=>({...x,c:shifted[j].c}));
          const text=formatQ(display,simpleNames);
          return {mp,text,preferred:form!==null,cost:text===undefined?0:multiCost([text],0,d)};
        }).filter(x=>x.text!==undefined));
        if(variants.some(row=>!row.length))continue;
        const rowMin=variants.map(row=>Math.min(...row.map(x=>x.cost))),remaining=Array(d+1).fill(0);
        for(let i=d-1;i>=0;i--)remaining[i]=remaining[i+1]+rowMin[i];
        // A fallback needs at least one row that fails conventional(). Force
        // each possible failing row in turn; other rows attain independent minima.
        // This is its exact eligible minimum, without a Cartesian enumeration.
        let minimum:number|undefined=remaining[0];
        if(!preferredOnly)for(const [i,row] of variants.entries())for(const v of row)if(!v.preferred){
          const cost=remaining[0]-rowMin[i]+v.cost;
          if(minimum===undefined||cost<minimum)minimum=cost;
        }
        if(minimum===undefined)continue;
        families.push({basis,q0,variants,remaining,minimum});
      }
      // The family key contains every exact input used to build these rows.
      // Retain all eligible families in their original order, independently of
      // the current bound; applying a later bound reuses work without changing
      // the candidate minimum or stable traversal order. Keep memory local.
      if(familyCache.size>=4)familyCache.delete(familyCache.keys().next().value!);
      familyCache.set(familyKey,families);
      }
      const limit=minimumOnly?undefined:costLimit;
      for(const {basis,q0,variants,remaining,minimum} of families){
        if(limit!==undefined&&minimum>limit)continue;
        // Row costs are nonnegative admissible lower bounds on multiCost.
        // An established winner can only improve: a branch above that bound
        // cannot reach the minimum. Keep equality and the traversal/tie order;
        // include the bound and search mode in the cache key for later searches.
        const visit=(i:number,mp:QMat,text:string[],cost:number,preferred:boolean)=>{
          if(limit!==undefined&&cost+remaining[i]>limit)return;
          if(i===d){
            const exactCost=multiCost(text,0);
            if(limit===undefined||exactCost<=limit)prepared!.push({basis,q0,form:text,conventional:preferred,mp,cost:exactCost});
            return;
          }
          for(const row of variants[i])visit(i+1,[...mp,row.mp],[...text,row.text!],cost+row.cost,preferred&&row.preferred);
        };
        visit(0,[],[],0,true);
      }
      // Descriptor sparsity is independent of phases and physical setting.
      prepared.sort((a,b)=>a.cost-b.cost||compareQForms(a.form,b.form));
      if(minimumOnly&&prepared.length)prepared=prepared.filter(x=>x.cost===prepared![0].cost);
      if(basisCache.size>=64)basisCache.delete(basisCache.keys().next().value!);basisCache.set(cacheKey,prepared);
    }
    // Rank descriptors before constructing symbols. A worse score cannot win
    // after an exact valid candidate has established the current best level.
    const score=(x:Prepared,phases=0):number[]=>[x.cost,phases,Number(!x.conventional),settingCost(x.form,bIndex)];
    const compare=(a:number[],b:number[])=>{for(let i=0;i<a.length;i++)if(a[i]!==b[i])return a[i]-b[i];return 0;};
    let applyInternal:ReturnType<typeof internalTransformer>|null|undefined;
    const zeroPhaseBases=new Map<string,number>();
    for(const item of prepared) {
      const zeroKey=Number(item.conventional)+'|'+item.form.join(';');
      const zeroBasis=zeroPhaseBases.get(zeroKey);
      // Within this affine basic setting and identical displayed form, a valid
      // zero-phase symbol already has the minimum phase content and full text.
      // A strictly worse HNF-relative basis cost cannot win a later tie.
      if(zeroBasis!==undefined&&basisCost(item.basis)>zeroBasis)continue;
      // Zero nontrivial phases is an admissible lower bound. Once achieved,
      // inferior settings can be skipped without hiding any phase improvement.
      if(winning&&item.cost>winning[0])break;
      if(winning&&compare(score(item),winning)>0)continue;
      const {basis,q0,form,mp,conventional:preferred}=item;
      // No raw-text prefix pruning: the semantic block order and complexity
      // rank must see all equal-score candidates, with their phases attached.
      const E=imul(basis,U),M=mp.map(row=>toM(row,P));
      if(applyInternal===undefined){const base=transformOps(ops,P,p,U,Array.from({length:d},()=>[Z,Z,Z]),Array(d).fill(Z));applyInternal=base?internalTransformer(base,p,transformedCache):null;}
      const newOps=applyInternal?applyInternal(basis,mp):transformOps(ops,P,p,E,M,Array(d).fill(Z));
      if(!newOps||newOps.some(o=>isIdW(o.W)&&o.w.slice(3).some(x=>!x.isZero())))continue;
      let operators=operatorKeys.get(newOps);
      if(operators===undefined){operators=newOps.map(opKey).sort().join(' ');operatorKeys.set(newOps,operators);}
      const sk=form.join(';')+'|'+operators;
      if(!symbols.has(sk))symbols.set(sk,ssgSymbolInAxes(newOps,form,system));
      const symbol=symbols.get(sk);
      if(symbol){
        const level=score(item,nontrivialPhaseCount(symbol));
        if(winning&&compare(level,winning)>0)continue;
        if(level[1]===0){
          zeroPhaseBases.set(zeroKey,Math.min(zeroBasis??Infinity,basisCost(item.basis)));
          if(independent&&preferred)independentZeroBases.set(zeroKeyForSetting,Math.min(independentZeroBases.get(zeroKeyForSetting)??Infinity,basisCost(item.basis)));
        }
        if(winning&&compare(level,winning)<0)list.length=0;
        winning=level;list.push({symbol,bIndex,P,p,E,M,basis,conventional:preferred,operators,qActual:q0.map((row,i)=>formatQ(row.map((x,j)=>({...x,c:x.c.add(mp[i][j])})),names))});
      }
    }

  }
  }
  if(list.length)break;
  }
  if(list.length)break;
  }
  return {sgNumber:id.number,list};
}

// Fitted algebraic rule: independent coordinate parameters, or an orbit of one
// such seed under the basic point group. Validate coefficient relations exactly.
// Parameter names follow seed axes, indexed only when independent q blocks coexist.
function conventional(q:Lin[][], names:string[], system:CrystalSystem, ops:SOp[]):string[]|null {
  const d=q.length, renamed=names.map(()=>''), owners=names.map(()=>new Set<number>());
  for(const [i,row] of q.entries()) for(const x of row) x.coef.forEach((v,k)=>{if(!v.isZero())owners[k].add(i)});
  const shared=owners.some(o=>o.size>1);
  if(shared) {
    const seed=q[0], introduced=new Set<number>();
    for(const [j,x] of seed.entries()) {
      const nz=x.coef.map((v,k)=>!v.isZero()?k:-1).filter(k=>k>=0);
      if(nz.length>1 || nz.some(k=>!x.coef[k].eq(ONE)) || (nz.length&&!x.c.isZero()))return null;
      for(const k of nz) {if(!introduced.has(k))renamed[k]=['a','b','g'][j];introduced.add(k);}
    }
    if(introduced.size!==names.length)return null;
    const support=seed.filter(x=>x.coef.some(v=>!v.isZero())).length;
    for(const row of q.slice(1)) {
      const orbit=ops.some(o=>row.every((x,j)=>names.every((_,k)=>x.coef[k].eq(seed.reduce((v,y,h)=>v.add(y.coef[k].mulInt(o.W[h][j])),Z)))));
      if(!orbit)return null;
      const first=row.flatMap(x=>x.coef).find(v=>!v.isZero());
      if(!first || (support>1 && first.sign()>0))return null;
    }
  } else {
    const supports:number[][]=[];
    for(const [i,row] of q.entries()) {
      const support:number[]=[];
      const coeff:QMat=[];
      for(const [j,x] of row.entries())if(x.coef.some(v=>!v.isZero())){support.push(j);coeff.push(x.coef);}
      if(!support.length)return null;
      // Independent row parameters may be any nonsingular linear mixture after
      // an allowed shear of physical axes. Compare their span, not their spelling.
      let rank=0;const A=coeff.map(r=>[...r]);
      for(let k=0;k<names.length&&rank<A.length;k++){
        const at=A.findIndex((r,j)=>j>=rank&&!r[k].isZero());if(at<0)continue;
        [A[rank],A[at]]=[A[at],A[rank]];const pivot=A[rank][k];A[rank]=A[rank].map(x=>x.div(pivot));
        for(let j=rank+1;j<A.length;j++){const f=A[j][k];A[j]=A[j].map((x,h)=>x.sub(f.mul(A[rank][h])));}rank++;
      }
      if(rank!==support.length)return null;supports.push(support);
    }
    for(let i=1;i<d;i++)if(supports[i-1].join(',')>supports[i].join(','))return null;
    // Conventional independent axial modulations occupy c, or b then c when
    // their external directions differ. Independent planar rows occupy ab.
    const axial=supports.filter(s=>s.length===1).map(s=>s[0]);
    const directions=[...new Set(axial)].sort();
    if(directions.length===1&&directions[0]!==2)return null;
    if(directions.length===2&&directions.join(',')!=='1,2')return null;
    if(supports.some(s=>s.length===2&&s.join(',')!=='0,1'))return null;
  }
  for(const row of q) for(const x of row) if(x.c.cmp(Q.of(-1,2))<0||x.c.cmp(ONE)>0||![1n,2n,3n].includes(x.c.d))return null;
  if(!shared)return q.map((row,i)=>row.map((x,j)=>x.coef.some(v=>!v.isZero())?['a','b','g'][j]+String(i+1):x.c.toString()).join(','));
  return q.map(row=>formatQ(row,renamed));
}

function standardiseMulti(ops:SOp[],q:string[]):Standardised | Fail {
  const c=candidatesMulti(ops,q); if('ok' in c)return c;
  if(!c.list.length)return {ok:false,reason:'no valid transformed superspace lattice'};
  const form=(x:MultiCandidate)=>(x.symbol.q.match(/\([^)]*\)/g)??[]).map(row=>row.slice(1,-1));
  const cost=(x:MultiCandidate)=>multiCost(form(x),0);
  const phases=(x:MultiCandidate)=>nontrivialPhaseCount(x.symbol);
  // Total order: plain/star eligibility, occurrences, offsets, phase letters,
  // conventional preference, sign/normalisation and basic-axis setting cost,
  // semantic q-block order, tied-block phase order, full symbol in Unicode codepoint order; HNF-relative
  // basis complexity and entries; actual q; then exact serialized affine map.
  // Never rank by the input-relative E: identity depends on the supplied setting.
  const cmp=(a:string,b:string)=>a<b?-1:a>b?1:0;
  const basisCost=(E:number[][])=>E.reduce((s,r,i)=>s+r.reduce((v,x,j)=>v+(i!==j&&x!==0?100:0)+(i===j&&x<=0?10:0),0),0);
  const key=(x:MultiCandidate)=>[mkey(x.basis),x.operators,x.qActual.join(';')].join('|');
  const affineKey=(x:MultiCandidate)=>[mkey(x.P),vkey(x.p),mkey(x.E),mkey(x.M)].join('|');
  c.list.sort((a,b)=>cost(a)-cost(b)||phases(a)-phases(b)||Number(b.conventional)-Number(a.conventional)||settingCost(form(a),a.bIndex)-settingCost(form(b),b.bIndex)||compareQForms(form(a),form(b))||phaseOrderCost(a.symbol)-phaseOrderCost(b.symbol)||cmp(a.symbol.full,b.symbol.full)||basisCost(a.basis)-basisCost(b.basis)||cmp(key(a),key(b))||cmp(affineKey(a),affineKey(b)));
  const best=c.list[0];
  return {ok:true,sgNumber:c.sgNumber,symbol:best.symbol,qActual:best.qActual,
    P:best.P.map(r=>r.map(String)),p:best.p.map(String),E:best.E,M:best.M.map(r=>r.map(String)),p4:q.map(()=>"0"),alternates:[]};
}


const typeFamilyCache=new Map<string,string[]>();
const typeNameCache=new Map<string,Standardised|Fail>();
export const _namingStats={fallbacks:0};

/** Full affine solution of Q R - E Q = M. Translation phases do not restrict
 * Q: the actual modulation is a specialisation of this operation-defined type.
 * Free parameters are independent, even when the supplied modulation is not. */
export function typeQFamily(ops:SOp[],d:number):string[] {
  if(d<2||d>3||!ops.length||ops.some(o=>o.W.length!==3+d))throw new Error('expected d=2 or 3 operations');
  const matrices=[...new Map(ops.map(o=>[mkey(o.W),o.W])).entries()].sort(([a],[b])=>a<b?-1:a>b?1:0);
  const key=matrices.map(([k])=>k).join('|'),cached=typeFamilyCache.get(key);if(cached)return [...cached];
  const size=3*d,rows:QMat=[];
  for(const [,W] of matrices)for(let i=0;i<d;i++)for(let j=0;j<3;j++){
    const row=Array.from({length:size},(_,k)=>{
      const l=Math.floor(k/3),h=k%3;return qi((i===l?W[h][j]:0)-(h===j?W[3+i][3+l]:0));
    });
    const rhs=qi(W[3+i][j]);if(row.some(x=>!x.isZero())||!rhs.isZero())rows.push([...row,rhs]);
  }
  const pivots:number[]=[];let rank=0;
  for(let j=0;j<size&&rank<rows.length;j++){
    const at=rows.findIndex((r,i)=>i>=rank&&!r[j].isZero());if(at<0)continue;
    [rows[rank],rows[at]]=[rows[at],rows[rank]];const pivot=rows[rank][j];rows[rank]=rows[rank].map(x=>x.div(pivot));
    for(let i=0;i<rows.length;i++)if(i!==rank&&!rows[i][j].isZero()){
      const f=rows[i][j];rows[i]=rows[i].map((x,k)=>x.sub(f.mul(rows[rank][k])));
    }
    pivots.push(j);rank++;
  }
  if(rows.slice(rank).some(r=>!r[size].isZero()))throw new Error('inconsistent operation q-family');
  const free=Array.from({length:size},(_,i)=>i).filter(i=>!pivots.includes(i)),names=free.map((_,i)=>'p'+(i+1));
  const forms=Array.from({length:size},(_,i):Lin=>{
    const pivot=pivots.indexOf(i);return pivot<0?{c:Z,coef:free.map(j=>i===j?ONE:Z)}:
      {c:rows[pivot][size],coef:free.map(j=>rows[pivot][j].neg())};
  });
  const family=Array.from({length:d},(_,i)=>formatQ(forms.slice(i*3,i*3+3),names));
  if(typeFamilyCache.size>=1024)typeFamilyCache.delete(typeFamilyCache.keys().next().value!);
  typeFamilyCache.set(key,family);return [...family];
}

/** Exact affine inclusion, allowing a proper specialisation of the type. */
export function isQSpecialisation(actual:string[],family:string[]):boolean {
  if(actual.length!==family.length)return false;
  const an=paramsOf(actual),fn=paramsOf(family),parse=(q:string[],names:string[])=>q.flatMap(r=>r.split(',').map(x=>parseLinear(x,names)));
  const a=parse(actual,an),f=parse(family,fn);if(a.length!==f.length)return false;
  const coordinates=coordinateSolver(fn.map((_,k)=>f.map(x=>x.coef[k])));
  return an.every((_,k)=>coordinates(a.map(x=>x.coef[k]))!==null)&&coordinates(a.map((x,i)=>x.c.sub(f[i].c)))!==null;
}

/** Apply precisely the emitted coordinate map, without renaming or enlarging
 * the subgroup's actual modulation family. */
export function mapQFamily(q:string[],P:QMat,E:number[][],M:QMat):string[] {
  const names=paramsOf(q),a=q.map(r=>r.split(',').map(x=>parseLinear(x,names))),Pi=matInv(P);
  return E.map((row,i)=>formatQ([0,1,2].map(j=>({
    coef:names.map((_,k)=>row.reduce((s,e,l)=>s.add(a[l].reduce((t,x,h)=>t.add(x.coef[k].mul(Pi[h][j]).mulInt(e)),Z)),Z)),
    c:row.reduce((s,e,l)=>s.add(a[l].reduce((t,x,h)=>t.add(x.c.mul(Pi[h][j]).mulInt(e)),Z)),Z).add(M[i].reduce((s,x,h)=>s.add(x.mul(Pi[h][j])),Z))
  })),names));
}

function standardiseTypeFamily(ops:SOp[],actual:string[]):Standardised|Fail {
  let family:string[];
  try{family=typeQFamily(ops,actual.length);}catch(e){return {ok:false,reason:String(e)};}
  if(!isQSpecialisation(actual,family))return {ok:false,reason:'actual q is not a specialisation of the operation-defined type family'};
  const finish=(st:Standardised):Standardised|Fail=>{
    const qActual=mapQFamily(actual,st.P.map(r=>r.map(Q.parse)),st.E,st.M.map(r=>r.map(Q.parse)));
    const form=[...st.symbol.q.matchAll(/\(([^()]*)\)/g)].map(m=>m[1]);
    return isQSpecialisation(qActual,form)?{...st,qActual}:{ok:false,reason:'actual q is outside the named type family'};
  };
  // The type depends on the operations, not on which exact member was pasted.
  // Keep operation order in the key and copy all mutable result fields.
  const nameKey=ops.map(opKey).join(' ');
  const cached=typeNameCache.get(nameKey);
  const clone=(st:Standardised):Standardised=>({...st,symbol:{...st.symbol,shifts:[...st.symbol.shifts]},P:st.P.map(r=>[...r]),p:[...st.p],E:st.E.map(r=>[...r]),M:st.M.map(r=>[...r]),p4:[...st.p4],qActual:[...st.qActual],alternates:[...st.alternates]});
  const attempt=():Standardised|Fail=>{try{return standardiseMulti(ops,family);}catch(e){return {ok:false,reason:String(e)};}};
  const named=cached?(cached.ok?clone(cached):cached):attempt();
  if(!cached){if(typeNameCache.size>=512)typeNameCache.delete(typeNameCache.keys().next().value!);typeNameCache.set(nameKey,named.ok?clone(named):named);}
  if(named.ok){const result=finish(named);if(result.ok)return result;}
  // Preserve the exact historical class name whenever the preferred writer has
  // no admissible form. This code snapshot contains rules, never catalogue data.
  for(const q of [actual,family]){
    const legacy=legacyStandardise(ops,q);if(!legacy.ok)continue;
    const result=finish(legacy);if(result.ok){_namingStats.fallbacks++;return {...result,namingFallback:'d244096'};}
  }
  return named.ok?{ok:false,reason:'no exact member of an admissible named family'}:named;
}
