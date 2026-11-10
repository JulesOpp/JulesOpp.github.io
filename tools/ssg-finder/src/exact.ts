// Exact t0 classes of commensurate sections of a (3+d)D superspace group.
// TypeScript port of reference/exact.py (PLAN.md §4, CONTRACTS.md §3.1).
//
//   op g = (W, w), W = [[R, 0], [M, E]], w = (tE, tI);  M + E Q = Q R.
//   A = E − 1, c = tI − Q tE. g survives at t iff A t + c ∈ Λ = Z^d + Q Z^3,
//   and its matched n is the coset with Q n ≡ A t + c (mod Z^d).
//   t = B s with B a basis of Λ (columns), s ∈ R^d/Z^d.  Ã = B⁻¹AB (integer), c̃ = B⁻¹c.

import { Q, QMat, QVec, denLcm, matInv, matMul, matVec, transpose, gcd } from "./rational";
import { IMat, hnf, hnfLower, snf, imatRank, imatDet } from "./intmat";
import { GroupRecord, SOp, buildOps, evalQ, splitOp, formatOp } from "./ssg";
import type { Op3 } from "./identify";

export const P_GENERIC = 10007n;
export const GENERIC_SEEDS = [3001n, 4723n, 7177n];
export const GENERIC_STEP = 1009n;
export const EXAMPLE_MAX_DEN = 24;

// ------------------------------------------------------------------ helpers

const Z = Q.ZERO;
const fq = (n: bigint | number, d: bigint | number = 1n) => Q.of(BigInt(n), BigInt(d));
const qcmpTuple = (a: Q[], b: Q[]) => {
  for (let i = 0; i < Math.min(a.length, b.length); i++) {
    const c = a[i].cmp(b[i]);
    if (c) return c;
  }
  return a.length - b.length;
};
const tupleKey = (a: Q[]) => a.map(String).join(",");
const toBig = (A: number[][]): IMat => A.map((r) => r.map(BigInt));
const qOfBigMat = (A: IMat): QMat => A.map((r) => r.map((x) => fq(x)));
function intOfQMat(A: QMat, what: string): IMat {
  return A.map((r) => r.map((x) => {
    if (!x.isInt()) throw new Error(`${what}: non-integer entry ${x}`);
    return x.n;
  }));
}
const absB = (x: bigint) => (x < 0n ? -x : x);

// ------------------------------------------------------------------ flats

/** Connected affine subtorus {s ∈ R^d/Z^d : C s ≡ f (mod 1)}, C saturated row HNF. */
export class Flat {
  readonly C: IMat;
  readonly f: Q[];
  readonly d: number;
  readonly key: string;
  readonly pivots: number[];

  constructor(C: IMat, f: Q[], d: number) {
    this.C = C;
    this.f = f;
    this.d = d;
    this.key = C.map((r) => r.join(",")).join(";") + "|" + f.map(String).join(",");
    this.pivots = C.map((r) => r.findIndex((x) => x !== 0n));
  }
  get r() { return this.C.length; }
  get dim() { return this.d - this.C.length; }

  contains(s: Q[]): boolean {
    return this.C.every((row, i) => row.reduce((acc, c, j) => acc.add(s[j].mulInt(c)), Z).sub(this.f[i]).mod1().isZero());
  }

  /** Point with the non-pivot coordinates = free (in order), pivots solved bottom-up. */
  point(free: Q[]): Q[] {
    const s: (Q | null)[] = Array(this.d).fill(null);
    let k = 0;
    for (let j = 0; j < this.d; j++) if (!this.pivots.includes(j)) s[j] = free[k++];
    for (let i = this.r - 1; i >= 0; i--) {
      const p = this.pivots[i];
      const row = this.C[i];
      let acc = this.f[i];
      for (let j = p + 1; j < this.d; j++) acc = acc.sub((s[j] as Q).mulInt(row[j]));
      s[p] = acc.div(fq(row[p]));
    }
    return s as Q[];
  }

  toJSON() {
    return { C: this.C.map((r) => r.map(Number)), f: this.f.map(String) };
  }
}

export function canonicalFlat(C: IMat, f: Q[], d: number): Flat {
  if (!C.length) return new Flat([], [], d);
  const { H, U } = hnf(C);
  if (H.length !== C.length) throw new Error("C not of full row rank");
  const f2 = U.slice(0, H.length).map((Ur) => Ur.reduce((acc, u, j) => acc.add(f[j].mulInt(u)), Z).mod1());
  return new Flat(H, f2, d);
}

/** Connected components of {s : G s ≡ h (mod 1)}. */
export function solve(G: IMat, h: Q[], d: number): Flat[] {
  if (!G.length) return [new Flat([], [], d)];
  const { U, V, diag } = snf(G);
  const k = G.length;
  const Uh = U.map((Ur) => Ur.reduce((acc, u, j) => acc.add(h[j].mulInt(u)), Z));
  const rho = diag.length;
  for (let i = rho; i < k; i++) if (!Uh[i].isInt()) return [];
  // V^{-1} rows: invert the unimodular V exactly
  const Vi = intOfQMat(matInv(qOfBigMat(V)), "V^-1");
  const C = Vi.slice(0, rho);
  const out: Flat[] = [];
  const js = Array(rho).fill(0n);
  for (;;) {
    out.push(canonicalFlat(C, js.map((jv, i) => Uh[i].add(fq(jv)).div(fq(diag[i]))), d));
    let i = 0;
    while (i < rho) {
      js[i]++;
      if (js[i] < diag[i]) break;
      js[i] = 0n;
      i++;
    }
    if (i === rho) break;
  }
  return out;
}

export const intersect = (F1: Flat, F2: Flat) => solve([...F1.C, ...F2.C], [...F1.f, ...F2.f], F1.d);

/** G ⊆ F ? */
export function subset(G: Flat, F: Flat): boolean {
  if (G.dim > F.dim) return false;
  if (F.r && imatRank([...G.C, ...F.C]) !== G.r) return false;
  return F.contains(G.point(Array(G.dim).fill(Z)));
}

// ------------------------------------------------------------------ result types

export interface OpInfo {
  R: number[][]; M: number[][]; E: number[][];
  tE: Q[]; tI: Q[]; A: number[][]; c: Q[];
  At: IMat; Et: IMat; EtInv: IMat; ct: Q[];
  flats: Flat[] | null; always: boolean;
}

export interface TClass {
  flat: Flat;
  dim: number;
  genericS: Q[];
  genericT: Q[];
  exampleT: Q[];
  surviving: [number, number[]][];
  ops3: Op3[];
  orbit: number;
  locus: string;
  itLabel?: string;
}

// ------------------------------------------------------------------ main analysis

export class Exact {
  readonly rec: GroupRecord;
  readonly d: number;
  readonly params: Record<string, Q>;
  readonly Q: QMat;
  readonly ops: SOp[];
  N!: number[];
  B!: QMat;
  Binv!: QMat;
  cells!: number[][];
  hnfSuper!: IMat;
  index!: number;
  nlookup!: Map<string, number[]>;
  lamReps!: Q[][];
  opinfo!: OpInfo[];
  flats!: Map<string, Flat>;
  sub!: Map<string, Flat[]>;
  classes!: TClass[];
  itN?: number;

  constructor(rec: GroupRecord, params: Record<string, Q | string>) {
    this.rec = rec;
    this.d = rec.d;
    this.params = Object.fromEntries(Object.entries(params).map(([k, v]) => [k, typeof v === "string" ? Q.parse(v) : v]));
    this.Q = evalQ(rec, this.params);
    this.ops = buildOps(rec);
    this.setupLattices();
    this.setupOps();
    this.closure();
    this.buildClasses();
  }

  // ---------------------------------------------------------------- lattices
  private setupLattices() {
    const { d, Q: Qm } = this;
    this.N = [0, 1, 2].map((j) => Number(denLcm(Qm.map((r) => r[j]))));
    const Dd = denLcm(Qm.flat());
    const gens: IMat = [
      ...Array.from({ length: d }, (_, i) => Array.from({ length: d }, (_, j) => (i === j ? Dd : 0n))),
      ...[0, 1, 2].map((j) => Array.from({ length: d }, (_, i) => Qm[i][j].mul(fq(Dd)).n)),
    ];
    const H = hnfLower(gens);
    const Lrows = H.map((r) => r.map((x) => fq(x, Dd)));
    this.B = transpose(Lrows);
    this.Binv = matInv(this.B);
    // supercell L_super = {n : Q n ∈ Z^d}, box points in lexicographic order
    const box: number[][] = [];
    for (let a = 0; a < this.N[0]; a++)
      for (let b = 0; b < this.N[1]; b++)
        for (let c = 0; c < this.N[2]; c++) box.push([a, b, c]);
    const Qn = (n: number[]) => Qm.map((r) => r.reduce((acc, x, j) => acc.add(x.mulInt(n[j])), Z));
    this.cells = box.filter((n) => Qn(n).every((x) => x.isInt()));
    const sgens: IMat = [
      ...this.cells.filter((n) => n.some((x) => x)).map((n) => n.map(BigInt)),
      ...[0, 1, 2].map((j) => [0, 1, 2].map((i) => (i === j ? BigInt(this.N[j]) : 0n))),
    ];
    this.hnfSuper = hnfLower(sgens);
    this.index = Number(absB(imatDet(this.hnfSuper)));
    this.nlookup = new Map();
    for (const n of box) {
      const k = tupleKey(Qn(n).map((x) => x.mod1()));
      if (!this.nlookup.has(k)) this.nlookup.set(k, n);
    }
    if (this.nlookup.size !== this.index) throw new Error("internal: |Λ/Z^d| ≠ supercell index");
    this.lamReps = [...this.nlookup.values()].map((n) => Qn(n).map((x) => x.mod1()));
  }

  // ---------------------------------------------------------------- operators
  private setupOps() {
    const { d, Q: Qm } = this;
    this.opinfo = this.ops.map((op) => {
      const { R, M, E, tE, tI } = splitOp(op, d);
      for (let i = 0; i < d; i++)
        for (let j = 0; j < 3; j++) {
          let lhs = fq(M[i][j]);
          for (let k = 0; k < d; k++) lhs = lhs.add(Qm[k][j].mulInt(E[i][k]));
          let rhs = Z;
          for (let k = 0; k < 3; k++) rhs = rhs.add(Qm[i][k].mulInt(R[k][j]));
          if (!lhs.eq(rhs)) throw new Error(`operator inconsistent with q: ${formatOp(op)}`);
        }
      const A = E.map((r, i) => r.map((x, j) => x - (i === j ? 1 : 0)));
      const c = tI.map((ti, i) => Qm[i].reduce((acc, x, j) => acc.sub(x.mul(tE[j])), ti));
      const conj = (X: number[][]) => matMul(matMul(this.Binv, X.map((r) => r.map((x) => fq(x)))), this.B);
      const At = intOfQMat(conj(A), "Λ not A-invariant");
      const Et = intOfQMat(conj(E), "Λ not E-invariant");
      const EtInv = intOfQMat(matInv(qOfBigMat(Et)), "Ẽ^-1");
      const ct = matVec(this.Binv, c);
      const nonzero = At.some((r) => r.some((x) => x !== 0n));
      return {
        R, M, E, tE, tI, A, c, At, Et, EtInv, ct,
        flats: nonzero ? solve(At, ct.map((x) => x.neg()), d) : null,
        always: nonzero ? false : ct.every((x) => x.isInt()),
      };
    });
  }

  tOfS(s: Q[]): Q[] { return matVec(this.B, s).map((x) => x.mod1()); }
  sOfT(t: Q[]): Q[] { return matVec(this.Binv, t).map((x) => x.mod1()); }

  /** [(op index, n)] for exact rational t (basic setting). */
  survivingAtT(t: Q[]): [number, number[]][] {
    const out: [number, number[]][] = [];
    this.opinfo.forEach((o, k) => {
      const v = o.A.map((row, i) => row.reduce((acc, a, j) => acc.add(t[j].mulInt(a)), o.c[i]).mod1());
      const n = this.nlookup.get(tupleKey(v));
      if (n) out.push([k, [...n]]);
    });
    return out;
  }

  // ---------------------------------------------------------------- closure
  private closure() {
    const torus = new Flat([], [], this.d);
    this.flats = new Map([[torus.key, torus]]);
    const queue: Flat[] = [];
    for (const o of this.opinfo) if (o.flats) queue.push(...o.flats);
    while (queue.length) {
      const F = queue.shift()!;
      if (this.flats.has(F.key)) continue;
      const existing = [...this.flats.values()];
      this.flats.set(F.key, F);
      for (const G of existing) for (const H of intersect(F, G)) if (!this.flats.has(H.key)) queue.push(H);
    }
    const fl = [...this.flats.values()];
    this.sub = new Map(fl.map((F) => [F.key, fl.filter((G) => G.key !== F.key && G.dim < F.dim && subset(G, F))]));
  }

  inOpenPart(F: Flat, s: Q[]): boolean {
    return F.contains(s) && !this.sub.get(F.key)!.some((G) => G.contains(s));
  }

  classifyS(s: Q[]): Flat {
    const cont = [...this.flats.values()].filter((F) => F.contains(s));
    const m = Math.min(...cont.map((F) => F.dim));
    const best = cont.filter((F) => F.dim === m);
    if (best.length !== 1) throw new Error("internal: point lies in several minimal flats");
    return best[0];
  }

  // ---------------------------------------------------------------- points
  genericPoint(F: Flat): Q[] {
    for (let a = 0n; a < 1000n; a++) {
      const free: Q[] = [];
      for (let m = 0; m < F.dim; m++) {
        let k = (GENERIC_SEEDS[m] + GENERIC_STEP * a) % P_GENERIC;
        if (k === 0n) k = 1n;
        free.push(fq(k, P_GENERIC));
      }
      const s = F.point(free);
      if (this.inOpenPart(F, s)) return s.map((x) => x.mod1());
    }
    throw new Error("no generic point found");
  }

  examplePoint(F: Flat): Q[] | null {
    const keyOf = (t: Q[]) => ({ den: denLcm(t), nz: t.filter((x) => !x.isZero()).length, t });
    const better = (a: ReturnType<typeof keyOf>, b: ReturnType<typeof keyOf>) =>
      a.den !== b.den ? a.den < b.den : a.nz !== b.nz ? a.nz < b.nz : qcmpTuple(a.t, b.t) < 0;
    for (let m = 1; m <= EXAMPLE_MAX_DEN; m++) {
      let best: ReturnType<typeof keyOf> | null = null;
      const ks = Array(F.dim).fill(0);
      for (;;) {
        const s = F.point(ks.map((k) => fq(k, m)));
        if (this.inOpenPart(F, s)) {
          const t0 = matVec(this.B, s);
          for (const lam of this.lamReps) {
            const cand = keyOf(t0.map((x, i) => x.add(lam[i]).mod1()));
            if (!best || better(cand, best)) best = cand;
          }
        }
        let i = F.dim - 1;
        while (i >= 0) {
          ks[i]++;
          if (ks[i] < m) break;
          ks[i] = 0;
          i--;
        }
        if (i < 0) break;
      }
      if (best) return best.t;
    }
    return null;
  }

  // ---------------------------------------------------------------- 3D group
  ops3(surv: [number, number[]][]): Op3[] {
    const out: Op3[] = [];
    const seen = new Set<string>();
    const N = this.N;
    for (const [k, n] of surv) {
      const o = this.opinfo[k];
      const Rs = o.R.map((r, i) => r.map((x, j) => {
        const v = fq(x * N[j], N[i]);
        if (!v.isInt()) throw new Error("internal: supercell not invariant");
        return Number(v.n);
      }));
      const tau = o.tE.map((x, i) => x.add(fq(n[i])));
      for (const cell of this.cells) {
        const t = tau.map((x, i) => x.add(fq(cell[i])).div(fq(N[i])).mod1());
        const key = Rs.flat().join(",") + "|" + tupleKey(t);
        if (!seen.has(key)) {
          seen.add(key);
          out.push({ R: Rs, t });
        }
      }
    }
    return out;
  }

  /** Flat image under s ↦ Ẽ s + c̃. */
  image(F: Flat, o: OpInfo): Flat {
    if (!F.C.length) return F;
    const Ci = F.C.map((row) => o.EtInv[0].map((_, j) => row.reduce((acc, x, k) => acc + x * o.EtInv[k][j], 0n)));
    const f = F.f.map((fi, i) => Ci[i].reduce((acc, a, j) => acc.add(o.ct[j].mulInt(a)), fi));
    return canonicalFlat(Ci, f, this.d);
  }

  // ---------------------------------------------------------------- classes
  private buildClasses() {
    const d = this.d;
    const classes: TClass[] = [];
    for (const F of this.flats.values()) {
      const g = this.genericPoint(F);
      const t = this.tOfS(g);
      const surv = this.survivingAtT(t);
      classes.push({
        flat: F, dim: F.dim, genericS: g, genericT: t, exampleT: [], surviving: surv,
        ops3: this.ops3(surv), orbit: -1, locus: "",
      });
    }
    // order: (-dim, generic_t); the reference also sorts by sg number and setting,
    // which needs identification — done by the caller (finder.ts) after identify().
    classes.sort((a, b) => b.dim - a.dim || qcmpTuple(a.genericT, b.genericT));
    // SSG orbits (union-find)
    const parent = new Map(classes.map((c) => [c.flat.key, c.flat.key]));
    const find = (x: string): string => {
      while (parent.get(x) !== x) {
        parent.set(x, parent.get(parent.get(x)!)!);
        x = parent.get(x)!;
      }
      return x;
    };
    for (const c of classes)
      for (const o of this.opinfo) {
        const img = this.image(c.flat, o);
        if (!parent.has(img.key)) throw new Error("internal: closure not invariant under the SSG");
        const a = find(c.flat.key), b = find(img.key);
        if (a !== b) parent.set(b, a);
      }
    const roots = new Map<string, number>();
    for (const c of classes) {
      const r = find(c.flat.key);
      if (!roots.has(r)) roots.set(r, roots.size);
      c.orbit = roots.get(r)!;
    }
    for (const c of classes) {
      const ex = this.examplePoint(c.flat);
      c.exampleT = ex ?? c.genericT;
      if (this.classifyS(this.sOfT(c.exampleT)) !== c.flat) throw new Error("internal: example point in wrong class");
    }
    for (const oid of new Set(classes.map((c) => c.orbit))) {
      const members = classes.filter((c) => c.orbit === oid);
      const txt = this.locusText(members);
      for (const c of members) c.locus = txt;
    }
    if (d === 1) {
      this.itN = Number(Q.ONE.div(this.B[0][0]).n);
      for (const c of classes) c.itLabel = this.itLabel(c);
    }
    this.classes = classes;
  }

  // ---------------------------------------------------------------- loci text
  tvar(j: number) { return this.d === 1 ? "t" : `t${j + 1}`; }

  linform(c: bigint[]): string {
    let s = "";
    c.forEach((x, j) => {
      if (x === 0n) return;
      const sign = !s ? (x < 0n ? "-" : "") : x < 0n ? " - " : " + ";
      const mag = absB(x) === 1n ? "" : String(absB(x));
      s += `${sign}${mag}${this.tvar(j)}`;
    });
    return s;
  }

  locusText(members: TClass[]): string {
    const d = this.d;
    if (members[0].dim === d) return "general";
    const groups = new Map<string, { Wb: QMat; R: Map<string, Q[]> }>();
    for (const c of members) {
      const K = matMul(qOfBigMat(c.flat.C), this.Binv);
      const den = denLcm(K.flat());
      const Wb = hnf(K.map((r) => r.map((x) => x.mul(fq(den)).n))).H.map((r) => r.map((x) => fq(x, den)));
      const key = Wb.map(tupleKey).join(";");
      const phi = Wb.map((r) => r.reduce((acc, w, j) => acc.add(w.mul(c.genericT[j])), Z).mod1());
      if (!groups.has(key)) groups.set(key, { Wb, R: new Map() });
      groups.get(key)!.R.set(tupleKey(phi), phi);
    }
    // Python sorts the group keys as tuples of tuples of Fractions
    const sorted = [...groups.values()].sort((a, b) => {
      for (let i = 0; i < Math.min(a.Wb.length, b.Wb.length); i++) {
        const c = qcmpTuple(a.Wb[i], b.Wb[i]);
        if (c) return c;
      }
      return a.Wb.length - b.Wb.length;
    });
    const texts = sorted.map((g) => this.groupText(g.Wb, [...g.R.values()]));
    if (texts.length === 1) return texts[0];
    return [...texts].sort().map((t) => (t.includes(",") && !t.startsWith("(") ? `[${t}]` : t)).join(" or ");
  }

  private groupText(Wb: QMat, Rlist: Q[][]): string {
    const r = Wb.length;
    const R = [...Rlist].sort(qcmpTuple);
    const Rset = new Set(R.map(tupleKey));
    const x0 = R[0];
    const H: Q[][] = [];
    for (const x of R) {
      const h = x.map((a, i) => a.sub(x0[i]).mod1());
      const shifted = new Set(R.map((y) => tupleKey(y.map((a, i) => a.add(h[i]).mod1()))));
      if (shifted.size === Rset.size && [...shifted].every((k) => Rset.has(k))) H.push(h);
    }
    const den = denLcm(H.flat());
    const gens: IMat = [
      ...Array.from({ length: r }, (_, i) => Array.from({ length: r }, (_, j) => (i === j ? den : 0n))),
      ...H.map((h) => h.map((x) => x.mul(fq(den)).n)),
    ];
    const L = hnfLower(gens).map((row) => row.map((x) => fq(x, den)));
    const dual = transpose(matInv(L));
    const Zm = hnf(intOfQMat(dual, "dual lattice")).H;
    const rows = Zm.map((z) => Array.from({ length: this.d }, (_, j) => z.reduce((acc, zk, k) => acc.add(Wb[k][j].mulInt(zk)), Z)));
    const valsMap = new Map<string, Q[]>();
    for (const x of R) {
      const v = Zm.map((z) => z.reduce((acc, zk, k) => acc.add(x[k].mulInt(zk)), Z).mod1());
      valsMap.set(tupleKey(v), v);
    }
    const vals = [...valsMap.values()].sort(qcmpTuple);
    const pmod = (a: Q, m: Q) => a.modQ(m);
    const conds = rows.map((w, i) => {
      const m0 = denLcm(w);
      const v = w.map((x) => x.mul(fq(m0)).n);
      const g = v.reduce((acc, x) => gcd(acc, x), 0n);
      let cvec = v.map((x) => x / g);
      const mu = fq(m0, g);
      const sign = cvec.find((x) => x !== 0n)! > 0n ? 1 : -1;
      cvec = cvec.map((x) => x * BigInt(sign));
      return { cvec, mu, sign, i };
    });
    const out: string[] = [];
    const varying = conds.map((c) => c.i).filter((i) => new Set(vals.map((v) => String(v[i]))).size > 1);
    for (const c of conds)
      if (!varying.includes(c.i)) {
        const val = pmod(vals[0][c.i].mulInt(c.sign).mul(c.mu), c.mu);
        out.push(`${this.linform(c.cvec)} ≡ ${val} (mod ${c.mu})`);
      }
    if (varying.length) {
      const cs = varying.map((i) => conds[i]);
      const tupMap = new Map<string, Q[]>();
      for (const v of vals) {
        const tp = cs.map((c) => pmod(v[c.i].mulInt(c.sign).mul(c.mu), c.mu));
        tupMap.set(tupleKey(tp), tp);
      }
      const tups = [...tupMap.values()].sort(qcmpTuple);
      if (cs.length === 1) {
        out.push(`${this.linform(cs[0].cvec)} ≡ ${tups.map((x) => String(x[0])).join(" or ")} (mod ${cs[0].mu})`);
      } else {
        const lhs = "(" + cs.map((c) => this.linform(c.cvec)).join(", ") + ")";
        const alts = tups.map((tp) => "(" + tp.map(String).join(", ") + ")").join(" or ");
        const mods = cs.map((c) => String(c.mu)).join(", ");
        out.push(`${lhs} ≡ ${alts} (mod ${mods})`);
      }
    }
    return out.join(", ");
  }

  itLabel(c: TClass): string {
    if (c.dim === 1) return "general";
    const a = c.genericT[0].mulInt(this.itN!).mod1();
    const names: Record<string, string> = { "0": "n/N", "1/4": "1/4N+n/N", "1/2": "1/2N+n/N", "3/4": "3/4N+n/N" };
    return names[String(a)] ?? `${a}/N+n/N`;
  }
}

export type { QVec };
