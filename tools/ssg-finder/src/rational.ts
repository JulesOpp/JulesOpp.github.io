// Exact rational arithmetic on BigInt. Values are immutable and always reduced,
// with a positive denominator.

export class Q {
  readonly n: bigint;
  readonly d: bigint;

  private constructor(n: bigint, d: bigint) {
    this.n = n;
    this.d = d;
  }

  static of(n: bigint | number, d: bigint | number = 1n): Q {
    let nn = BigInt(n);
    let dd = BigInt(d);
    if (dd === 0n) throw new RangeError("zero denominator");
    if (dd < 0n) {
      nn = -nn;
      dd = -dd;
    }
    const g = gcd(nn < 0n ? -nn : nn, dd);
    return new Q(nn / g, dd / g);
  }

  /** Parse "3", "-1/4", "+2/6", "0.25" (finite decimals only). */
  static parse(s: string): Q {
    const t = s.trim().replace(/^\+/, "");
    let m = /^(-?\d+)(?:\/(\d+))?$/.exec(t);
    if (m) return Q.of(BigInt(m[1]), m[2] ? BigInt(m[2]) : 1n);
    m = /^(-?)(\d*)\.(\d+)$/.exec(t);
    if (m) {
      const den = 10n ** BigInt(m[3].length);
      const num = BigInt((m[2] || "0") + m[3]);
      return Q.of(m[1] === "-" ? -num : num, den);
    }
    throw new SyntaxError(`not a rational number: "${s}"`);
  }

  static readonly ZERO = Q.of(0n);
  static readonly ONE = Q.of(1n);

  add(o: Q): Q { return Q.of(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o: Q): Q { return Q.of(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o: Q): Q { return Q.of(this.n * o.n, this.d * o.d); }
  div(o: Q): Q {
    if (o.n === 0n) throw new RangeError("division by zero");
    return Q.of(this.n * o.d, this.d * o.n);
  }
  neg(): Q { return new Q(-this.n, this.d); }
  mulInt(k: bigint | number): Q { return Q.of(this.n * BigInt(k), this.d); }

  isZero(): boolean { return this.n === 0n; }
  isInt(): boolean { return this.d === 1n; }
  eq(o: Q): boolean { return this.n === o.n && this.d === o.d; }
  cmp(o: Q): number {
    const l = this.n * o.d, r = o.n * this.d;
    return l < r ? -1 : l > r ? 1 : 0;
  }
  sign(): number { return this.n < 0n ? -1 : this.n > 0n ? 1 : 0; }

  floor(): bigint {
    const q = this.n / this.d;
    return this.n < 0n && q * this.d !== this.n ? q - 1n : q;
  }
  /** Fractional part in [0,1). */
  mod1(): Q { return Q.of(this.n - this.floor() * this.d, this.d); }
  /** Reduce into [0, m) for a positive rational modulus m. */
  modQ(m: Q): Q {
    const k = this.div(m).floor();
    return this.sub(m.mulInt(k));
  }

  toNumber(): number { return Number(this.n) / Number(this.d); }
  toString(): string { return this.d === 1n ? this.n.toString() : `${this.n}/${this.d}`; }
  toJSON(): string { return this.toString(); }
}

export function gcd(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) [a, b] = [b, a % b];
  return a;
}

export function lcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const r = (a / gcd(a, b)) * b;
  return r < 0n ? -r : r;
}

export const q = (s: string | number | bigint): Q =>
  typeof s === "string" ? Q.parse(s) : Q.of(BigInt(s));

export type QVec = Q[];
export type QMat = Q[][];

export const qvec = (xs: (string | number)[]): QVec => xs.map((x) => q(x as string));
export const qmat = (rows: (string | number)[][]): QMat => rows.map(qvec);

export function denLcm(xs: Iterable<Q>): bigint {
  let l = 1n;
  for (const x of xs) l = lcm(l, x.d);
  return l;
}

export function vecAdd(a: QVec, b: QVec): QVec { return a.map((x, i) => x.add(b[i])); }
export function vecSub(a: QVec, b: QVec): QVec { return a.map((x, i) => x.sub(b[i])); }
export function vecMod1(a: QVec): QVec { return a.map((x) => x.mod1()); }
export function vecIsInt(a: QVec): boolean { return a.every((x) => x.isInt()); }
export function vecEq(a: QVec, b: QVec): boolean { return a.length === b.length && a.every((x, i) => x.eq(b[i])); }
export function vecCmp(a: QVec, b: QVec): number {
  for (let i = 0; i < a.length; i++) {
    const c = a[i].cmp(b[i]);
    if (c) return c;
  }
  return a.length - b.length;
}

export function matVec(A: QMat, v: QVec): QVec {
  return A.map((row) => row.reduce((s, a, j) => s.add(a.mul(v[j])), Q.ZERO));
}
export function matMul(A: QMat, B: QMat): QMat {
  const m = B[0]?.length ?? 0;
  return A.map((row) =>
    Array.from({ length: m }, (_, j) => row.reduce((s, a, k) => s.add(a.mul(B[k][j])), Q.ZERO)),
  );
}
export function matEq(A: QMat, B: QMat): boolean {
  return A.length === B.length && A.every((r, i) => vecEq(r, B[i]));
}
export function transpose<T>(A: T[][]): T[][] {
  return A.length ? A[0].map((_, j) => A.map((r) => r[j])) : [];
}
export function identity(n: number): QMat {
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? Q.ONE : Q.ZERO)));
}

/** Exact inverse by Gauss–Jordan; throws if singular. */
export function matInv(A: QMat): QMat {
  const n = A.length;
  const M = A.map((r, i) => [...r, ...identity(n)[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    while (p < n && M[p][c].isZero()) p++;
    if (p === n) throw new RangeError("singular matrix");
    [M[c], M[p]] = [M[p], M[c]];
    const inv = Q.ONE.div(M[c][c]);
    M[c] = M[c].map((x) => x.mul(inv));
    for (let r = 0; r < n; r++) {
      if (r !== c && !M[r][c].isZero()) {
        const f = M[r][c];
        M[r] = M[r].map((x, j) => x.sub(f.mul(M[c][j])));
      }
    }
  }
  return M.map((r) => r.slice(n));
}
