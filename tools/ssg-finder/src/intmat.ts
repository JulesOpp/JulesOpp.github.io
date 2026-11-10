// Exact integer matrix algorithms on BigInt: row Hermite normal form, Smith normal
// form with unimodular transforms, saturation. Matrices are arrays of rows.

export type IMat = bigint[][];
export type IVec = bigint[];

const abs = (x: bigint) => (x < 0n ? -x : x);

/** Floor division for bigint (rounds toward -infinity). */
export function fdiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return (a % b !== 0n && (a < 0n) !== (b < 0n)) ? q - 1n : q;
}

export function imatIdentity(n: number): IMat {
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1n : 0n)));
}
export const imatClone = (A: IMat): IMat => A.map((r) => [...r]);
export function imatMul(A: IMat, B: IMat): IMat {
  const m = B[0]?.length ?? 0;
  return A.map((row) => Array.from({ length: m }, (_, j) => row.reduce((s, a, k) => s + a * B[k][j], 0n)));
}
export function imatEq(A: IMat, B: IMat): boolean {
  return A.length === B.length && A.every((r, i) => r.length === B[i].length && r.every((x, j) => x === B[i][j]));
}
export function imatTranspose(A: IMat): IMat {
  return A.length ? A[0].map((_, j) => A.map((r) => r[j])) : [];
}

/** Determinant (Bareiss, exact). */
export function imatDet(A: IMat): bigint {
  const n = A.length;
  if (n === 0) return 1n;
  const M = imatClone(A);
  let sign = 1n, prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (M[k][k] === 0n) {
      let p = k + 1;
      while (p < n && M[p][k] === 0n) p++;
      if (p === n) return 0n;
      [M[k], M[p]] = [M[p], M[k]];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++)
      for (let j = k + 1; j < n; j++) M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev;
    prev = M[k][k];
  }
  return sign * M[n - 1][n - 1];
}

/**
 * Row-style Hermite normal form: returns H (zero rows removed) and unimodular U with
 * U·A = [H; 0]. H is in row echelon form, pivots positive, entries above each pivot
 * reduced into [0, pivot). Unique for a given row lattice.
 */
export function hnf(A: IMat): { H: IMat; U: IMat; pivots: number[] } {
  const m = A.length, n = A[0]?.length ?? 0;
  const M = imatClone(A);
  const U = imatIdentity(m);
  const swap = (i: number, j: number) => { [M[i], M[j]] = [M[j], M[i]]; [U[i], U[j]] = [U[j], U[i]]; };
  const addRow = (dst: number, src: number, k: bigint) => {
    if (k === 0n) return;
    for (let c = 0; c < n; c++) M[dst][c] += k * M[src][c];
    for (let c = 0; c < m; c++) U[dst][c] += k * U[src][c];
  };
  const negRow = (i: number) => { M[i] = M[i].map((x) => -x); U[i] = U[i].map((x) => -x); };
  const pivots: number[] = [];
  let r = 0;
  for (let c = 0; c < n && r < m; c++) {
    // Euclid on column c among rows r..m-1
    for (;;) {
      let best = -1;
      for (let i = r; i < m; i++) if (M[i][c] !== 0n && (best < 0 || abs(M[i][c]) < abs(M[best][c]))) best = i;
      if (best < 0) break;
      swap(r, best);
      let done = true;
      for (let i = r + 1; i < m; i++) {
        if (M[i][c] !== 0n) {
          addRow(i, r, -(M[i][c] / M[r][c]));
          if (M[i][c] !== 0n) done = false;
        }
      }
      if (done) break;
    }
    if (M[r][c] === 0n) continue;
    if (M[r][c] < 0n) negRow(r);
    for (let i = 0; i < r; i++) addRow(i, r, -fdiv(M[i][c], M[r][c]));
    pivots.push(c);
    r++;
  }
  return { H: M.slice(0, r), U, pivots };
}

/**
 * Smith normal form: U·A·V = D with U (m×m), V (n×n) unimodular, D diagonal,
 * d_1 | d_2 | … | d_r > 0, zeros after. Returns diag = [d_1..d_r] (r = rank).
 * Alternates row- and column-HNF (which keeps entries reduced), then repairs the
 * divisibility chain by adding columns and repeating.
 */
export function snf(A: IMat): { D: IMat; U: IMat; V: IMat; diag: bigint[] } {
  const m = A.length, n = A[0]?.length ?? 0;
  let D = imatClone(A);
  let U = imatIdentity(m), V = imatIdentity(n);
  const pad = (H: IMat, rows: number, cols: number): IMat =>
    Array.from({ length: rows }, (_, i) => (i < H.length ? [...H[i]] : Array(cols).fill(0n)));
  const isDiag = (M: IMat) => M.every((r, i) => r.every((x, j) => i === j || x === 0n));
  const diagonalise = () => {
    for (let guard = 0; !isDiag(D); guard++) {
      if (guard > 10000) throw new Error("snf: no convergence");
      const r = hnf(D);
      D = pad(r.H, m, n);
      U = imatMul(r.U, U);
      if (isDiag(D)) break;
      const c = hnf(imatTranspose(D));
      D = imatTranspose(pad(c.H, n, m));
      V = imatMul(V, imatTranspose(c.U));
    }
  };
  if (m && n) diagonalise();
  const k = Math.min(m, n);
  for (let guard = 0; ; guard++) {
    if (guard > 10000) throw new Error("snf: divisibility repair did not converge");
    let fixed = true;
    outer: for (let i = 0; i < k; i++)
      for (let j = i + 1; j < k; j++) {
        const di = D[i][i], dj = D[j][j];
        if (dj !== 0n && (di === 0n || dj % di !== 0n)) {
          // add column j to column i, then re-diagonalise (yields gcd, lcm)
          for (const row of D) row[i] += row[j];
          for (const row of V) row[i] += row[j];
          fixed = false;
          break outer;
        }
      }
    if (fixed) break;
    diagonalise();
  }
  for (let i = 0; i < k; i++)
    if (D[i][i] < 0n) {
      D[i] = D[i].map((x) => -x);
      U[i] = U[i].map((x) => -x);
    }
  const diag: bigint[] = [];
  for (let i = 0; i < k && D[i][i] !== 0n; i++) diag.push(D[i][i]);
  return { D, U, V, diag };
}

/**
 * Lower-triangular HNF basis of the row lattice of A (full column rank expected):
 * pivots on the diagonal, entries below each pivot in [0, pivot). Reverse the
 * columns, take the row HNF, reverse columns and rows back (as reference/intmat.py).
 */
export function hnfLower(A: IMat): IMat {
  const R = A.map((r) => [...r].reverse());
  const H = hnf(R).H.map((r) => [...r].reverse());
  return H.reverse();
}

/** Rank of an integer matrix. */
export function imatRank(A: IMat): number {
  return A.length ? hnf(A).H.length : 0;
}

/** Inverse of a unimodular matrix (exact). */
export function unimodularInverse(A: IMat): IMat {
  const n = A.length;
  const aug = A.map((r, i) => [...r, ...imatIdentity(n)[i]]);
  const { H } = hnf(aug);
  // For unimodular A, the HNF of [A | I] is [I | A^{-1}].
  if (H.length !== n || !H.every((r, i) => r.slice(0, n).every((x, j) => x === (i === j ? 1n : 0n))))
    throw new RangeError("matrix is not unimodular");
  return H.map((r) => r.slice(n));
}

/**
 * Saturation: a basis (in HNF) of span_Q(rows) ∩ Z^n — the primitive lattice
 * containing the row lattice of A with the same rank.
 */
export function saturate(A: IMat): IMat {
  if (A.length === 0) return [];
  const { V, diag } = snf(A);
  const Vinv = unimodularInverse(V);
  return hnf(Vinv.slice(0, diag.length)).H;
}

/** Integer kernel basis: rows k with A·k = 0 (k as a row vector), in HNF. */
export function kernel(A: IMat, n: number): IMat {
  if (A.length === 0) return imatIdentity(n);
  const { V, diag } = snf(A);
  const cols = imatTranspose(V).slice(diag.length);
  return cols.length ? hnf(cols).H : [];
}
