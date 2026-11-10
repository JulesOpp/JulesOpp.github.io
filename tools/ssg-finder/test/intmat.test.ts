import { describe, expect, it } from "vitest";
import fc from "fast-check";
import { hnf, snf, imatMul, imatDet, saturate, kernel, unimodularInverse, imatEq, IMat } from "../src/intmat";

const big = (A: number[][]): IMat => A.map((r) => r.map(BigInt));
const matArb = (m: number, n: number) =>
  fc.array(fc.array(fc.integer({ min: -12, max: 12 }), { minLength: n, maxLength: n }), { minLength: m, maxLength: m }).map(big);
const shapeArb = fc.tuple(fc.integer({ min: 1, max: 6 }), fc.integer({ min: 1, max: 6 })).chain(([m, n]) => matArb(m, n));

describe("snf", () => {
  it("U·A·V = D, unimodular, divisibility chain", () => {
    fc.assert(
      fc.property(shapeArb, (A) => {
        const { D, U, V, diag } = snf(A);
        expect(imatEq(imatMul(imatMul(U, A), V), D)).toBe(true);
        expect([1n, -1n]).toContain(imatDet(U));
        expect([1n, -1n]).toContain(imatDet(V));
        D.forEach((r, i) => r.forEach((x, j) => { if (i !== j || i >= diag.length) expect(x).toBe(0n); }));
        diag.forEach((d, i) => { expect(d > 0n).toBe(true); if (i) expect(d % diag[i - 1]).toBe(0n); });
      }),
      { numRuns: 400, endOnFailure: true },
    );
  });
  it("known example", () => {
    expect(snf(big([[2, 4, 4], [-6, 6, 12], [10, -4, -16]])).diag).toEqual([2n, 6n, 12n]);
  });
});

describe("hnf", () => {
  it("U·A = [H;0], echelon, reduced, and invariant under unimodular row ops", () => {
    fc.assert(
      fc.property(shapeArb, (A) => {
        const { H, U, pivots } = hnf(A);
        const UA = imatMul(U, A);
        expect(imatEq(UA.slice(0, H.length), H)).toBe(true);
        UA.slice(H.length).forEach((r) => r.forEach((x) => expect(x).toBe(0n)));
        expect([1n, -1n]).toContain(imatDet(U));
        pivots.forEach((c, i) => {
          expect(H[i][c] > 0n).toBe(true);
          for (let k = 0; k < i; k++) expect(H[k][c] >= 0n && H[k][c] < H[i][c]).toBe(true);
          for (let k = i + 1; k < H.length; k++) expect(H[k][c]).toBe(0n);
        });
        // uniqueness: scrambling rows with a random unimodular matrix gives the same HNF
        const m = A.length;
        const S = big(Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : j === i + 1 ? 3 : 0))));
        expect(imatEq(hnf(imatMul(S, A)).H, H)).toBe(true);
      }),
      { numRuns: 400, endOnFailure: true },
    );
  });
});

describe("saturate / kernel / inverse", () => {
  it("saturation of 2·(1,1,0),(0,0,4) is (1,1,0),(0,0,1)", () => {
    expect(saturate(big([[2, 2, 0], [0, 0, 4]]))).toEqual(big([[1, 1, 0], [0, 0, 1]]));
  });
  it("kernel vectors are annihilated", () => {
    fc.assert(
      fc.property(shapeArb, (A) => {
        const n = A[0].length;
        for (const k of kernel(A, n)) A.forEach((r) => expect(r.reduce((s, a, j) => s + a * k[j], 0n)).toBe(0n));
        const rank = snf(A).diag.length;
        expect(kernel(A, n).length).toBe(n - rank);
      }),
      { numRuns: 300, endOnFailure: true },
    );
  });
  it("unimodular inverse", () => {
    const A = big([[2, 3, 1], [1, 2, 1], [1, 1, 1]]);
    expect(imatMul(A, unimodularInverse(A))).toEqual(big([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
  });
});
