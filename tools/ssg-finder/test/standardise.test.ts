import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { Q } from "../src/rational";
import { buildOps, GroupRecord, paramsOf, parseLinear, splitOp, type SOp } from "../src/ssg";
import { loadMoyo } from "../src/moyo";
import { tSubgroups } from "../src/subgroups";
import { crystalSystemOf } from "../src/symbol";
import { formatQ, standardise, transformOps, type Standardised } from "../src/standardise";

// P-3c1(0,0,g)000, IT A No. 165 generators in hexagonal axes with the 4th coordinate.
const P3c1: GroupRecord = {
  sym: "P-3c1(0,0,g)000", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
  gen: ["y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"],
};

const isOk = (r: ReturnType<typeof standardise>): r is Standardised => !("ok" in r && r.ok === false);

beforeAll(async () => {
  await loadMoyo();
});

describe("standardise (public)", () => {
  it("a group already in its standard setting keeps its symbol", () => {
    const r = standardise(buildOps(P3c1), P3c1.q);
    expect(isOk(r)).toBe(true);
    if (isOk(r)) {
      expect(r.sgNumber).toBe(165);
      expect(r.symbol.full).toBe("P-3c1(0,0,g)000");
      expect(r.qActual).toBe("0,0,g");
    }
  });
  it("P-3c1: its C2/c (index 3) subgroup class is monoclinic No. 15", () => {
    const cls = tSubgroups(P3c1).find((c) => c.index === 3)!;
    const r = standardise(cls.ops, P3c1.q);
    expect(isOk(r)).toBe(true);
    if (isOk(r)) {
      expect(r.sgNumber).toBe(15);
      expect(crystalSystemOf(r.sgNumber)).toBe("monoclinic");
      expect(r.symbol.basic).toMatch(/^[PABCI]2(_1)?\/[abcmn]$/);
    }
  });
  it("the returned map takes the operations into the reported setting", () => {
    const cls = tSubgroups(P3c1).find((c) => c.index === 3)!;
    const r = standardise(cls.ops, P3c1.q);
    if (!isOk(r)) throw new Error("not ok");
    const S = transformOps(cls.ops, r.P.map((row) => row.map(Q.parse)), r.p.map(Q.parse), r.eta, r.m.map((x) => Q.of(x)), Q.parse(r.p4));
    expect(S).not.toBeNull();
    const again = standardise(S!, [r.qActual]);
    expect(isOk(again) && again.symbol).toEqual(r.symbol);
  });
  it("refuses d != 1", () => {
    const ops = buildOps({ ...P3c1, d: 2, gen: ["y,-x+y,-z,-t,-u"] });
    const r = standardise(ops, ["0,0,g", "0,0,g"]);
    expect(r).toMatchObject({ ok: false });
    expect((r as { reason: string }).reason).toMatch(/d = 1/);
  });
});

// ------------------------------------------------------------------ local-only (BYU-derived oracle)

const LOCAL = join(__dirname, "..", "local");
const DATA = join(LOCAL, "data");
const SUBS = join(LOCAL, "fixtures", "subgroups");
const have = existsSync(DATA) && existsSync(SUBS);

function allD1(): (GroupRecord & { label: string })[] {
  const out: (GroupRecord & { label: string })[] = [];
  for (const f of readdirSync(DATA).filter((f) => /^sg\d{3}\.json$/.test(f)).sort()) {
    const groups: Record<string, GroupRecord> = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups;
    for (const [label, r] of Object.entries(groups)) if (r.d === 1) out.push({ ...r, label });
  }
  return out;
}

/** Deterministic PRNG (mulberry32). */
function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const mul3 = (A: number[][], B: number[][]) => A.map((r) => [0, 1, 2].map((j) => r.reduce((s, x, k) => s + x * B[k][j], 0)));
const det3 = (A: number[][]) =>
  A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
const I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
/** All 24 proper signed permutation matrices. */
const SIGNED_PERMS: number[][][] = (() => {
  const out: number[][][] = [];
  for (const p of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]])
    for (let s = 0; s < 8; s++) {
      const M = p.map((j, i) => [0, 1, 2].map((k) => (k === j ? ((s >> i) & 1 ? -1 : 1) : 0)));
      if (det3(M) === 1) out.push(M);
    }
  return out;
})();
const HEX6 = [[1, -1, 0], [1, 0, 0], [0, 0, 1]]; // 6-fold about c, hexagonal axes
const HEXSWAP = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]; // a <-> b, c -> -c

/** A random proper setting change allowed by the crystal system of the basic group. */
function randomLinear(sg: number, rnd: () => number): number[][] {
  const pick = <T>(xs: T[]) => xs[Math.floor(rnd() * xs.length)];
  const sys = crystalSystemOf(sg);
  if (sys === "triclinic" || sys === "monoclinic") {
    // integer shear(s) keeping the origin, then (monoclinic) an axis permutation
    let M = I3;
    for (let k = 0; k < 2; k++) {
      const E = I3.map((r) => [...r]);
      const [i, j] = pick(sys === "monoclinic" ? [[0, 1], [1, 0]] : [[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]]);
      E[i][j] = pick([-1, 1]);
      M = mul3(E, M);
    }
    return mul3(pick(SIGNED_PERMS), M);
  }
  if (sys === "orthorhombic") return pick(SIGNED_PERMS);
  if (sys === "tetragonal") return pick([I3, [[0, -1, 0], [1, 0, 0], [0, 0, 1]], HEXSWAP, [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]]);
  // trigonal / hexagonal: powers of the 6-fold, optionally with the a <-> b swap
  let M = I3;
  for (let k = Math.floor(rnd() * 6); k > 0; k--) M = mul3(HEX6, M);
  return rnd() < 0.5 ? M : mul3(HEXSWAP, M);
}

/** Apply x' = Px + p, x4' = m.x + eta x4 + p4 to a record: (ops, q text). */
function randomSetting(r: GroupRecord, sg: number, rnd: () => number): { ops: SOp[]; q: string } {
  const P = randomLinear(sg, rnd).map((row) => row.map((x) => Q.of(x)));
  const frac = () => Q.of(Math.floor(rnd() * 12), 12);
  const p = [frac(), frac(), frac()];
  const eta: 1 | -1 = rnd() < 0.5 ? 1 : -1;
  const m = [0, 1, 2].map(() => Q.of(Math.floor(rnd() * 3) - 1));
  const ops = transformOps(buildOps(r), P, p, eta, m, frac());
  if (!ops) throw new Error(`setting change failed for ${r.sym}`);
  // q' = (eta q + m) P^-1
  const names = paramsOf(r.q);
  const q = r.q[0].split(",").map((c) => parseLinear(c, names));
  const Pinv = inv3(P);
  const qn = [0, 1, 2].map((j) => ({
    coef: names.map((_, k) => [0, 1, 2].reduce((s, i) => s.add(q[i].coef[k].mulInt(eta).mul(Pinv[i][j])), Q.ZERO)),
    c: [0, 1, 2].reduce((s, i) => s.add(q[i].c.mulInt(eta).add(m[i]).mul(Pinv[i][j])), Q.ZERO),
  }));
  return { ops, q: formatQ(qn, names) };
}

function inv3(A: Q[][]): Q[][] {
  const M = A.map((r, i) => [...r, ...I3[i].map((x) => Q.of(x))]);
  for (let c = 0; c < 3; c++) {
    let p = c;
    while (M[p][c].isZero()) p++;
    [M[c], M[p]] = [M[p], M[c]];
    const f = Q.ONE.div(M[c][c]);
    M[c] = M[c].map((x) => x.mul(f));
    for (let r = 0; r < 3; r++) if (r !== c && !M[r][c].isZero()) {
      const g = M[r][c];
      M[r] = M[r].map((x, j) => x.sub(g.mul(M[c][j])));
    }
  }
  return M.map((r) => r.slice(3));
}

function report(bad: string[], total: number): string {
  const groups = new Map<string, number>();
  for (const b of bad) {
    const k = b.split("\t")[0];
    groups.set(k, (groups.get(k) ?? 0) + 1);
  }
  return `${bad.length}/${total} mismatches; by space group: ${[...groups].map(([k, n]) => `${k}:${n}`).join(" ")}\n` + bad.slice(0, 30).join("\n");
}

describe.skipIf(!have)("standardise vs the BYU (3+1)D tables (local-only)", () => {
  it("(a) every record, after 3 seeded random setting changes, standardises to its own symbol", () => {
    const recs = allD1();
    expect(recs).toHaveLength(775);
    const rnd = rng(20261110);
    const bad: string[] = [];
    let n = 0;
    for (const r of recs) {
      const sg = Number(r.label.split(".")[0]);
      for (let k = 0; k < 3; k++) {
        n++;
        const { ops, q } = randomSetting(r, sg, rnd);
        let got: string;
        try {
          const s = standardise(ops, [q]);
          got = isOk(s) ? s.symbol.full : `fail ${s.reason}`;
        } catch (e) {
          got = `throw ${(e as Error).message}`;
        }
        if (got !== r.sym) bad.push(`${sg}\t${r.label}\t${r.sym}\t(q ${q})\tgot ${got}`);
      }
    }
    expect(bad, report(bad, n)).toEqual([]);
  });

  it("(b) every k-index-1 subgroup class of every record standardises to its BYU type", () => {
    const recs = new Map(allD1().map((r) => [r.label, r]));
    const files = readdirSync(SUBS).filter((f) => f.endsWith(".json") && !f.startsWith("_"));
    expect(files).toHaveLength(775);
    const bad: string[] = [];
    let n = 0;
    for (const f of files) {
      const fx = JSON.parse(readFileSync(join(SUBS, f), "utf8"));
      const parent = recs.get(fx.label)!;
      const ops = buildOps(parent);
      for (const cls of fx.classes) {
        n++;
        const keep = new Set((cls.pg_ops as number[]).map((i) => splitOp(ops[i], 1).R.join(";")));
        const H = ops.filter((o) => keep.has(splitOp(o, 1).R.join(";")));
        let got: string;
        try {
          const s = standardise(H, parent.q);
          got = isOk(s) ? s.symbol.full : `fail ${s.reason}`;
        } catch (e) {
          got = `throw ${(e as Error).message}`;
        }
        if (got !== cls.byu_sym) bad.push(`${cls.byu.split(".")[0]}\t${fx.label} ${fx.sym} [${cls.pg_ops}]\t${cls.byu_sym}\tgot ${got}`);
      }
    }
    expect(n).toBe(8204);
    expect(bad, report(bad, n)).toEqual([]);
  });
});
