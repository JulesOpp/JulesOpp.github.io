import { describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { buildOps, GroupRecord, opMul, formatOp } from "../src/ssg";
import { pointGroup, tSubgroups } from "../src/subgroups";

// P-3c1(0,0,g)000, IT A No. 165 generators in hexagonal axes with the 4th coordinate.
// q = (0,0,g): R^T q = eps q gives eps = -1 for -3 (z -> -z) and for the 2-fold along a (z -> -z).
const P3c1: GroupRecord = {
  sym: "P-3c1(0,0,g)000", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
  gen: ["y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"],
};

// IT A generators: unique-axis-b P2/m, axial mirrors Pmmm, and P4/mmm.
const handRecords = [P3c1,
  { ...P3c1, sym: "P2/m", gen: ["-x,y,-z,-t", "x,-y,z,t"] },
  { ...P3c1, sym: "Pmmm", gen: ["-x,y,z,t", "x,-y,z,t", "x,y,-z,-t"] },
  { ...P3c1, sym: "P4/mmm", gen: ["-y,x,z,t", "x,-y,z,t", "x,y,-z,-t"] },
];

// Independently compare each class to every conjugate of every larger class.
// Matrix products, rather than subgroup implementation indices, determine containment.
function checkMaximal(rec: GroupRecord) {
  const K = pointGroup(buildOps(rec), rec.d).map(p => p.R);
  const mul = (a: number[][], b: number[][]) => a.map(row => [0,1,2].map(j => row.reduce((s,x,k) => s+x*b[k][j],0)));
  const rk = (a: number[][]) => JSON.stringify(a);
  const id = rk([[1,0,0],[0,1,0],[0,0,1]]);
  const classes = tSubgroups(rec);
  for (const p of classes) {
    const keys = p.pointOps.map(o => rk(o.R));
    const intermediate = classes.some(q => q.pointOps.length > keys.length && q.pointOps.length < K.length &&
      K.some(g => {
        const inv = K.find(h => rk(mul(g,h)) === id)!;
        const conjugate = new Set(q.pointOps.map(h => rk(mul(mul(g,h.R),inv))));
        return keys.every(k => conjugate.has(k));
      }));
    expect(p.maximal, `${rec.sym}: index ${p.index}`).toBe(keys.length < K.length && !intermediate);
  }
}

it("maximal, all classes: four hand groups", () => handRecords.forEach(checkMaximal));

const key = (o: { W: number[][]; w: { toString(): string }[] }) => JSON.stringify(o.W) + o.w.map(String).join(",");

function checkClosedAndIndex(rec: GroupRecord) {
  const ops = buildOps(rec);
  const classes = tSubgroups(rec);
  for (const c of classes) {
    expect(ops.length).toBe(c.index * c.ops.length);
    const have = new Set(c.ops.map(key));
    expect(have.size).toBe(c.ops.length);
    for (const a of c.ops) for (const b of c.ops) expect(have.has(key(opMul(a, b)))).toBe(true);
  }
  return classes;
}

describe("tSubgroups (hand-written P-3c1(0,0,g)000)", () => {
  it("has 10 classes with the expected indices and sizes", () => {
    const cl = tSubgroups(P3c1);
    expect(cl.map((c) => c.index)).toEqual([1, 2, 2, 2, 3, 4, 6, 6, 6, 12]);
    expect(cl.reduce((s, c) => s + c.classSize, 0)).toBe(16);
    const three = cl.filter((c) => c.classSize === 3).map((c) => c.pointOps.length).sort();
    expect(three).toEqual([2, 2, 4]);
    expect(cl.filter((c) => c.classSize === 1)).toHaveLength(7);
    expect(cl[0].pointOps).toHaveLength(12);
    expect(cl[0].maximal).toBe(false);
    expect(cl[0].pointOps[0].E[0][0]).toBe(1);
  });
  it("closed under opMul, order index times smaller", () => {
    checkClosedAndIndex(P3c1);
  });
  it("pointGroup returns identity first with E = [[+-1]]", () => {
    const pg = pointGroup(buildOps(P3c1), 1);
    expect(pg).toHaveLength(12);
    expect(pg[0].R).toEqual([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    expect(pg.every((p) => p.E.length === 1 && p.E[0].length === 1 && Math.abs(p.E[0][0]) === 1)).toBe(true);
  });
  it("marks maximal subgroups", () => {
    const cl = tSubgroups(P3c1);
    const byOrder = (n: number) => cl.filter((c) => c.pointOps.length === n);
    // -3m1 (12) is the parent; 3m1, 321, -3 (order 6) are maximal; trivial is not.
    expect(byOrder(6).every((c) => c.maximal)).toBe(true);
    expect(byOrder(1)[0].maximal).toBe(false);
  });
  it("refuses d > 3 before reading operators", () => {
    expect(() => tSubgroups({ ...P3c1, d: 4 })).toThrow("d > 3 is not supported");
    expect(() => pointGroup([], 4)).toThrow("d > 3 is not supported");
  });
});

const LOCAL = join(__dirname, "..", "local");
const FIX = join(LOCAL, "fixtures", "subgroups", "d1");
const DATA = join(LOCAL, "data");
const have = existsSync(FIX) && existsSync(DATA);

function allRecords(): GroupRecord[] {
  const out: GroupRecord[] = [];
  for (const f of readdirSync(DATA).filter((f) => /^sg\d+\.json$/.test(f))) {
    const groups = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups as Record<string, GroupRecord>;
    for (const [label, r] of Object.entries(groups)) if (r.d === 1) out.push({ ...r, label });
  }
  return out;
}

describe.skipIf(!have)("tSubgroups vs oracle fixtures", () => {
  it("all d = 1 records: class count, index and class_size multisets", () => {
    const recs = allRecords();
    expect(recs.length).toBe(775);
    for (const rec of recs) {
      const fx = JSON.parse(readFileSync(join(FIX, `${rec.label}.json`), "utf8"));
      const cl = tSubgroups(rec);
      const ms = (a: number[]) => [...a].sort((x, y) => x - y);
      expect(cl.length, rec.label).toBe(fx.classes.length);
      expect(ms(cl.map((c) => c.index)), rec.label).toEqual(ms(fx.classes.map((c: any) => c.index)));
      expect(ms(cl.map((c) => c.classSize)), rec.label).toEqual(ms(fx.classes.map((c: any) => c.class_size)));
    }
  });
  it("maximal, all classes: all 775 local parents", () => {
    const recs = allRecords();
    expect(recs).toHaveLength(775);
    recs.forEach(checkMaximal);
  });
  it("15.1.7.3 keeps its centring translation in every class", () => {
    const rec = allRecords().find((r) => r.label === "15.1.7.3")!;
    expect(rec.cen).toHaveLength(1);
    const cen = rec.cen[0];
    const classes = checkClosedAndIndex(rec);
    for (const c of classes) expect(c.ops.map(formatOp)).toContain(centringStr(cen));
  });
});

function centringStr(cen: string): string {
  // "1/2,0,1/2,0" -> "x+1/2,y,z+1/2,t"
  return cen.split(",").map((c, i) => "xyzt"[i] + (c === "0" ? "" : "+" + c)).join(",");
}


// IT A No. 191, hexagonal axes: 6_z, m_(vertical), m_(horizontal).
// The two independent q parameters are parallel to c; their E blocks are ±I2.
// These written IT A operations contain no local-catalogue data.
const P6mmm2: GroupRecord = {
  sym: "P6/mmm(0,0,g1)0000(0,0,g2)0000", cls: "", d: 2,
  q: ["0,0,g1", "0,0,g2"], params: ["g1", "g2"], cen: [],
  gen: ["x-y,x,z,t,u", "y,x,z,t,u", "x,y,-z,-t,-u"],
};

// Signed axial permutations: cubic m-3m with the three q arms along a,b,c.
// E has off-diagonal entries, so scalar/sign-only bookkeeping cannot pass.
const cubic3: GroupRecord = {
  sym: "hand cubic three-arm star", cls: "", d: 3,
  q: ["a,0,0", "0,a,0", "0,0,a"], params: ["a"], cen: [],
  gen: ["-y,x,z,-u,t,v", "y,z,x,u,v,t", "-x,y,z,-t,u,v"],
};

// Independent exhaustive method: enumerate every set of up to three generators,
// rather than grow the lattice by joining cyclic subgroups as production does.
// Crystallographic point groups in 3D need at most three generators.
const bruteCountCache = new Map<string, number>();
function brutePointSubgroupCount(rec: GroupRecord): number {
  const matrices = buildOps(rec).map(op => op.W.slice(0, 3).map(row => row.slice(0, 3)));
  const unique = new Map(matrices.map(R => [JSON.stringify(R), R]));
  const cacheKey = [...unique.keys()].sort().join("|");
  const cached = bruteCountCache.get(cacheKey);
  if (cached !== undefined) return cached;
  const K = [...unique.values()];
  const pos = new Map(K.map((R, i) => [JSON.stringify(R), i]));
  const mul = (A: number[][], B: number[][]) => A.map(row => [0, 1, 2].map(j => row.reduce((sum, x, k) => sum + x * B[k][j], 0)));
  const table = K.map(A => K.map(B => pos.get(JSON.stringify(mul(A, B)))!));
  const id = pos.get(JSON.stringify([[1,0,0], [0,1,0], [0,0,1]]))!;
  const groups = new Set<bigint>();
  for (let a = 0; a < K.length; a++) for (let b = a; b < K.length; b++) for (let c = b; c < K.length; c++) {
    const generators = [a, b, c];
    const queue = [id];
    let mask = 1n << BigInt(id);
    for (let cursor = 0; cursor < queue.length; cursor++) for (const g of generators) {
      const x = table[queue[cursor]][g];
      const bit = 1n << BigInt(x);
      if (!(mask & bit)) { mask |= bit; queue.push(x); }
    }
    groups.add(mask);
  }
  bruteCountCache.set(cacheKey, groups.size);
  return groups.size;
}

describe("tSubgroups matrix internal actions", () => {
  it("hand P6/mmm with parallel q has 32 conjugacy classes and 54 subgroups", () => {
    const classes = checkClosedAndIndex(P6mmm2);
    // D6h ≅ S3 × C2². By projection onto S3, subgroup conjugacy-class
    // counts for 1, C2, C3, S3 are 5,11,5,11; the C2 case has 3 conjugates.
    expect(classes).toHaveLength(5 + 11 + 5 + 11);
    expect(classes.reduce((n, c) => n + c.classSize, 0)).toBe(5 + 3 * 11 + 5 + 11);
    const pg = pointGroup(buildOps(P6mmm2), 2);
    expect(pg).toHaveLength(24);
    expect(pg[0].E).toEqual([[1,0], [0,1]]);
    expect(pg.some(p => JSON.stringify(p.E) === "[[-1,0],[0,-1]]")).toBe(true);
    checkMaximal(P6mmm2);
  });
  it("retains full off-diagonal E matrices for a cubic three-arm star", () => {
    const pg = pointGroup(buildOps(cubic3), 3);
    expect(pg).toHaveLength(48);
    expect(pg[0].E).toEqual([[1,0,0], [0,1,0], [0,0,1]]);
    for (const p of pg) expect(p.E).toEqual(p.R);
    const classes = checkClosedAndIndex(cubic3);
    expect(classes.reduce((n, c) => n + c.classSize, 0)).toBe(brutePointSubgroupCount(cubic3));
  });
});

interface EnumerationFixture {
  label: string;
  classes: { index: number; class_size: number }[];
}

function dimensionRecords(d: number): GroupRecord[] {
  const out: GroupRecord[] = [];
  for (const f of readdirSync(DATA).filter(f => /^sg\d+\.json$/.test(f)).sort()) {
    const groups = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups as Record<string, GroupRecord>;
    for (const [label, r] of Object.entries(groups)) if (r.d === d) out.push({ ...r, label });
  }
  return out;
}

const SUBGROUP_ROOT = join(LOCAL, "fixtures", "subgroups");
describe.skipIf(!existsSync(DATA))("tSubgroups vs dimension fixtures", () => {
  for (const d of [2, 3]) it(`d = ${d}: fixture counts, indices, class sizes and cubic subgroup totals`, async () => {
    const records = dimensionRecords(d);
    const dir = join(SUBGROUP_ROOT, `d${d}`);
    expect(existsSync(dir), `missing d${d} fixtures; generate Task 5 fixtures first`).toBe(true);
    const labels = readdirSync(dir).filter(f => /^\d+\..*\.json$/.test(f)).map(f => f.slice(0, -5));
    const byLabel = new Map(records.map(r => [r.label!, r]));
    const selected = labels.map(label => {
      const rec = byLabel.get(label);
      expect(rec, `unknown fixture ${label}`).toBeDefined();
      return rec!;
    });
    expect(selected.length).toBeGreaterThan(0);
    if (d === 2) {
      expect(records).toHaveLength(3338);
      expect([...labels].sort()).toEqual(records.map(r => r.label!).sort());
    } else {
      expect(new Set(selected.map(r => r.cls))).toEqual(new Set(records.map(r => r.cls)));
      expect(new Set(selected.map(r => r.label!.split(".")[0]))).toEqual(new Set(records.map(r => r.label!.split(".")[0])));
    }
    const ms = (a: number[]) => [...a].sort((x, y) => x - y);
    for (const rec of selected) {
      const fixture: EnumerationFixture = JSON.parse(readFileSync(join(dir, `${rec.label}.json`), "utf8"));
      expect(fixture.label).toBe(rec.label);
      const classes = tSubgroups(rec);
      expect(classes, rec.label).toHaveLength(fixture.classes.length);
      expect(ms(classes.map(c => c.index)), rec.label).toEqual(ms(fixture.classes.map(c => c.index)));
      expect(ms(classes.map(c => c.classSize)), rec.label).toEqual(ms(fixture.classes.map(c => c.class_size)));
      if (Number(rec.label!.split(".")[0]) >= 195) {
        expect(classes.reduce((n, c) => n + c.classSize, 0), rec.label).toBe(brutePointSubgroupCount(rec));
      }
      // Keep Phase A's event-loop fairness during large local-corpus validation.
      await new Promise<void>(resolve => setImmediate(resolve));
    }
  }, 180000);
});
