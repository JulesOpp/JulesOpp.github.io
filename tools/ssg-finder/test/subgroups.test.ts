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
    expect(cl[0].pointOps[0].eps).toBe(1);
  });
  it("closed under opMul, order index times smaller", () => {
    checkClosedAndIndex(P3c1);
  });
  it("pointGroup returns identity first with eps = +-1", () => {
    const pg = pointGroup(buildOps(P3c1), 1);
    expect(pg).toHaveLength(12);
    expect(pg[0].R).toEqual([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    expect(pg.every((p) => Math.abs(p.eps) === 1)).toBe(true);
  });
  it("marks maximal subgroups", () => {
    const cl = tSubgroups(P3c1);
    const byOrder = (n: number) => cl.filter((c) => c.pointOps.length === n);
    // -3m1 (12) is the parent; 3m1, 321, -3 (order 6) are maximal; trivial is not.
    expect(byOrder(6).every((c) => c.maximal)).toBe(true);
    expect(byOrder(1)[0].maximal).toBe(false);
  });
  it("refuses d != 1", () => {
    expect(() => tSubgroups({ ...P3c1, d: 2 })).toThrow(/only d = 1/);
  });
});

const LOCAL = join(__dirname, "..", "local");
const FIX = join(LOCAL, "fixtures", "subgroups");
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
