import { beforeAll, describe, expect, it, vi } from "vitest";
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { loadMoyo } from "../src/moyo";
import * as subfinder from "../src/subfinder";
import * as standardiser from "../src/standardise";
import { keep, runSubgroups, SubgroupError, type SubRow } from "../src/subfinder";

beforeAll(async () => {
  await loadMoyo();
});

// Hand-written BYU-style page text (no BYU data): P-3c1(0,0,g)000, IT A No. 165.
const page = (q: string[], ops: string[]) =>
  `Superspace group: 165.1.24.1 P-3c1(0,0,g)000 Independent modulation vectors: ` +
  q.map((v, i) => `q${i + 1} = (${v})`).join(" ") +
  ` Centering: (${"0,".repeat(q.length + 2)}0) Non-lattice generators: Non-lattice operators: ` +
  ops.map((o) => `(${o})`).join(" ");
const P3C1 = page(["0,0,g"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);

describe("runSubgroups (P-3c1(0,0,g)000)", () => {
  let all: ReturnType<typeof runSubgroups>;
  beforeAll(() => {
    all = runSubgroups(P3C1, { kIndex: 1 });
  });
  it("lists 10 classes", () => {
    expect(all.rows).toHaveLength(10);
    expect(all.rows.map((r) => r.sgNumber).sort((a, b) => a - b)).toEqual([1, 2, 5, 9, 15, 143, 147, 150, 158, 165]);
  });
  it("row for No. 150 is P321(0,0,g)000 with index 2", () => {
    const r = all.rows.find((x) => x.sgNumber === 150)!;
    expect(r.type!.symbol.full).toBe("P321(0,0,g)000");
    expect(r.index).toBe(2);
    expect(r.tIndex).toBe(2);
    expect(r.kIndex).toBe(1);
    expect(r.pointGroup).toBe("321");
    expect(r.inParentAxes).not.toBeNull();
    expect(r.generators.length).toBeGreaterThan(0);
  });
  it("system filter trigonal gives 5 rows", () => {
    expect(runSubgroups(P3C1, { kIndex: 1, system: "trigonal" }).rows).toHaveLength(5);
  });
  it("pointGroup filter 321 gives 1 row", () => {
    expect(runSubgroups(P3C1, { kIndex: 1, pointGroup: "321" }).rows).toHaveLength(1);
  });
  it("centring filter", () => {
    const rows = runSubgroups(P3C1, { kIndex: 1, centring: "B" }).rows;
    expect(rows).toHaveLength(3);
    expect(rows.every((r) => r.type!.symbol.basic[0] === "B")).toBe(true);
  });
  it("monoclinic subgroup drops out of the parent axes", () => {
    expect(all.rows.find((x) => x.sgNumber === 15)!.inParentAxes).toBeNull();
  });
  it("refuses k-index != 1 and d != 1", () => {
    expect(() => runSubgroups(P3C1, { kIndex: 2 })).toThrow(SubgroupError);
    expect(() => runSubgroups(P3C1, { kIndex: 2 })).toThrow("only k-index 1 is supported for now");
    const d2 = page(["0,0,g", "0,0,h"], ["x,y,z,t,u", "y,-x+y,-z,-t,-u"]);
    expect(() => runSubgroups(d2, { kIndex: 1 })).toThrow("only d = 1 is supported for now");
  });
});

it("caches the parent once for repeated canonical operators and q text", () => {
  const raw = page(["0,0,g+99"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);
  const spy = vi.spyOn(standardiser, "standardise");
  const hits = vi.spyOn(subfinder, "_cacheHits");
  const before = subfinder._cacheHits();
  const a = runSubgroups(raw, { kIndex: 1 });
  const first = spy.mock.calls.length;
  runSubgroups(raw, { kIndex: 1 });
  expect(first).toBe(a.rows.length + 1);
  expect(spy.mock.calls.length - first).toBe(a.rows.length);
  expect(subfinder._cacheHits()).toBe(before + 1);
  expect(hits).toHaveBeenCalledTimes(2);
  spy.mockRestore(); hits.mockRestore();
});

it("parent cache uses canonical operations and q text, and evicts beyond eight entries", () => {
  const original = page(["0,0,g+101"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);
  const reordered = page(["0,0,g+101"], ["x-y,-y,-z+1/2,-t", "y,-x+y,-z,-t", "x,y,z,t"]);
  runSubgroups(original, { kIndex: 1 });
  const before = subfinder._cacheHits();
  runSubgroups(reordered, { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(before + 1);
  // Identity-only groups keep this capacity check cheap; the q text must be part
  // of the key even though every one has the same basic space-group number.
  const raws = Array.from({ length: 9 }, (_, i) => page([`0,0,g+${i+201}`], ["x,y,z,t"]));
  const start = subfinder._cacheHits();
  for (const raw of raws) runSubgroups(raw, { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start);
  runSubgroups(raws[8], { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start + 1);
  runSubgroups(raws[0], { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start + 1);
});

describe("keep", () => {
  it("rows with an unidentified type pass every filter", () => {
    const r: SubRow = {
      index: 2, tIndex: 2, kIndex: 1, classSize: 1, maximal: true, pointGroup: "?", sgNumber: 0, sgSymbol: "?",
      inParentAxes: null, type: null, generators: [],
    };
    expect(keep(r, { kIndex: 1, pointGroup: "321", system: "cubic", centring: "F" })).toBe(true);
  });
});

// Local-only: 98.1.21.1 I4_122(0,0,g)000 has two classes with the same BYU type B2(a,b,0)0
// (checked in local/fixtures/subgroups); they must stay two rows with equal type symbols.
const REC = join(__dirname, "..", "local", "data", "sg098.json");
describe.skipIf(!existsSync(REC))("duplicate types (local)", () => {
  it("98.1.21.1 keeps two rows with equal type.symbol.full", () => {
    const g = JSON.parse(readFileSync(REC, "utf8")).groups["98.1.21.1"];
    const cen = g.cen.map((c: string) => `(${c})`).join(" ");
    const raw =
      `Superspace group: 98.1.21.1 ${g.sym} Independent modulation vectors: q1 = (${g.q[0]}) ` +
      `Centering: (0,0,0,0) ${cen} Non-lattice generators: Non-lattice operators: ` +
      ["x,y,z,t", ...g.gen].map((o: string) => `(${o})`).join(" ");
    const rows = runSubgroups(raw, { kIndex: 1 }).rows;
    const by = new Map<string, number>();
    for (const r of rows) if (r.type) by.set(r.type.symbol.full, (by.get(r.type.symbol.full) ?? 0) + 1);
    expect(Math.max(...by.values())).toBeGreaterThanOrEqual(2);
  });
});
