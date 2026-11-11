import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { loadMoyo } from "../src/moyo";
import { runSubgroups, SubgroupError } from "../src/subfinder";

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
