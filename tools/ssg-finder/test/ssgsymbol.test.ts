import { describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { buildOps, GroupRecord } from "../src/ssg";
import { crystalSystemOf } from "../src/symbol";
import { ssgSymbolInAxes } from "../src/ssgsymbol";

// P-3c1(0,0,g)000, IT A No. 165 generators in hexagonal axes with the 4th coordinate.
const P3c1: GroupRecord = {
  sym: "P-3c1(0,0,g)000", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
  gen: ["y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"],
};

describe("ssgSymbolInAxes (public)", () => {
  it("hand-written P-3c1(0,0,g)000", () => {
    const s = ssgSymbolInAxes(buildOps(P3c1), P3c1.q, "trigonal");
    expect(s).toEqual({ basic: "P-3c1", q: "(0,0,g)", shifts: "000", full: "P-3c1(0,0,g)000" });
  });
  it("phase letters: P6_2 with t' = t + 1/6 gives h", () => {
    const r: GroupRecord = { sym: "", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [], gen: ["x-y,x,z+1/3,t+1/6"] };
    expect(ssgSymbolInAxes(buildOps(r), r.q, "hexagonal")?.full).toBe("P6_2(0,0,g)h");
  });
  it("refuses d != 1", () => {
    const ops = buildOps({ ...P3c1, d: 2, gen: ["y,-x+y,-z,-t,-u"] });
    expect(() => ssgSymbolInAxes(ops, ["0,0,g", "0,0,g"], "trigonal")).toThrow(/d = 1/);
  });
});

// The oracle: every (3+1)D BYU record's own symbol. BYU-derived, so local-only.
const DATA = join(__dirname, "..", "local", "data");
const have = existsSync(DATA);

function allD1(): (GroupRecord & { label: string })[] {
  const out: (GroupRecord & { label: string })[] = [];
  for (const f of readdirSync(DATA).filter((f) => /^sg\d{3}\.json$/.test(f)).sort()) {
    const groups: Record<string, GroupRecord> = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups;
    for (const [label, r] of Object.entries(groups)) if (r.d === 1) out.push({ ...r, label });
  }
  return out;
}

describe.skipIf(!have)("ssgSymbolInAxes vs all (3+1)D BYU symbols", () => {
  it("reproduces every record's own symbol", () => {
    const recs = allD1();
    expect(recs).toHaveLength(775);
    const bad: string[] = [];
    for (const r of recs) {
      const sg = Number(r.label.split(".")[0]);
      let got: string;
      try {
        got = ssgSymbolInAxes(buildOps(r), r.q, crystalSystemOf(sg))?.full ?? "null";
      } catch (e) {
        got = `throw ${(e as Error).message}`;
      }
      if (got !== r.sym) bad.push(`${sg}\t${r.label}\t${r.sym}\tgot ${got}`);
    }
    if (bad.length) {
      const bySg = new Map<string, string[]>();
      for (const b of bad) {
        const k = b.split("\t")[0];
        if (!bySg.has(k)) bySg.set(k, []);
        bySg.get(k)!.push(b);
      }
      console.log(`${bad.length} mismatches of ${recs.length} in ${bySg.size} space groups; first 20:\n` + bad.slice(0, 20).join("\n"));
    }
    expect(bad.length).toBe(0);
  });
});
