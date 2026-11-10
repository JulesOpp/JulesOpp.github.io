import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { loadMoyo } from "../src/moyo";
import { runFinder, FinderError } from "../src/finder";
import { parseByuPage } from "../src/paste";
import type { GroupRecord } from "../src/ssg";

const CACHE = join(__dirname, "..", "local", "cache");
const page = (label: string) => parseByuPage(readFileSync(join(CACHE, label + ".html"), "utf8"));
const haveCache = existsSync(join(CACHE, "62.1.9.2.html"));

beforeAll(async () => { await loadMoyo(); });

// A hand-written (3+1)D group (no BYU data): Pnma(00g)0s0-like generators.
const PNMA: GroupRecord = {
  label: "62.1.9.2", sym: "test", cls: "1.9", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
  gen: ["-x+1/2,-y,z+1/2,t", "-x,y+1/2,-z,-t+1/2", "-x,-y,-z,-t"],
};

describe("runFinder", () => {
  it("validates parameters", () => {
    expect(() => runFinder(PNMA, { g: "" })).toThrow(FinderError);
    expect(() => runFinder(PNMA, { g: "one third" })).toThrow(/not a rational/);
  });
  it("refuses huge supercells", () => {
    expect(() => runFinder(PNMA, { g: "1/3000001" })).toThrow(/stops at/);
  });
  it("returns rows ordered general → points with identified groups", () => {
    const r = runFinder(PNMA, { g: "1/3" });
    expect(r.rows[0].kind).toBe("general");
    expect(r.classes.every((c) => c.id.ok)).toBe(true);
    expect(r.supercell.conv).toEqual([1, 1, 3]);
  });
  it.skipIf(!haveCache)("IT finder: Pnma(00g)0s0, g = 1/3 → P21, P21/a, P212121", () => {
    const r = runFinder(page("62.1.9.2"), { g: "1/3" });
    expect(new Set(r.classes.map((c) => c.itShort))).toEqual(new Set(["P21", "P21/a", "P212121"]));
    expect(r.warnings).toEqual([]);
  });
  it.skipIf(!haveCache)("In8IrSe9: 55.3.122.96 with a1=b1=g2=1/4", () => {
    const r = runFinder(page("55.3.122.96"), { a1: "1/4", b1: "1/4", g2: "1/4" });
    expect(r.supercell.index).toBe(32);
    const nums = new Set(r.classes.map((c) => (c.id.ok ? c.id.number : -1)));
    // P1, Pc, P2, P2/c (README §4, A-centred a-unique cell reduced to primitive)
    expect(nums).toEqual(new Set([1, 7, 3, 13]));
    expect(r.rows.some((row) => row.locus.includes("t3 ≡ 1/8 (mod 1/4)"))).toBe(true);
  });
});
