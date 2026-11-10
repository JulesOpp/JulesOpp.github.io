// IUCr IT Superspace Group Finder, 165.1.24.1 P-3c1(0,0,g)000, commensurate table
// (https://it.iucr.org/resources/finder/?search=P63, retrieved 2026-11-10).
// First non-orthogonal (hexagonal-axes) check. Needs the gitignored BYU page in local/cache;
// skipped when it is absent.
import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { loadMoyo } from "../src/moyo";
import { runFinder } from "../src/finder";
import { parseByuPage } from "../src/paste";

const FILE = join(__dirname, "..", "local", "cache", "165.1.24.1.html");
const have = existsSync(FILE);

beforeAll(async () => { await loadMoyo(); });

type Col = "general" | "n/N" | "1/4N+n/N" | "1/2N+n/N" | "3/4N+n/N";
// IT columns → [space-group number, IT short symbol in the parental hexagonal axes]; absent = empty.
const D_ODD: Partial<Record<Col, [number, string]>> = {
  general: [143, "P3"], "n/N": [147, "P-3"], "1/4N+n/N": [150, "P321"],
  "1/2N+n/N": [147, "P-3"], "3/4N+n/N": [150, "P321"],
};
const D_EVEN_N_ODD: Partial<Record<Col, [number, string]>> = {
  general: [158, "P3c1"], "n/N": [165, "P-3c1"], "1/2N+n/N": [165, "P-3c1"],
};

const CASES: [string, Partial<Record<Col, [number, string]>>][] = [
  ...["1/3", "1/5", "3/4", "1/2", "1/4"].map((g) => [g, D_ODD] as [string, typeof D_ODD]),
  ...["2/3", "2/5", "4/7"].map((g) => [g, D_EVEN_N_ODD] as [string, typeof D_ODD]),
];

describe.skipIf(!have)("IT finder: P-3c1(0,0,g)000 (165.1.24.1)", () => {
  const rec = have ? parseByuPage(readFileSync(FILE, "utf8")) : null;
  it.each(CASES)("g = %s", (g, expected) => {
    const r = runFinder(rec!, { g });
    expect(r.warnings).toEqual([]);
    const got: Record<string, [number, string]> = {};
    for (const c of r.classes) {
      expect(c.id.ok).toBe(true);
      const v: [number, string] = [c.id.ok ? c.id.number : -1, c.itShort ?? ""];
      const col = c.itLabel!;
      if (col in got) expect(got[col]).toEqual(v); // one entry per column
      got[col] = v;
    }
    expect(got).toEqual(expected);
  });
});
