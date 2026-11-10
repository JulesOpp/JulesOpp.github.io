// The TS port (src/exact.ts + identify) must reproduce the golden fixtures written by
// reference/make_golden.py (CONTRACTS.md §3). BYU-derived, so they live in local/.
import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { Exact } from "../src/exact";
import { formatOp, GroupRecord } from "../src/ssg";
import { loadMoyo } from "../src/moyo";
import { identify } from "../src/identify";

const LOCAL = join(__dirname, "..", "local");
const GOLD = join(LOCAL, "fixtures", "golden");
const DATA = join(LOCAL, "data");
const files = existsSync(GOLD) ? readdirSync(GOLD).filter((f) => f.endsWith(".json") && !f.startsWith("_")) : [];
const have = files.length > 0 && existsSync(DATA);

const recCache = new Map<string, Record<string, GroupRecord>>();
function record(label: string): GroupRecord {
  const sg = label.split(".")[0].padStart(3, "0");
  if (!recCache.has(sg)) recCache.set(sg, JSON.parse(readFileSync(join(DATA, `sg${sg}.json`), "utf8")).groups);
  const r = recCache.get(sg)![label];
  if (!r) throw new Error(`no record for ${label}`);
  return { ...r, label };
}

beforeAll(async () => {
  if (have) await loadMoyo();
});

describe.skipIf(!have)("golden fixtures (reference/exact.py)", () => {
  it.each(files)("%s", (file) => {
    const g = JSON.parse(readFileSync(join(GOLD, file), "utf8"));
    const ex = new Exact(record(g.label), g.params);
    expect(ex.Q.map((r) => r.map(String))).toEqual(g.Q);
    expect(ex.ops.map(formatOp)).toEqual(g.ops);
    expect(ex.N).toEqual(g.supercell.conv);
    expect(ex.index).toBe(g.supercell.index);
    expect(ex.hnfSuper.map((r) => r.map(Number))).toEqual(g.supercell.hnf);
    expect(ex.cells.filter((n) => n.some((x) => x)).map((n) => n.map((x, i) => (x === 0 ? "0" : simplify(x, ex.N[i])))))
      .toEqual(g.supercell.centring);
    expect(ex.B.map((r) => r.map(String))).toEqual(g.Lambda);
    if (g.d === 1) expect(ex.itN).toBe(g.it_N);

    // classes: compare by canonical flat (class order depends on setting strings)
    const byFlat = new Map(ex.classes.map((c) => [JSON.stringify(c.flat.toJSON()), c]));
    expect(ex.classes.length).toBe(g.classes.length);
    const orbitPairs: [number, number][] = [];
    for (const gc of g.classes) {
      const c = byFlat.get(JSON.stringify(gc.flat));
      expect(c, `missing flat ${JSON.stringify(gc.flat)}`).toBeDefined();
      expect(c!.dim).toBe(gc.dim);
      expect(c!.genericT.map(String)).toEqual(gc.generic_t);
      expect(c!.exampleT.map(String)).toEqual(gc.example_t);
      expect(c!.surviving).toEqual(gc.surviving);
      expect(c!.locus).toBe(gc.locus);
      if (g.d === 1) expect(c!.itLabel).toBe(gc.it_label);
      const id = identify(c!.ops3, [[ex.N[0], 0, 0], [0, ex.N[1], 0], [0, 0, ex.N[2]]]);
      expect(id.ok ? id.number : id.reason).toBe(gc.sg_number);
      orbitPairs.push([c!.orbit, gc.orbit]);
    }
    // orbit partitions agree (ids may be numbered differently)
    const fwd = new Map<number, number>(), back = new Map<number, number>();
    for (const [a, b] of orbitPairs) {
      expect(fwd.get(a) ?? b).toBe(b);
      expect(back.get(b) ?? a).toBe(a);
      fwd.set(a, b);
      back.set(b, a);
    }
  });
});

function simplify(n: number, d: number): string {
  const g = (a: number, b: number): number => (b ? g(b, a % b) : a);
  const k = g(n, d);
  return d / k === 1 ? String(n / k) : `${n / k}/${d / k}`;
}
