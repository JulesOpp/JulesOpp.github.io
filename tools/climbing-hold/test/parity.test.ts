// The ported holds against the desktop app's own OpenSCAD renders
// (scripts/parity_reference.py -> parity-reference.json), defaults, no screws.
import { beforeAll, describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { initKernel } from "../src/kernel";
import { HOLDS } from "../src/holds";
import { buildHold, DEFAULT_OPTIONS } from "../src/pipeline";
import { defaults } from "../src/spec";

const ref: Record<string, { volume: number; extents: number[] }> = JSON.parse(
  readFileSync(new URL("./parity-reference.json", import.meta.url), "utf8"));

beforeAll(async () => { await initKernel(); });

describe("matches the desktop app", () => {
  it.each(Object.keys(ref))("%s", (id) => {
    const h = HOLDS.find((x) => x.id === id)!;
    const r = buildHold(h, defaults(h), { ...DEFAULT_OPTIONS, screws: false, segments: 48 });
    const want = ref[id];
    // both are 48-segment tessellations of the same curved solids, and the
    // superellipse/loft templates are re-implementations: allow 4 %
    expect(Math.abs(r.volume - want.volume) / want.volume).toBeLessThan(0.04);
    r.bbox.forEach((e, i) => expect(Math.abs(e - want.extents[i])).toBeLessThan(1.5));
  });
});
