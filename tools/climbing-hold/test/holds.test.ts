import { beforeAll, describe, expect, it } from "vitest";
import { initKernel } from "../src/kernel";
import { HOLDS } from "../src/holds";
import { buildHold, DEFAULT_OPTIONS } from "../src/pipeline";
import { defaults, type Values } from "../src/spec";

beforeAll(async () => { await initKernel(); });

const corners = (h: (typeof HOLDS)[number]): [string, Values][] => [
  ["default", defaults(h)],
  ["min", Object.fromEntries(h.params.map((p) => [p.key, p.lo]))],
  ["max", Object.fromEntries(h.params.map((p) => [p.key, p.hi]))],
];

/** Every vertex at x <= 0 must be on the wall plane, and some must be. */
function checkFlatBack(pos: Float32Array) {
  let minX = Infinity, onWall = 0;
  for (let i = 0; i < pos.length; i += 3) { minX = Math.min(minX, pos[i]); if (Math.abs(pos[i]) < 1e-3) onWall++; }
  expect(minX).toBeGreaterThan(-1e-3);
  expect(onWall).toBeGreaterThan(2);
}

describe.each(HOLDS.map((h) => [h.id, h] as const))("%s", (_id, h) => {
  it.each(corners(h))("%s settings build a closed solid with a flat back", (_name, v) => {
    const r = buildHold(h, v, { ...DEFAULT_OPTIONS, segments: 24 });
    expect(r.volume).toBeGreaterThan(100);
    expect(r.triangles).toBeGreaterThan(20);
    checkFlatBack(r.mesh.pos);
  }, 60_000);
});
