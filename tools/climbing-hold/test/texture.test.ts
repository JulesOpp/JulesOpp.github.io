import { beforeAll, expect, it } from "vitest";
import { initKernel } from "../src/kernel";
import { HOLDS } from "../src/holds";
import { buildHold, DEFAULT_OPTIONS } from "../src/pipeline";
import { defaults } from "../src/spec";
import { TEX_MODES } from "../src/texture";

beforeAll(async () => { await initKernel(); });

const pinch = HOLDS.find((h) => h.id === "pinch")!;
const plain = () => buildHold(pinch, defaults(pinch), { ...DEFAULT_OPTIONS, segments: 32 });

it.each(TEX_MODES.filter((m) => m !== "none"))("%s keeps a closed solid with a flat back and its screw holes", (mode) => {
  const base = plain();
  const r = buildHold(pinch, defaults(pinch), { ...DEFAULT_OPTIONS, segments: 32, texMode: mode, texHeight: 1.5, texSize: 6 });
  expect(r.genus).toBe(base.genus); // the two through-holes survive
  expect(Math.abs(r.volume - base.volume) / base.volume).toBeLessThan(0.05);
  expect(r.triangles).toBeGreaterThan(base.triangles);
  let minX = Infinity;
  for (let i = 0; i < r.mesh.pos.length; i += 3) minX = Math.min(minX, r.mesh.pos[i]);
  expect(minX).toBeGreaterThan(-1e-3);
}, 60_000);
