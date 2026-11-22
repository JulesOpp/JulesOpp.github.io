// Values in, finished mesh out: body -> texture -> screws -> rotation.

import { arena, setSegments, toTriMesh, union, type Manifold, type TriMesh } from "./kernel";
import type { HoldDef, Values } from "./spec";
import { applyTexture, type TexMode } from "./texture";

export interface BuildOptions {
  screws: boolean;
  texMode: TexMode;
  texHeight: number;
  texSize: number;
  density: number;
  rotation: number; // degrees about the wall normal
  segments: number; // circular resolution
}

export const DEFAULT_OPTIONS: BuildOptions = {
  screws: true, texMode: "none", texHeight: 0.8, texSize: 6, density: 1, rotation: 0, segments: 48,
};

export interface BuildResult {
  mesh: TriMesh;
  volume: number; // mm^3
  bbox: [number, number, number]; // off-wall x along-wall x up, mm
  triangles: number;
  genus: number;
}

/** Turn a finished Manifold into the result record (inside an arena). */
export function finish(m: Manifold): BuildResult {
  const b = m.boundingBox();
  return {
    mesh: toTriMesh(m),
    volume: m.volume(),
    bbox: [b.max[0] - b.min[0], b.max[1] - b.min[1], b.max[2] - b.min[2]],
    triangles: m.numTri(),
    genus: m.genus(),
  };
}

export function buildHold(def: HoldDef, values: Values, opt: BuildOptions = DEFAULT_OPTIONS): BuildResult {
  setSegments(opt.segments);
  return arena(() => {
    const { body, screws } = def.build(values, opt.screws);
    let m = applyTexture(body, opt.texMode, opt.texHeight, opt.texSize, opt.density);
    if (screws.length) m = m.subtract(union(screws));
    if (opt.rotation) m = m.rotate([opt.rotation, 0, 0]);
    if (m.isEmpty()) throw new Error("these settings leave no solid: try other values");
    return finish(m);
  });
}
