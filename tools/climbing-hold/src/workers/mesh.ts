// Module worker: builds holds off the main thread so the sliders stay smooth.

import { arena, fromTriMesh, initKernel, setSegments } from "../kernel";
import { HOLDS } from "../holds";
import { buildHold, finish, type BuildOptions, type BuildResult } from "../pipeline";
import { overhangFraction } from "../printcheck";
import type { Values } from "../spec";
import { fromBinarySTL } from "../stl";
import { applyTexture } from "../texture";
import type { TriMesh } from "../kernel";

export type MeshRequest =
  | { kind: "build"; id: number; hold: string; values: Values; opt: BuildOptions }
  | { kind: "scan"; id: number; url: string; opt: BuildOptions };

export type MeshReply =
  | { id: number; ok: true; result: BuildResult; overhang: number }
  | { id: number; ok: false; error: string };

const ready = initKernel((f) => new URL(f, import.meta.url).href);
const scans = new Map<string, Promise<TriMesh>>();

function scanMesh(url: string): Promise<TriMesh> {
  let p = scans.get(url);
  if (!p) {
    p = fetch(url).then((r) => {
      if (!r.ok) throw new Error(`could not load ${url} (${r.status})`);
      return r.arrayBuffer();
    }).then(fromBinarySTL);
    scans.set(url, p);
  }
  return p;
}

self.onmessage = async (ev: MessageEvent<MeshRequest>) => {
  const req = ev.data;
  try {
    await ready;
    let result: BuildResult;
    if (req.kind === "build") {
      const def = HOLDS.find((h) => h.id === req.hold);
      if (!def) throw new Error(`unknown hold ${req.hold}`);
      result = buildHold(def, req.values, req.opt);
    } else {
      const tm = await scanMesh(req.url);
      setSegments(req.opt.segments);
      result = arena(() => {
        let m = applyTexture(fromTriMesh({ pos: tm.pos.slice(), idx: tm.idx.slice() }), req.opt.texMode, req.opt.texHeight, req.opt.texSize, req.opt.density);
        if (req.opt.rotation) m = m.rotate([req.opt.rotation, 0, 0]);
        return finish(m);
      });
    }
    const reply: MeshReply = { id: req.id, ok: true, result, overhang: overhangFraction(result.mesh) };
    (self as unknown as Worker).postMessage(reply, [result.mesh.pos.buffer, result.mesh.idx.buffer]);
  } catch (err) {
    const reply: MeshReply = { id: req.id, ok: false, error: (err as Error).message ?? String(err) };
    (self as unknown as Worker).postMessage(reply);
  }
};
