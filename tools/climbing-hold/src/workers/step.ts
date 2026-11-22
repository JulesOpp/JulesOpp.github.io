// Module worker: OpenCascade (replicad) for the machinable edge's STEP.
// Only loaded when someone asks for a STEP file -- the wasm is ~23 MB.

import opencascade from "replicad-opencascadejs";
import { setOC } from "replicad";
import { machinableEdgeSolid } from "../brep/machinable-edge";
import type { EdgeParams } from "../edge/params";

export interface StepRequest { id: number; params: EdgeParams }
export type StepReply =
  | { id: number; ok: true; step: ArrayBuffer; volume: number; area: number; faces: number; flatArea: number; bbox: [number, number, number]; mesh: { pos: Float32Array; idx: Uint32Array } }
  | { id: number; ok: false; error: string };

const ready = (opencascade as unknown as (o: object) => Promise<unknown>)({
  locateFile: (f: string) => new URL(f.endsWith(".wasm") ? "replicad_single.wasm" : f, import.meta.url).href,
}).then((oc) => setOC(oc as never));

self.onmessage = async (ev: MessageEvent<StepRequest>) => {
  const { id, params } = ev.data;
  try {
    await ready;
    const r = machinableEdgeSolid(params);
    const step = await r.solid.blobSTEP().arrayBuffer();
    const m = r.solid.mesh({ tolerance: 0.05, angularTolerance: 0.2 });
    const mesh = { pos: new Float32Array(m.vertices), idx: new Uint32Array(m.triangles) };
    r.solid.delete();
    const reply: StepReply = { id, ok: true, step, volume: r.volume, area: r.area, faces: r.faces, flatArea: r.flatArea, bbox: r.bbox, mesh };
    (self as unknown as Worker).postMessage(reply, [step, mesh.pos.buffer, mesh.idx.buffer]);
  } catch (err) {
    (self as unknown as Worker).postMessage({ id, ok: false, error: (err as Error).message ?? String(err) } satisfies StepReply);
  }
};
