import opencascade from "replicad-opencascadejs";
import { setOC } from "replicad";
import { createRequire } from "node:module";

let ready: Promise<void> | null = null;
/** Load OpenCascade once per test file. */
export function loadOC(): Promise<void> {
  if (!ready) {
    const require = createRequire(import.meta.url);
    const wasm = require.resolve("replicad-opencascadejs/wasm");
    ready = (opencascade as unknown as (o: object) => Promise<unknown>)({ locateFile: () => wasm }).then((oc) => setOC(oc as never));
  }
  return ready;
}
