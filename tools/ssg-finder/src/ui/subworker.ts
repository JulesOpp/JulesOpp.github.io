// Module worker: computes only the viewer's pasted data, never fetches group tables.
import { loadMoyo } from "../moyo";
import { runSubgroups, type SubFilter, type SubResult } from "../subfinder";

export type SubRequest = { type: "run"; raw: string; filter: SubFilter };
export type SubResponse =
  | { type: "progress"; done: number; total: number }
  | { type: "rows"; result: SubResult; complete: boolean }
  | { type: "error"; message: string };

const ready = loadMoyo(new URL("./moyo_wasm_bg.wasm", import.meta.url));
ready.catch(() => {});
const send = (message: SubResponse) => self.postMessage(message);
self.onmessage = async (event: MessageEvent<SubRequest>) => {
  if (event.data.type !== "run") return;
  try {
    await ready;
    const result = runSubgroups(event.data.raw, event.data.filter, (done, total, partial) => {
      send({ type: "progress", done, total });
      send({ type: "rows", result: partial, complete: false });
    });
    send({ type: "rows", result, complete: true });
  } catch (error) {
    send({ type: "error", message: error instanceof Error ? error.message : String(error) });
  }
};
