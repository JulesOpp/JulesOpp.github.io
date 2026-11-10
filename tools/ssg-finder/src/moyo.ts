// Loader for @spglib/moyo-wasm, usable both in the browser and in Node (tests).
//
// Browser: `await loadMoyo(urlOfWasm)` fetches and instantiates the wasm (the default
// wasm-bindgen `init`). Node: `await loadMoyo()` reads the wasm next to the package's JS
// with `initSync({module: bytes})`. Either way the module is initialised once and cached;
// `getMoyo()` then returns it synchronously (identify() needs a synchronous handle).

import init, * as moyoNs from "@spglib/moyo-wasm";

export type Moyo = typeof moyoNs;

let ready: Promise<Moyo> | null = null;
let loaded: Moyo | null = null;

function isNode(): boolean {
  return typeof process !== "undefined" && !!process.versions?.node && typeof window === "undefined";
}

async function loadNode(wasmPath?: string | URL): Promise<void> {
  // Specifiers in variables so a browser bundler does not try to resolve them.
  const fsName = "node:fs";
  const modName = "node:module";
  const fs: typeof import("node:fs") = await import(/* @vite-ignore */ fsName);
  let path = wasmPath;
  if (path === undefined) {
    const { createRequire }: typeof import("node:module") = await import(/* @vite-ignore */ modName);
    const require = createRequire(import.meta.url);
    path = require.resolve("@spglib/moyo-wasm/moyo_wasm_bg.wasm");
  }
  moyoNs.initSync({ module: fs.readFileSync(path) });
}

/**
 * Initialise moyo-wasm (idempotent). `wasmUrl` is required in the browser (the URL of
 * moyo_wasm_bg.wasm as deployed); in Node it may be omitted (resolved from node_modules)
 * or be a file path / file: URL.
 */
export async function loadMoyo(wasmUrl?: string | URL): Promise<Moyo> {
  if (!ready) {
    ready = (async () => {
      const fileLike = wasmUrl instanceof URL ? wasmUrl.protocol === "file:" : wasmUrl !== undefined && !/^[a-z]+:\/\//i.test(wasmUrl);
      if (isNode() && (wasmUrl === undefined || fileLike)) {
        await loadNode(wasmUrl);
      } else {
        await init({ module_or_path: wasmUrl ?? new URL("moyo_wasm_bg.wasm", import.meta.url) });
      }
      loaded = moyoNs;
      return moyoNs;
    })();
    ready.catch(() => {
      ready = null;
    });
  }
  return ready;
}

/** The initialised module; throws if loadMoyo() has not completed. */
export function getMoyo(): Moyo {
  if (!loaded) throw new Error("moyo-wasm not loaded: await loadMoyo() first");
  return loaded;
}

/** moyo returns 3x3 matrices as flat column-major arrays; convert to row-major number[][]. */
export function fromColMajor(a: ArrayLike<number>): number[][] {
  return [0, 1, 2].map((i) => [0, 1, 2].map((j) => a[3 * j + i]));
}

/** Row-major number[][] to moyo's flat column-major layout. */
export function toColMajor(m: number[][]): [number, number, number, number, number, number, number, number, number] {
  return [m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]];
}
