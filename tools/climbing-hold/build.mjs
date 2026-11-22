// Build the Climbing Hold Generator page scripts.
//   src/ui.ts           -> ../../assets/blog/climbing-hold/chg.js     (page; ES module)
//   src/workers/mesh.ts -> ../../assets/blog/climbing-hold/chgw.js    (module worker: manifold-3d)
//   src/workers/step.ts -> ../../assets/blog/climbing-hold/chstepw.js (module worker: OpenCascade, lazy)
//   manifold-3d/manifold.wasm, replicad-opencascadejs/replicad_single.wasm copied alongside,
//   with each dependency's licence text.
import { build } from "esbuild";
import { copyFileSync, mkdirSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const out = join(here, "..", "..", "assets", "blog", "climbing-hold");
mkdirSync(out, { recursive: true });
const pkgDir = (name) => join(here, "node_modules", name);
const version = (name) => JSON.parse(readFileSync(join(pkgDir(name), "package.json"), "utf8")).version;

// Licence texts are taken from the installed packages; pin the versions
// so an upgrade forces a look at whether the terms changed.
const PINNED = { "manifold-3d": "3.5.5", replicad: "1.2.1", "replicad-opencascadejs": "1.2.1" };
for (const [name, v] of Object.entries(PINNED))
  if (version(name) !== v) throw new Error(`${name} is ${version(name)}, expected ${v}: re-check its licence and update PINNED`);

const common = {
  bundle: true,
  format: "esm",
  platform: "browser",
  target: ["es2020"],
  minify: true,
  sourcemap: true,
  legalComments: "eof",
  logLevel: "warning",
  // emscripten glue mentions node built-ins behind runtime checks
  external: ["node:*", "fs", "path", "module", "url", "crypto", "worker_threads"],
};
const note = "Licences of bundled libraries: see LICENSES.txt next to this file.";
await build({ ...common, entryPoints: [join(here, "src", "ui.ts")], outfile: join(out, "chg.js"), banner: { js: `/* Climbing Hold Generator. ${note} */` } });
await build({ ...common, entryPoints: [join(here, "src", "workers", "mesh.ts")], outfile: join(out, "chgw.js"), banner: { js: `/* Climbing Hold Generator mesh worker (manifold-3d, Apache-2.0). ${note} */` } });
await build({ ...common, entryPoints: [join(here, "src", "workers", "step.ts")], outfile: join(out, "chstepw.js"), banner: { js: `/* Climbing Hold Generator STEP worker (replicad MIT; OpenCascade LGPL-2.1). ${note} */` } });

copyFileSync(join(pkgDir("manifold-3d"), "manifold.wasm"), join(out, "manifold.wasm"));
copyFileSync(join(pkgDir("replicad-opencascadejs"), "dist", "replicad_single.wasm"), join(out, "replicad_single.wasm"));

const lic = (name) => readFileSync(join(pkgDir(name), "LICENSE"), "utf8");
writeFileSync(join(out, "LICENSES.txt"), [
  "Third-party software in this folder",
  "===================================",
  "",
  `manifold-3d ${PINNED["manifold-3d"]} (https://github.com/elalish/manifold): manifold.wasm and the`,
  "bindings bundled in chgw.js. Apache License 2.0:",
  "",
  lic("manifold-3d"),
  "",
  `replicad ${PINNED.replicad} (https://github.com/sgenoud/replicad): bundled in chstepw.js. MIT:`,
  "",
  lic("replicad"),
  "",
  `replicad-opencascadejs ${PINNED["replicad-opencascadejs"]}: replicad_single.wasm and its bindings in chstepw.js.`,
  "A WebAssembly build of Open CASCADE Technology (https://dev.opencascade.org) made with",
  "opencascade.js (https://github.com/donalffons/opencascade.js), distributed unmodified from",
  "the npm package. Source code: https://github.com/sgenoud/replicad (packages/replicad-opencascadejs),",
  "https://github.com/donalffons/opencascade.js and https://github.com/Open-Cascade-SAS/OCCT.",
  "GNU Lesser General Public License 2.1 (with the Open CASCADE exception):",
  "",
  lic("replicad-opencascadejs"),
].join("\n"));

const kb = (f) => (statSync(join(out, f)).size / 1024).toFixed(0) + " KB";
console.log(["chg.js", "chgw.js", "chstepw.js", "manifold.wasm", "replicad_single.wasm"].map((f) => `${f} ${kb(f)}`).join(", "));
