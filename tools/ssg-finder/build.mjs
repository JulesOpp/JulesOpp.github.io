// Build the Superspace Group Finder page script.
//   src/ui/main.ts      → ../../assets/blog/ssg-finder/ssgf.js   (finder page; ES module, minified, sourcemap)
//   src/ui/subgroups.ts → ../../assets/blog/ssg-finder/ssgsub.js (subgroup page; same options)
//   @spglib/moyo-wasm/moyo_wasm_bg.wasm → ../../assets/blog/ssg-finder/
// Ships no BYU data (decision D1): only our code, moyo's JS glue and moyo's wasm.
import { build } from "esbuild";
import { copyFileSync, existsSync, mkdirSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { createRequire } from "node:module";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const out = join(here, "..", "..", "assets", "blog", "ssg-finder");
mkdirSync(out, { recursive: true });

const common = {
  bundle: true,
  format: "esm",
  platform: "browser",
  target: ["es2020"],
  minify: true,
  sourcemap: true,
  legalComments: "eof",
  logLevel: "warning",
};
const moyoNote = "Includes moyo-wasm (MIT, (c) 2023 Kohei Shinohara): see moyo_LICENSE.txt";
await build({
  ...common,
  entryPoints: [join(here, "src", "ui", "main.ts")],
  outfile: join(out, "ssgf.js"),
  banner: { js: `/* Superspace Group Finder. ${moyoNote} */` },
});
await build({
  ...common,
  entryPoints: [join(here, "src", "ui", "subgroups.ts")],
  outfile: join(out, "ssgsub.js"),
  banner: { js: `/* Superspace subgroups (k-index 1). ${moyoNote} */` },
});

const require = createRequire(import.meta.url);
const pkgDir = dirname(require.resolve("@spglib/moyo-wasm/package.json"));
copyFileSync(join(pkgDir, "moyo_wasm_bg.wasm"), join(out, "moyo_wasm_bg.wasm"));

// Licence: moyo is MIT OR Apache-2.0; we distribute it under MIT. The npm package
// ships no LICENSE file, so the MIT text (with its copyright line) is vendored in
// licenses/moyo-LICENSE-MIT.txt, taken from the spglib/moyo repository at the tag
// matching the installed version, and shipped next to the wasm.
const pkg = JSON.parse(readFileSync(join(pkgDir, "package.json"), "utf8"));
const MOYO_LICENSED_VERSION = "0.22.0";
if (pkg.version !== MOYO_LICENSED_VERSION)
  throw new Error(`@spglib/moyo-wasm is ${pkg.version}; refresh licenses/moyo-LICENSE-MIT.txt from the v${pkg.version} tag and update MOYO_LICENSED_VERSION`);
const mit = readFileSync(join(here, "licenses", "moyo-LICENSE-MIT.txt"), "utf8");
writeFileSync(
  join(out, "moyo_LICENSE.txt"),
  `moyo_wasm_bg.wasm and the moyo-wasm JavaScript bindings bundled in ssgf.js and ssgsub.js are\n` +
    `${pkg.name} ${pkg.version} (https://github.com/spglib/moyo), dual-licensed\n` +
    `MIT OR Apache-2.0, distributed here under the MIT licence:\n\n` + mit,
);

const kb = (f) => (statSync(join(out, f)).size / 1024).toFixed(1) + " KB";
console.log(`ssgf.js ${kb("ssgf.js")}, ssgsub.js ${kb("ssgsub.js")}, moyo_wasm_bg.wasm ${kb("moyo_wasm_bg.wasm")} → ${out}`);
