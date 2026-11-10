// Build the Superspace Group Finder page script.
//   src/ui/main.ts  → ../../assets/blog/ssg-finder/ssgf.js (ES module, minified, sourcemap)
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

await build({
  entryPoints: [join(here, "src", "ui", "main.ts")],
  outfile: join(out, "ssgf.js"),
  bundle: true,
  format: "esm",
  platform: "browser",
  target: ["es2020"],
  minify: true,
  sourcemap: true,
  legalComments: "eof",
  logLevel: "warning",
});

const require = createRequire(import.meta.url);
const pkgDir = dirname(require.resolve("@spglib/moyo-wasm/package.json"));
copyFileSync(join(pkgDir, "moyo_wasm_bg.wasm"), join(out, "moyo_wasm_bg.wasm"));

// Licence: copy the package's LICENSE file(s) if it ships any; otherwise write a notice
// from its package.json (moyo-wasm 0.22.0 ships none; it is "MIT OR Apache-2.0").
const pkg = JSON.parse(readFileSync(join(pkgDir, "package.json"), "utf8"));
const licFiles = ["LICENSE", "LICENSE.md", "LICENSE.txt", "LICENSE-MIT", "LICENSE-APACHE", "COPYING"].filter((f) => existsSync(join(pkgDir, f)));
for (const f of licFiles) copyFileSync(join(pkgDir, f), join(out, `moyo_wasm_${f}${/\.\w+$/.test(f) ? "" : ".txt"}`));
if (!licFiles.length) {
  const repo = typeof pkg.repository === "string" ? pkg.repository : pkg.repository?.url ?? "";
  writeFileSync(
    join(out, "moyo_wasm_NOTICE.txt"),
    `moyo_wasm_bg.wasm and the moyo-wasm JavaScript bindings bundled in ssgf.js are\n` +
      `${pkg.name} ${pkg.version}, part of the spglib project's moyo library.\n` +
      `License: ${pkg.license}\nSource: ${repo.replace(/^git\+/, "").replace(/\.git$/, "")}\n` +
      `The npm package ships no separate LICENSE file; see the repository for the licence texts.\n`,
  );
}

const kb = (f) => (statSync(join(out, f)).size / 1024).toFixed(1) + " KB";
console.log(`ssgf.js ${kb("ssgf.js")}, moyo_wasm_bg.wasm ${kb("moyo_wasm_bg.wasm")} → ${out}`);
