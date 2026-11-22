// Page logic for blog/climbing-hold-generator.html.

import { HOLDS } from "./holds";
import { edgeParams } from "./holds/machinable";
import { validateEdge } from "./edge/params";
import { DEFAULT_OPTIONS, type BuildOptions, type BuildResult } from "./pipeline";
import { ALUMINIUM_DENSITY, MAX_OVERHANG_FRACTION, printEstimate } from "./printcheck";
import { defaults, quantize, type HoldDef, type Param, type Values } from "./spec";
import { toBinarySTL } from "./stl";
import { TEX_MODES, type TexMode } from "./texture";
import { createViewer } from "./viewer";
import type { MeshReply, MeshRequest } from "./workers/mesh";
import type { StepReply } from "./workers/step";

interface Scan {
  id: string; file: string; name: string; source: string; author: string; authorUrl: string; license: string; licenseUrl: string;
  /** set on load: which asset folder it came from */
  dir?: string;
  /** default "rotation on the wall", degrees */
  rotation?: number;
}
/** Fixed-shape STLs: asset folder, menu group, whether exports must stay share-alike. */
const LIBRARIES = [
  { dir: "lafaye", group: "Scanned (CC BY-SA)", shareAlike: true },
  { dir: "wupinstick", group: "Fixed shapes (CC BY)", shareAlike: false },
];

const here = new URL(".", import.meta.url);
const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

const canvas = $<HTMLCanvasElement>("ch-gl");
const msg = $("ch-msg");
const status = $("ch-status");
const viewer = createViewer(canvas);
if (!viewer) { msg.hidden = false; msg.textContent = "WebGL is not available in this browser, so there is no preview. The downloads still work."; }

const holdSel = $<HTMLSelectElement>("ch-hold");
const blurb = $("ch-blurb");
const credit = $("ch-credit");
const paramsBox = $("ch-params");
const screwsBox = $<HTMLInputElement>("ch-screws");
const texMode = $<HTMLSelectElement>("ch-tex");
const texH = $<HTMLInputElement>("ch-texh");
const texS = $<HTMLInputElement>("ch-texs");
const texD = $<HTMLInputElement>("ch-texd");
const rot = $<HTMLInputElement>("ch-rot");
const quality = $<HTMLSelectElement>("ch-quality");
const out = $("ch-out");
const warn = $("ch-warn");
const stlBtn = $<HTMLButtonElement>("ch-stl");
const jsonBtn = $<HTMLButtonElement>("ch-json");
const stepBtn = $<HTMLButtonElement>("ch-step");
const stepBox = $("ch-stepbox");
const stepOut = $("ch-stepout");

let scans: Scan[] = [];
const FIRST = HOLDS.find((h) => h.id === "system-edge") ?? HOLDS[0];
/** The holds starred in the desktop app (the licensed ones that made it here). */
const STARRED = new Set(["system-edge", "gradient-edge", "wupinstick-big-incut"]);
let hold: HoldDef | Scan = FIRST;
let values: Values = defaults(FIRST);
let last: { result: BuildResult; overhang: number } | null = null;
const isScan = (h: HoldDef | Scan): h is Scan => !("build" in h);

// ------------------------------------------------------------------ workers

const meshWorker = new Worker(new URL("chgw.js", here), { type: "module" });
let reqId = 0, inFlight = false, queued: MeshRequest | null = null;
const pending = new Map<number, (r: MeshReply) => void>();
meshWorker.onmessage = (e: MessageEvent<MeshReply>) => {
  const cb = pending.get(e.data.id);
  pending.delete(e.data.id);
  cb?.(e.data);
};

function options(final = false): BuildOptions {
  return {
    ...DEFAULT_OPTIONS,
    screws: screwsBox.checked,
    texMode: texMode.value as TexMode,
    texHeight: +texH.value,
    texSize: +texS.value,
    density: +texD.value,
    rotation: +rot.value,
    segments: final ? 96 : +quality.value,
  };
}

function request(final = false): MeshRequest {
  const id = ++reqId;
  return isScan(hold)
    ? { kind: "scan", id, url: new URL(`${hold.dir}/${hold.file}`, here).href, opt: options(final) }
    : { kind: "build", id, hold: hold.id, values: { ...values }, opt: options(final) };
}

function send(req: MeshRequest): Promise<MeshReply> {
  return new Promise((res) => { pending.set(req.id, res); meshWorker.postMessage(req); });
}

/** Rebuild the preview; requests made while one is running coalesce. */
function rebuild() {
  const req = request();
  if (inFlight) { queued = req; return; }
  run(req);
}

async function run(req: MeshRequest) {
  inFlight = true;
  status.textContent = "building…";
  const t0 = performance.now();
  const r = await send(req);
  inFlight = false;
  if (queued) { const q = queued; queued = null; run(q); return; }
  if (!r.ok) { status.textContent = ""; showError(r.error); return; }
  last = { result: r.result, overhang: r.overhang };
  viewer?.setMesh(r.result.mesh.pos, r.result.mesh.idx);
  status.textContent = `${r.result.triangles.toLocaleString()} triangles · ${Math.round(performance.now() - t0)} ms`;
  report();
}

function showError(e: string) {
  warn.replaceChildren(li(`Could not build: ${e}`));
}

function li(text: string) { const el = document.createElement("li"); el.textContent = text; return el; }

function report() {
  if (!last) return;
  const { result: r, overhang } = last;
  const { grams, minutes } = printEstimate(r.volume);
  const rows: [string, string][] = [
    ["Size", `${r.bbox.map((x) => x.toFixed(0)).join(" × ")} mm (off wall × along × up)`],
    ["Volume", `${(r.volume / 1000).toFixed(1)} cm³`],
    ["Print (PLA, 30 %)", `≈ ${grams.toFixed(0)} g, ${minutes >= 90 ? (minutes / 60).toFixed(1) + " h" : minutes.toFixed(0) + " min"}`],
    ["Needs support", `${(100 * overhang).toFixed(0)} % of the surface`],
  ];
  if (!isScan(hold) && hold.step) rows.push(["In 6061 aluminium", `≈ ${((r.volume / 1000) * ALUMINIUM_DENSITY).toFixed(0)} g`]);
  out.replaceChildren(...rows.flatMap(([k, v]) => { const dt = document.createElement("dt"); dt.textContent = k; const dd = document.createElement("dd"); dd.textContent = v; return [dt, dd]; }));
  const w: string[] = [];
  if (overhang > MAX_OVERHANG_FRACTION) w.push(`More than ${100 * MAX_OVERHANG_FRACTION} % of the surface overhangs when printed back-down: expect support marks on the grip.`);
  if (r.genus < 0) w.push("The mesh came out with separate pieces: check the settings.");
  if (!isScan(hold) && hold.step && texMode.value !== "none") w.push("The texture is in the preview and the STL only. STEP files are always smooth: the bead blast gives the metal its texture.");
  warn.replaceChildren(...w.map(li));
}

// ---------------------------------------------------------------- controls

function rowFor(p: Param): HTMLElement {
  const row = document.createElement("div");
  row.className = "ch-row";
  const id = `ch-p-${p.key}`;
  const label = document.createElement("label");
  label.htmlFor = id;
  label.textContent = p.label;
  if (p.tip) label.title = p.tip;
  const range = document.createElement("input");
  range.type = "range"; range.id = id; range.min = String(p.lo); range.max = String(p.hi); range.step = String(p.step ?? 1);
  range.value = String(values[p.key]);
  const num = document.createElement("input");
  num.type = "number"; num.min = range.min; num.max = range.max; num.step = range.step; num.value = range.value;
  num.setAttribute("aria-label", `${p.label} value`);
  const unit = document.createElement("span");
  unit.className = "unit";
  unit.textContent = (p.unit ?? " mm").trim();
  const set = (raw: number, from: HTMLInputElement) => {
    if (!Number.isFinite(raw)) return;
    const v = quantize(p, raw);
    values[p.key] = v;
    if (from !== range) range.value = String(v);
    if (from !== num) num.value = String(v);
    rebuild();
    touched();
  };
  range.addEventListener("input", () => set(+range.value, range));
  num.addEventListener("change", () => set(+num.value, num));
  row.append(label, range, num, unit);
  if (p.tip) { const n = document.createElement("p"); n.className = "ch-tip"; n.textContent = p.tip; row.append(n); }
  return row;
}

/** `preset` keeps the given slider values; `keepRot` keeps the rotation box (undo/redo). */
function selectHold(id: string, preset?: Values, keepRot = !!preset) {
  const def = HOLDS.find((h) => h.id === id);
  const scan = scans.find((s) => s.id === id);
  const same = hold.id === (def ?? scan ?? HOLDS[0]).id;
  hold = def ?? scan ?? HOLDS[0];
  values = def ? { ...defaults(def), ...preset } : {};
  holdSel.value = hold.id;
  if (!keepRot) rot.value = String(scan?.rotation ?? 0);
  paramsBox.replaceChildren(...(def ? def.params.map(rowFor) : []));
  $("ch-shape").hidden = !def;
  screwsBox.disabled = !def;
  if (def) {
    blurb.textContent = def.blurb;
    credit.hidden = true;
  } else if (scan) {
    const lib = LIBRARIES.find((l) => l.dir === scan.dir)!;
    blurb.textContent = lib.shareAlike
      ? "A fixed scanned shape: texture, rotation and STL export apply, sliders do not. No screw holes are added; drill or design your own fixing."
      : "A fixed shape: texture, rotation and STL export apply, sliders do not. It keeps the designer's own bolt and screw holes (M10 bolt, 5 mm wood screws).";
    credit.hidden = false;
    credit.innerHTML = "";
    credit.append("“", link(scan.name, scan.source), "” by ", link(scan.author, scan.authorUrl), ", ", link(scan.license, scan.licenseUrl),
      lib.shareAlike
        ? ". Decimated and re-oriented here; anything you export from it is shared under the same licence."
        : ". Re-oriented here; credit the designer if you share what you export.");
  }
  const step = !!def?.step;
  stepBox.hidden = !step;
  stepOut.replaceChildren();
  if (step && !preset) texMode.value = "none";
  last = null;
  if (!same) viewer?.resetView();
  rebuild();
  touched();
}

function link(text: string, href: string) {
  const a = document.createElement("a");
  a.href = href; a.textContent = text; a.rel = "noopener";
  return a;
}

function fillHoldMenu() {
  holdSel.replaceChildren();
  const groups = new Map<string, HTMLOptGroupElement>();
  const add = (group: string, id: string, name: string) => {
    let g = groups.get(group);
    if (!g) { g = document.createElement("optgroup"); g.label = group; groups.set(group, g); holdSel.append(g); }
    const o = document.createElement("option");
    o.value = id; o.textContent = STARRED.has(id) ? `★ ${name}` : name;
    g.append(o);
  };
  for (const h of HOLDS) add(h.group, h.id, h.name);
  for (const s of scans) add(LIBRARIES.find((l) => l.dir === s.dir)!.group, s.id, s.name);
  holdSel.value = hold.id;
}

for (const m of TEX_MODES) { const o = document.createElement("option"); o.value = m; o.textContent = m; texMode.append(o); }
texMode.value = "none";
holdSel.addEventListener("change", () => selectHold(holdSel.value));
for (const el of [screwsBox, texMode, quality]) el.addEventListener("change", () => { report(); rebuild(); touched(); });
for (const el of [texH, texS, texD, rot]) el.addEventListener("change", () => { rebuild(); touched(); });
texMode.addEventListener("change", () => { texD.disabled = !["divots", "grit"].includes(texMode.value); });
texD.disabled = true;
$("ch-view").addEventListener("click", () => viewer?.resetView());

// view only: none of these rebuild the hold
const wallBox = $("ch-wall") as HTMLInputElement, angleBox = $("ch-angle") as HTMLInputElement;
const projBtn = $("ch-proj") as HTMLButtonElement;
const syncScene = () => {
  const a = Number(angleBox.value);
  viewer?.setScene({ showWall: wallBox.checked, wallAngle: Number.isFinite(a) ? Math.max(-30, Math.min(90, a)) : 0 });
};
wallBox.addEventListener("change", syncScene);
angleBox.addEventListener("input", syncScene);
projBtn.addEventListener("click", () => {
  const on = projBtn.getAttribute("aria-pressed") !== "true";
  projBtn.setAttribute("aria-pressed", String(on));
  viewer?.setScene({ ortho: on });
});
$("ch-front").addEventListener("click", () => viewer?.view("front"));
$("ch-side").addEventListener("click", () => viewer?.view("side"));


// ------------------------------------------------------- undo, redo, random

const undoBtn = $<HTMLButtonElement>("ch-undo");
const redoBtn = $<HTMLButtonElement>("ch-redo");
const randBtn = $<HTMLButtonElement>("ch-rand");

interface Snap { hold: string; values: Values; screws: boolean; tex: string; texH: string; texS: string; texD: string; rot: string }
const snap = (): Snap => ({
  hold: hold.id, values: { ...values }, screws: screwsBox.checked,
  tex: texMode.value, texH: texH.value, texS: texS.value, texD: texD.value, rot: rot.value,
});
const history: Snap[] = [];
let at = -1, restoring = false, settle = 0;

/** Record the state once it has stopped changing (a slider drag is one step). */
function touched() {
  if (restoring) return;
  clearTimeout(settle);
  settle = window.setTimeout(record, 350);
}
function record() {
  const s = snap();
  if (at >= 0 && JSON.stringify(history[at]) === JSON.stringify(s)) return;
  history.splice(at + 1);
  history.push(s);
  if (history.length > 200) history.shift();
  at = history.length - 1;
  syncHistory();
}
function syncHistory() { undoBtn.disabled = at <= 0; redoBtn.disabled = at >= history.length - 1; }

function restore(s: Snap) {
  restoring = true;
  screwsBox.checked = s.screws;
  texMode.value = s.tex; texH.value = s.texH; texS.value = s.texS; texD.value = s.texD; rot.value = s.rot;
  texD.disabled = !["divots", "grit"].includes(texMode.value);
  selectHold(s.hold, s.values);
  restoring = false;
  syncHistory();
}
function move(by: number) {
  clearTimeout(settle);
  record(); // keep a change made in the last moment
  const to = at + by;
  if (to < 0 || to >= history.length) return;
  at = to;
  restore(history[at]);
}
undoBtn.addEventListener("click", () => move(-1));
redoBtn.addEventListener("click", () => move(1));
document.addEventListener("keydown", (e) => {
  if (!(e.metaKey || e.ctrlKey) || e.altKey) return;
  const t = e.target as HTMLElement;
  // leave typing fields their own text undo
  if (t instanceof HTMLInputElement && t.type === "number") return;
  const k = e.key.toLowerCase();
  if (k === "z" && !e.shiftKey) { e.preventDefault(); move(-1); }
  else if ((k === "z" && e.shiftKey) || k === "y") { e.preventDefault(); move(1); }
});

/** A random hold, and every shape slider somewhere in its range (retried until valid). */
randBtn.addEventListener("click", () => {
  const all: (HoldDef | Scan)[] = [...HOLDS, ...scans];
  const others = all.filter((h) => h.id !== hold.id);
  const pick = others[Math.floor(Math.random() * others.length)];
  if (isScan(pick)) { selectHold(pick.id); return; }
  const def = pick;
  let v: Values = {};
  for (let tries = 0; tries < 50; tries++) {
    v = {};
    for (const p of def.params) v[p.key] = quantize(p, p.lo + Math.random() * (p.hi - p.lo));
    if (!def.step) break;
    try { validateEdge(edgeParams(v, screwsBox.checked)); break; } catch { /* draw again */ }
  }
  selectHold(def.id, v, false);
});

// ---------------------------------------------------------------- downloads

function save(data: BlobPart, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([data], { type }));
  const a = document.createElement("a");
  a.href = url; a.download = name;
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

const slug = () => (isScan(hold) ? hold.id : hold.id);

stlBtn.addEventListener("click", async () => {
  stlBtn.disabled = true;
  status.textContent = "building the export at full resolution…";
  const r = await send(request(true));
  stlBtn.disabled = false;
  if (!r.ok) { showError(r.error); status.textContent = ""; return; }
  const note = isScan(hold) ? `${hold.name} by ${hold.author}, ${hold.license}` : `climbing hold: ${hold.name}`;
  save(toBinarySTL(r.result.mesh, note), `${slug()}.stl`, "model/stl");
  status.textContent = `saved ${slug()}.stl (${r.result.triangles.toLocaleString()} triangles)`;
});

jsonBtn.addEventListener("click", () => {
  const o = options();
  const data = { hold: hold.id, values, screws: o.screws, texture: { mode: o.texMode, height: o.texHeight, size: o.texSize, density: o.density }, rotation: o.rotation };
  save(JSON.stringify(data, null, 2), `${slug()}.json`, "application/json");
});

let stepWorker: Worker | null = null;
let stepId = 0;
stepBtn.addEventListener("click", () => {
  if (isScan(hold) || !hold.step) return;
  stepBtn.disabled = true;
  const first = !stepWorker;
  stepOut.replaceChildren(li(first ? "Loading the CAD kernel (about 7 MB, once)…" : "Building the solid…"));
  stepWorker ??= new Worker(new URL("chstepw.js", here), { type: "module" });
  const id = ++stepId;
  const params = edgeParams(values, screwsBox.checked);
  stepWorker.onmessage = (e: MessageEvent<StepReply>) => {
    if (e.data.id !== id) return;
    stepBtn.disabled = false;
    const r = e.data;
    if (!r.ok) { stepOut.replaceChildren(li(`OpenCascade could not build this one: ${r.error}`)); return; }
    save(r.step, `${slug()}.step`, "application/step");
    // show the exact B-rep (fillets included) in the viewer
    viewer?.setMesh(r.mesh.pos, r.mesh.idx);
    status.textContent = "showing the STEP solid";
    stepOut.replaceChildren(
      li(`Saved ${slug()}.step: one solid, ${r.faces} exact faces (planes, cylinders, cones, spheres and fillet blends), millimetres.`),
      li(`${r.bbox.map((x) => x.toFixed(1)).join(" × ")} mm · ${(r.volume / 1000).toFixed(2)} cm³ · ≈ ${((r.volume / 1000) * ALUMINIUM_DENSITY).toFixed(0)} g in 6061`),
      ...(r.flatArea > 0 ? [li(`Vice flats: ${r.flatArea.toFixed(0)} mm² each, ${r.bbox[1].toFixed(1)} mm across the jaws.`)] : [li("No vice flats: the shop will need soft jaws or a fixture to hold the rounded ends.")]),
    );
  };
  stepWorker.onerror = (e) => { stepBtn.disabled = false; stepOut.replaceChildren(li(`The CAD kernel failed to load: ${e.message}`)); };
  stepWorker.postMessage({ id, params });
});

// -------------------------------------------------------------------- start

fillHoldMenu();
selectHold(FIRST.id);
// browsers restore typed form fields on reload (and on back/forward), at
// different moments; make the view match them whenever that may happen
syncScene();
window.addEventListener("load", syncScene);
window.addEventListener("pageshow", syncScene);
Promise.all(LIBRARIES.map((lib) =>
  fetch(new URL(`${lib.dir}/manifest.json`, here))
    .then((r) => (r.ok ? r.json() : []))
    .then((list: Scan[]) => list.map((s) => ({ ...s, dir: lib.dir })))
    .catch(() => [] as Scan[]), // fixed shapes are optional
)).then((lists) => {
  scans = lists.flat();
  fillHoldMenu();
});
