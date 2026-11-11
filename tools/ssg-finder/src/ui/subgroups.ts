// Superspace subgroups (k-index 1): page script (vanilla DOM). Bundled by build.mjs into
// assets/blog/ssg-finder/ssgsub.js and loaded by blog/superspace-subgroups.html.
//
// Paste mode only (decision D1): the page ships no BYU data and never fetches BYU. Subgroup
// types are named algorithmically (standardise.ts). "Send to finder" hands the subgroup,
// computed from the viewer's own paste, to the finder page through sessionStorage.

import { loadMoyo } from "../moyo";
import { parseByuPage, PasteError } from "../paste";
import type { GroupRecord } from "../ssg";
import { keep, runSubgroups, SubgroupError, type SubFilter, type SubRow } from "../subfinder";
import type { CrystalSystem } from "../symbol";
import { LABEL_RE, esc, symbolHtml } from "./flow";
import { SEND_KEY, affineText, subgroupRecord } from "./subflow";

const $ = <T extends HTMLElement = HTMLElement>(id: string): T => {
  const el = document.getElementById(id);
  if (!el) throw new Error(`#${id} missing`);
  return el as T;
};

const moyoReady = loadMoyo(new URL("./moyo_wasm_bg.wasm", import.meta.url));
moyoReady.catch(() => {});

const elLabel = $<HTMLInputElement>("ssgs-label");
const elByu = $<HTMLFormElement>("ssgs-byu");
const elByuLabel = $<HTMLInputElement>("ssgs-byu-label");
const elByuNote = $("ssgs-byu-note");
const elPaste = $<HTMLTextAreaElement>("ssgs-paste");
const elPasteStatus = $("ssgs-paste-status");
const elForm = $<HTMLFormElement>("ssgs-form");
const elK = $<HTMLSelectElement>("ssgs-k");
const elOnlyK = $<HTMLInputElement>("ssgs-onlyk");
const elPg = $<HTMLSelectElement>("ssgs-pg");
const elSys = $<HTMLSelectElement>("ssgs-sys");
const elCen = $<HTMLSelectElement>("ssgs-cen");
const elRun = $<HTMLButtonElement>("ssgs-run");
const elRunStatus = $("ssgs-run-status");
const elResults = $("ssgs-results");

interface State {
  parent: GroupRecord;
  name: string;
  rows: SubRow[];
  warnings: string[];
  ms: number;
}
let state: State | null = null;
let runId = 0;

// ------------------------------------------------------------------ step 1: label + paste

function syncLabel(): void {
  const v = elLabel.value.trim();
  const ok = LABEL_RE.test(v);
  elByuLabel.value = ok ? v : "";
  elByu.hidden = !ok;
  elByuNote.hidden = ok;
  elLabel.setAttribute("aria-invalid", String(!!v && !ok));
}
elLabel.addEventListener("input", syncLabel);
elByu.addEventListener("submit", (e) => {
  if (!LABEL_RE.test(elByuLabel.value)) e.preventDefault();
});

function onPaste(): void {
  state = null;
  elResults.innerHTML = "";
  elRunStatus.textContent = "";
  const text = elPaste.value;
  if (!text.trim()) {
    elPasteStatus.textContent = "";
    elPasteStatus.className = "ssgf-msg";
    elRun.disabled = true;
    return;
  }
  try {
    const rec = parseByuPage(text);
    if (rec.label && elLabel.value.trim() !== rec.label) {
      elLabel.value = rec.label;
      syncLabel();
    }
    const dNote = rec.d === 1 ? "" : ` This group has d = ${rec.d}; only d = 1 is supported for now.`;
    elPasteStatus.innerHTML = `Read ${rec.label ? esc(rec.label) + " " : ""}<span class="mono">${esc(rec.sym)}</span>, ${rec.nops ?? "?"} operators.${esc(dNote)}`;
    elPasteStatus.className = rec.d === 1 ? "ssgf-msg" : "ssgf-msg err";
    elRun.disabled = rec.d !== 1;
  } catch (e) {
    elPasteStatus.textContent = e instanceof PasteError ? e.message : "Could not read this text as a BYU group page: " + (e instanceof Error ? e.message : String(e));
    elPasteStatus.className = "ssgf-msg err";
    elRun.disabled = true;
  }
}
elPaste.addEventListener("input", onPaste);

// ------------------------------------------------------------------ step 2: compute (all rows once; filters are applied here)

const nextFrame = () =>
  new Promise<void>((r) => {
    let done = false;
    const go = () => {
      if (!done) {
        done = true;
        r();
      }
    };
    requestAnimationFrame(() => setTimeout(go, 0));
    setTimeout(go, 50);
  });

async function run(): Promise<void> {
  const id = ++runId;
  elRun.disabled = true;
  elRunStatus.innerHTML = '<span class="ssgf-spin" aria-hidden="true"></span> Computing…';
  elResults.setAttribute("aria-busy", "true");
  const done = (msg: string) => {
    elRun.disabled = false;
    elResults.removeAttribute("aria-busy");
    elRunStatus.innerHTML = msg;
  };
  try {
    await moyoReady;
  } catch (e) {
    done(`<span class="err">Could not load the symmetry library (moyo-wasm): ${esc(e instanceof Error ? e.message : String(e))}</span>`);
    return;
  }
  await nextFrame();
  if (id !== runId) return;
  const t0 = performance.now();
  try {
    const parent = parseByuPage(elPaste.value);
    const out = runSubgroups(elPaste.value, { kIndex: Number(elK.value) });
    state = { parent, name: out.parent, rows: out.rows, warnings: out.warnings, ms: performance.now() - t0 };
  } catch (e) {
    state = null;
    elResults.innerHTML = "";
    const msg = e instanceof SubgroupError || e instanceof PasteError ? e.message : "Internal error: " + (e instanceof Error ? e.message : String(e));
    done(`<span class="err">${esc(msg)}</span>`);
    return;
  }
  fillPointGroups(state.rows);
  done("");
  render();
}

elForm.addEventListener("submit", (e) => {
  e.preventDefault();
  if (!elRun.disabled) void run();
});
for (const el of [elOnlyK, elPg, elSys, elCen]) el.addEventListener("change", () => state && render());

function fillPointGroups(rows: SubRow[]): void {
  const old = elPg.value;
  const pgs = [...new Set(rows.filter((r) => r.type).map((r) => r.pointGroup))];
  elPg.innerHTML = `<option value="">any</option>` + pgs.map((p) => `<option value="${esc(p)}">${esc(p)}</option>`).join("");
  elPg.value = pgs.includes(old) ? old : "";
}

function filter(): SubFilter {
  return {
    kIndex: 1,
    onlyK: elOnlyK.checked,
    pointGroup: elPg.value || undefined,
    system: (elSys.value || undefined) as CrystalSystem | undefined,
    centring: elCen.value || undefined,
  };
}

// ------------------------------------------------------------------ results

const sym = (s: string) => `<span class="sym">${symbolHtml(s)}</span>`;

function typeCell(r: SubRow): string {
  if (!r.type) return `<span class="err">not identified</span>`;
  let h = sym(r.type.symbol.full);
  if (r.type.alternates.length)
    h += `<div class="alt">≡ ${r.type.alternates.map(sym).join(" ≡ ")}</div>`;
  return h;
}

function render(): void {
  if (!state) return;
  const s = state;
  const f = filter();
  const rows = s.rows.filter((r) => keep(r, f));
  let h = "";
  if (s.warnings.length)
    h += `<div class="ssgf-warn" role="alert"><strong>Warnings</strong><ul>${s.warnings.map((w) => `<li>${esc(w)}</li>`).join("")}</ul></div>`;
  const time = s.ms < 1000 ? Math.round(s.ms) + " ms" : (s.ms / 1000).toFixed(1) + " s";
  h += `<p class="ssgf-head">Subgroups of ${s.parent.label ? esc(s.parent.label) + " " : ""}${sym(s.name)} with k-index 1: ` +
    `<strong>${rows.length}</strong> of ${s.rows.length} conjugacy class${s.rows.length === 1 ? "" : "es"} shown (${time}).</p>`;
  h += `<div class="table-scroll ssgs-table"><table><caption class="sr-only">Subgroup classes of k-index 1</caption><thead><tr>` +
    `<th scope="col">Index</th><th scope="col">t-index</th><th scope="col">k-index</th><th scope="col">Point group</th>` +
    `<th scope="col">Basic SG</th><th scope="col">Superspace group (type)</th><th scope="col">q</th>` +
    `<th scope="col">Conjugates</th><th scope="col">Maximal</th>` +
    `<th scope="col"><abbr title="New coordinates (x', y', z'; t') of the type's setting in terms of the parent's (x, y, z; t)">Transformation</abbr></th>` +
    `</tr></thead><tbody>`;
  rows.forEach((r, i) => {
    const did = `ssgs-det-${i}`;
    h += `<tr>` +
      `<td><button type="button" class="ssgs-x" aria-expanded="false" aria-controls="${did}" title="Show generators">` +
      `<span aria-hidden="true">▸</span> ${r.index}</button></td>` +
      `<td>${r.tIndex}</td><td>${r.kIndex}</td>` +
      `<td>${r.type ? sym(r.pointGroup) : "?"}</td>` +
      `<td>${r.type ? `${sym(r.sgSymbol)} <span class="muted">(${r.sgNumber})</span>` : "?"}</td>` +
      `<td class="ty">${typeCell(r)}</td>` +
      `<td class="mono">${r.type ? "(" + esc(r.type.qActual) + ")" : "?"}</td>` +
      `<td>${r.classSize}</td><td>${r.maximal ? "yes" : "no"}</td>` +
      `<td class="mono">${r.type ? esc(affineText(r.type)) : "?"}</td></tr>`;
    h += `<tr id="${did}" class="ssgs-det" hidden><td colspan="10"><div class="ssgs-detbox">` +
      `<p><span class="muted">In the parent's axes:</span> ${r.inParentAxes ? sym(r.inParentAxes) : `<span class="muted">no symbol (different crystal system)</span>`}</p>` +
      `<p><span class="muted">Generators (parent axes, besides the lattice${s.parent.cen.length ? " and centring" : ""} translations):</span></p>` +
      `<pre>${r.generators.length ? esc(r.generators.join("\n")) : "none (only translations)"}</pre>` +
      `<p><a href="superspace-group-finder.html" data-send="${s.rows.indexOf(r)}">Send to finder →</a> ` +
      `<span class="muted">opens the Superspace Group Finder with this subgroup's operators</span></p>` +
      `</div></td></tr>`;
  });
  h += `</tbody></table></div>`;
  if (!rows.length) h += `<p class="muted">No subgroup class matches these filters.</p>`;
  elResults.innerHTML = h;
}

elResults.addEventListener("click", (e) => {
  const t = e.target as HTMLElement;
  const btn = t.closest<HTMLButtonElement>("button.ssgs-x");
  if (btn) {
    const det = document.getElementById(btn.getAttribute("aria-controls")!)!;
    const open = det.hidden;
    det.hidden = !open;
    btn.setAttribute("aria-expanded", String(open));
    btn.querySelector("span")!.textContent = open ? "▾" : "▸";
    return;
  }
  const a = t.closest<HTMLAnchorElement>("a[data-send]");
  if (a && state) {
    const row = state.rows[Number(a.dataset.send)];
    try {
      sessionStorage.setItem(SEND_KEY, JSON.stringify(subgroupRecord(state.parent, row)));
    } catch {
      e.preventDefault();
      const note = a.nextElementSibling as HTMLElement | null;
      if (note) note.textContent = "Could not hand the subgroup over (browser storage is blocked).";
    }
  }
});

syncLabel();
elRun.disabled = true;
if (elPaste.value) onPaste(); // text kept by the browser on reload
