// Superspace subgroups: page script (vanilla DOM). Bundled by build.mjs into
// assets/blog/ssg-finder/ssgsub.js and loaded by blog/superspace-subgroups.html.
//
// Paste mode only (decision D1): the page ships no BYU data and never fetches BYU. Subgroup
// types are named algorithmically (standardise.ts). "Send to finder" hands the subgroup,
// computed from the viewer's own paste, to the finder page through sessionStorage.

import { loadMoyo } from "../moyo";
import { parseByuPage, PasteError } from "../paste";
import type { GroupRecord } from "../ssg";
import { maxK, keep, subgroupSteps, SubgroupError, type SubFilter, type SubRow, type SubResult } from "../subfinder";
import type { SubRequest, SubResponse } from "./subworker";
import type { CrystalSystem } from "../symbol";
import { LABEL_RE, esc, symbolHtml } from "./flow";
import { SEND_KEY, affineText, displayMap, prepareSend } from "./subflow";

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
const elKind = $<HTMLSelectElement>("ssgs-kind");
const elKNote = $("ssgs-k-note");
const elOnlyK = $<HTMLInputElement>("ssgs-onlyk");
const elPg = $<HTMLSelectElement>("ssgs-pg");
const elSys = $<HTMLSelectElement>("ssgs-sys");
const elCen = $<HTMLSelectElement>("ssgs-cen");
const elRun = $<HTMLButtonElement>("ssgs-run");
const elCancel = $<HTMLButtonElement>("ssgs-cancel");
const elRunStatus = $("ssgs-run-status");
const elResults = $("ssgs-results");

interface State {
  parent: GroupRecord;
  parentSgNumber: number;
  name: string;
  rows: SubRow[];
  warnings: string[];
  ms: number;
}
let state: State | null = null;
let runId = 0;
let pasted: GroupRecord | null = null;
let activeWorker: Worker | null = null;
const supported = (rec: GroupRecord) => rec.d >= 1 && rec.d <= 3;
function stopRun(): void {
  ++runId;
  activeWorker?.terminate();
  activeWorker = null;
  elCancel.hidden = true;
}
const expanded = new Set<string>();
const rowKey = (row: SubRow) => JSON.stringify([row.kIndex,row.lattice?.basis,[...row.generators].sort()]);

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
  stopRun();
  pasted = null;
  state = null;
  expanded.clear();
  elResults.removeAttribute("aria-busy");
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
    pasted = rec;
    const cap = maxK(rec.d);
    for(const option of Array.from(elK.options)) option.disabled=Number(option.value)>cap;
    if(Number(elK.value)>cap)elK.value=String(cap);
    elKNote.textContent=`d = ${rec.d}: k ≤ ${cap}. Uncheck “Only this k-index” to include smaller k.`;
    if (rec.label && elLabel.value.trim() !== rec.label) {
      elLabel.value = rec.label;
      syncLabel();
    }
    const dNote = rec.d === 1 ? "" : ` d = ${rec.d}: type names use this tool’s own canonical setting.${rec.d === 3 ? " Cubic groups can take a minute or more." : ""}`;
    elPasteStatus.innerHTML = `Read ${rec.label ? esc(rec.label) + " " : ""}<span class="mono">${esc(rec.sym)}</span>, ${rec.nops ?? "?"} operators.${esc(dNote)}`;
    elPasteStatus.className = "ssgf-msg";
    elRun.disabled = !supported(rec);
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
  const parent = pasted;
  if (!parent || !supported(parent)) return;
  stopRun();
  const id = runId;
  const raw = elPaste.value;
  const f = { kIndex: Number(elK.value) };
  state = null;
  expanded.clear();
  elResults.innerHTML = "";
  elRun.disabled = true;
  elRunStatus.classList.remove("err");
  elCancel.hidden = false;
  elRunStatus.innerHTML = '<span class="ssgf-spin" aria-hidden="true"></span> Computing…';
  elResults.setAttribute("aria-busy", "true");
  const t0 = performance.now();
  const finish = (message: string) => {
    if (id !== runId) return;
    activeWorker?.terminate();
    activeWorker = null;
    elCancel.hidden = true;
    elRun.disabled = !(pasted && supported(pasted));
    elResults.removeAttribute("aria-busy");
    elRunStatus.textContent = message;
  };
  const accept = (out: SubResult) => {
    if (id !== runId) return;
    state = { parent, parentSgNumber: out.parentSgNumber, name: out.parent, rows: out.rows, warnings: out.warnings, ms: performance.now() - t0 };
    fillPointGroups(out.rows);
    render();
  };
  let worker: Worker | null = null;
  try {
    worker = new Worker(new URL("./ssgsubw.js", import.meta.url), { type: "module" });
  } catch {
    // Construction can fail in browsers or policies without module-worker support.
  }
  if (worker) {
    activeWorker = worker;
    worker.onmessage = (event: MessageEvent<SubResponse>) => {
      if (id !== runId) return;
      const message = event.data;
      if (message.type === "progress") {
        elRunStatus.textContent = `Named ${message.done} of ${message.total} classes…`;
      } else if (message.type === "rows") {
        accept(message.result);
        if (message.complete) finish("");
      } else {
        finish(message.message);
        elRunStatus.classList.add("err");
      }
    };
    worker.onerror = (event) => {
      if (id !== runId) return;
      event.preventDefault();
      finish(event.message || "Subgroup worker failed.");
      elRunStatus.classList.add("err");
    };
    elRunStatus.classList.remove("err");
    worker.postMessage({ type: "run", raw, filter: f } satisfies SubRequest);
    return;
  }
  elRunStatus.textContent = "Worker unavailable; computing here (Cancel is checked between classes)…";
  try {
    await moyoReady;
    await nextFrame();
    if (id !== runId) return;
    for (const step of subgroupSteps(parent, f)) {
      if (id !== runId) return;
      accept(step.result);
      elRunStatus.textContent = `Named ${step.done} of ${step.total} classes (main thread)…`;
      await nextFrame();
      if (id !== runId) return;
    }
    finish("");
  } catch (e) {
    if (id !== runId) return;
    const message = e instanceof SubgroupError || e instanceof PasteError ? e.message : "Internal error: " + (e instanceof Error ? e.message : String(e));
    finish(message);
    elRunStatus.classList.add("err");
  }
}

elCancel.addEventListener("click", () => {
  stopRun();
  elResults.removeAttribute("aria-busy");
  elRun.disabled = !(pasted && supported(pasted));
  elRunStatus.textContent = "Cancelled.";
});

elForm.addEventListener("submit", (e) => {
  e.preventDefault();
  if (!elRun.disabled) void run();
});
for (const el of [elOnlyK, elKind, elPg, elSys, elCen]) el.addEventListener("change", () => state && render());

elK.addEventListener("change",()=>{
  stopRun(); state=null; expanded.clear(); elResults.innerHTML="";
  elRun.disabled=!(pasted&&supported(pasted));
  elRunStatus.textContent="Press List subgroups to compute this k-index range.";
});

function fillPointGroups(rows: SubRow[]): void {
  const old = elPg.value;
  const pgs = [...new Set(rows.filter((r) => r.type).map((r) => r.pointGroup))];
  elPg.innerHTML = `<option value="">any</option>` + pgs.map((p) => `<option value="${esc(p)}">${esc(p)}</option>`).join("");
  elPg.value = pgs.includes(old) ? old : "";
}

function filter(): SubFilter {
  return {
    kIndex: Number(elK.value),
    kinds: elKind.value ? [elKind.value as "3d" | "internal" | "mixed"] : undefined,
    onlyK: elOnlyK.checked,
    pointGroup: elPg.value || undefined,
    system: (elSys.value || undefined) as CrystalSystem | undefined,
    centring: elCen.value || undefined,
  };
}

// ------------------------------------------------------------------ results

const sym = (s: string) => `<span class="sym">${symbolHtml(s)}</span>`;

function typeCell(r: SubRow): string {
  if (!r.type) return `<span class="err">unidentified</span>`;
  let h = sym(r.type.symbol.full);
  if (r.type.alternates.length)
    h += `<div class="alt">≡ ${r.type.alternates.map(sym).join(" ≡ ")}</div>`;
  if (r.type.E.length >= 2)
    h += `<div class="alt">symbol given in our setting; BYU may list the same group in an axis- or q-permuted setting.</div>`;
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
  h += `<p class="ssgf-head">Subgroups of ${s.parent.label ? esc(s.parent.label) + " " : ""}${sym(s.name)} with k-index ≤ ${elK.value}: ` +
    `<strong>${rows.length}</strong> of ${s.rows.length} conjugacy class${s.rows.length === 1 ? "" : "es"} shown (${time}).</p>`;
  h += `<div class="table-scroll ssgs-table"><table><caption class="sr-only">Subgroup conjugacy classes</caption><thead><tr>` +
    `<th scope="col">Index</th><th scope="col">t-index</th><th scope="col">k-index</th><th scope="col">Point group</th>` +
    `<th scope="col">Basic SG</th><th scope="col">Superspace group (type)</th><th scope="col">q</th>` +
    `<th scope="col">Lattice</th><th scope="col">Conjugates</th><th scope="col">Maximal</th>` +
    `<th scope="col"><abbr title="New external and internal coordinates of the type's setting in terms of the parent's coordinates">Transformation</abbr></th>` +
    `</tr></thead><tbody>`;
  rows.forEach((r, i) => {
    const did = `ssgs-det-${i}`;
    const open = expanded.has(rowKey(r));
    h += `<tr>` +
      `<td><button type="button" class="ssgs-x" aria-expanded="${open}" data-row="${s.rows.indexOf(r)}" aria-controls="${did}" title="Show generators">` +
      `<span aria-hidden="true">${open ? "▾" : "▸"}</span> ${r.index}</button></td>` +
      `<td>${r.tIndex}</td><td>${r.kIndex}</td>` +
      `<td>${r.type ? sym(r.pointGroup) : "?"}</td>` +
      `<td>${r.type ? `${sym(r.sgSymbol)} <span class="muted">(${r.sgNumber})</span>` : "?"}</td>` +
      `<td class="ty">${typeCell(r)}</td>` +
      `<td class="mono">${r.type ? "(" + esc(r.type.qActual.join(")(")) + ")" : "?"}</td>` +
      `<td>${r.lattice ? `${esc(r.kIndex === 1 ? "unchanged" : r.lattice.kind === "3d" ? "3D supercell" : r.lattice.kind)}<div class="mono">${esc(r.lattice.basis)}</div>` : "?"}</td>` +
      `<td>${r.classSize}</td><td>${r.maximal ? "yes" : "no"}</td>` +
      `<td class="mono">${r.type ? esc(affineText(displayMap(r.transformation ?? r.type))) : "?"}</td></tr>`;
    h += `<tr id="${did}" class="ssgs-det"${open ? "" : " hidden"}><td colspan="11"><div class="ssgs-detbox">` +
      `<p><span class="muted">In the parent's axes:</span> ${r.inParentAxes ? sym(r.inParentAxes) : `<span class="muted">no symbol in the parent’s axes</span>`}</p>` +
      `<p><span class="muted">Operators generating H (parent axes, besides the displayed subgroup lattice translations):</span></p>` +
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
    const row = state?.rows[Number(btn.dataset.row)];
    if (row) {
      if (open) expanded.add(rowKey(row));
      else expanded.delete(rowKey(row));
    }
    return;
  }
  sendToFinder(e);
});
elResults.addEventListener("auxclick", (e) => {
  if (e.button === 1) sendToFinder(e);
});

function sendToFinder(e: MouseEvent): void {
  const a = (e.target as HTMLElement).closest<HTMLAnchorElement>("a[data-send]");
  if (!a || !state) return;
  const row = state.rows[Number(a.dataset.send)];
  const out = prepareSend(state.parent, row, state.parentSgNumber);
  const note = a.nextElementSibling as HTMLElement | null;
  if (!out.ok) {
    e.preventDefault();
    if (note) note.textContent = `Could not build the subgroup record: ${out.reason}`;
    return;
  }
  const middle = e.type === "auxclick" && e.button === 1;
  if (middle) e.preventDefault();
  try {
    sessionStorage.setItem(SEND_KEY, out.json);
  } catch {
    e.preventDefault();
    if (note) note.textContent = "Could not hand the subgroup over (browser storage is blocked).";
    return;
  }
  if (middle) {
    // Some browsers do not clone sessionStorage for middle-clicked links. Initialize a
    // same-origin blank tab explicitly before navigating it to the finder.
    let child: Window | null = null;
    try {
      child = window.open("about:blank", "_blank");
    } catch {
      // Treat an exception the same as the popup blocker's null result below.
    }
    if (!child) {
      if (note) note.textContent = "Could not open the finder in a new tab (browser popups are blocked).";
      return;
    }
    try {
      child.sessionStorage.setItem(SEND_KEY, out.json);
    } catch {
      child.close();
      if (note) note.textContent = "Could not hand the subgroup over (browser storage is blocked).";
      return;
    }
    child.location.href = a.href;
  }
  if (note) note.textContent = "opens the Superspace Group Finder with this subgroup's operators";
}

syncLabel();
elRun.disabled = true;
if (elPaste.value) onPaste(); // text kept by the browser on reload
