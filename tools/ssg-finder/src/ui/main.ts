// Superspace Group Finder: page script (vanilla DOM). Bundled by build.mjs into
// assets/blog/ssg-finder/ssgf.js and loaded by blog/superspace-group-finder.html.
//
// Paste mode only (decision D1): the page ships no BYU data and never fetches BYU.
// The "Open BYU page" button is a plain form post from the viewer's own browser.

import { loadMoyo } from "../moyo";
import type { GroupRecord } from "../ssg";
import type { FinderResult, OrbitRow, ClassResult } from "../finder";
import {
  LABEL_RE, parsePaste, validateParam, encodeHash, decodeHash, compute, esc, symbolHtml, classSymbol,
  summaryColumns, subDigits, tText, transformText, supercellText, resultJson, exportName, storageKey,
} from "./flow";
import { SEND_KEY, decodeSent } from "./subflow";

const $ = <T extends HTMLElement = HTMLElement>(id: string): T => {
  const el = document.getElementById(id);
  if (!el) throw new Error(`#${id} missing`);
  return el as T;
};

// ------------------------------------------------------------------ storage (per viewer; may be unavailable)

function store(label: string, text: string): void {
  try {
    localStorage.setItem(storageKey(label), text);
  } catch {
    /* private mode, quota, blocked storage: the paste just isn't remembered */
  }
}
function recall(label: string): string | null {
  try {
    return localStorage.getItem(storageKey(label));
  } catch {
    return null;
  }
}

// ------------------------------------------------------------------ state

let rec: GroupRecord | null = null;
let pending: Record<string, string> = {}; // parameter values from a permalink, applied once inputs exist
let result: FinderResult | null = null;
let runId = 0;
const moyoReady = loadMoyo(new URL("./moyo_wasm_bg.wasm", import.meta.url));
moyoReady.catch(() => {});

// ------------------------------------------------------------------ elements

const elLabel = $<HTMLInputElement>("ssgf-label");
const elByu = $<HTMLFormElement>("ssgf-byu");
const elByuLabel = $<HTMLInputElement>("ssgf-byu-label");
const elByuNote = $("ssgf-byu-note");
const elPaste = $<HTMLTextAreaElement>("ssgf-paste");
const elPasteStatus = $("ssgf-paste-status");
const elGroup = $("ssgf-group");
const elPForm = $<HTMLFormElement>("ssgf-pform");
const elParams = $("ssgf-params");
const elRun = $<HTMLButtonElement>("ssgf-run");
const elRunStatus = $("ssgf-run-status");
const elPermalink = $<HTMLAnchorElement>("ssgf-permalink");
const elResults = $("ssgf-results");

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
  const out = parsePaste(elPaste.value);
  if (!out.ok) {
    rec = null;
    elGroup.hidden = true;
    elGroup.innerHTML = "";
    elPasteStatus.textContent = out.error;
    elPasteStatus.className = out.error ? "ssgf-msg err" : "ssgf-msg";
    renderParams();
    return;
  }
  rec = out.rec;
  elPasteStatus.textContent = "";
  elPasteStatus.className = "ssgf-msg";
  if (rec.label) {
    store(rec.label, elPaste.value);
    if (elLabel.value.trim() !== rec.label) {
      elLabel.value = rec.label;
      syncLabel();
    }
  }
  showGroup(rec);
}
elPaste.addEventListener("input", onPaste);

function showGroup(rec: GroupRecord): void {
  const qs = rec.q.map((q, i) => `<span class="nw">q${subDigits(i + 1)} = (${esc(q)})</span>`).join(", ");
  elGroup.innerHTML =
    `<div><dt>Group</dt><dd>${rec.label ? esc(rec.label) + " " : ""}<span class="mono">${esc(rec.sym)}</span></dd></div>` +
    `<div><dt>d</dt><dd>${rec.d}</dd></div>` +
    `<div><dt>q-vectors (basic setting)</dt><dd class="mono">${qs}</dd></div>` +
    `<div><dt>Operators</dt><dd>${rec.nops ?? "?"} (modulo lattice translations, incl. centrings)</dd></div>`;
  elGroup.hidden = false;
  renderParams();
}

/** A subgroup sent from the subgroup page (sessionStorage, read once). Its operators are
 *  in the parent's axes and were computed from the viewer's own paste there. */
function receiveSubgroup(): boolean {
  let text: string | null = null;
  try {
    text = sessionStorage.getItem(SEND_KEY);
    sessionStorage.removeItem(SEND_KEY);
  } catch {
    return false;
  }
  const sent = decodeSent(text);
  if (!sent) return false;
  rec = sent;
  elPaste.value = "";
  elPasteStatus.className = "ssgf-msg";
  elPasteStatus.textContent = "Loaded a subgroup from the subgroup page (operators in the parent group's axes). Paste a BYU page to replace it.";
  showGroup(rec);
  return true;
}

// ------------------------------------------------------------------ step 2: parameters

function currentParams(): Record<string, string> {
  const out: Record<string, string> = {};
  for (const inp of Array.from(elParams.querySelectorAll<HTMLInputElement>("input[data-param]"))) out[inp.dataset.param!] = inp.value.trim();
  return out;
}

function renderParams(): void {
  const names = rec?.params ?? [];
  const old = currentParams();
  const have = Array.from(elParams.querySelectorAll<HTMLInputElement>("input[data-param]")).map((i) => i.dataset.param);
  if (have.join(",") !== names.join(",")) {
    elParams.innerHTML = "";
    for (const p of names) {
      const id = "ssgf-p-" + p;
      const row = document.createElement("div");
      row.className = "ssgf-prow";
      const where = rec!.q.map((q, i) => (new RegExp(`\\b${p}\\b`).test(q) ? `q${subDigits(i + 1)}` : "")).filter(Boolean).join(", ");
      row.innerHTML =
        `<label for="${id}"><span class="mono">${esc(p)}</span>${where ? ` <span class="muted">(in ${where})</span>` : ""}</label>` +
        `<input id="${id}" type="text" inputmode="decimal" autocomplete="off" spellcheck="false" data-param="${esc(p)}" aria-describedby="${id}-err" placeholder="e.g. 1/4">` +
        `<span class="ssgf-perr" id="${id}-err"></span>`;
      const inp = row.querySelector("input")!;
      inp.value = pending[p] ?? old[p] ?? "";
      inp.addEventListener("input", () => validateInput(inp, false));
      elParams.appendChild(row);
      if (inp.value) validateInput(inp, false);
    }
  }
  elPForm.hidden = !rec;
  elRun.disabled = !rec;
  if (!rec) elRunStatus.textContent = "";
}

function validateInput(inp: HTMLInputElement, strict: boolean): boolean {
  const err = inp.value.trim() || strict ? validateParam(inp.value) : null;
  const msg = document.getElementById(inp.id + "-err")!;
  msg.textContent = err ?? "";
  inp.setAttribute("aria-invalid", String(!!err));
  return !err;
}

// ------------------------------------------------------------------ step 3: compute

// Yield so the spinner paints. rAF is paused in hidden tabs, so a timeout also resolves it.
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
  if (!rec) return;
  const inputs = Array.from(elParams.querySelectorAll<HTMLInputElement>("input[data-param]"));
  const bad = inputs.filter((i) => !validateInput(i, true));
  if (bad.length) {
    bad[0].focus();
    return;
  }
  const params = currentParams();
  pending = {};
  const id = ++runId;
  elRun.disabled = true;
  elRunStatus.innerHTML = '<span class="ssgf-spin" aria-hidden="true"></span> Computing…';
  elResults.setAttribute("aria-busy", "true");
  try {
    await moyoReady;
  } catch (e) {
    elRunStatus.textContent = "Could not load the symmetry library (moyo-wasm): " + (e instanceof Error ? e.message : String(e));
    elRun.disabled = false;
    elResults.removeAttribute("aria-busy");
    return;
  }
  await nextFrame(); // let the spinner paint before the synchronous computation
  if (id !== runId) return;
  const t0 = performance.now();
  const out = compute(rec, params);
  const ms = performance.now() - t0;
  elRun.disabled = false;
  elResults.removeAttribute("aria-busy");
  const hash = encodeHash(rec.label, params);
  if (rec.label) {
    try {
      history.replaceState(null, "", hash);
    } catch {
      /* ignore */
    }
    elPermalink.href = hash;
    elPermalink.hidden = false;
  }
  if (!out.ok) {
    result = null;
    elRunStatus.innerHTML = `<span class="err">${esc(out.error)}</span>`;
    elResults.innerHTML = "";
    return;
  }
  result = out.result;
  elRunStatus.textContent = `${result.rows.length} row${result.rows.length === 1 ? "" : "s"}, ${result.classes.length} class${result.classes.length === 1 ? "" : "es"} in ${ms < 1000 ? Math.round(ms) + " ms" : (ms / 1000).toFixed(1) + " s"}.`;
  renderResults(result);
}

elPForm.addEventListener("submit", (e) => {
  e.preventDefault();
  void run();
});

// ------------------------------------------------------------------ results

function renderResults(r: FinderResult): void {
  const sc = supercellText(r);
  const cols = summaryColumns(r);
  const tvar = "t₀";
  const params = Object.entries(r.params).map(([k, v]) => `${esc(k)} = ${esc(v)}`).join(", ");
  const qrows = r.Q.map((row, i) => `q${subDigits(i + 1)} = (${row.map(esc).join(", ")})`).join("; ");

  let h = "";
  if (r.warnings.length) {
    h += `<div class="ssgf-warn" role="alert"><strong>Warnings</strong><ul>${r.warnings.map((w) => `<li>${esc(w)}</li>`).join("")}</ul></div>`;
  }
  h += `<p class="ssgf-head"><span class="mono">${esc(r.symbol)}</span>${r.label ? ` (${esc(r.label)})` : ""}, ${params || "no free parameters"}: <span class="mono nw">${qrows}</span>. Supercell ${esc(sc.conv)}, index ${sc.index}.</p>`;

  // (a) IT-finder-style summary: one column per SSG orbit of t₀ classes
  const isD1 = r.d === 1;
  h += `<div class="table-scroll ssgf-sum"><table><caption class="sr-only">3D space groups by locus of ${tvar}</caption><tbody>`;
  h += `<tr><th scope="row">Row</th>${cols.map((c, i) => `<td><a href="#ssgf-row-${c.orbit}">${i + 1}</a></td>`).join("")}</tr>`;
  h += `<tr><th scope="row">Kind</th>${cols.map((c) => `<td>${esc(c.kind)}</td>`).join("")}</tr>`;
  h += `<tr><th scope="row">Locus of ${tvar}</th>${cols.map((c) => `<td class="mono lw">${esc(c.locus)}</td>`).join("")}</tr>`;
  if (isD1) h += `<tr><th scope="row">IT column</th>${cols.map((c) => `<td class="mono">${esc(c.itLabels.join(", "))}</td>`).join("")}</tr>`;
  h += `<tr><th scope="row">Space group</th>${cols.map((c) => `<td class="sym">${c.symbol.html}</td>`).join("")}</tr>`;
  h += `<tr><th scope="row">No.</th>${cols.map((c) => `<td>${esc(c.number)}</td>`).join("")}</tr>`;
  h += `</tbody></table></div>`;

  h += `<div class="ssgf-actions"><button type="button" id="ssgf-export">Export JSON</button></div>`;

  // (b) one card per row
  h += r.rows.map((row, i) => card(r, row, i)).join("");
  elResults.innerHTML = h;

  $("ssgf-export").addEventListener("click", () => {
    const blob = new Blob([resultJson(r)], { type: "application/json" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = exportName(r);
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  });
  for (const b of Array.from(elResults.querySelectorAll<HTMLButtonElement>("button[data-copy]"))) {
    b.addEventListener("click", () => copy(b, document.getElementById(b.dataset.copy!)!.textContent ?? ""));
  }
}

function card(r: FinderResult, row: OrbitRow, i: number): string {
  const c: ClassResult = row.rep;
  const sym = classSymbol(c);
  const sc = supercellText(r);
  const id = c.id;
  const dl: [string, string][] = [];
  dl.push(["Locus", `<span class="mono">${esc(row.locus)}</span>${c.itLabel && r.d === 1 ? ` <span class="muted">(IT: ${esc([...new Set(row.members.map((m) => m.itLabel))].join(", "))})</span>` : ""}`]);
  dl.push(["Example t₀", `<span class="mono">${esc(tText(c.exampleT))}</span>`]);
  if (row.members.length > 1) {
    dl.push(["Classes in this orbit", `${row.members.length}, at t₀ = <span class="mono">${row.members.map((m) => esc(tText(m.exampleT))).join("; ")}</span> (related by superspace operations, so the same crystal)`]);
  }
  dl.push(["Symbol in SSG axes", c.setting ? `<span class="sym">${symbolHtml(c.setting)}</span>` : `<span class="muted">not available for this axis choice</span>`]);
  if (id.ok) {
    dl.push(["Standard symbol", `<span class="sym">${symbolHtml(id.hmShort.replace(/\s+/g, ""))}</span> (No. ${id.number}; full <span class="sym">${symbolHtml(id.hmFull)}</span>)`]);
    dl.push(["Hall number", String(id.hallNumber)]);
  } else {
    dl.push(["Identification", `<span class="err">failed: ${esc(id.reason)}</span>`]);
  }
  dl.push(["Supercell", `conventional ${esc(sc.conv)}, index ${sc.index}` +
    (sc.centring.length ? `; centring vectors <span class="mono">${sc.centring.map(esc).join(", ")}</span>` : "; primitive")]);

  const cid = `ssgf-ops3-${row.orbit}`;
  let h = `<article class="ssgf-card" id="ssgf-row-${row.orbit}" aria-labelledby="ssgf-row-${row.orbit}-h">`;
  h += `<h3 id="ssgf-row-${row.orbit}-h"><span class="muted">${i + 1}.</span> <span class="sym">${sym.html}</span>${id.ok ? ` <span class="muted">No. ${id.number}</span>` : ""} <span class="tag">${esc(row.kind)}</span></h3>`;
  if (row.inconsistent) h += `<p class="err">The classes of this orbit gave different space-group types.</p>`;
  h += `<dl class="ssgf-dl">${dl.map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`).join("")}</dl>`;
  if (id.ok) {
    h += `<details><summary>Transformation to the standard setting</summary>` +
      `<p class="note">Fractional coordinates of the conventional supercell (x_conv) in terms of moyo's standard setting (x_std), at t₀ = ${esc(tText(c.exampleT))}.</p>` +
      `<pre>${esc(transformText(c).join("\n"))}</pre></details>`;
  }
  h += `<details><summary>Surviving superspace operations (${c.surviving.length})</summary>` +
    `<div class="table-scroll"><table><thead><tr><th scope="col">Operation (basic setting)</th><th scope="col">matched n</th></tr></thead><tbody>` +
    c.surviving.map((s) => `<tr><td class="mono">${esc(s.op)}</td><td class="mono">(${s.n.join(", ")})</td></tr>`).join("") +
    `</tbody></table></div></details>`;
  h += `<details><summary>3D operations (${c.ops3.length}, conventional supercell)</summary>` +
    `<pre id="${cid}">${esc(c.ops3.join("\n"))}</pre>` +
    `<button type="button" data-copy="${cid}">Copy operations</button> <span class="muted ssgf-copied" aria-live="polite"></span></details>`;
  h += `</article>`;
  return h;
}

async function copy(btn: HTMLButtonElement, text: string): Promise<void> {
  const note = btn.nextElementSibling as HTMLElement | null;
  let ok = false;
  try {
    await navigator.clipboard.writeText(text);
    ok = true;
  } catch {
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    try {
      ok = document.execCommand("copy");
    } catch {
      ok = false;
    }
    ta.remove();
  }
  if (note) {
    note.textContent = ok ? "Copied." : "Copy failed; select the text instead.";
    setTimeout(() => (note.textContent = ""), 2500);
  }
}

// ------------------------------------------------------------------ permalink restore

function restore(): void {
  const st = decodeHash(location.hash);
  if (!st.label) return;
  pending = st.params;
  elLabel.value = st.label;
  syncLabel();
  const saved = recall(st.label);
  if (saved) {
    elPaste.value = saved;
    onPaste();
    if (rec && rec.label === st.label && rec.params.every((p) => st.params[p])) void run();
  } else {
    elPasteStatus.className = "ssgf-msg";
    elPasteStatus.textContent = `This link is for group ${st.label}. Open its BYU page with the button above, copy the whole page and paste it here; the parameter values from the link will then be filled in.`;
  }
}

syncLabel();
renderParams();
if (!receiveSubgroup()) restore();
window.addEventListener("hashchange", () => {
  const st = decodeHash(location.hash);
  if (st.label && (!rec || rec.label !== st.label || Object.entries(st.params).some(([k, v]) => currentParams()[k] !== v))) {
    rec = null;
    elPaste.value = "";
    elParams.innerHTML = "";
    restore();
  }
});
