// UI logic without the DOM: paste parsing, parameter validation, permalinks, result
// formatting and JSON export. main.ts wires these to the page; the Node check calls
// the same functions, so the code path is shared.
//
// Licence rule (D1): nothing here contains BYU data. The permalink holds the label and
// parameter values only; operators come from the viewer's own paste.

import { Q } from "../rational";
import { GroupRecord } from "../ssg";
import { parseByuPage, PasteError } from "../paste";
import { runFinder, FinderResult, ClassResult, OrbitRow, FinderError } from "../finder";
import { itShortForm } from "../symbol";
import { numToQ } from "../identify";

export const LABEL_RE = /^\d+\.\d+\.\d+\.\d+$/;
export const BYU_FORM_URL = "https://iso.byu.edu/ssginfoform.php";

// ------------------------------------------------------------------ paste

export type PasteOutcome = { ok: true; rec: GroupRecord } | { ok: false; error: string };

export function parsePaste(text: string): PasteOutcome {
  if (!text.trim()) return { ok: false, error: "" };
  try {
    return { ok: true, rec: parseByuPage(text) };
  } catch (e) {
    if (e instanceof PasteError) return { ok: false, error: e.message };
    return { ok: false, error: "Could not read this text as a BYU group page: " + (e instanceof Error ? e.message : String(e)) };
  }
}

// ------------------------------------------------------------------ parameters

/** Inline validation of one parameter value; null when valid. */
export function validateParam(raw: string): string | null {
  const v = raw.trim();
  if (!v) return "Enter a rational value, e.g. 1/4.";
  try {
    Q.parse(v);
    return null;
  } catch {
    return "Not a rational number. Use e.g. 1/4, -3/8 or 0.125.";
  }
}

// ------------------------------------------------------------------ permalink (#ssg=<label>&a1=1/4&...)

export interface HashState {
  label: string | null;
  params: Record<string, string>;
}

export function encodeHash(label: string | undefined, params: Record<string, string>): string {
  const parts: string[] = [];
  if (label) parts.push("ssg=" + encodeURIComponent(label));
  for (const [k, v] of Object.entries(params)) parts.push(encodeURIComponent(k) + "=" + encodeURIComponent(v.trim()).replace(/%2F/gi, "/"));
  return "#" + parts.join("&");
}

export function decodeHash(hash: string): HashState {
  const out: HashState = { label: null, params: {} };
  const h = hash.replace(/^#/, "");
  if (!h) return out;
  for (const part of h.split("&")) {
    const i = part.indexOf("=");
    if (i <= 0) continue;
    let k: string, v: string;
    try {
      k = decodeURIComponent(part.slice(0, i));
      v = decodeURIComponent(part.slice(i + 1));
    } catch {
      continue;
    }
    if (k === "ssg") {
      if (LABEL_RE.test(v)) out.label = v;
    } else if (/^[a-z]\d?$/i.test(k) && v.length < 40) {
      out.params[k] = v;
    }
  }
  return out;
}

// ------------------------------------------------------------------ compute

export type ComputeOutcome = { ok: true; result: FinderResult } | { ok: false; error: string };

/** runFinder with errors turned into messages. loadMoyo() must have completed. */
export function compute(rec: GroupRecord, params: Record<string, string>): ComputeOutcome {
  try {
    return { ok: true, result: runFinder(rec, params) };
  } catch (e) {
    if (e instanceof FinderError) return { ok: false, error: e.message };
    return { ok: false, error: "Internal error: " + (e instanceof Error ? e.message : String(e)) };
  }
}

// ------------------------------------------------------------------ symbols

export function esc(s: string): string {
  return s.replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]!);
}

/** HTML for an HM symbol: "-3" → 3 with overbar, "2_1" (or "2§1") → 2 with subscript 1. */
export function symbolHtml(s: string): string {
  return esc(s)
    .replace(/-(\d)/g, '<span class="ovl">$1</span>')
    .replace(/(\d)[_§](\d)/g, "$1<sub>$2</sub>");
}

/** IT-finder short symbol of a class, as HTML (screw subscripts kept), with its plain-text form. */
export function classSymbol(c: ClassResult): { text: string; html: string } {
  if (c.setting && c.itShort) {
    let html: string;
    try {
      html = symbolHtml(itShortForm(c.setting.replace(/(\d)_(\d)/g, "$1§$2")));
    } catch {
      html = esc(c.itShort);
    }
    return { text: c.itShort, html };
  }
  if (c.id.ok) {
    const t = c.id.hmShort.replace(/\s+/g, "");
    return { text: t, html: symbolHtml(t) };
  }
  return { text: "?", html: "?" };
}

// ------------------------------------------------------------------ summary table

export interface SummaryColumn {
  orbit: number;
  kind: string;
  locus: string;
  itLabels: string[];
  symbol: { text: string; html: string };
  number: string;
  nClasses: number;
}

export function summaryColumns(r: FinderResult): SummaryColumn[] {
  return r.rows.map((row: OrbitRow) => ({
    orbit: row.orbit,
    kind: row.kind,
    locus: row.locus,
    itLabels: [...new Set(row.members.map((m) => m.itLabel).filter((x): x is string => !!x))],
    symbol: classSymbol(row.rep),
    number: row.rep.id.ok ? String(row.rep.id.number) : "?",
    nClasses: row.members.length,
  }));
}

// ------------------------------------------------------------------ formatting helpers

const SUB = "₀₁₂₃₄₅₆₇₈₉";
export const subDigits = (n: number): string => String(n).replace(/\d/g, (c) => SUB[Number(c)]);

export function tText(t: string[]): string {
  return t.length === 1 ? t[0] : "(" + t.join(", ") + ")";
}

export function fracText(x: number): string {
  try {
    return numToQ(x).toString();
  } catch {
    return String(x);
  }
}

/** Transformation conventional supercell → moyo standard setting, as lines of text. */
export function transformText(c: ClassResult): string[] {
  if (!c.id.ok) return [];
  const P = c.id.stdLinear;
  const p = c.id.stdOriginShift.map(fracText);
  const w = Math.max(...P.flat().map((x) => fracText(x).length));
  const pw = Math.max(...p.map((x) => x.length));
  return P.map((row, i) => {
    const l = i === 1 ? "x_conv = " : "         ";
    const op = i === 1 ? " · x_std + " : " ".repeat(11);
    return `${l}( ${row.map((x) => fracText(x).padStart(w)).join("  ")} )${op}( ${p[i].padStart(pw)} )`;
  });
}

export function supercellText(r: FinderResult): { conv: string; centring: string[]; index: number } {
  const [a, b, c] = r.supercell.conv;
  return {
    conv: `${a}a × ${b}b × ${c}c`,
    centring: r.supercell.centring.map((v) => "(" + v.join(", ") + ")"),
    index: r.supercell.index,
  };
}

// ------------------------------------------------------------------ export

export function resultJson(r: FinderResult): string {
  return JSON.stringify(
    { tool: "Superspace Group Finder", generated: new Date().toISOString(), ...r },
    (_k, v) => (v instanceof Q ? v.toString() : typeof v === "bigint" ? v.toString() : v),
    1,
  );
}

export function exportName(r: FinderResult): string {
  const p = Object.entries(r.params).map(([k, v]) => `${k}=${v.replace(/\//g, "_")}`).join(",");
  return `ssg-finder_${r.label ?? "group"}${p ? "__" + p : ""}.json`;
}

// ------------------------------------------------------------------ storage keys

export const storageKey = (label: string) => "ssgf:paste:" + label;
