// Paste mode: turn the text (or HTML source) of a BYU ISO(3+d)D group page
// (ssginfoform.php, "Detailed information") into a GroupRecord (CONTRACTS.md §2).
// Mirrors reference/ssg_descent.py::parse_ssg, which reads the BASIC setting block.

import { GroupRecord, paramsOf, parseOp, formatOp, buildOps } from "./ssg";
import { Q } from "./rational";

export class PasteError extends Error {}

function normalise(raw: string): string {
  let t = raw;
  if (/<[a-zA-Z][^>]*>/.test(t)) {
    t = t.replace(/<[^>]+>/g, "\n");
    t = t.replace(/&nbsp;/g, " ").replace(/&amp;/g, "&").replace(/&lt;/g, "<").replace(/&gt;/g, ">")
      .replace(/&#(\d+);/g, (_, n) => String.fromCharCode(Number(n)));
  }
  t = t.replace(/[   ]/g, " ").replace(/[−–]/g, "-");
  // alternative coordinate notations x1..x6 / xs1..xs6 -> x,y,z,t,u,v
  t = t.replace(/\bxs?([1-6])\b/g, (_, k) => "xyztuv"[Number(k) - 1]);
  return t.replace(/\s+/g, " ");
}

export function parseByuPage(raw: string): GroupRecord {
  const t = normalise(raw);
  const head = /Superspace group:\s*(\d+\.\d+\.\d+\.\d+)\s+(\S+)/.exec(t);
  const i0 = t.indexOf("Independent modulation vectors:");
  if (i0 < 0) throw new PasteError('Could not find "Independent modulation vectors:". Paste the whole BYU group page (basic space group setting).');
  let i1 = t.indexOf("SUPERCENTERED SETTING", i0);
  if (i1 < 0) i1 = t.indexOf("Reflection conditions:", i0);
  if (i1 < 0) i1 = t.length;
  const blk = t.slice(i0, i1);
  const q = [...blk.matchAll(/q\d+\s*=\s*\(([^)]*)\)/g)].map((m) => m[1].replace(/\s+/g, ""));
  const d = q.length;
  if (d < 1 || d > 3) throw new PasteError(`Found ${d} modulation vectors; expected 1–3.`);
  const nv = 3 + d;
  const cenM = /Centering:\s*(.*?)\s*Non-lattice generators/.exec(blk);
  if (!cenM) throw new PasteError('Could not find the "Centering:" line of the basic setting.');
  const cenAll = [...cenM[1].matchAll(/\(([^)]*)\)/g)].map((m) => m[1].replace(/\s+/g, ""));
  const opsM = /Non-lattice operators:\s*(.*)$/.exec(blk);
  if (!opsM) throw new PasteError('Could not find "Non-lattice operators:".');
  const opStrs = [...opsM[1].matchAll(/\(([^)]*)\)/g)].map((m) => m[1].replace(/\s+/g, ""));
  if (!opStrs.length) throw new PasteError("No operators found.");

  const gen: string[] = [];
  for (const s of opStrs) {
    const op = parseOp(s, nv);
    const norm = { W: op.W, w: op.w.map((x) => x.mod1()) };
    const isId = norm.W.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0))) && norm.w.every((x) => x.isZero());
    if (!isId) gen.push(formatOp(norm));
  }
  const cen: string[] = [];
  for (const c of cenAll) {
    const v = c.split(",").map((x) => Q.parse(x).mod1());
    if (v.length !== nv) throw new PasteError(`Centring vector (${c}) has the wrong length.`);
    if (v.some((x) => !x.isZero())) cen.push(v.map(String).join(","));
  }
  const label = head?.[1];
  const m = label ? /^\d+\.(\d+\.\d+)\.\d+$/.exec(label) : null;
  const rec: GroupRecord = { label, sym: head?.[2] ?? "?", cls: m ? m[1] : "", d, q, params: paramsOf(q), cen, gen };
  rec.nops = buildOps(rec).length;
  return rec;
}
