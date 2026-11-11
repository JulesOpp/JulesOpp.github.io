// Orchestration: (group record, rational parameters) → everything the UI shows.
// Pure (no DOM). Requires loadMoyo() to have completed before runFinder().

import { Q } from "./rational";
import { GroupRecord, evalQ, formatOp } from "./ssg";
import { Exact, TClass } from "./exact";
import { identify, Identification, Op3 } from "./identify";
import { crystalSystemOf, itShortForm, settingSymbol } from "./symbol";

export const MAX_BOX = 2_000_000; // ∏ N_j: conventional-supercell volume in basic cells

export class FinderError extends Error {}

export interface ClassResult {
  dim: number;
  kind: "general" | "plane" | "line" | "point";
  locus: string;
  orbit: number;
  exampleT: string[];
  genericT: string[];
  itLabel?: string;
  /** surviving superspace operations with their matched basic translation n */
  surviving: { op: string; n: number[] }[];
  /** 3D operations of the section, conventional supercell fractions (x,y,z) */
  ops3: string[];
  id: Identification;
  setting: string | null;
  itShort: string | null;
}

export interface OrbitRow {
  orbit: number;
  dim: number;
  kind: ClassResult["kind"];
  locus: string;
  members: ClassResult[];
  /** representative (first member) */
  rep: ClassResult;
  /** true if members disagree on the space-group type (should not happen) */
  inconsistent: boolean;
}

export interface FinderResult {
  label?: string;
  symbol: string;
  d: number;
  params: Record<string, string>;
  Q: string[][];
  nOps: number;
  supercell: { conv: number[]; centring: string[][]; index: number; hnf: number[][] };
  Lambda: string[][];
  itN?: number;
  classes: ClassResult[];
  rows: OrbitRow[];
  warnings: string[];
}

export function kindOf(dim: number, d: number): ClassResult["kind"] {
  if (dim === d) return "general";
  return dim === 0 ? "point" : dim === 1 ? "line" : "plane";
}

const VAR3 = "xyz";
export function formatOp3(o: Op3): string {
  return o.R.map((row, i) => {
    let s = "";
    row.forEach((c, j) => {
      if (!c) return;
      s += (c < 0 ? "-" : s ? "+" : "") + (Math.abs(c) === 1 ? "" : String(Math.abs(c))) + VAR3[j];
    });
    const t = o.t[i].mod1();
    if (!t.isZero()) s += (s ? "+" : "") + t.toString();
    return s || "0";
  }).join(",");
}

/** Parse user parameter values; throws FinderError naming the bad field. */
export function parseParams(rec: GroupRecord, raw: Record<string, string>): Record<string, Q> {
  const out: Record<string, Q> = {};
  for (const p of rec.params) {
    const v = (raw[p] ?? "").trim();
    if (!v) throw new FinderError(`Enter a value for ${p}.`);
    try {
      out[p] = Q.parse(v);
    } catch {
      throw new FinderError(`${p} = "${v}" is not a rational number (use e.g. 1/4, -3/8 or 0.125).`);
    }
  }
  return out;
}

export function runFinder(rec: GroupRecord, rawParams: Record<string, string>): FinderResult {
  const params = parseParams(rec, rawParams);
  const warnings: string[] = [];
  // cheap pre-check of the supercell size before the exact analysis
  const Qm = evalQ(rec, params);
  const N = [0, 1, 2].map((j) => Qm.reduce((l, r) => lcmNum(l, Number(r[j].d)), 1));
  const box = N[0] * N[1] * N[2];
  if (box > MAX_BOX)
    throw new FinderError(`The conventional supercell would be ${N.join("×")} basic cells; this tool stops at ${MAX_BOX.toLocaleString()}.`);
  for (let i = 0; i < Qm.length; i++) {
    if (Qm[i].every((x) => x.isInt())) warnings.push(`q${i + 1} is a lattice vector (≡ 0 mod Z³) for these values.`);
    for (let j = i + 1; j < Qm.length; j++) {
      const same = Qm[i].every((x, k) => x.sub(Qm[j][k]).isInt());
      const opp = Qm[i].every((x, k) => x.add(Qm[j][k]).isInt());
      if (same || opp) warnings.push(`q${i + 1} ≡ ${same ? "" : "−"}q${j + 1} (mod Z³) for these values; the section is still defined, but the q's are not independent.`);
    }
  }

  let ex: Exact;
  try {
    ex = new Exact(rec, params);
  } catch (e) {
    throw new FinderError(e instanceof Error ? e.message : String(e));
  }
  const system = crystalSystemOf(Number((rec.label ?? "").split(".")[0]) || rec.sgNumber || guessSg(rec));
  const basis = [[ex.N[0], 0, 0], [0, ex.N[1], 0], [0, 0, ex.N[2]]];

  const classes: ClassResult[] = ex.classes.map((c: TClass) => {
    const id = identify(c.ops3, basis);
    let setting: string | null = null;
    try {
      setting = settingSymbol(c.ops3, system);
    } catch {
      setting = null;
    }
    return {
      dim: c.dim,
      kind: kindOf(c.dim, ex.d),
      locus: c.locus,
      orbit: c.orbit,
      exampleT: c.exampleT.map(String),
      genericT: c.genericT.map(String),
      itLabel: c.itLabel,
      surviving: c.surviving.map(([k, n]) => ({ op: formatOp(ex.ops[k]), n })),
      ops3: c.ops3.map(formatOp3),
      id,
      setting,
      itShort: setting ? itShortForm(setting) : null,
    };
  });
  const num = (c: ClassResult) => (c.id.ok ? c.id.number : 999);
  classes.sort((a, b) => b.dim - a.dim || num(a) - num(b) || (a.setting ?? "").localeCompare(b.setting ?? "") || 0);

  const rowsMap = new Map<number, ClassResult[]>();
  for (const c of classes) {
    if (!rowsMap.has(c.orbit)) rowsMap.set(c.orbit, []);
    rowsMap.get(c.orbit)!.push(c);
  }
  const rows: OrbitRow[] = [...rowsMap.entries()].map(([orbit, members]) => ({
    orbit,
    dim: members[0].dim,
    kind: members[0].kind,
    locus: members[0].locus,
    members,
    rep: members[0],
    inconsistent: new Set(members.map(num)).size > 1,
  }));
  for (const r of rows) if (r.inconsistent) warnings.push(`Classes in the orbit "${r.locus}" gave different space-group types; please report this.`);
  for (const c of classes) if (!c.id.ok) warnings.push(`Identification failed at t₀ = (${c.exampleT.join(", ")}): ${c.id.reason}`);

  return {
    label: rec.label,
    symbol: rec.sym,
    d: ex.d,
    params: Object.fromEntries(Object.entries(params).map(([k, v]) => [k, v.toString()])),
    Q: ex.Q.map((r) => r.map(String)),
    nOps: ex.ops.length,
    supercell: {
      conv: ex.N,
      centring: ex.cells.filter((n) => n.some((x) => x)).map((n) => n.map((x, i) => Q.of(BigInt(x), BigInt(ex.N[i])).toString())),
      index: ex.index,
      hnf: ex.hnfSuper.map((r) => r.map(Number)),
    },
    Lambda: ex.B.map((r) => r.map(String)),
    itN: ex.itN,
    classes,
    rows,
    warnings,
  };
}

function lcmNum(a: number, b: number): number {
  const g = (x: number, y: number): number => (y ? g(y, x % y) : x);
  return (a / g(a, b)) * b;
}

/** Fall back to the symbol's lattice letter when no label was pasted (triclinic as a last resort). */
function guessSg(_rec: GroupRecord): number {
  return 1;
}
