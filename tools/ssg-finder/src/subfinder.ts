// Orchestration for the subgroup tool: pasted BYU page + filters -> rows of finite-index
// subgroup classes, each named in this tool's deterministic canonical setting. Pure (no DOM). Requires loadMoyo() first.
//
// pointGroup notation: the HM point-group symbol of the subgroup's basic space group, obtained
// from the standardised basic symbol by dropping the centring letter, screw subscripts
// (2_1 -> 2) and mapping every glide letter (a b c d e n) to m. So P-3c1 -> "-3m1",
// P321 -> "321", C2/c -> "2/m". It is the point group of the subgroup, in the type's own
// setting (its axis order is not necessarily the parent's).

import { GroupRecord, SOp, buildOps, formatOp } from "./ssg";
import { parseByuPage } from "./paste";
import { kSubgroups, rebase, type LatticeKind } from "./ksubgroups";
import { Q, matMul, type QMat } from "./rational";
import { isFail, standardise, type Standardised } from "./standardise";
import { ssgSymbolInAxes } from "./ssgsymbol";
import { crystalSystemOf, type CrystalSystem } from "./symbol";

export class SubgroupError extends Error {}

export interface AffineMap { P: string[][]; p: string[]; E: (string | number)[][]; M: string[][]; p4: string[] }
export const maxK = (d: number): number => d === 1 ? 3 : d === 2 ? 2 : d === 3 ? 1 : 0;

/** Compose type coordinates with the subgroup's primitive rebase. */
function composeMap(type: Standardised, S: QMat): AffineMap {
  const d = type.E.length;
  const T = [...type.P.map(row => [...row.map(Q.parse), ...Array.from({length:d},()=>Q.ZERO)]),
    ...type.M.map((row,i)=>[...row.map(Q.parse), ...type.E[i].map(x=>Q.of(x))])];
  const C = matMul(T,S);
  return {P:C.slice(0,3).map(row=>row.slice(0,3).map(String)),p:type.p,
    E:C.slice(3).map(row=>row.slice(3).map(String)),M:C.slice(3).map(row=>row.slice(0,3).map(String)),p4:type.p4};
}
function basisText(B: QMat): string {
  const names = ["a","b","c","t","u","v"];
  const vectors = B[0].map((_,j)=>B.map((row,i)=>{
    const x=row[j]; if(x.isZero())return "";
    const a=x.sign()<0?x.neg():x;
    return (x.sign()<0?"-":"+")+(a.eq(Q.ONE)?"":a.toString())+names[i];
  }).join("").replace(/^\+/,"") || "0");
  return `(${vectors.slice(0,3).join(", ")}; ${vectors.slice(3).join(", ")})`;
}

export interface SubFilter {
  kIndex: number;
  /** When false, include every k up to kIndex. */
  onlyK?: boolean;
  kinds?: LatticeKind[];
  pointGroup?: string;
  system?: CrystalSystem;
  /** Lattice letter of the subgroup type's basic symbol (P, A, B, C, I, F, R). */
  centring?: string;
}

export interface SubRow {
  index: number;
  tIndex: number;
  kIndex: number;
  lattice?: { kind: LatticeKind; basis: string };
  rebased?: GroupRecord;
  transformation?: AffineMap;
  classSize: number;
  maximal: boolean;
  pointGroup: string;
  sgNumber: number;
  sgSymbol: string;
  /** Symbol in the parent's axes; null when the crystal system differs from the parent's or none exists. */
  inParentAxes: string | null;
  /** Canonical type of the subgroup; null only if standardisation failed (a warning is then added). */
  type: Standardised | null; // when null, pointGroup "?", sgNumber 0, sgSymbol "?"
  /** xyzt operators of H in the parent axes, centrings excluded. */
  generators: string[];
}

export function pointGroupOf(basic: string): string {
  return basic
    .replace(/^[A-Z]/, "")
    .replace(/_\d/g, "")
    .replace(/[abcden]/g, "m")
    .replace(/\s+/g, "");
}

const isId = (W: number[][]) => W.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0)));
const key = (o: SOp) => JSON.stringify(o.W) + o.w.map((x) => x.mod1().toString()).join(",");

/** A small generating set of H (greedy), as operator strings; centring translations are excluded. */
function generatorsOf(rec: GroupRecord, ops: SOp[]): string[] {
  const gens: string[] = [];
  let have = new Set(buildOps({ ...rec, gen: [] }).map(key));
  for (const o of ops) {
    if (isId(o.W) || have.has(key(o))) continue;
    gens.push(formatOp(o));
    have = new Set(buildOps({ ...rec, gen: gens }).map(key));
  }
  return gens;
}

export interface SubResult {
  parent: string;
  parentSgNumber: number;
  rows: SubRow[];
  warnings: string[];
}

/** The third argument is the current filtered snapshot, for incremental displays. */
export type SubProgress = (done: number, total: number, partial: SubResult) => void;

export function runSubgroups(raw: string | GroupRecord, f: SubFilter, onProgress?: SubProgress): SubResult {
  let result!: SubResult;
  for (const step of subgroupSteps(raw, f)) {
    result = step.result;
    onProgress?.(step.done, step.total, result);
  }
  return result;
}

/** Shared by the synchronous API and the yielding main-thread fallback. */
export function* subgroupSteps(raw: string | GroupRecord, f: SubFilter): Generator<{ done: number; total: number; result: SubResult }> {
  const rec = typeof raw === "string" ? parseByuPage(raw) : raw;
  if (!Number.isInteger(f.kIndex) || f.kIndex < 1 || f.kIndex > maxK(rec.d))
    throw new SubgroupError(`d=${rec.d} supports k-index 1 to ${maxK(rec.d)}`);

  const warnings: string[] = [];
  const parentSgNumber = parentNumber(rec);
  const parentSystem = crystalSystemOf(parentSgNumber);
  const rows: SubRow[] = [];
  // Enumerate all k <= the selected cap so optional filters can change without recomputing.
  const classes = Array.from({length:f.kIndex},(_,i)=>kSubgroups(rec,i+1)).flat();
  let done = 0;
  const snapshot = () => ({ done, total: classes.length, result: { parent: rec.sym, parentSgNumber, rows: rows.filter(r => keep(r, f)), warnings: [...warnings] } });
  yield snapshot();
  for (const c of classes) {
    const idx = `index ${c.index}, k-index ${c.k}, class ${done+1}`;
    const row: SubRow = {
      index:c.index,tIndex:c.tIndex,kIndex:c.k,classSize:c.classSize,maximal:c.maximal,
      pointGroup:"?",sgNumber:0,sgSymbol:"?",inParentAxes:null,type:null,
      generators:c.k===1?generatorsOf(rec,c.ops):c.ops.map(formatOp),
      lattice:{kind:c.kind,basis:"unavailable"},
    };
    try {
      const rb = rebase(rec,c.ops,c.sigma);
      row.lattice = {kind:c.kind,basis:basisText(rb.B)};
      // Keep k=1 input unchanged, including conventional centring and axis conventions.
      const ops = c.k===1 ? c.ops : rb.ops;
      const q = c.k===1 ? rec.q : rb.record.q;
      if(c.k>1) row.rebased = rb.record;
      const st = standardise(ops,q);
      if(isFail(st)) warnings.push(`Subgroup (${idx}) is unidentified: ${st.reason}`);
      else {
        row.type=st;row.sgNumber=st.sgNumber;row.sgSymbol=st.symbol.basic;row.pointGroup=pointGroupOf(st.symbol.basic);
        row.transformation = c.k===1 ? st : composeMap(st,rb.S);
        if(c.k>1) row.rebased!.sgNumber=st.sgNumber;
        if(c.k===1 && crystalSystemOf(st.sgNumber)===parentSystem) {
          try { row.inParentAxes=ssgSymbolInAxes(c.ops,rec.q,parentSystem)?.full??null; }
          catch { /* A missing parent-axes symbol does not invalidate the identified type. */ }
        }
      }
    } catch(e) {
      row.type=null; row.transformation=undefined; row.sgNumber=0; row.sgSymbol="?"; row.pointGroup="?";
      warnings.push(`Subgroup (${idx}) is unidentified: ${e instanceof Error?e.message:String(e)}`);
    }
    rows.push(row);
    ++done;
    yield snapshot();
  }
}

/** Filter test. Rows whose type could not be identified always pass (they carry a warning). */
export function keep(r: SubRow, f: SubFilter): boolean {
  if (r.type === null) return true;
  if (r.kIndex > f.kIndex || (f.onlyK && r.kIndex !== f.kIndex)) return false;
  if (f.kinds && r.lattice && !f.kinds.includes(r.lattice.kind)) return false;
  if (f.pointGroup && r.pointGroup !== f.pointGroup) return false;
  if (f.system && crystalSystemOf(r.sgNumber) !== f.system) return false;
  if (f.centring && r.type?.symbol.basic[0] !== f.centring) return false;
  return true;
}

const parentCache = new Map<string, number>();
let cacheHits = 0;
/** Test-only cumulative parent cache hit counter. */
export function _cacheHits(): number { return cacheHits; }

/** Space-group number of the parent's basic group, from its identity-section standardisation. */
function parentNumber(rec: GroupRecord): number {
  const ops = buildOps(rec);
  const cacheKey = JSON.stringify([ops.map(formatOp).sort(), rec.q]);
  const hit = parentCache.get(cacheKey);
  if (hit !== undefined) { cacheHits++; return hit; }
  const st = standardise(ops, rec.q);
  if (isFail(st)) throw new SubgroupError(`could not identify the parent group: ${st.reason}`);
  if (parentCache.size === 8) parentCache.delete(parentCache.keys().next().value!);
  parentCache.set(cacheKey, st.sgNumber);
  return st.sgNumber;
}
