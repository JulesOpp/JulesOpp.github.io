// Orchestration for the subgroup tool: pasted BYU page + filters -> rows of k-index-1
// subgroup classes, each identified as a BYU type. Pure (no DOM). Requires loadMoyo() first.
//
// pointGroup notation: the HM point-group symbol of the subgroup's basic space group, obtained
// from the standardised basic symbol by dropping the centring letter, screw subscripts
// (2_1 -> 2) and mapping every glide letter (a b c d e n) to m. So P-3c1 -> "-3m1",
// P321 -> "321", C2/c -> "2/m". It is the point group of the subgroup, in the type's own
// setting (its axis order is not necessarily the parent's).

import { GroupRecord, SOp, buildOps, formatOp } from "./ssg";
import { parseByuPage } from "./paste";
import { tSubgroups } from "./subgroups";
import { isFail, standardise, type Standardised } from "./standardise";
import { ssgSymbolInAxes } from "./ssgsymbol";
import { crystalSystemOf, type CrystalSystem } from "./symbol";

export class SubgroupError extends Error {}

export interface SubFilter {
  kIndex: number;
  /** Reserved for the UI ("only k-index 1"); all rows are k-index 1, so it changes nothing. */
  onlyK?: boolean;
  pointGroup?: string;
  system?: CrystalSystem;
  /** Lattice letter of the subgroup type's basic symbol (P, A, B, C, I, F, R). */
  centring?: string;
}

export interface SubRow {
  index: number;
  tIndex: number;
  kIndex: 1;
  classSize: number;
  maximal: boolean;
  pointGroup: string;
  sgNumber: number;
  sgSymbol: string;
  /** Symbol in the parent's axes; null when the crystal system differs from the parent's or none exists. */
  inParentAxes: string | null;
  /** BYU type of the subgroup; null only if standardisation failed (a warning is then added). */
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

export function runSubgroups(raw: string | GroupRecord, f: SubFilter): { parent: string; parentSgNumber: number; rows: SubRow[]; warnings: string[] } {
  const rec = typeof raw === "string" ? parseByuPage(raw) : raw;
  if (rec.d !== 1) throw new SubgroupError("only d = 1 is supported for now");
  if (f.kIndex !== 1) throw new SubgroupError("only k-index 1 is supported for now");

  const warnings: string[] = [];
  const parentSgNumber = parentNumber(rec);
  const parentSystem = crystalSystemOf(parentSgNumber);
  const rows: SubRow[] = [];
  for (const c of tSubgroups(rec)) {
    const st = standardise(c.ops, rec.q);
    const idx = `index ${c.index}, ${c.pointOps.length} point operations`;
    if (isFail(st)) {
      warnings.push(`Subgroup (${idx}) could not be identified: ${st.reason}`);
      rows.push({
        index: c.index, tIndex: c.index, kIndex: 1, classSize: c.classSize, maximal: c.maximal,
        pointGroup: "?", sgNumber: 0, sgSymbol: "?", inParentAxes: null, type: null,
        generators: generatorsOf(rec, c.ops),
      });
      continue;
    }
    const type = st;
    const sys = crystalSystemOf(type.sgNumber);
    const inParent = sys === parentSystem ? ssgSymbolInAxes(c.ops, rec.q, parentSystem)?.full ?? null : null;
    rows.push({
      index: c.index,
      tIndex: c.index,
      kIndex: 1,
      classSize: c.classSize,
      maximal: c.maximal,
      pointGroup: pointGroupOf(type.symbol.basic),
      sgNumber: type.sgNumber,
      sgSymbol: type.symbol.basic,
      inParentAxes: inParent,
      type,
      generators: generatorsOf(rec, c.ops),
    });
  }
  return { parent: rec.sym, parentSgNumber, rows: rows.filter((r) => keep(r, f)), warnings };
}

/** Filter test. Rows whose type could not be identified always pass (they carry a warning). */
export function keep(r: SubRow, f: SubFilter): boolean {
  if (r.type === null) return true;
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
