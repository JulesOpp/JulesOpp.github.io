import { beforeAll, describe, expect, it, vi } from "vitest";
import { existsSync, readFileSync, readdirSync } from "node:fs";
import { join } from "node:path";
import { Q } from "../src/rational";
import { buildOps, type GroupRecord } from "../src/ssg";
import { affineText, prepareSend, subgroupRecord } from "../src/ui/subflow";
import { loadMoyo } from "../src/moyo";
import * as subfinder from "../src/subfinder";
import * as standardiser from "../src/standardise";
import { keep, runSubgroups, SubgroupError, type SubRow } from "../src/subfinder";

beforeAll(async () => {
  await loadMoyo();
});

// Hand-written BYU-style page text (no BYU data): P-3c1(0,0,g)000, IT A No. 165.
const page = (q: string[], ops: string[]) =>
  `Superspace group: 165.1.24.1 P-3c1(0,0,g)000 Independent modulation vectors: ` +
  q.map((v, i) => `q${i + 1} = (${v})`).join(" ") +
  ` Centering: (${"0,".repeat(q.length + 2)}0) Non-lattice generators: Non-lattice operators: ` +
  ops.map((o) => `(${o})`).join(" ");
const P3C1 = page(["0,0,g"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);

describe("runSubgroups (P-3c1(0,0,g)000)", () => {
  let all: ReturnType<typeof runSubgroups>;
  beforeAll(() => {
    all = runSubgroups(P3C1, { kIndex: 1 });
  });
  it("lists 10 classes", () => {
    expect(all.rows).toHaveLength(10);
    expect(all.rows.map((r) => r.sgNumber).sort((a, b) => a - b)).toEqual([1, 2, 5, 9, 15, 143, 147, 150, 158, 165]);
  });
  it("row for No. 150 is P321(0,0,g)000 with index 2", () => {
    const r = all.rows.find((x) => x.sgNumber === 150)!;
    expect(r.type!.symbol.full).toBe("P321(0,0,g)000");
    expect(r.index).toBe(2);
    expect(r.tIndex).toBe(2);
    expect(r.kIndex).toBe(1);
    expect(r.pointGroup).toBe("321");
    expect(r.inParentAxes).not.toBeNull();
    expect(r.generators.length).toBeGreaterThan(0);
  });
  it("system filter trigonal gives 5 rows", () => {
    expect(runSubgroups(P3C1, { kIndex: 1, system: "trigonal" }).rows).toHaveLength(5);
  });
  it("pointGroup filter 321 gives 1 row", () => {
    expect(runSubgroups(P3C1, { kIndex: 1, pointGroup: "321" }).rows).toHaveLength(1);
  });
  it("centring filter", () => {
    const rows = runSubgroups(P3C1, { kIndex: 1, centring: "B" }).rows;
    expect(rows).toHaveLength(3);
    expect(rows.every((r) => r.type!.symbol.basic[0] === "B")).toBe(true);
  });
  it("monoclinic subgroup drops out of the parent axes", () => {
    expect(all.rows.find((x) => x.sgNumber === 15)!.inParentAxes).toBeNull();
  });
  it("refuses values beyond the dimensional cap", () => {
    expect(() => runSubgroups(P3C1, { kIndex: 4 })).toThrow(SubgroupError);
  });
});

it("caches the parent once for repeated canonical operators and q text", () => {
  const raw = page(["0,0,g+99"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);
  const spy = vi.spyOn(standardiser, "standardise");
  const hits = vi.spyOn(subfinder, "_cacheHits");
  const before = subfinder._cacheHits();
  const a = runSubgroups(raw, { kIndex: 1 });
  const first = spy.mock.calls.length;
  runSubgroups(raw, { kIndex: 1 });
  expect(first).toBe(a.rows.length + 1);
  expect(spy.mock.calls.length - first).toBe(a.rows.length);
  expect(subfinder._cacheHits()).toBe(before + 1);
  expect(hits).toHaveBeenCalledTimes(2);
  spy.mockRestore(); hits.mockRestore();
});

it("parent cache uses canonical operations and q text, and evicts beyond eight entries", () => {
  const original = page(["0,0,g+101"], ["x,y,z,t", "y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"]);
  const reordered = page(["0,0,g+101"], ["x-y,-y,-z+1/2,-t", "y,-x+y,-z,-t", "x,y,z,t"]);
  runSubgroups(original, { kIndex: 1 });
  const before = subfinder._cacheHits();
  runSubgroups(reordered, { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(before + 1);
  // Identity-only groups keep this capacity check cheap; the q text must be part
  // of the key even though every one has the same basic space-group number.
  const raws = Array.from({ length: 9 }, (_, i) => page([`0,0,g+${i+201}`], ["x,y,z,t"]));
  const start = subfinder._cacheHits();
  for (const raw of raws) runSubgroups(raw, { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start);
  runSubgroups(raws[8], { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start + 1);
  runSubgroups(raws[0], { kIndex: 1 });
  expect(subfinder._cacheHits()).toBe(start + 1);
});

describe("keep", () => {
  it("rows with an unidentified type pass every filter", () => {
    const r: SubRow = {
      index: 2, tIndex: 2, kIndex: 1, classSize: 1, maximal: true, pointGroup: "?", sgNumber: 0, sgSymbol: "?",
      inParentAxes: null, type: null, generators: [],
    };
    expect(keep(r, { kIndex: 1, pointGroup: "321", system: "cubic", centring: "F" })).toBe(true);
  });
});

// Local-only: 98.1.21.1 I4_122(0,0,g)000 has two classes with the same BYU type B2(a,b,0)0
// (checked in local/fixtures/subgroups); they must stay two rows with equal type symbols.
const REC = join(__dirname, "..", "local", "data", "sg098.json");
describe.skipIf(!existsSync(REC))("duplicate types (local)", () => {
  it("98.1.21.1 keeps two rows with equal type.symbol.full", () => {
    const g = JSON.parse(readFileSync(REC, "utf8")).groups["98.1.21.1"];
    const cen = g.cen.map((c: string) => `(${c})`).join(" ");
    const raw =
      `Superspace group: 98.1.21.1 ${g.sym} Independent modulation vectors: q1 = (${g.q[0]}) ` +
      `Centering: (0,0,0,0) ${cen} Non-lattice generators: Non-lattice operators: ` +
      ["x,y,z,t", ...g.gen].map((o: string) => `(${o})`).join(" ");
    const rows = runSubgroups(raw, { kIndex: 1 }).rows;
    const by = new Map<string, number>();
    for (const r of rows) if (r.type) by.set(r.type.symbol.full, (by.get(r.type.symbol.full) ?? 0) + 1);
    expect(Math.max(...by.values())).toBeGreaterThanOrEqual(2);
  });
});

// Hand-built operations and q vectors; no catalogue record or expected catalogue name.
const HEX_D2: GroupRecord = {
  sym: "hand P6/mmm d=2", cls: "hexagonal", d: 2,
  q: ["0,0,g1", "0,0,g2"], params: ["g1", "g2"], cen: [],
  gen: ["x-y,x,z,t,u", "y,x,z,t,u", "x,y,-z,-t,-u"],
};

it("names the hand P6/mmm d=2 classes invariantly in a random setting", () => {
  const progress: number[] = [];
  const out = runSubgroups(HEX_D2, { kIndex: 1 }, (done, total) => {
    expect(total).toBeGreaterThan(0);
    progress.push(done);
  });
  expect(out.rows.length).toBeGreaterThan(0);
  expect(out.warnings).toEqual([]);
  expect(progress).toEqual(Array.from({ length: out.rows.length + 1 }, (_, i) => i));
  const row = out.rows.find(r => r.sgNumber === 177);
  expect(row).toBeDefined();
  expect(row!.type).not.toBeNull();
  const rec = subgroupRecord(HEX_D2, row!, out.parentSgNumber);
  const sign = () => Math.random() < 0.5 ? -1 : 1;
  const swap = Math.random() < 0.5;
  const P = (swap ? [[0,sign(),0],[sign(),0,0],[0,0,1]] : [[sign(),0,0],[0,sign(),0],[0,0,1]]).map(r=>r.map(x=>Q.of(x)));
  const p = Array.from({length:3},()=>Q.parse(`${Math.floor(Math.random()*5)+1}/7`));
  const ops = standardiser.transformOps(buildOps(rec), P, p, [[1,0],[0,1]], Array.from({length:2},()=>[Q.of(0),Q.of(0),Q.of(0)]), [Q.of(0),Q.of(0)]);
  expect(ops).not.toBeNull();
  const named = standardiser.standardise(ops!, rec.q);
  expect(standardiser.isFail(named)).toBe(false);
  if (standardiser.isFail(named)) throw Error(named.reason);
  expect(row!.type!.symbol.full).toBe(named.symbol.full);
});

const CUBIC_DATA = join(__dirname, "..", "local", "data");
describe.skipIf(!existsSync(join(CUBIC_DATA, "sg230.json")))("cubic d=3 timing (local)", () => {
  // Exact-output naming optimisation: 53.171 s for all 33 classes of 225.3.209.1.
  it("names the largest stored cubic d=3 parent within 120 s", () => {
    const parents: GroupRecord[] = [];
    for (const file of readdirSync(CUBIC_DATA).filter(f=>/^sg\d+\.json$/.test(f) && Number(f.slice(2,5)) >= 195)) {
      const groups = JSON.parse(readFileSync(join(CUBIC_DATA,file),"utf8")).groups as Record<string,GroupRecord>;
      parents.push(...Object.values(groups).filter(r=>r.d===3));
    }
    const parent = parents.sort((a,b)=>(b.nops??0)-(a.nops??0))[0];
    expect(parent).toBeDefined();
    const start = performance.now();
    const out = runSubgroups(parent, {kIndex:1});
    expect(out.rows.length).toBeGreaterThan(0);
    expect(out.rows.every(r=>r.type!==null)).toBe(true);
    expect(out.warnings).toEqual([]);
    expect((performance.now()-start)/1000).toBeLessThanOrEqual(120);
  }, 600_000);
});

const P1: GroupRecord = {sym:"hand P1",cls:"triclinic",d:1,q:["a,b,g"],params:["a","b","g"],cen:[],gen:[]};
it("includes smaller k by default and filters exact k and lattice kinds",()=>{
 const all=runSubgroups(P1,{kIndex:2});
 expect(all.rows).toHaveLength(16);
 expect(runSubgroups(P1,{kIndex:2,onlyK:true}).rows).toHaveLength(15);
 expect(runSubgroups(P1,{kIndex:3,onlyK:true}).rows).toHaveLength(40);
 expect(all.rows.every(r=>r.type!==null)).toBe(true);
 expect(all.rows.filter(r=>r.kIndex===2).every(r=>r.rebased?.cen.length===0)).toBe(true);
 expect(all.rows.every(r=>r.lattice?.basis)).toBe(true);
 const physical=runSubgroups(P1,{kIndex:2,onlyK:true,kinds:["3d"]});
 expect(physical.rows.length).toBeGreaterThan(0);
 expect(physical.rows.every(r=>r.lattice?.kind==="3d")).toBe(true);
});
it("retains every unidentified class after returned failures or thrown naming errors",()=>{
 runSubgroups(P1,{kIndex:1}); // cache the parent before simulating class failures
 const spy=vi.spyOn(standardiser,"standardise");
 let calls=0;
 spy.mockImplementation(()=>{if(calls++%2)throw new Error("bounded search failure");return {ok:false,reason:"no setting"};});
 try{
  const out=runSubgroups(P1,{kIndex:2,onlyK:true,pointGroup:"321",centring:"F",kinds:["mixed"]});
  expect(out.rows).toHaveLength(16);
  expect(out.rows.every(r=>r.type===null)).toBe(true);
  expect(out.warnings).toHaveLength(16);
  expect(out.warnings.every(w=>w.includes("unidentified"))).toBe(true);
 }finally{spy.mockRestore();}
});
it("enforces d=2 and d=3 caps",()=>{
 expect(()=>runSubgroups({...P1,d:2,q:["a,b,g","a2,b2,g2"]},{kIndex:3})).toThrow("d=2 supports k-index 1 to 2");
 expect(()=>runSubgroups({...P1,d:3,q:["a,b,g","a2,b2,g2","a3,b3,g3"]},{kIndex:2})).toThrow("d=3 supports k-index 1 to 1");
});

it("composes a fractional internal rebase and sends the subgroup lattice to the finder",()=>{
 const out=runSubgroups(P1,{kIndex:2,onlyK:true});
 const internal=out.rows.find(r=>r.lattice?.basis==="(a, b, c; 2t)")!;
 expect(internal).toBeDefined();
 expect(affineText(internal.transformation!)).toBe("(x, y, z; 1/2t)");
 const record=subgroupRecord(P1,internal,out.parentSgNumber);
 expect(record.q).toEqual(["1/2a,1/2b,1/2g"]);
 expect(record.cen).toEqual([]);
 expect(buildOps(record)).toHaveLength(1);
 expect(prepareSend(P1,{...internal,rebased:undefined,type:null},out.parentSgNumber).ok).toBe(false);
});
