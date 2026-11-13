import { Q } from "../src/rational";
import { describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { buildOps, GroupRecord } from "../src/ssg";
import { crystalSystemOf } from "../src/symbol";
import { ssgSymbolInAxes, _stats } from "../src/ssgsymbol";

// P-3c1(0,0,g)000, IT A No. 165 generators in hexagonal axes with the 4th coordinate.
const P3c1: GroupRecord = {
  sym: "P-3c1(0,0,g)000", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
  gen: ["y,-x+y,-z,-t", "x-y,-y,-z+1/2,-t"],
};

describe("ssgSymbolInAxes (public)", () => {
  it("hand-written P-3c1(0,0,g)000", () => {
    const s = ssgSymbolInAxes(buildOps(P3c1), P3c1.q, "trigonal");
    expect(s).toEqual({ basic: "P-3c1", q: "(0,0,g)", shifts: ["000"], full: "P-3c1(0,0,g)000" });
  });
  it("phase letters: P6_2 with t' = t + 1/6 gives h", () => {
    const r: GroupRecord = { sym: "", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [], gen: ["x-y,x,z+1/3,t+1/6"] };
    expect(ssgSymbolInAxes(buildOps(r), r.q, "hexagonal")?.full).toBe("P6_2(0,0,g)h");
  });
  it("keeps the minimum letter when nearest copies disagree (synthetic)", () => {
    // Hand axial IT A P4/mmm generators, plus an internal half translation. This
    // deliberately supplies coincident physical copies carrying both 0 and s.
    const rec = { ...P3c1, gen: ["-y,x,z,t", "x,-y,z,t", "x,y,-z,-t", "x,y,z,t+1/2"] };
    const before = _stats.disagree;
    expect(ssgSymbolInAxes(buildOps(rec), rec.q, "tetragonal")?.full).toBe("P4/mmm(0,0,g)0000");
    expect(_stats.disagree).toBeGreaterThan(before);
  });
  it("writes one phase string per q in hand P6/mmm", () => {
    const rec: GroupRecord = {
      sym: "", cls: "", d: 2, q: ["0,0,g1", "0,0,g2"],
      params: ["g1", "g2"], cen: [],
      gen: ["x-y,x,z,t,u", "x,x-y,z,t,u", "x,y,-z,-t,-u"],
    };
    expect(ssgSymbolInAxes(buildOps(rec), rec.q, "hexagonal")).toEqual({
      basic: "P6/mmm", q: "(0,0,g1)(0,0,g2)", shifts: ["0000", "0000"],
      full: "P6/mmm(0,0,g1)0000(0,0,g2)0000",
    });
  });
  it("preserves the signed sixfold phase in hand P622", () => {
    const rec: GroupRecord = {
      sym: "", cls: "", d: 2, q: ["0,0,g1", "0,0,g2"],
      params: ["g1", "g2"], cen: [],
      gen: ["x-y,x,z,t-1/6,u", "-x,-x+y,-z,-t,-u"],
    };
    const result = ssgSymbolInAxes(buildOps(rec), rec.q, "hexagonal");
    expect(result?.shifts).toEqual(["-h00", "000"]);
    expect(result?.full).toBe("P622(0,0,g1)-h00(0,0,g2)000");
  });
  it("changes the signed phase under reversal of the first internal axis", () => {
    const rec: GroupRecord = {
      sym: "", cls: "", d: 2, q: ["0,0,g1", "0,0,g2"],
      params: ["g1", "g2"], cen: [],
      gen: ["x-y,x,z,t+1/6,u", "-x,-x+y,-z,-t,-u"],
    };
    expect(ssgSymbolInAxes(buildOps(rec), rec.q, "hexagonal")?.shifts).toEqual(["h00", "000"]);
  });
  it("rejects inconsistent dimensions", () => {
    expect(() => ssgSymbolInAxes(buildOps(P3c1), ["0,0,g1", "0,0,g2"], "trigonal"))
      .toThrow(/dimension/);
  });
});

// The oracle: every (3+1)D BYU record's own symbol. BYU-derived, so local-only.
const DATA = join(__dirname, "..", "local", "data");
const have = existsSync(DATA);

function allD1(): (GroupRecord & { label: string })[] {
  const out: (GroupRecord & { label: string })[] = [];
  for (const f of readdirSync(DATA).filter((f) => /^sg\d{3}\.json$/.test(f)).sort()) {
    const groups: Record<string, GroupRecord> = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups;
    for (const [label, r] of Object.entries(groups)) if (r.d === 1) out.push({ ...r, label });
  }
  return out;
}

describe.skipIf(!have)("ssgSymbolInAxes vs all (3+1)D BYU symbols", () => {
  it("reproduces every record's own symbol", () => {
    const recs = allD1();
    expect(recs).toHaveLength(775);
    _stats.disagree = 0;
    const affected: string[] = [];
    const bad: string[] = [];
    for (const r of recs) {
      const sg = Number(r.label.split(".")[0]);
      const before = _stats.disagree;
      let got: string;
      try {
        got = ssgSymbolInAxes(buildOps(r), r.q, crystalSystemOf(sg))?.full ?? "null";
      } catch (e) {
        got = `throw ${(e as Error).message}`;
      }
      if (_stats.disagree > before) affected.push(`${sg} (${r.q.join(";")})`);
      if (got !== r.sym) bad.push(`${sg}\t${r.label}\t${r.sym}\tgot ${got}`);
    }
    console.log(`nearest-copy disagreements: ${_stats.disagree}; parents: ${[...new Set(affected)].join(", ")}`);
    expect(_stats.disagree).toBe(2);
    expect([...new Set(affected)].sort()).toEqual(["126 (1/2,1/2,g)", "133 (1/2,1/2,g)"]);
    if (bad.length) {
      const bySg = new Map<string, string[]>();
      for (const b of bad) {
        const k = b.split("\t")[0];
        if (!bySg.has(k)) bySg.set(k, []);
        bySg.get(k)!.push(b);
      }
      console.log(`${bad.length} mismatches of ${recs.length} in ${bySg.size} space groups; first 20:\n` + bad.slice(0, 20).join("\n"));
    }
    expect(bad.length).toBe(0);
  });
});

// Local oracle reads remain inside skipIf; nothing from local/data is embedded here.
type LabeledRecord = GroupRecord & { label: string };

function allDimension(d: number): LabeledRecord[] {
  const out: LabeledRecord[] = [];
  for (const f of readdirSync(DATA).filter(f => /^sg\d{3}\.json$/.test(f)).sort()) {
    const groups = JSON.parse(readFileSync(join(DATA, f), "utf8")).groups as Record<string, GroupRecord>;
    for (const [label, r] of Object.entries(groups)) if (r.d === d) out.push({ ...r, label });
  }
  return out;
}

// Match reference/check_equiv_d23.py: first record per Bravais class and basic SG,
// preserving file/JSON insertion order, then lexical label order.
function stratified(recs: LabeledRecord[]): LabeledRecord[] {
  const byClass = new Map<string, LabeledRecord>(), bySg = new Map<string, LabeledRecord>();
  for (const r of recs) {
    if (!byClass.has(r.cls)) byClass.set(r.cls, r);
    const sg = r.label.split(".")[0];
    if (!bySg.has(sg)) bySg.set(sg, r);
  }
  return [...new Map([...byClass.values(), ...bySg.values()].map(r => [r.label, r])).values()]
    .sort((a, b) => a.label < b.label ? -1 : a.label > b.label ? 1 : 0);
}

// Keep parameter sharing but ignore cosmetic parameter names when grouping failures.
function qPattern(q: string[]): string {
  const names = new Map<string, string>();
  return q.join(";").replace(/[a-z]+\d*/g, name => {
    if (!names.has(name)) names.set(name, `p${names.size + 1}`);
    return names.get(name)!;
  });
}

// Local TSV only: label, expected symbol, actual symbol. No oracle rows live in source.
function knownSymbolExceptions(d: number): Map<string, { expected: string; actual: string }> {
  const file = join(__dirname, "..", "local", "fixtures", `ssgsymbol-known-d${d}.txt`);
  const exceptions = new Map<string, { expected: string; actual: string }>();
  if (!existsSync(file)) return exceptions;
  for (const line of readFileSync(file, "utf8").split(/\r?\n/).filter(line => line.trim())) {
    const fields = line.split("\t");
    if (fields.length !== 3 || fields.some(field => !field) || fields[1] === fields[2])
      throw new Error(`invalid local d=${d} exception: expected label, expected symbol, different actual symbol`);
    const [label, expected, actual] = fields;
    if (exceptions.has(label)) throw new Error(`duplicate local d=${d} exception: ${label}`);
    exceptions.set(label, { expected, actual });
  }
  return exceptions;
}

describe.skipIf(!have)("ssgSymbolInAxes higher-dimensional local symbols", () => {
  for (const d of [2, 3]) it(`reproduces d=${d} own symbols`, async () => {
    const all = allDimension(d);
    expect(all).toHaveLength(d === 2 ? 3338 : 12584);
    const recs = d === 3 && process.env.SSG_SWEEP !== "d3" ? stratified(all) : all;
    if (d === 3 && process.env.SSG_SWEEP !== "d3") expect(recs).toHaveLength(430);
    const exceptions = d === 2 ? knownSymbolExceptions(2) : new Map<string, { expected: string; actual: string }>();
    const remaining = new Set(exceptions.keys());
    const exceptionErrors: string[] = [];
    const bad = new Map<string, string[]>(), disagreements = new Map<string, number>();
    _stats.disagree = 0;
    let nbad = 0, accepted = 0;
    for (const r of recs) {
      const sg = Number(r.label.split(".")[0]), bucket = `${sg} (${qPattern(r.q)})`;
      const before = _stats.disagree;
      let got: string;
      try { got = ssgSymbolInAxes(buildOps(r), r.q, crystalSystemOf(sg))?.full ?? "null"; }
      catch (error) { got = `throw ${(error as Error).message}`; }
      if (_stats.disagree !== before)
        disagreements.set(bucket, (disagreements.get(bucket) ?? 0) + _stats.disagree - before);
      const exception = exceptions.get(r.label);
      if (exception) {
        remaining.delete(r.label);
        if (got === r.sym) exceptionErrors.push(`${r.label}: fixed: remove from list`);
        else if (r.sym !== exception.expected || got !== exception.actual)
          exceptionErrors.push(`${r.label}: exception changed: expected ${r.sym}; got ${got}`);
        else accepted++;
      }
      if (got !== r.sym) {
        if (!bad.has(bucket)) bad.set(bucket, []);
        bad.get(bucket)!.push(`${r.label}: expected ${r.sym}; got ${got}`);
        nbad++;
      }
      // A parent may be quick, but a complete sweep must keep Vitest RPC responsive.
      await new Promise<void>(resolve => setImmediate(resolve));
    }
    console.log(`d=${d}: ${recs.length - nbad}/${recs.length}; nearest-copy disagreements ${_stats.disagree}`);
    console.log([...disagreements].map(([k, n]) => `${k}: ${n}`).join("\n"));
    if (nbad) {
      let left = 20;
      const lines: string[] = [];
      for (const [bucket, failures] of bad) {
        if (!left) break;
        lines.push(`${bucket} (${failures.length} mismatches):`, ...failures.slice(0, left));
        left -= Math.min(left, failures.length);
      }
      console.log(`${nbad} mismatches in ${bad.size} (basic SG, q-pattern) buckets; first 20:\n${lines.join("\n")}`);
    }
    for (const label of remaining) exceptionErrors.push(`${label}: listed record absent from sweep`);
    expect(exceptionErrors).toEqual([]);
    expect(nbad - accepted, "unexpected symbol mismatches outside local exception list").toBe(0);
  }, 600_000);
});

describe.skipIf(!have)("ssgSymbolInAxes signed local symbols", () => {
  it("preserves every signed symbol, including d3 records outside the sample", async () => {
    const signed = [...allDimension(2), ...allDimension(3)].filter(r => /-[thq]/.test(r.sym));
    expect(signed).toHaveLength(28);
    const exceptions = knownSymbolExceptions(3);
    const remaining = new Set(exceptions.keys());
    const bad: string[] = [];
    for (const r of signed) {
      const sg = Number(r.label.split(".")[0]);
      let got: string;
      try { got = ssgSymbolInAxes(buildOps(r), r.q, crystalSystemOf(sg))?.full ?? "null"; }
      catch (error) { got = `throw ${(error as Error).message}`; }
      const exception = r.d === 3 ? exceptions.get(r.label) : undefined;
      if (exception) {
        remaining.delete(r.label);
        if (got === r.sym) bad.push(`${r.label}: fixed: remove from list`);
        else if (r.sym !== exception.expected || got !== exception.actual)
          bad.push(`${r.label}: exception changed: expected ${r.sym}; got ${got}`);
      } else if (got !== r.sym) bad.push(`${sg} (${qPattern(r.q)}): expected ${r.sym}; got ${got}`);
      await new Promise<void>(resolve => setImmediate(resolve));
    }
    for (const label of remaining) bad.push(`${label}: listed record absent from signed sweep`);
    if (bad.length) console.log(`${bad.length}/${signed.length} signed mismatches:\n${bad.slice(0, 20).join("\n")}`);
    expect(bad).toHaveLength(0);
  });
});

describe("ssgSymbolInAxes matrix star (public)", () => {
  // Hand construction: the positive physical90-degree rotation sends
  // (q1,q2) to(-q2,q1), and reflection normal to x swaps q1 and q2.
  // The half translation on that swap is removable by an INTEGER internal
  // translate before projection, but the diagonal mirror retains(s,0).
  const diagonalStar: GroupRecord = {
    sym: "", cls: "", d: 2, q: ["a,a,0", "-a,a,0"], params: ["a"], cen: [],
    gen: ["-y,x,z,-u,t", "-x,y,z,u+1/2,t+1/2"],
  };
  it("keeps conventional square internal axes and coupled phase candidates", () => {
    expect(ssgSymbolInAxes(buildOps(diagonalStar), diagonalStar.q, "tetragonal")?.full)
      .toBe("P4mm(a,a,0)00s(-a,a,0)000");
  });
  it("is invariant under a nontrivial two-component internal origin shift", () => {
    // new coordinates t_new=t_old-p; w_new=w+(E-I)p.
    const p = [Q.of(1, 7), Q.of(2, 7)];
    const shifted = buildOps(diagonalStar).map(op => ({
      W: op.W,
      w: op.w.map((x, i) => i < 3 ? x : op.W[i].slice(3).reduce((s, a, j) => s.add(p[j].mulInt(a)), x.sub(p[i - 3])).mod1()),
    }));
    expect(ssgSymbolInAxes(shifted, diagonalStar.q, "tetragonal")?.full)
      .toBe("P4mm(a,a,0)00s(-a,a,0)000");
  });
});
