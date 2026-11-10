import { describe, expect, it } from "vitest";
import { existsSync, readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { Q, q } from "../src/rational";
import { parseLinear, parseOp, formatOp, buildOps, evalQ, checkConsistency, paramsOf, GroupRecord } from "../src/ssg";
import { parseByuPage, PasteError } from "../src/paste";

const LOCAL = join(__dirname, "..", "local");
const hasLocal = existsSync(join(LOCAL, "cache")) && existsSync(join(LOCAL, "fixtures", "records"));

describe("rationals", () => {
  it("parse, reduce, mod1, decimals", () => {
    expect(Q.parse("6/-4" .replace("-", "")).toString()).toBe("3/2");
    expect(Q.parse("-3/4").mod1().toString()).toBe("1/4");
    expect(Q.parse("0.125").toString()).toBe("1/8");
    expect(q("7/3").floor()).toBe(2n);
    expect(q("-7/3").floor()).toBe(-3n);
    expect(() => Q.parse("1/0")).toThrow();
    expect(() => Q.parse("abc")).toThrow();
  });
});

describe("operators", () => {
  it("parse/format round trip", () => {
    for (const s of ["x,y,z,t", "-x+1/2,-y,z+1/2,-t,-u,v", "x,y,-z,t+1/2,-u+v,v", "-y,x-y,z+1/3,-t+u,t"]) {
      const nv = s.split(",").length;
      expect(formatOp(parseOp(s, nv))).toBe(s);
    }
  });
  it("linear expressions", () => {
    const { coef, c } = parseLinear("2a1-b1+1/4", ["a1", "b1"]);
    expect(coef.map(String)).toEqual(["2", "-1"]);
    expect(c.toString()).toBe("1/4");
    expect(() => parseLinear("a+", ["a"])).toThrow();
    expect(() => parseLinear("w", ["a"])).toThrow();
  });
  it("closure of Pnma(00g)0s0-like generators has identity first", () => {
    const rec: GroupRecord = {
      sym: "test", cls: "1.9", d: 1, q: ["0,0,g"], params: ["g"], cen: [],
      gen: ["-x+1/2,-y,z+1/2,t", "-x,y+1/2,-z,-t+1/2", "-x,-y,-z,-t"],
    };
    const ops = buildOps(rec);
    expect(ops.length).toBe(8);
    expect(formatOp(ops[0])).toBe("x,y,z,t");
    expect(checkConsistency(ops, evalQ(rec, { g: q("1/3") }), 1)).toBeNull();
  });
  it("params in first-appearance order", () => {
    expect(paramsOf(["a1,b1,0", "-a1,b1,0", "0,0,g2"])).toEqual(["a1", "b1", "g2"]);
  });
});

describe("paste mode", () => {
  it("rejects junk with a helpful error", () => {
    expect(() => parseByuPage("hello")).toThrow(PasteError);
  });
  it.skipIf(!hasLocal)("matches the Python reference on every cached BYU page (HTML and page text)", () => {
    const files = readdirSync(join(LOCAL, "cache")).filter((f) => f.endsWith(".html"));
    expect(files.length).toBeGreaterThan(0);
    for (const f of files) {
      const label = f.slice(0, -5);
      const html = readFileSync(join(LOCAL, "cache", f), "utf8");
      const ref = JSON.parse(readFileSync(join(LOCAL, "fixtures", "records", label + ".json"), "utf8"));
      // as pasted HTML source
      const rec = parseByuPage(html);
      expect(rec.label).toBe(label);
      expect(rec.q).toEqual(ref.record.q);
      expect(rec.params).toEqual(ref.record.params);
      expect(rec.cen).toEqual(ref.record.cen);
      expect(rec.gen).toEqual(ref.record.gen);
      expect(buildOps(rec).map(formatOp)).toEqual(ref.ops);
      // as copied page text (tags gone, line breaks kept)
      const text = html.replace(/<[^>]+>/g, "\n").replace(/&nbsp;/g, " ");
      expect(buildOps(parseByuPage(text)).map(formatOp)).toEqual(ref.ops);
    }
  });
});
