import { beforeAll, describe, expect, it } from "vitest";
import { loadMoyo } from "../src/moyo";
import { parseByuPage } from "../src/paste";
import { runSubgroups } from "../src/subfinder";
import { affineText, decodeSent, subgroupRecord } from "../src/ui/subflow";
import { compute } from "../src/ui/flow";

beforeAll(async () => {
  await loadMoyo();
});

// Hand-written BYU-style page text (no BYU data): P-3c1(0,0,g)000, IT A No. 165.
const P3C1 =
  "Superspace group: 165.1.24.1 P-3c1(0,0,g)000 Independent modulation vectors: q1 = (0,0,g) " +
  "Centering: (0,0,0,0) Non-lattice generators: Non-lattice operators: (x,y,z,t) (y,-x+y,-z,-t) (x-y,-y,-z+1/2,-t)";

describe("affineText", () => {
  it("writes x_std = P x + p and t' = m·x + η t + p4 as coordinates", () => {
    const t = { P: [["0", "-1", "0"], ["1", "-1", "0"], ["0", "0", "-1"]], p: ["0", "0", "-1/4"], eta: -1 as const, m: [0, 0, 0], p4: "0" };
    expect(affineText(t)).toBe("(-y, x-y, -z-1/4; -t)");
    const u = { P: [["1/2", "0", "0"], ["0", "0", "1"], ["-1/2", "1", "0"]], p: ["1/2", "0", "0"], eta: 1 as const, m: [1, 0, -2], p4: "1/4" };
    expect(affineText(u)).toBe("(1/2x+1/2, z, -1/2x+y; x-2z+t+1/4)");
  });
});

describe("send to finder", () => {
  let parent: ReturnType<typeof parseByuPage>;
  let parentSg: number;
  let rows: ReturnType<typeof runSubgroups>["rows"];
  beforeAll(() => {
    parent = parseByuPage(P3C1);
    const out = runSubgroups(P3C1, { kIndex: 1 });
    rows = out.rows;
    parentSg = out.parentSgNumber;
  });
  it("round-trips the P321 subgroup through JSON and runs in the finder", () => {
    const row = rows.find((r) => r.sgNumber === 150)!;
    const rec = subgroupRecord(parent, row, parentSg);
    expect(rec.label).toBeUndefined();
    expect(parentSg).toBe(165);
    expect(rec.sgNumber).toBe(165); // the parent's axes, not the P321 type's
    expect(rec.nops).toBe(6);
    const back = decodeSent(JSON.stringify(rec));
    expect(back).toEqual(rec);
    const out = compute(back!, { g: "1/3" });
    expect(out.ok).toBe(true);
  });
  it("the P-3c1 parent row keeps the identity map", () => {
    const parentRow = rows.find((r) => r.index === 1)!;
    expect(parentRow.sgNumber).toBe(165);
    expect(affineText(parentRow.type!)).toBe("(x, y, z; t)");
  });
  it("P1 subgroup (no generators) is valid", () => {
    const rec = subgroupRecord(parent, rows.find((r) => r.sgNumber === 1)!, parentSg);
    expect(rec.gen).toEqual([]);
    expect(decodeSent(JSON.stringify(rec))!.nops).toBe(1);
  });
  it("rejects missing or malformed records", () => {
    expect(decodeSent(null)).toBeNull();
    expect(decodeSent("not json")).toBeNull();
    expect(decodeSent(JSON.stringify({ sym: "x" }))).toBeNull();
    expect(decodeSent(JSON.stringify({ sym: "x", cls: "", d: 1, q: ["0,0,g"], params: ["g"], cen: [], gen: ["x,y"] }))).toBeNull();
  });
});
