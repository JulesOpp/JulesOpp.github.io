import { beforeAll, describe, expect, it } from "vitest";
import { existsSync, readFileSync, readdirSync } from "node:fs";
import { join } from "node:path";
import { indexKSublattices, isInvariant, reduce } from "../src/sublattices";
import { kSubgroups, lattice, rebase } from "../src/ksubgroups";
import { tSubgroups } from "../src/subgroups";
import { buildOps, type GroupRecord } from "../src/ssg";
import { loadMoyo } from "../src/moyo";
import { standardise } from "../src/standardise";
import { imatIdentity } from "../src/intmat";
import { Q, matInv, matMul, matVec } from "../src/rational";
beforeAll(() => loadMoyo());
const hand = (d = 1, inversion = false): GroupRecord => ({ sym: "hand", cls: "triclinic", d, q: d === 1 ? ["a,b,g"] : ["a1,b1,g1", "a2,b2,g2"], params: d === 1 ? ["a", "b", "g"] : ["a1", "b1", "g1", "a2", "b2", "g2"], cen: [], gen: inversion ? ["-x,-y,-z,-t"] : [] });
it("enumerates all lower HNFs at the supported small indices", () => {
    expect(indexKSublattices(4, 2)).toHaveLength(15);
    expect(indexKSublattices(4, 3)).toHaveLength(40);
    expect(indexKSublattices(6, 2)).toHaveLength(63);
    expect(isInvariant(indexKSublattices(4, 2)[0], [imatIdentity(4)])).toBe(true);
});
it("counts translation-only and inversion complements without merging classes", () => {
    expect(kSubgroups(hand(), 2)).toHaveLength(15);
    expect(kSubgroups(hand(), 3)).toHaveLength(40);
    const rows = kSubgroups(hand(1, true), 2);
    expect(rows.filter(r => r.pointOps.length === 2)).toHaveLength(30);
    expect(rows).toHaveLength(45); // 30 inversion classes + 15 translations-only classes.
    expect(rows.every(r => r.classSize === 1)).toBe(true);
    expect(rows.filter(r => r.maximal)).toHaveLength(30);
    expect(kSubgroups(hand(2), 2)).toHaveLength(31);
});
it("k=1 preserves the t-subgroup algorithm", () => {
    const rec = hand(1, true);
    expect(kSubgroups(rec, 1).map(({ k, tIndex, sigma, kind, ...c }) => c)).toEqual(tSubgroups(rec));
});
it("rebases a doubled physical axis and removes its rational q component", () => {
    const rec: GroupRecord = { sym: "hand mirror", cls: "monoclinic", d: 1, q: ["1/2,0,g"], params: ["g"], cen: [], gen: ["-x,-y,z,t"] };
    const H = imatIdentity(4);
    H[0][0] = 2n;
    const rb = rebase(rec, buildOps(rec), H);
    expect(rb.record.q).toEqual(["1,0,g"]);
    const named = standardise(rb.ops, rb.record.q);
    expect(named.ok).toBe(true);
    if (named.ok)
        expect(named.symbol.q).toContain("0,0,g");
    expect(buildOps(rb.record)).toHaveLength(buildOps(rec).length);
    // Independently apply the inverse coordinate map to every rebased operation.
    for (let i = 0; i < rb.ops.length; i++) {
        expect(matMul(matMul(rb.B, rb.ops[i].W.map(r => r.map(x => Q.of(x)))), matInv(rb.B)).map(r => r.map(String))).toEqual(buildOps(rec)[i].W.map(r => r.map(String)));
        expect(matVec(rb.B, rb.ops[i].w).map(String)).toEqual(buildOps(rec)[i].w.map(String));
    }
});
const data = join(__dirname, "../local/data"), fixtures = join(__dirname, "../local/fixtures/ksubgroups/d1/k2");
describe.skipIf(!existsSync(fixtures))("oracle fixture comparisons (local)", () => {
    it("matches counts, kinds, and indices for primitive and centred parents", () => {
        const db: Record<string, GroupRecord> = {};
        for (const file of readdirSync(data).filter(f => /^sg\d+\.json$/.test(f)))
            Object.assign(db, JSON.parse(readFileSync(join(data, file), "utf8")).groups);
        const labels = [...new Set([1, 2, 4, 5, 7, 9, 16, 19].map(sg => Object.keys(db).find(label => label.startsWith(sg + ".") && db[label].d === 1)!))];
        for (const label of labels) {
            const fixture = JSON.parse(readFileSync(join(fixtures, label + ".json"), "utf8"));
            const rows = kSubgroups(db[label], 2);
            expect(rows.map(r => `${r.index}/${r.kind}/${r.classSize}`).sort(), label).toEqual(fixture.classes.map((r: any) => `${r.index}/${r.kind}/${r.class_size}`).sort());
            // Compare actual complement translations, not just class counts.
            const Li = matInv(lattice(db[label])), ops = buildOps(db[label]);
            const key = (H: bigint[][], aff: {
                W: number[][];
                w: Q[];
            }[]) => JSON.stringify([H.map(r => r.map(String)), aff.map(o => JSON.stringify([o.W, reduce(matVec(Li, o.w), H).map(String)])).sort()]);
            expect(rows.map(r => key(r.sigma, r.ops)).sort(), label).toEqual(fixture.classes.map((r: any) => key(r.sigma.map((row: any[]) => row.map(BigInt)), r.pg_ops.map((j: number, i: number) => ({ W: ops[j].W, w: r.tw[i].map((x: any) => Q.parse(String(x))) })))).sort());
        }
    }, 120000);
});
