// Finite-index complements in G/Σ′, conjugated in G. All lattice arithmetic is exact.
import { Q, denLcm, identity, matInv, matMul, matVec, toSafeInt, transpose, vecCmp, type QMat } from "./rational";
import { hnf, hnfLower, imatDet, imatIdentity, type IMat } from "./intmat";
import { buildOps, formatOp, parseLinear, type GroupRecord, type SOp } from "./ssg";
import { pointGroup, tSubgroups, type SubgroupClass } from "./subgroups";
import { formatQ } from "./standardise";
import { indexKSublattices, reduce } from "./sublattices";
export type LatticeKind = "3d" | "internal" | "mixed";
export interface KSubgroupClass extends SubgroupClass {
    k: number;
    tIndex: number;
    sigma: IMat;
    kind: LatticeKind;
}
interface Aff {
    W: QMat;
    w: Q[];
}
interface Pair {
    H: IMat;
    ops: Aff[];
}
const qm = (A: IMat | number[][]): QMat => A.map(r => r.map(x => Q.of(x)));
const ints = (A: QMat): IMat => A.map(r => r.map(x => { toSafeInt(x, "lattice rotation"); return x.n; }));
const wk = (W: QMat) => W.map(r => r.join(",")).join(";");
const ak = (a: Aff) => wk(a.W) + "|" + a.w.join(",");
const compare = (a: Aff, b: Aff) => vecCmp(a.W.flat(), b.W.flat()) || vecCmp(a.w, b.w);
const pairKey = (p: Pair) => p.H.map(r => r.join(",")).join(";") + "/" + p.ops.map(ak).join("/");
const canonical = (H: IMat, ops: Aff[]): Pair => ({ H, ops: ops.map(o => ({ W: o.W, w: reduce(o.w, H) })).sort(compare) });
function comparePair(a: Pair, b: Pair): number {
    const c = vecCmp(a.H.flat().map(x => Q.of(x)), b.H.flat().map(x => Q.of(x)));
    if (c)
        return c;
    for (let i = 0; i < Math.min(a.ops.length, b.ops.length); i++) {
        const c = compare(a.ops[i], b.ops[i]);
        if (c)
            return c;
    }
    return a.ops.length - b.ops.length;
}
/** Column basis of the full translation lattice, including centrings. */
export function lattice(rec: GroupRecord): QMat {
    const rows = [...identity(3 + rec.d), ...rec.cen.map(c => c.split(",").map(Q.parse))];
    return transpose(rowBasis(rows));
}
function rowBasis(rows: QMat): QMat {
    const D = denLcm(rows.flat());
    return hnf(rows.map(r => r.map(x => x.n * (D / x.d)))).H.map(r => r.map(x => Q.of(x, D)));
}
function multiply(a: Aff, b: Aff, H: IMat): Aff {
    const v = matVec(a.W, b.w);
    return { W: matMul(a.W, b.W), w: reduce(a.w.map((x, i) => x.add(v[i])), H) };
}
function close(gens: Aff[], H: IMat, limit: number): Aff[] | null {
    const e = { W: identity(H.length), w: H.map(() => Q.ZERO) };
    const seen = new Map([[ak(e), e]]), todo = [e];
    while (todo.length) {
        const a = todo.pop()!;
        for (const g of gens) {
            const b = multiply(a, g, H), key = ak(b);
            if (!seen.has(key)) {
                seen.set(key, b);
                if (seen.size > limit)
                    return null;
                todo.push(b);
            }
        }
    }
    return [...seen.values()];
}
function conjugate(p: Pair, g: Aff): Pair {
    const H = hnfLower(ints(transpose(matMul(g.W, transpose(qm(p.H))))));
    const Vi = matInv(g.W);
    return canonical(H, p.ops.map(o => {
        const W = matMul(matMul(g.W, o.W), Vi), Vw = matVec(g.W, o.w), Wv = matVec(W, g.w);
        return { W, w: g.w.map((v, i) => v.add(Vw[i]).sub(Wv[i])) };
    }));
}
function orbit(p: Pair, gens: Aff[]): Pair[] {
    const seen = new Map([[pairKey(p), p]]), todo = [p];
    while (todo.length) {
        const a = todo.pop()!;
        for (const g of gens) {
            const b = conjugate(a, g), key = pairKey(b);
            if (!seen.has(key)) {
                seen.set(key, b);
                todo.push(b);
            }
        }
    }
    return [...seen.values()];
}
function pointSubs(table: number[][]): {
    members: number[];
    gens: number[];
}[] {
    const closure = (gens: number[]) => { const seen = new Set([0]), todo = [0]; while (todo.length) {
        const a = todo.pop()!;
        for (const g of gens) {
            const b = table[a][g];
            if (!seen.has(b)) {
                seen.add(b);
                todo.push(b);
            }
        }
    } return [...seen].sort((a, b) => a - b); };
    const initial = { members: [0], gens: [] as number[] }, seen = new Map([["0", initial]]), todo = [initial];
    while (todo.length) {
        const a = todo.pop()!;
        for (let g = 1; g < table.length; g++) {
            if (a.members.includes(g))
                continue;
            const gens = [...a.gens, g], members = closure(gens), key = members.join(",");
            if (!seen.has(key)) {
                const b = { members, gens };
                seen.set(key, b);
                todo.push(b);
            }
        }
    }
    return [...seen.values()];
}
function choices<T>(items: T[], n: number): T[][] { if (!n)
    return [[]]; return items.flatMap(x => choices(items, n - 1).map(r => [x, ...r])); }
function representatives(H: IMat): Q[][] {
    let out: Q[][] = [[]];
    for (let i = 0; i < H.length; i++)
        out = out.flatMap(r => Array.from({ length: toSafeInt(H[i][i], "index") }, (_, v) => [...r, Q.of(v)]));
    return out;
}
export function latticeKind(rec: GroupRecord, H: IMat): LatticeKind {
    const L = lattice(rec), parent = rowBasis(transpose(L)), child = rowBasis(transpose(matMul(L, transpose(qm(H)))));
    const eq = (A: QMat, B: QMat) => A.length === B.length && A.every((r, i) => !vecCmp(r, B[i]));
    if (eq(parent.slice(3), child.slice(3)))
        return "3d";
    return eq(rowBasis(parent.map(r => r.slice(0, 3))), rowBasis(child.map(r => r.slice(0, 3)))) ? "internal" : "mixed";
}
export function kSubgroups(rec: GroupRecord, k: number): KSubgroupClass[] {
    if (![1, 2, 3].includes(rec.d))
        throw new RangeError("d must be 1, 2 or 3");
    if (![1, 2, 3].includes(k))
        throw new RangeError("k must be 1, 2 or 3");
    const n = 3 + rec.d;
    if (k === 1)
        return tSubgroups(rec).map(c => ({ ...c, k, tIndex: c.index, sigma: imatIdentity(n), kind: "3d" }));
    const L = lattice(rec), Li = matInv(L), first = new Map<string, SOp>();
    for (const o of buildOps(rec))
        if (!first.has(JSON.stringify(o.W)))
            first.set(JSON.stringify(o.W), o);
    const lifts = [...first.values()].map(o => ({ W: matMul(matMul(Li, qm(o.W)), L), w: matVec(Li, o.w) }));
    const pos = new Map(lifts.map((o, i) => [wk(o.W), i]));
    const table = lifts.map(a => lifts.map(b => pos.get(wk(matMul(a.W, b.W)))!));
    const subs = pointSubs(table), seen = new Set<string>(), out: KSubgroupClass[] = [];
    const conjugators = [...lifts, ...identity(n).map(w => ({ W: identity(n), w }))];
    for (const H of indexKSublattices(n, k)) {
        const reps = representatives(H);
        const stabiliser = new Set(lifts.flatMap((g, i) => H.every(h => reduce(matVec(g.W, h.map(x => Q.of(x))), H).every(x => x.isZero())) ? [i] : []));
        for (const P of subs) {
            if (!P.gens.every(g => stabiliser.has(g)))
                continue;
            for (const shifts of choices(reps, P.gens.length)) {
                const gens = P.gens.map((g, i) => ({ W: lifts[g].W, w: lifts[g].w.map((v, j) => v.add(shifts[i][j])) }));
                const ops = close(gens, H, P.members.length);
                if (!ops)
                    continue;
                const p = canonical(H, ops);
                if (seen.has(pairKey(p)))
                    continue;
                const orb = orbit(p, conjugators);
                for (const o of orb)
                    seen.add(pairKey(o));
                const best = orb.sort(comparePair)[0];
                const physical: SOp[] = best.ops.map(o => ({ W: matMul(matMul(L, o.W), Li).map(r => r.map(x => toSafeInt(x, "physical rotation"))), w: matVec(L, o.w) }));
                const tIndex = lifts.length / P.members.length;
                out.push({ ops: physical, pointOps: pointGroup(physical, rec.d), index: k * tIndex, tIndex, k, sigma: best.H, kind: latticeKind(rec, best.H), classSize: orb.length, maximal: tIndex === 1 });
            }
        }
    }
    return out.sort((a, b) => a.index - b.index || vecCmp(a.sigma.flat().map(x => Q.of(x)), b.sigma.flat().map(x => Q.of(x))) || a.ops.map(formatOp).join(";").localeCompare(b.ops.map(formatOp).join(";")));
}
/** Block triangular primitive basis and coordinate map for Σ′. */
export function rebase(rec: GroupRecord, ops: SOp[], sigma: IMat): {
    record: GroupRecord;
    ops: SOp[];
    S: QMat;
    B: QMat;
} {
    const B = transpose(rowBasis(transpose(matMul(lattice(rec), transpose(qm(sigma))))));
    const scale = denLcm(B.flat());
    if (imatDet(B.slice(0, 3).map(r => r.slice(0, 3).map(x => x.n * (scale / x.d)))) < 0n)
        for (const row of B)
            row[0] = row[0].neg();
    const S = matInv(B);
    const rebased = ops.map(o => ({ W: matMul(matMul(S, qm(o.W)), B).map(r => r.map(x => toSafeInt(x, "rebased rotation"))), w: matVec(S, o.w).map(x => x.mod1()) }));
    const A = B.slice(0, 3).map(r => r.slice(0, 3)), C = B.slice(3).map(r => r.slice(0, 3)), Di = matInv(B.slice(3).map(r => r.slice(3)));
    const parsed = rec.q.map(row => row.split(",").map(x => parseLinear(x, rec.params)));
    const coefficients = Array.from({ length: rec.params.length + 1 }, (_, p) => {
        const Qm = parsed.map(row => row.map(x => p === rec.params.length ? x.c : x.coef[p]));
        const QA = matMul(Qm, A);
        return matMul(Di, QA.map((r, i) => r.map((x, j) => p === rec.params.length ? x.sub(C[i][j]) : x)));
    });
    const q = rec.q.map((_, i) => formatQ([0, 1, 2].map(j => ({ coef: coefficients.slice(0, -1).map(m => m[i][j]), c: coefficients.at(-1)![i][j] })), rec.params));
    const record: GroupRecord = { sym: "rebased subgroup", cls: rec.cls, d: rec.d, q, params: [...rec.params], cen: [], gen: rebased.map(formatOp) };
    record.nops = buildOps(record).length;
    return { record, ops: rebased, S, B };
}
