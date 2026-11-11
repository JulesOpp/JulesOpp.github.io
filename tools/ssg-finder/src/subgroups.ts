// k-index-1 (translationengleiche) subgroups of a (3+1)D superspace group.
//
// For incommensurate q the point group K_s of G is isomorphic to the basic point group K
// ((R, eps) -> R is injective), and a subgroup H with the same superspace lattice is the full
// preimage of a subgroup P of K. Classes up to conjugacy in G are K-conjugacy classes of
// subgroups of K. Mirrors reference/subgroups.py (without the BYU identification).

import { GroupRecord, SOp, buildOps, splitOp } from "./ssg";

export interface PointOp {
  R: number[][];
  eps: number; // +1 or -1
}

export interface SubgroupClass {
  pointOps: PointOp[]; // P, identity first
  ops: SOp[]; // canonical list of H, centrings included
  index: number; // |K| / |P|
  classSize: number; // number of conjugates in K
  maximal: boolean; // no subgroup class strictly between P and K
}

const rkey = (R: number[][]) => R.map((r) => r.join(",")).join(";");

function mul3(a: number[][], b: number[][]): number[][] {
  return a.map((row) => [0, 1, 2].map((j) => row.reduce((s, x, k) => s + x * b[k][j], 0)));
}

/** Point operations of a (3+1)D op list, one per distinct rotation part, in order of first appearance. */
export function pointGroup(ops: SOp[], d: number): PointOp[] {
  if (d !== 1) throw new Error("only d = 1 is supported for now");
  const seen = new Set<string>();
  const out: PointOp[] = [];
  for (const op of ops) {
    const { R, E } = splitOp(op, d);
    const k = rkey(R);
    if (seen.has(k)) continue;
    seen.add(k);
    out.push({ R, eps: E[0][0] });
  }
  return out;
}

export function tSubgroups(rec: GroupRecord): SubgroupClass[] {
  if (rec.d !== 1) throw new Error("only d = 1 is supported for now");
  const ops = buildOps(rec);
  const K = pointGroup(ops, 1);
  const n = K.length;
  const pos = new Map(K.map((p, i) => [rkey(p.R), i]));
  const table = K.map((a) => K.map((b) => pos.get(rkey(mul3(a.R, b.R)))!));
  const inv = K.map((_, i) => table[i].indexOf(0));

  // Subgroups = closure of {1} under joining with cyclic subgroups.
  const cyclic = K.map((_, g) => {
    const s = new Set([0]);
    for (let x = g; !s.has(x); x = table[x][g]) s.add(x);
    return [...s];
  });
  const closure = (gens: number[]): number[] => {
    const G = new Set([0]);
    let frontier = [0];
    while (frontier.length) {
      const next: number[] = [];
      for (const a of frontier)
        for (const g of gens) {
          const b = table[a][g];
          if (!G.has(b)) {
            G.add(b);
            next.push(b);
          }
        }
      frontier = next;
    }
    return [...G].sort((x, y) => x - y);
  };
  const subs = new Map<string, number[]>();
  const add = (S: number[]) => subs.set(S.join(","), S);
  add([0]);
  let todo = [[0]];
  while (todo.length) {
    const next: number[][] = [];
    for (const S of todo)
      for (let g = 1; g < n; g++) {
        const T = closure([...S, ...cyclic[g]]);
        if (!subs.has(T.join(","))) {
          add(T);
          next.push(T);
        }
      }
    todo = next;
  }
  const all = [...subs.values()];
  const cmpIdx = (a: number[], b: number[]) => {
    for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i] - b[i];
    return a.length - b.length;
  };
  const conj = (S: number[], g: number) => S.map((h) => table[table[g][h]][inv[g]]).sort((x, y) => x - y);

  const classes: { rep: number[]; size: number }[] = [];
  const seen = new Set<string>();
  for (const S of [...all].sort((a, b) => a.length - b.length || cmpIdx(a, b))) {
    if (seen.has(S.join(","))) continue;
    const orbit = new Map<string, number[]>();
    for (let g = 0; g < n; g++) {
      const c = conj(S, g);
      orbit.set(c.join(","), c);
    }
    for (const k of orbit.keys()) seen.add(k);
    classes.push({ rep: [...orbit.values()].sort(cmpIdx)[0], size: orbit.size });
  }

  const isSub = (a: number[], b: number[]) => a.every((x) => b.includes(x));
  const out: SubgroupClass[] = classes.map(({ rep, size }) => {
    const inP = new Set(rep.map((i) => rkey(K[i].R)));
    // A proper subgroup with nothing strictly between it and K (P = K itself is not maximal).
    const maximal = rep.length < n && !all.some((Q) => Q.length > rep.length && Q.length < n && isSub(rep, Q));
    return {
      pointOps: rep.map((i) => K[i]),
      ops: ops.filter((op) => inP.has(rkey(splitOp(op, 1).R))),
      index: n / rep.length,
      classSize: size,
      maximal,
      _idx: rep,
    } as SubgroupClass & { _idx: number[] };
  });
  out.sort((a, b) => a.index - b.index || cmpIdx((a as any)._idx, (b as any)._idx));
  for (const c of out) delete (c as any)._idx;
  return out;
}
