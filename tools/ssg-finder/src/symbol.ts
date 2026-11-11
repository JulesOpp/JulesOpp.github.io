// Hermann–Mauguin-like symbol of a 3D space group in the GIVEN (parental) axes of the
// conventional supercell (plan §3, "The symbol in SSG (parental) axes"). Generalises the
// prototype's hm_in_setting (orthogonal axes only) to all crystal systems, using the
// symmetry directions of IT A Table 2.1.3.1 for the crystal system of the parent.
//
// Output: "<centring> <pos1> <pos2> ...", tokens separated by single spaces, e.g.
// "P 2_1/n 2_1/m 2_1/a", "C 1 1 2/e", "P 4_2/m 2_1/n 2/c", "R -3 2/m 1", "B -1".
//  * centring: P A B C I F R (R = obverse, "R(rev)" = reverse in hexagonal axes), or the explicit
//    non-zero centring vectors "{0,1/2,1/2;1/3,0,0}" if none fits ("{Npts}" when more than 8 points).
//  * screws are relative to the shortest lattice vector along the axis and glide letters to the
//    shortest lattice vectors along a, b, c (= the axes themselves for every conventional cell).
//  * mirror/glide priority m > e (two axial glides in one plane) > axial (cyclic: b,c / c,a / a,b)
//    > n > d > "g?"; rotation priority: pure rotation over screw (lowest screw index).
//  * one token per symmetry-direction position of the parent system; a position whose
//    equivalent directions carry different elements in the subgroup is written
//    "[e1;e2;...]" in the order of the direction list.
//  * if no position carries an element: "<centring> 1" or "<centring> -1".
// Returns null if some rotation/mirror axis is not one of the parental directions.

import { Q, matVec, type QMat, type QVec } from "./rational";
import { numToQ, type Op3 } from "./identify";

export type CrystalSystem =
  | "triclinic"
  | "monoclinic"
  | "orthorhombic"
  | "tetragonal"
  | "trigonal"
  | "hexagonal"
  | "cubic";

type V3 = [number, number, number];

/** Symmetry directions (IT A Table 2.1.3.1), one array of equivalent directions per position. */
export const SYMMETRY_DIRECTIONS: Record<CrystalSystem, V3[][]> = {
  triclinic: [],
  // full monoclinic symbols list all three axes (e.g. "P 1 2_1/c 1"), like orthorhombic
  monoclinic: [[[1, 0, 0]], [[0, 1, 0]], [[0, 0, 1]]],
  orthorhombic: [[[1, 0, 0]], [[0, 1, 0]], [[0, 0, 1]]],
  tetragonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0]], [[1, -1, 0], [1, 1, 0]]],
  trigonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0], [-1, -1, 0]], [[1, -1, 0], [1, 2, 0], [-2, -1, 0]]],
  hexagonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0], [-1, -1, 0]], [[1, -1, 0], [1, 2, 0], [-2, -1, 0]]],
  cubic: [
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
    [[1, -1, 0], [1, 1, 0], [0, 1, -1], [0, 1, 1], [-1, 0, 1], [1, 0, 1]],
  ],
};

/** Power k of the cubic 3-fold C: (u,v,w) -> (w,u,v) that maps a direction onto the first
 *  member family of its position ([001]; [1-10]/[110]); used to relabel a/b/c glides. */
const CUBIC_SHIFT: Record<string, number> = {
  "0,0,1": 0, "1,0,0": 2, "0,1,0": 1,
  "1,-1,0": 0, "1,1,0": 0, "0,1,-1": 2, "0,1,1": 2, "-1,0,1": 1, "1,0,1": 1,
};

/** Parent crystal system from a space-group number. */
export function crystalSystemOf(sgNumber: number): CrystalSystem {
  if (sgNumber <= 2) return "triclinic";
  if (sgNumber <= 15) return "monoclinic";
  if (sgNumber <= 74) return "orthorhombic";
  if (sgNumber <= 142) return "tetragonal";
  if (sgNumber <= 167) return "trigonal";
  if (sgNumber <= 194) return "hexagonal";
  return "cubic";
}

// ------------------------------------------------------------------ helpers

const h = Q.of(1n, 2n);
const third = Q.of(1n, 3n);
const twoThirds = Q.of(2n, 3n);
const Z = Q.ZERO;
const vkey = (v: QVec) => v.map(String).join(",");

interface Elem {
  W: QMat;
  ws: QVec[]; // all translations (mod 1) occurring with this W
  proper: boolean;
  order: number; // order of the proper part (W or -W)
  axis: QVec | null; // kernel of W - det(W) I, null for 1 and -1
}

function trace(W: QMat): Q {
  return W[0][0].add(W[1][1]).add(W[2][2]);
}
function det(W: QMat): Q {
  return W[0][0].mul(W[1][1].mul(W[2][2]).sub(W[1][2].mul(W[2][1])))
    .sub(W[0][1].mul(W[1][0].mul(W[2][2]).sub(W[1][2].mul(W[2][0]))))
    .add(W[0][2].mul(W[1][0].mul(W[2][1]).sub(W[1][1].mul(W[2][0]))));
}
function cross(a: QVec, b: QVec): QVec {
  return [a[1].mul(b[2]).sub(a[2].mul(b[1])), a[2].mul(b[0]).sub(a[0].mul(b[2])), a[0].mul(b[1]).sub(a[1].mul(b[0]))];
}
const isZeroV = (v: QVec) => v.every((x) => x.isZero());
function dot(a: QVec, b: QVec): Q {
  return a.reduce((s, x, i) => s.add(x.mul(b[i])), Z);
}

/** Axis (kernel of W - s I) for a rank-2 matrix W - s I. */
function axisOf(W: QMat, s: Q): QVec | null {
  const M = W.map((r, i) => r.map((x, j) => (i === j ? x.sub(s) : x)));
  for (const [i, j] of [[0, 1], [0, 2], [1, 2]]) {
    const c = cross(M[i], M[j]);
    if (!isZeroV(c)) return c;
  }
  return null;
}

/** u parallel to v (either sense)? */
function parallel(u: QVec, v: QVec): boolean {
  return isZeroV(cross(u, v));
}

function classify(W: QMat): { proper: boolean; order: number } {
  const d = det(W).toNumber();
  const tr = trace(W).toNumber() * d; // trace of the proper part
  const order = { 3: 1, 2: 6, 1: 4, 0: 3, [-1]: 2 }[Math.round(tr)];
  if (order === undefined) throw new Error(`not a crystallographic rotation (trace ${tr})`);
  return { proper: d > 0, order };
}

function frac(x: Q): string {
  return x.toString();
}

// ------------------------------------------------------------------ centring

function centringSymbol(cents: QVec[]): string {
  const set = new Set(cents.map(vkey));
  const S = (...vs: QVec[]) => vs.length === set.size && vs.every((v) => set.has(vkey(v)));
  if (set.size === 0) return "P";
  if (S([Z, h, h])) return "A";
  if (S([h, Z, h])) return "B";
  if (S([h, h, Z])) return "C";
  if (S([h, h, h])) return "I";
  if (S([Z, h, h], [h, Z, h], [h, h, Z])) return "F";
  if (S([twoThirds, third, third], [third, twoThirds, twoThirds])) return "R";
  if (S([third, twoThirds, third], [twoThirds, third, twoThirds])) return "R(rev)";
  if (cents.length > 7) return `{${cents.length + 1}pts}`;
  const sorted = [...cents].sort((a, b) => vkey(a).localeCompare(vkey(b)));
  return "{" + sorted.map((v) => v.map(frac).join(",")).join(";") + "}";
}

// ------------------------------------------------------------------ per-direction symbols

function sumPowers(W: QMat, n: number, w: QVec): QVec {
  // (1/n) sum_{k<n} W^k w : intrinsic (screw/glide) translation
  let acc: QVec = [Z, Z, Z];
  let v = w;
  for (let k = 0; k < n; k++) {
    acc = acc.map((x, i) => x.add(v[i]));
    v = matVec(W, v);
  }
  const inv = Q.of(1n, BigInt(n));
  return acc.map((x) => x.mul(inv));
}

/** Lattice context: centring vectors (incl. 0, mod 1) and shortest lattice vectors along axes. */
interface Ctx {
  cents: QVec[];
  axisLen: Q[]; // shortest lattice translation along a, b, c (fraction of the axis; 1 for the usual cells)
}

/** Smallest lambda in (0,1] with lambda*d in the lattice Z^3 + cents (d an integer direction). */
function shortestAlong(d: QVec, cents: QVec[]): Q {
  let best = Q.ONE;
  const k = d.findIndex((x) => !x.isZero());
  const dk = d[k].n < 0n ? -d[k].n : d[k].n;
  for (const c of cents) {
    for (let m = 0n; m < dk; m++) {
      // lambda * d_k = c_k + m (mod 1)  =>  lambda = (c_k + m)/d_k, up to sign of d_k
      let lam = c[k].add(Q.of(m)).div(d[k]).mod1();
      if (lam.isZero()) continue;
      if (d.every((x, i) => lam.mul(x).sub(c[i]).isInt()) && lam.cmp(best) < 0) best = lam;
    }
  }
  return best;
}

/** Integer shifts {0,1}^3 (with the centrings these generate every parallel element mod 2L). */
const SHIFTS: QVec[] = [0, 1].flatMap((i) => [0, 1].flatMap((j) => [0, 1].map((k) => [i, j, k].map((x) => Q.of(BigInt(x))))));

/** Screw index p (0..n-1) of translation part w along u = (shortest lattice vector along the axis). */
function screwIndex(W: QMat, n: number, w: QVec, u: QVec): number {
  const s = sumPowers(W, n, w);
  const k = u.findIndex((x) => !x.isZero());
  const pn = s[k].div(u[k]).mulInt(n);
  if (!pn.isInt()) return -1; // not a multiple of u/n: inconsistent operation
  return Number(((pn.n % BigInt(n)) + BigInt(n)) % BigInt(n));
}

/** Glide name for intrinsic glide vector g (components in units of the axis lattice lengths). */
function glideName(g: QVec): string {
  const r = g.map((x) => x.mod1());
  const nz = r.map((x, i) => (x.isZero() ? -1 : i)).filter((i) => i >= 0);
  if (!nz.length) return "m";
  if (nz.every((i) => r[i].eq(h))) return nz.length === 1 ? "abc"[nz[0]] : "n";
  if (nz.length >= 2 && nz.every((i) => r[i].eq(Q.of(1n, 4n)) || r[i].eq(Q.of(3n, 4n)))) return "d";
  return "g?";
}

/** Rotation sense: +1 if W turns counter-clockwise about u (right-handed axes assumed). */
function sense(W: QMat, u: QVec): number {
  const basisVecs: QVec[] = [
    [Q.ONE, Z, Z],
    [Z, Q.ONE, Z],
    [Z, Z, Q.ONE],
  ];
  for (const v of basisVecs) {
    if (parallel(u, v)) continue;
    return dot(u, cross(v, matVec(W, v))).sign();
  }
  return 0;
}

/** Translations of all parallel copies of an element: coset reps (incl. centrings) + {0,1}^3. */
function allTranslations(ws: QVec[]): QVec[] {
  const out = new Map<string, QVec>();
  for (const w of ws)
    for (const t of SHIFTS) {
      const v = w.map((x, i) => x.add(t[i]));
      out.set(vkey(v), v);
    }
  return [...out.values()];
}

/** Axial-glide preference for a plane with normal d: cyclic (b,c | c,a | a,b), as in Pbca/Ibca/Ia-3d. */
function axialOrder(d: QVec): string[] {
  const nz = d.map((x, i) => (x.isZero() ? -1 : i)).filter((i) => i >= 0);
  if (nz.length === 1) {
    const i = nz[0];
    return ["abc"[(i + 1) % 3], "abc"[(i + 2) % 3]];
  }
  return ["a", "b", "c"];
}

function mirrorSymbol(m: Elem, d: QVec, ctx: Ctx): string {
  // group the glide names by plane location (the anti-invariant part of the translation)
  const byLoc = new Map<string, Set<string>>();
  for (const w of allTranslations(m.ws)) {
    const g = sumPowers(m.W, 2, w);
    const loc = vkey(w.map((x, i) => x.sub(g[i]).mod1()));
    const name = glideName(g.map((x, i) => x.div(ctx.axisLen[i])));
    if (!byLoc.has(loc)) byLoc.set(loc, new Set());
    byLoc.get(loc)!.add(name);
  }
  const all = new Set([...byLoc.values()].flatMap((s) => [...s]));
  if (all.has("m")) return "m";
  for (const s of byLoc.values()) if (["a", "b", "c"].filter((x) => s.has(x)).length >= 2) return "e";
  for (const x of axialOrder(d)) if (all.has(x)) return x;
  if (all.has("n")) return "n";
  if (all.has("d")) return "d";
  return "g?";
}

function dirSymbol(elems: Elem[], d: QVec, ctx: Ctx): string {
  const along = elems.filter((e) => e.axis && parallel(e.axis, d));
  const proper = along.filter((e) => e.proper);
  const improper = along.filter((e) => !e.proper);
  let n = 1;
  for (const e of proper) n = Math.max(n, e.order);
  const has = (ord: number) => improper.some((e) => e.order === ord);
  // rotation token
  let rot = "";
  if (n > 1) {
    const lam = shortestAlong(d, ctx.cents);
    const u = d.map((x) => x.mul(lam));
    const cands = proper.filter((e) => e.order === n && (n === 2 || sense(e.W, d) > 0));
    let best = Infinity;
    for (const e of cands)
      for (const w of allTranslations(e.ws)) {
        const p = screwIndex(e.W, n, w, u);
        if (p >= 0 && p < best) best = p;
        if (best === 0) break;
      }
    rot = best === Infinity ? `${n}(?)` : best === 0 ? `${n}` : `${n}_${best}`;
  }
  // mirror/glide perpendicular to d (improper part of order 2)
  const m = improper.find((e) => e.order === 2);
  const mir = m ? mirrorSymbol(m, d, ctx) : "";
  if (n === 4 || n === 6) return mir ? `${rot}/${mir}` : rot;
  if (n === 3) {
    if (has(6)) return "-6";
    if (has(3)) return "-3";
    return rot;
  }
  if (n === 2) {
    if (has(4)) return "-4";
    return mir ? `${rot}/${mir}` : rot;
  }
  return mir || "1";
}

// ------------------------------------------------------------------ main entry

/**
 * Setting symbol of the group `ops` (conventional-supercell axes = parental axes scaled)
 * for a parent of crystal system `parentSystem`. Null if the group's symmetry axes are not
 * parental symmetry directions.
 */
export function settingSymbol(ops: Op3[], parentSystem: CrystalSystem): string | null {
  // group operations by W, translations mod 1
  const byW = new Map<string, { W: QMat; ws: Map<string, QVec> }>();
  for (const o of ops) {
    const W = o.R.map((r) => r.map(numToQ));
    const k = W.map(vkey).join(";");
    const w = o.t.map((x) => x.mod1());
    let g = byW.get(k);
    if (!g) byW.set(k, (g = { W, ws: new Map() }));
    g.ws.set(vkey(w), w);
  }
  const I = "1,0,0;0,1,0;0,0,1";
  const cents = byW.has(I) ? [...byW.get(I)!.ws.values()].filter((v) => !isZeroV(v)) : [];
  const cen = centringSymbol(cents);
  const zero: QVec = [Z, Z, Z];
  const all = [zero, ...cents];
  const ctx: Ctx = {
    cents: all,
    axisLen: [0, 1, 2].map((i) => shortestAlong([0, 1, 2].map((j) => (i === j ? Q.ONE : Z)), all)),
  };

  const elems: Elem[] = [];
  let inversion = false;
  for (const [k, g] of byW) {
    if (k === I) continue;
    const { proper, order } = classify(g.W);
    if (!proper && order === 1) {
      inversion = true;
      continue;
    }
    const axis = axisOf(g.W, proper ? Q.ONE : Q.ONE.neg());
    if (!axis) throw new Error("rotation without axis");
    elems.push({ W: g.W, ws: [...g.ws.values()], proper, order, axis });
  }

  const positions = SYMMETRY_DIRECTIONS[parentSystem];
  // every element axis must be a parental direction
  const allDirs = positions.flat().map((d) => d.map((x) => Q.of(BigInt(x))));
  for (const e of elems) if (!allDirs.some((d) => parallel(d, e.axis!))) return null;

  const tokens = positions.map((dirs) => {
    const syms = dirs.map((d) => dirSymbol(elems, d.map((x) => Q.of(BigInt(x))), ctx));
    // cubic: members are related by the 3-fold (x->y->z->x), which relabels the axial glides
    const canon = dirs.map((d, j) => {
      const k = parentSystem === "cubic" ? CUBIC_SHIFT[d.join(",")] ?? 0 : 0;
      return syms[j].replace(/[abc]/g, (c) => "abc"[("abc".indexOf(c) + k) % 3]);
    });
    return canon.every((s) => s === canon[0]) ? syms[0] : `[${syms.join(";")}]`;
  });
  if (tokens.every((t) => t === "1")) return `${cen} ${inversion ? "-1" : "1"}`;
  return [cen, ...tokens].join(" ");
}

/**
 * Port of short() from reference/validate_it_finder.py: full setting symbol to IT short
 * form ("P 2_1/n 2/m 2_1/n" -> "Pnmn"; monoclinic keeps only the unique axis). Identical to
 * the Python for triclinic, monoclinic and orthorhombic symbols; generalised so that a 4- or
 * 6-fold primary axis keeps its rotation ("P 4_2/m 2_1/n 2/m" -> "P42/mnm"), cubic symbols
 * reduce every n/m ("F 4/m -3 2/m" -> "Fm-3m"), "-3 2/m" becomes "-3m", and "R ... 1" drops the trailing 1.
 */
export function itShortForm(setting: string): string {
  const { cen, tokens } = shortTokens(setting);
  return cen + tokens.map((t) => t.token.replace(/_/g, "")).join("");
}

/**
 * The short-form decision behind itShortForm, with screw underscores kept ("4_2/m") and the
 * symmetry-direction position (index into the space-separated tokens after the centring) of
 * every kept token, e.g. "P 1 1 2_1/b" -> [{pos 2, "2_1/b"}], "P 2_1/b 2_1/n 2/m" -> b, n, m.
 */
export function shortTokens(setting: string): { cen: string; tokens: { pos: number; token: string }[] } {
  const [cen, ...rest] = setting.split(/\s+/).filter(Boolean);
  let ax = rest.map((a, pos) => ({ pos, token: a }));
  const non1 = ax.filter((a) => a.token !== "1");
  if (non1.length === 1) return { cen, tokens: non1 };
  const cubic = ax.length === 3 && /^-?3$/.test(ax[1].token);
  const has = (a: { token: string }) => a.token.includes("/");
  const slashes = ax.filter(has).length;
  if (ax.every(has) || (slashes >= 2 && (cubic || /^[46]/.test(ax[0].token))) || (slashes >= 1 && ax[0].token === "-3")) {
    ax = ax.map((a, i) =>
      !has(a) || (i === 0 && !cubic && /^[46]/.test(a.token)) ? a : { pos: a.pos, token: a.token.split("/")[1] },
    );
  }
  if (cen.startsWith("R") && ax.length === 3 && ax[2].token === "1") ax = ax.slice(0, 2);
  return { cen, tokens: ax };
}

// re-exported for tests
export const _internal = { glideName, centringSymbol, sense, shortestAlong, parallel };
