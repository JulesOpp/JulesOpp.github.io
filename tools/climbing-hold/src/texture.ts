// Procedural surface textures, ported from the desktop app's
// geometry/texture.py. The mesh is refined so every edge is shorter than
// a fraction of the feature size, then each vertex is pushed along its
// normal; the push fades to zero over the first few millimetres off the
// wall so the back stays flat.

import { fromTriMesh, toTriMesh, type Manifold, type TriMesh } from "./kernel";

type V3 = [number, number, number];
export type TexMode = "none" | "random" | "rock" | "stripes" | "knurl" | "scales" | "divots" | "grit" | "pebbles";
export const TEX_MODES: TexMode[] = ["none", "random", "rock", "stripes", "knurl", "scales", "divots", "grit", "pebbles"];
export const TEX_FADE_MM = 3;

const fract = (x: number) => x - Math.floor(x);

function vnoise3(x: number, y: number, z: number, seed = 0): number {
  const ix = Math.floor(x), iy = Math.floor(y), iz = Math.floor(z);
  const fx = x - ix, fy = y - iy, fz = z - iz;
  const ux = fx * fx * (3 - 2 * fx), uy = fy * fy * (3 - 2 * fy), uz = fz * fz * (3 - 2 * fz);
  let n = 0;
  for (let dx = 0; dx < 2; dx++) for (let dy = 0; dy < 2; dy++) for (let dz = 0; dz < 2; dz++) {
    const s = Math.sin((ix + dx) * 12.9898 + (iy + dy) * 78.233 + (iz + dz) * 37.719 + seed * 113.1) * 43758.5453;
    const w = (dx ? ux : 1 - ux) * (dy ? uy : 1 - uy) * (dz ? uz : 1 - uz);
    n += w * (fract(s) * 2 - 1);
  }
  return n;
}

function fbm3(x: number, y: number, z: number, octaves = 4): number {
  let total = 0, amp = 1, freq = 1, norm = 0;
  for (let o = 0; o < octaves; o++) {
    total += amp * vnoise3(x * freq, y * freq, z * freq, 31 * o);
    norm += amp; amp *= 0.5; freq *= 2;
  }
  return total / norm;
}

function pebble(x: number, y: number, z: number, size: number, seed = 11): number {
  const q = [x / size, y / size, z / size];
  const b = q.map(Math.floor);
  const hash = (c: number[], k: number) => fract(Math.sin(c[0] * 127.1 + c[1] * 311.7 + c[2] * 74.7 + k * 269.5 + seed * 113.1) * 43758.5453);
  const r2 = 0.85 * 0.85;
  let g = 0;
  for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
    const c = [b[0] + dx, b[1] + dy, b[2] + dz];
    const cx = c[0] + 0.15 + 0.7 * hash(c, 0), cy = c[1] + 0.15 + 0.7 * hash(c, 1), cz = c[2] + 0.15 + 0.7 * hash(c, 2);
    const d2 = (q[0] - cx) ** 2 + (q[1] - cy) ** 2 + (q[2] - cz) ** 2;
    if (d2 < r2) g += (0.55 + 0.9 * hash(c, 3)) * (1 - d2 / r2) ** 2;
  }
  return g;
}

/** Small seeded PRNG (mulberry32) so divots/grit are repeatable. */
function rng(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Area-weighted vertex normals of an indexed mesh. */
export function vertexNormals(t: TriMesh): Float32Array {
  const { pos: p, idx } = t;
  const n = new Float32Array(p.length);
  for (let i = 0; i < idx.length; i += 3) {
    const a = 3 * idx[i], b = 3 * idx[i + 1], c = 3 * idx[i + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2];
    const vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    for (const k of [a, b, c]) { n[k] += nx; n[k + 1] += ny; n[k + 2] += nz; }
  }
  for (let i = 0; i < n.length; i += 3) {
    const l = Math.hypot(n[i], n[i + 1], n[i + 2]) || 1;
    n[i] /= l; n[i + 1] /= l; n[i + 2] /= l;
  }
  return n;
}

/**
 * Refine `m` to `maxEdge` and push each vertex along its normal by
 * field(p). Used by the textures below and by holds with relief baked
 * into their shape.
 */
export function displaceField(m: Manifold, maxEdge: number, field: (p: V3) => number, faceBudget = 1_500_000): Manifold {
  let r = m.refineToLength(maxEdge);
  if (r.numTri() > faceBudget) r = m.refineToLength(maxEdge * Math.sqrt(r.numTri() / faceBudget));
  const t = toTriMesh(r);
  const n = vertexNormals(t);
  const p = t.pos;
  for (let i = 0; i < p.length; i += 3) {
    const d = field([p[i], p[i + 1], p[i + 2]]);
    p[i] += d * n[i]; p[i + 1] += d * n[i + 1]; p[i + 2] += d * n[i + 2];
  }
  return fromTriMesh(t);
}

function maxEdgeFor(mode: TexMode, size: number): number {
  switch (mode) {
    case "divots": return Math.max(size / 5, 1);
    case "grit": return Math.max(size / 4, 0.8);
    case "rock": return Math.max(size / 6, 1);
    case "pebbles": return Math.max(size / 6, 0.9);
    case "knurl": case "scales": return Math.max(size / 3, 1.2);
    default: return Math.max(size / 2, 1.4);
  }
}

/** Divots (sign -1) or grit (sign +1): scattered quartic bumps. */
function scatterOffsets(t: TriMesh, height: number, size: number, maxEdge: number, density: number, sign: number, seed = 7): Float32Array {
  const { pos: v, idx } = t;
  const nv = v.length / 3, nt = idx.length / 3;
  const off = new Float32Array(nv);
  const area = new Float64Array(nt);
  let total = 0;
  for (let f = 0; f < nt; f++) {
    const a = 3 * idx[3 * f], b = 3 * idx[3 * f + 1], c = 3 * idx[3 * f + 2];
    const ux = v[b] - v[a], uy = v[b + 1] - v[a + 1], uz = v[b + 2] - v[a + 2];
    const wx = v[c] - v[a], wy = v[c + 1] - v[a + 1], wz = v[c + 2] - v[a + 2];
    area[f] = 0.5 * Math.hypot(uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx);
    total += area[f];
  }
  if (total <= 0) return off;
  const rand = rng(seed);
  const gauss = () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());
  const count = Math.round(Math.min(4000, Math.max(8, (density * total) / (2.4 * size * size))));
  const cdf = new Float64Array(nt);
  let acc = 0;
  for (let f = 0; f < nt; f++) { acc += area[f]; cdf[f] = acc; }
  const rLo = Math.max(0.45 * size, 1.6 * maxEdge), rHi = Math.max(1.1 * size, 1.5 * rLo);
  type F = { c: V3; r: number; d: number };
  const feats: F[] = [];
  for (let k = 0; k < count; k++) {
    const x = rand() * total;
    let lo = 0, hi = nt - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (cdf[mid] < x) lo = mid + 1; else hi = mid; }
    const s = Math.sqrt(rand()), u2 = rand();
    const w = [1 - s, s * (1 - u2), s * u2];
    const c: V3 = [0, 0, 0];
    for (let j = 0; j < 3; j++) for (let ax = 0; ax < 3; ax++) c[ax] += w[j] * v[3 * idx[3 * lo + j] + ax];
    const r = Math.min(rHi, Math.max(rLo, (size / 2) * Math.exp(0.4 * gauss())));
    const d = Math.min(height * Math.min(1, Math.max(0.35, r / (0.8 * size))), Math.min(0.3 * r, (0.07 * r * r) / maxEdge));
    feats.push({ c, r, d });
  }
  feats.sort((a, b) => b.r - a.r);
  const kept: F[] = [];
  for (const f of feats) {
    if (!kept.some((k) => (k.c[0] - f.c[0]) ** 2 + (k.c[1] - f.c[1]) ** 2 + (k.c[2] - f.c[2]) ** 2 < (0.8 * Math.max(k.r, f.r)) ** 2)) kept.push(f);
  }
  if (!kept.length) return off;
  const cell = Math.max(...kept.map((k) => k.r));
  const key = (x: number, y: number, z: number) => `${Math.floor(x / cell)},${Math.floor(y / cell)},${Math.floor(z / cell)}`;
  const buckets = new Map<string, number[]>();
  for (let i = 0; i < nv; i++) {
    const k = key(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
    let b = buckets.get(k);
    if (!b) buckets.set(k, (b = []));
    b.push(i);
  }
  for (const { c, r, d } of kept) {
    const k0 = c.map((x) => Math.floor((x - r) / cell)), k1 = c.map((x) => Math.floor((x + r) / cell));
    for (let kx = k0[0]; kx <= k1[0]; kx++) for (let ky = k0[1]; ky <= k1[1]; ky++) for (let kz = k0[2]; kz <= k1[2]; kz++) {
      const b = buckets.get(`${kx},${ky},${kz}`);
      if (!b) continue;
      for (const i of b) {
        const d2 = (v[3 * i] - c[0]) ** 2 + (v[3 * i + 1] - c[1]) ** 2 + (v[3 * i + 2] - c[2]) ** 2;
        if (d2 >= r * r) continue;
        const val = sign * d * (1 - d2 / (r * r)) ** 2;
        off[i] = sign < 0 ? Math.min(off[i], val) : Math.max(off[i], val);
      }
    }
  }
  // two Laplacian passes round off creases where features meet
  for (let pass = 0; pass < 2; pass++) {
    const sum = new Float32Array(nv), cnt = new Float32Array(nv);
    for (let i = 0; i < idx.length; i += 3)
      for (let e = 0; e < 3; e++) {
        const a = idx[i + e], b = idx[i + ((e + 1) % 3)];
        sum[a] += off[b]; sum[b] += off[a]; cnt[a]++; cnt[b]++;
      }
    for (let i = 0; i < nv; i++) off[i] = 0.5 * off[i] + (0.5 * sum[i]) / Math.max(cnt[i], 1);
  }
  for (let i = 0; i < nv; i++) off[i] = sign < 0 ? Math.min(off[i], 0) : Math.max(off[i], 0);
  return off;
}

/** Apply one of the app's texture modes to a hold body. */
export function applyTexture(m: Manifold, mode: TexMode, height: number, size: number, density = 1): Manifold {
  if (mode === "none" || height <= 0) return m;
  size = Math.max(size, 0.8);
  const maxEdge = maxEdgeFor(mode, size);
  let r = m.refineToLength(maxEdge);
  const budget = 1_200_000;
  if (r.numTri() > budget) r = m.refineToLength(maxEdge * Math.sqrt(r.numTri() / budget));
  const t = toTriMesh(r);
  const n = vertexNormals(t);
  const p = t.pos, nv = p.length / 3;
  const disp = new Float32Array(nv);
  if (mode === "divots" || mode === "grit") {
    disp.set(scatterOffsets(t, height, size, maxEdge, density, mode === "divots" ? -1 : 1));
  } else if (mode === "pebbles") {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < nv; i++) { const g = pebble(p[3 * i], p[3 * i + 1], p[3 * i + 2], size); disp[i] = g; lo = Math.min(lo, g); hi = Math.max(hi, g); }
    for (let i = 0; i < nv; i++) disp[i] = height * ((disp[i] - lo) / Math.max(hi - lo, 1e-9) - 0.5);
  } else {
    for (let i = 0; i < nv; i++) {
      const x = p[3 * i], y = p[3 * i + 1], z = p[3 * i + 2];
      let d: number;
      if (mode === "stripes") d = 0.5 * height * Math.sin((2 * Math.PI * x) / size);
      else if (mode === "knurl") {
        const u1 = (x + y) / Math.SQRT2, u2 = (x - y) / Math.SQRT2;
        d = 0.25 * height * (Math.sin((2 * Math.PI * u1) / size) + Math.sin((2 * Math.PI * u2) / size));
      } else if (mode === "scales") {
        const tt = (2 * Math.PI * x) / size;
        d = 0.5 * height * Math.sin(tt + 0.65 * Math.sin(tt));
      } else if (mode === "rock") d = 0.5 * height * fbm3(x / size, y / size, z / size);
      else d = 0.5 * height * vnoise3(x / size, y / size, z / size);
      disp[i] = d;
    }
  }
  for (let i = 0; i < nv; i++) {
    const fade = Math.min(1, Math.max(0, p[3 * i] / TEX_FADE_MM)) * disp[i];
    p[3 * i] += fade * n[3 * i]; p[3 * i + 1] += fade * n[3 * i + 1]; p[3 * i + 2] += fade * n[3 * i + 2];
  }
  return fromTriMesh(t);
}
