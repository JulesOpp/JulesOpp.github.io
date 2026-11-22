// Geometry kernel: manifold-3d with automatic memory management, plus the
// handful of primitives every hold is built from.
//
// manifold-3d objects live in WASM memory and must be .delete()d. Holds
// chain dozens of operations, so instead of tracking each one by hand
// every Manifold/CrossSection created inside `arena()` is recorded and
// freed when the arena closes. Builds return plain typed arrays, never
// live Manifolds, so nothing escapes an arena.
//
// House rule (from the desktop app): the mounting face is the x = 0
// plane, the hold stands off in +x, y runs along the wall, z is up.

import Module from "manifold-3d";
import type { CrossSection as CS, Manifold as M, ManifoldToplevel, Vec3 } from "manifold-3d";

export type Manifold = M;
export type CrossSection = CS;

let wasm: ManifoldToplevel | null = null;
let live: { delete(): void }[] | null = null;

/** Load manifold-3d once. `locateFile` points the browser at the .wasm. */
export async function initKernel(locateFile?: (f: string) => string): Promise<ManifoldToplevel> {
  if (wasm) return wasm;
  const m = await (Module as unknown as (o?: object) => Promise<ManifoldToplevel>)(locateFile ? { locateFile } : undefined);
  m.setup();
  patchForArena(m);
  wasm = m;
  return m;
}

export function K(): ManifoldToplevel {
  if (!wasm) throw new Error("geometry kernel not initialised");
  return wasm;
}

// Wrap every method that returns a Manifold/CrossSection so the result is
// recorded in the open arena.
function patchForArena(m: ManifoldToplevel) {
  for (const Cls of [m.Manifold, m.CrossSection] as unknown as { prototype: Record<string, unknown> }[]) {
    const wrap = (owner: Record<string, unknown>) => {
      for (const name of Object.getOwnPropertyNames(owner)) {
        if (name === "constructor" || name === "delete") continue;
        const desc = Object.getOwnPropertyDescriptor(owner, name);
        if (!desc || typeof desc.value !== "function") continue;
        const orig = desc.value as (...a: unknown[]) => unknown;
        owner[name] = function (this: unknown, ...args: unknown[]) {
          return track(orig.apply(this, args));
        };
      }
    };
    wrap(Cls.prototype);
    wrap(Cls as unknown as Record<string, unknown>);
  }
}

function track<T>(v: T): T {
  if (live && v && typeof v === "object") {
    if (Array.isArray(v)) v.forEach((x) => track(x));
    else if (typeof (v as { delete?: unknown }).delete === "function" &&
      (v instanceof K().Manifold || v instanceof K().CrossSection)) live.push(v as { delete(): void });
  }
  return v;
}

/** Run `fn`; free every Manifold/CrossSection it created. */
export function arena<T>(fn: () => T): T {
  const outer = live;
  live = [];
  try {
    return fn();
  } finally {
    const mine = live;
    live = outer;
    for (const o of mine) {
      try { o.delete(); } catch { /* already freed */ }
    }
  }
}

/** Record an object created by `new` (constructors are not patched). */
export const own = track;

// ---------------------------------------------------------------- primitives

export const DEG = Math.PI / 180;
let SEGS = 48;
/** Circular resolution used by every primitive (preview vs export). */
export function setSegments(n: number) { SEGS = Math.max(12, Math.round(n / 4) * 4); }
export function segments() { return SEGS; }

export function sphere(r: number, at: Vec3 = [0, 0, 0], segs = SEGS): Manifold {
  return K().Manifold.sphere(r, segs).translate(at);
}

/** Unit sphere scaled to semi-axes (a, b, c), then moved. */
export function ellipsoid(semi: Vec3, at: Vec3 = [0, 0, 0], segs = SEGS): Manifold {
  return K().Manifold.sphere(1, segs).scale(semi).translate(at);
}

export function hull(parts: (Manifold | Vec3)[]): Manifold {
  return K().Manifold.hull(parts);
}

export function union(parts: Manifold[]): Manifold {
  return parts.length === 1 ? parts[0] : K().Manifold.union(parts);
}

/** Keep x >= x0: the flat mounting back. */
export function cutWall(m: Manifold, x0 = 0): Manifold {
  return m.trimByPlane([1, 0, 0], x0);
}

/** Box of size (sx, sy, sz) centred on the origin, all edges rounded r. */
export function roundedBox(sx: number, sy: number, sz: number, r: number): Manifold {
  r = Math.min(r, sx / 2, sy / 2, sz / 2);
  const pts: Manifold[] = [];
  for (const ex of [-1, 1]) for (const ey of [-1, 1]) for (const ez of [-1, 1])
    pts.push(sphere(r, [ex * (sx / 2 - r), ey * (sy / 2 - r), ez * (sz / 2 - r)]));
  return hull(pts);
}

/** Polygon in OpenSCAD's rotate_extrude convention: (radius, height). */
export function revolve(profile: [number, number][], segs = SEGS): Manifold {
  return K().Manifold.revolve([ccw(profile)], segs);
}

/** Manifold's default fill rule drops clockwise contours: normalise. */
export function ccw(poly: [number, number][]): [number, number][] {
  let a = 0;
  for (let i = 0; i < poly.length; i++) {
    const p = poly[i], q = poly[(i + 1) % poly.length];
    a += p[0] * q[1] - q[0] * p[1];
  }
  return a < 0 ? [...poly].reverse() : poly;
}

/** Scale `m` so its bounding box is `size` (OpenSCAD resize()). */
export function resize(m: Manifold, size: Vec3): Manifold {
  const b = m.boundingBox();
  const s: Vec3 = [0, 1, 2].map((i) => {
    const ext = b.max[i] - b.min[i];
    return ext > 1e-9 && size[i] > 0 ? size[i] / ext : 1;
  }) as Vec3;
  return m.scale(s);
}

// --------------------------------------------------------------------- screws

/**
 * Removed metal/plastic for one flat-head wood screw (#8 / M4): a 13 mm
 * counterbore for the head and driver, a 90 deg countersink cone, and a
 * 5 mm shank running back through the wall face. `seatX` is where the
 * countersink's wide end (the counterbore floor) sits, measured forward
 * from the wall. The bore runs out through the front surface by `reach`.
 *
 * The same dimensions as the STEP builder (src/brep), so a printed and a
 * machined hold take the same screw.
 */
export const SCREW = { cbDia: 13, shankDia: 5, shankLen: 41.3 };

export function screwNegative(seatX: number, y: number, z: number, reach = 120): Manifold {
  const { cbDia, shankDia, shankLen } = SCREW;
  const cone = (cbDia - shankDia) / 2;
  // profile (radius, axial) revolved about the axial direction, then laid
  // along +x and moved into place
  const a0 = seatX - shankLen;
  const prof: [number, number][] = [
    [0, a0],
    [shankDia / 2, a0],
    [shankDia / 2, seatX - cone],
    [cbDia / 2, seatX],
    [cbDia / 2, seatX + reach],
    [0, seatX + reach],
  ];
  // revolve about z; then z -> x
  return revolve(prof, Math.max(24, SEGS)).rotate([0, 90, 0]).translate([0, y, z]);
}

/** Screw whose countersink starts `bury` mm below a surface at `surfaceX`. */
export function sunkScrew(surfaceX: number, y: number, z: number, bury = 3): Manifold {
  return screwNegative(Math.max(surfaceX - bury, 4), y, z);
}

// ----------------------------------------------------------------------- loft

/**
 * Skin a stack of closed rings (each the same vertex count, wound
 * counter-clockwise seen from +z of the stack) into a solid. The first
 * ring is capped flat; the last is capped flat, or closed to a point
 * when `apex` is given. This is the BOSL2 skin() the template holds used.
 */
export function loft(rings: Vec3[][], apex?: Vec3): Manifold {
  const n = rings[0].length;
  const L = rings.length;
  const verts: number[] = [];
  for (const ring of rings) {
    if (ring.length !== n) throw new Error("loft: rings differ in length");
    for (const p of ring) verts.push(p[0], p[1], p[2]);
  }
  const tris: number[] = [];
  for (let k = 0; k < L - 1; k++)
    for (let i = 0; i < n; i++) {
      const a = k * n + i, b = k * n + ((i + 1) % n), c = a + n, d = b + n;
      tris.push(a, b, d, a, d, c);
    }
  // bottom cap: fan from a centroid vertex (rings are star-shaped)
  const cap = (k: number, flip: boolean, tip?: Vec3) => {
    const ring = rings[k];
    const c = tip ?? ring.reduce<Vec3>((s, p) => [s[0] + p[0] / n, s[1] + p[1] / n, s[2] + p[2] / n], [0, 0, 0]);
    const ci = verts.length / 3;
    verts.push(c[0], c[1], c[2]);
    for (let i = 0; i < n; i++) {
      const a = k * n + i, b = k * n + ((i + 1) % n);
      if (flip) tris.push(ci, b, a); else tris.push(ci, a, b);
    }
  };
  cap(0, true);
  cap(L - 1, false, apex);
  const mesh = new (K().Mesh)({ numProp: 3, vertProperties: new Float32Array(verts), triVerts: new Uint32Array(tris) });
  mesh.merge();
  return own(new (K().Manifold)(mesh));
}

/** Cubic Bezier through 4 control points, `steps` segments (BOSL2 bezpath_curve, N=3). */
export function bezier(ctrl: [number, number][], steps: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i <= steps; i++) {
    const t = i / steps, u = 1 - t;
    const w = [u * u * u, 3 * u * u * t, 3 * u * t * t, t * t * t];
    out.push([0, 1].map((k) => w.reduce((s, wi, j) => s + wi * ctrl[j][k], 0)) as [number, number]);
  }
  return out;
}

/** Superellipse / squircle outline with `n` points, BOSL2-style squareness 0..1. */
export function squircle(sx: number, sy: number, squareness: number, n: number): [number, number][] {
  // BOSL2 squircle (fg type) approximated with a superellipse exponent
  // that runs from 2 (ellipse) at squareness 0 to ~10 near 1.
  const e = 2 + 8 * squareness * squareness;
  const out: [number, number][] = [];
  for (let i = 0; i < n; i++) {
    const t = (2 * Math.PI * i) / n;
    const c = Math.cos(t), s = Math.sin(t);
    out.push([(sx / 2) * Math.sign(c) * Math.abs(c) ** (2 / e), (sy / 2) * Math.sign(s) * Math.abs(s) ** (2 / e)]);
  }
  return out;
}

// ---------------------------------------------------------------- mesh output

export interface TriMesh {
  pos: Float32Array; // xyz per vertex
  idx: Uint32Array; // 3 per triangle
}

export function toTriMesh(m: Manifold): TriMesh {
  const mesh = m.getMesh();
  const np = mesh.numProp;
  const nv = mesh.vertProperties.length / np;
  const pos = new Float32Array(nv * 3);
  for (let i = 0; i < nv; i++) {
    pos[3 * i] = mesh.vertProperties[np * i];
    pos[3 * i + 1] = mesh.vertProperties[np * i + 1];
    pos[3 * i + 2] = mesh.vertProperties[np * i + 2];
  }
  return { pos, idx: new Uint32Array(mesh.triVerts) };
}

export function fromTriMesh(t: TriMesh): Manifold {
  const mesh = new (K().Mesh)({ numProp: 3, vertProperties: t.pos, triVerts: t.idx });
  mesh.merge();
  return own(new (K().Manifold)(mesh));
}
