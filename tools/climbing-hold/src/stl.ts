// Binary STL in and out.

import type { TriMesh } from "./kernel";

export function toBinarySTL(m: TriMesh, header = "climbing hold"): ArrayBuffer {
  const { pos: p, idx } = m;
  const n = idx.length / 3;
  const buf = new ArrayBuffer(84 + 50 * n);
  const dv = new DataView(buf);
  const h = new TextEncoder().encode(header.slice(0, 79));
  new Uint8Array(buf, 0, 80).set(h);
  dv.setUint32(80, n, true);
  let o = 84;
  for (let t = 0; t < n; t++) {
    const a = 3 * idx[3 * t], b = 3 * idx[3 * t + 1], c = 3 * idx[3 * t + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2];
    const vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const l = Math.hypot(nx, ny, nz) || 1;
    nx /= l; ny /= l; nz /= l;
    for (const val of [nx, ny, nz, p[a], p[a + 1], p[a + 2], p[b], p[b + 1], p[b + 2], p[c], p[c + 1], p[c + 2]]) {
      dv.setFloat32(o, val, true);
      o += 4;
    }
    dv.setUint16(o, 0, true);
    o += 2;
  }
  return buf;
}

/** Parse a binary STL and weld coincident vertices. */
export function fromBinarySTL(buf: ArrayBuffer): TriMesh {
  const dv = new DataView(buf);
  const n = dv.getUint32(80, true);
  if (84 + 50 * n !== buf.byteLength) throw new Error("not a binary STL");
  const map = new Map<string, number>();
  const pos: number[] = [];
  const idx = new Uint32Array(3 * n);
  for (let t = 0; t < n; t++) {
    for (let k = 0; k < 3; k++) {
      const o = 84 + 50 * t + 12 + 12 * k;
      const x = dv.getFloat32(o, true), y = dv.getFloat32(o + 4, true), z = dv.getFloat32(o + 8, true);
      const key = `${x},${y},${z}`;
      let i = map.get(key);
      if (i === undefined) { i = pos.length / 3; map.set(key, i); pos.push(x, y, z); }
      idx[3 * t + k] = i;
    }
  }
  return { pos: new Float32Array(pos), idx };
}
