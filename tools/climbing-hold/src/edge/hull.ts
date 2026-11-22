// Convex hull of a handful of points, returned one polygon per plane
// (coplanar triangles merged), each wound counter-clockwise about its
// outward normal. Brute force is fine: the machinable edge has <= 6
// points. Merging matters -- a triangulated hull would leave seams
// across flat panels, and offsetting across a seam asks OpenCascade
// for a zero-radius cylinder.

type P3 = [number, number, number];
const sub = (a: P3, b: P3): P3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a: P3, b: P3): P3 => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a: P3, b: P3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a: P3): P3 => {
  const l = Math.hypot(a[0], a[1], a[2]);
  return [a[0] / l, a[1] / l, a[2] / l];
};

export function coreHullFaces(pts: P3[], eps = 1e-7): P3[][] {
  const n = pts.length;
  if (n < 4) throw new Error("hull needs at least four points");
  const planes: { n: P3; d: number }[] = [];
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++)
      for (let k = j + 1; k < n; k++) {
        const c = cross(sub(pts[j], pts[i]), sub(pts[k], pts[i]));
        if (Math.hypot(...c) < 1e-9) continue;
        let nn = norm(c);
        let d = dot(nn, pts[i]);
        let pos = 0, neg = 0;
        for (const q of pts) {
          const s = dot(nn, q) - d;
          if (s > eps) pos++;
          else if (s < -eps) neg++;
        }
        if (pos && neg) continue;
        if (pos) { nn = [-nn[0], -nn[1], -nn[2]]; d = -d; }
        if (!planes.some((p) => dot(p.n, nn) > 1 - 1e-9 && Math.abs(p.d - d) < 1e-6)) planes.push({ n: nn, d });
      }
  if (planes.length < 4) throw new Error("degenerate hull: the points are coplanar");
  return planes.map(({ n: nn, d }) => {
    const on = pts.filter((q) => Math.abs(dot(nn, q) - d) < 1e-6);
    const c: P3 = [0, 1, 2].map((a) => on.reduce((s, q) => s + q[a], 0) / on.length) as P3;
    const seed: P3 = Math.abs(nn[0]) > 0.9 ? [0, 1, 0] : [1, 0, 0];
    const u = norm(cross(nn, seed));
    const v = cross(nn, u);
    return [...on].sort((a, b) => Math.atan2(dot(sub(a, c), v), dot(sub(a, c), u)) - Math.atan2(dot(sub(b, c), v), dot(sub(b, c), u)));
  });
}
