// Holds the desktop app drew as OpenSCAD/BOSL2 templates, rebuilt on
// manifold-3d. Each is modelled in the template's own frame -- wall on
// z = 0, hold standing up in +z, length along x -- and turned onto the
// wall at the end with the templates' own rotate([90, 0, 90]), which
// maps template (x, y, z) to hold (z, x, y).
//
// Only the templates written for the app itself are here; the ones
// carrying Jeremy Till's original geometry are not.

import {
  bezier, cutWall, ellipsoid, hull, K, loft, revolve, screwNegative, sphere, squircle, union,
  type Manifold,
} from "../kernel";
import { displaceField } from "../texture";
import { deg, frac, P, type HoldDef, type Values } from "../spec";

type V3 = [number, number, number];
const { sin, cos, tan, max, min, pow, sqrt, hypot, atan2, round, PI } = Math;
const sind = (d: number) => sin((d * PI) / 180);
const cosd = (d: number) => cos((d * PI) / 180);

/** Template frame -> wall frame, then the flat back. */
const onWall = (m: Manifold) => cutWall(m.rotate([90, 0, 90]));
/** A screw whose seat is at template point p (axis along template z). */
const tScrew = (p: V3) => screwNegative(p[2], p[0], p[1]);

/** hull(unit(i), unit(i+1)) for each segment, unioned. */
function chain(segs: number, unit: (i: number) => Manifold[]): Manifold {
  const parts: Manifold[] = [];
  let prev = unit(0);
  for (let i = 0; i < segs; i++) {
    const next = unit(i + 1);
    parts.push(hull([...prev, ...next]));
    prev = next;
  }
  return union(parts);
}

const SEGS = 24;

// -------------------------------------------------------------- gradient edge

const gradientEdge: HoldDef = {
  id: "gradient-edge",
  name: "Gradient edge",
  group: "Sculpted",
  blurb: "A long rail whose incut changes along its length, positive at one end and slopey at the other: one hold, many grades.",
  params: [
    P("length", "Length", 150, 500, 350),
    P("depth", "Depth", 14, 40, 25),
    P("drop", "Drop", 18, 50, 32),
    deg("angleStart", "Incut at start", -30, 35, 25),
    deg("angleEnd", "Incut at end", -30, 35, -15),
    P("edgeRadius", "Edge radius", 2, 8, 4),
    P("fasteners", "Screws", 2, 8, 3, 1, ""),
  ],
  build(v, screws) {
    const r = v.edgeRadius;
    const ex = (i: number) => -v.length / 2 + r + (v.length - 2 * r) * (i / SEGS);
    const ang = (t: number) => v.angleStart + (v.angleEnd - v.angleStart) * t;
    const body = chain(SEGS, (i) => {
      const t = i / SEGS, x = ex(i);
      return [sphere(r, [x, 0, v.depth - r], 24), sphere(r, [x, (v.depth - r) * tan((ang(t) * PI) / 180), 0], 24), sphere(r, [x, -(v.drop - r), 0], 24)];
    });
    const cuts: Manifold[] = [];
    if (screws)
      for (let i = 0; i < v.fasteners; i++) {
        const t = (2 * i + 1) / (2 * v.fasteners);
        cuts.push(tScrew([-v.length / 2 + r + (v.length - 2 * r) * t, -v.drop * 0.45, v.depth * 0.5]));
      }
    return { body: onWall(body), screws: cuts };
  },
};

// ---------------------------------------------------------------- system edge

// Port of system_edge.scad, which is in the hold frame already (no
// rotate): the hull of two half-balls at the top corners, a narrower
// dropped pair, and one ball stretched along the wall into the lip.
const systemEdge: HoldDef = {
  id: "system-edge",
  name: "System edge",
  group: "Sculpted",
  blurb: "A clean, uniform edge for symmetrical system-board training: two top corners, a narrower base and one long rounded lip. The Machinable edge is its CNC cousin.",
  params: [
    P("width", "Width", 60, 160, 100),
    P("depth", "Depth (lip)", 12, 35, 20),
    P("lipHeight", "Lip height", 0, 20, 7, 0.5),
    P("lipLength", "Lip length", 30, 120, 70),
    P("lowerWidth", "Base width", 30, 90, 50),
    P("drop", "Base drop", 15, 45, 28),
    P("edgeRadius", "Edge radius", 2, 6, 3, 0.5),
  ],
  build(v, screws) {
    const r = v.edgeRadius;
    const half = (y: number, z: number) => cutWall(sphere(r, [0, y, z], 32));
    const body = cutWall(hull([
      half(v.width / 2 - r, 0), half(-(v.width / 2 - r), 0),
      half(v.lowerWidth / 2 - r, -(v.drop - r)), half(-(v.lowerWidth / 2 - r), -(v.drop - r)),
      ellipsoid([r, v.lipLength / 2, r], [v.depth - r, 0, v.lipHeight], 48),
    ]));
    // the template's screws: the countersink's wide end sits 25.4 mm in
    // front of a plane (depth - 3) behind the wall; the centre one is
    // 3 mm deeper and further down
    const z0 = -0.18 * v.drop;
    const cuts = screws
      ? [
        screwNegative(28.4 - v.depth, v.width / 5, z0),
        screwNegative(28.4 - v.depth, -v.width / 5, z0),
        screwNegative(25.4 - v.depth, 0, z0 - 0.24 * v.drop),
      ]
      : [];
    return { body, screws: cuts };
  },
};

// ------------------------------------------------------------------ rock flake

const rockFlake: HoldDef = {
  id: "rock-flake",
  name: "Rock flake",
  group: "Sculpted",
  blurb: "A leaning, bowed blade of rock that fades out at both ends. Pull, pinch or side-pull it.",
  params: [
    P("length", "Length", 90, 260, 170),
    P("height", "Height", 30, 90, 55),
    P("thickness", "Thickness", 18, 50, 30),
    P("bow", "Bow", 0, 60, 30),
    deg("lean", "Lean", -30, 30, 15),
    frac("fade", "End fade", 0.15, 0.8, 0.35, 0.05),
  ],
  build(v, screws) {
    const bow = min(v.bow, v.length / 3);
    const h = (t: number) => v.height * (v.fade + (1 - v.fade) * sin(PI * t));
    const y = (t: number) => bow * sin(PI * t);
    const x = (t: number) => -v.length / 2 + (v.thickness * v.fade) / 2 + (v.length - v.thickness * v.fade) * t;
    const body = chain(SEGS, (i) => {
      const t = i / SEGS;
      return [K().Manifold.sphere(1, 24).scale([v.thickness / 4, v.thickness / 2, h(t)]).rotate([v.lean, 0, 0]).translate([x(t), y(t), 0])];
    });
    const cuts: Manifold[] = [];
    if (screws)
      for (const t of [0.25, 0.75]) {
        const hh = h(t);
        cuts.push(tScrew([x(t), y(t) - 0.8 * hh * sind(v.lean), hh * cosd(v.lean) - 6]));
      }
    return { body: onWall(body), screws: cuts };
  },
};

// -------------------------------------------------------------------- tufa rib

const tufaRib: HoldDef = {
  id: "tufa-rib",
  name: "Tufa rib",
  group: "Sculpted",
  blurb: "A snaking rib that swells in the middle and tapers at the ends, bent into a banana or an S.",
  params: [
    P("length", "Length", 100, 300, 200),
    P("width", "Width", 25, 70, 45),
    P("height", "Height", 20, 60, 40),
    P("bow", "Bow", 0, 60, 25),
    P("lobes", "Lobes (1 = banana, 2 = S)", 1, 2, 1, 1, ""),
    frac("taper", "Taper", 0.15, 0.8, 0.35, 0.05),
  ],
  build(v, screws) {
    const bow = min(v.bow, v.length / 3);
    const s = (t: number) => v.taper + (1 - v.taper) * sin(PI * t);
    const y = (t: number) => bow * sin(PI * v.lobes * t);
    const x = (t: number) => -v.length / 2 + (v.width * v.taper) / 2 + (v.length - v.width * v.taper) * t;
    const body = chain(SEGS, (i) => {
      const t = i / SEGS;
      return [K().Manifold.sphere(1, 24).scale([v.width / 4, (v.width / 2) * s(t), v.height * s(t)]).translate([x(t), y(t), 0])];
    });
    const cuts = screws ? [0.25, 0.75].map((t) => tScrew([x(t), y(t), v.height * s(t) - 6])) : [];
    return { body: onWall(body), screws: cuts };
  },
};

// ------------------------------------------------------------------- wave rail

const waveRail: HoldDef = {
  id: "wave-rail",
  name: "Wave rail",
  group: "Sculpted",
  blurb: "A rail whose top rolls up and down along its length: crests to pinch, troughs to crimp.",
  params: [
    P("length", "Length", 80, 260, 180),
    P("height", "Height", 18, 60, 35),
    P("depth", "Depth", 15, 50, 30),
    P("amp", "Wave height", 0, 20, 8),
    frac("waves", "Waves", 0.5, 4, 2, 0.5),
    P("phase", "Phase", 0, 360, 90, 5, "°"),
  ],
  build(v, screws) {
    const amp = max(0, min(v.amp, 0.45 * v.height));
    const crest = (t: number) => v.height / 2 + amp * sind(360 * v.waves * t + v.phase);
    const sx = (t: number) => -v.length / 2 + v.depth / 2 + (v.length - v.depth) * t;
    const body = chain(SEGS, (i) => {
      const t = i / SEGS;
      return [0, -(amp + 2)].map((dz) =>
        K().Manifold.sphere(1, 24).scale([v.depth / 4, v.depth / 2, v.height / 2]).translate([sx(t), 0, crest(t) + dz]));
    });
    const cuts: Manifold[] = [];
    if (screws) {
      const ts: number[] = [];
      for (let k = -3; k <= 2 * v.waves + 3; k++) {
        const t = (90 - v.phase + k * 180) / (360 * v.waves);
        if (t > 0.08 && t < 0.92) ts.push(t);
      }
      const use = ts.length >= 2 ? [ts[0], ts[ts.length - 1]] : [0.25, 0.75];
      for (const t of use) cuts.push(tScrew([sx(t), 0, crest(t) + v.height / 2 - 6]));
    }
    return { body: onWall(body), screws: cuts };
  },
};

// --------------------------------------------------------------- cobble sloper

const cobble: HoldDef = {
  id: "cobble",
  name: "Cobble sloper",
  group: "Sculpted",
  blurb: "A lofted superellipse dome with its crown drifted to one side. Pure open-hand friction, or add an incut scoop.",
  params: [
    P("w", "Width", 50, 160, 95),
    P("l", "Length", 40, 130, 70),
    P("d", "Depth", 20, 70, 38),
    frac("bulge", "Bulge", 1.2, 5, 2.6, 0.1),
    P("drift", "Crown drift", -40, 40, 16),
    frac("squareness", "Squareness", 0, 0.9, 0.25, 0.05),
    P("incut", "Incut scoop", 0, 20, 0),
    P("incutWidth", "Scoop width", 30, 100, 60),
    deg("incutTilt", "Scoop tilt", 0, 40, 25),
    P("fasteners", "Screws", 2, 4, 2, 1, ""),
  ],
  build(v, screws) {
    const bulge = v.bulge, drift = max(-v.w / 3, min(v.drift, v.w / 3));
    const incut = max(0, min(v.incut, v.d - 10));
    const f = (u: number) => pow(max(0, 1 - pow(u, bulge)), 1 / bulge);
    const n = 64, slices = 30;
    const rings: V3[][] = [];
    for (let i = 0; i <= slices; i++) {
      const u = i / slices, fu = max(f(u), 0.02);
      rings.push(squircle(v.w * fu, v.l * fu, v.squareness, n).map(([x, y]) => [x + drift * u * u, y, v.d * u] as V3));
    }
    let body = loft(rings);
    if (incut > 0) {
      const li = max(v.l * 0.55, incut + 8);
      body = body.subtract(
        K().Manifold.sphere(li / 2, 48).scale([max(20, min(v.incutWidth, v.w)) / li, 1, (2 * incut) / li]).rotate([v.incutTilt, 0, 0]).translate([drift, v.l * 0.06, v.d]),
      );
    }
    const cuts: Manifold[] = [];
    if (screws) {
      const us = pow(1 - pow(0.6, bulge), 1 / bulge);
      for (let k = 0; k < v.fasteners; k++) {
        const a = 180 / v.fasteners + (k * 360) / v.fasteners + (incut > 0 ? 90 : 0);
        const px = (0.6 * v.w) / 2 - 6 + drift * us * us;
        cuts.push(tScrew([px * cosd(a), px * sind(a), v.d * us]));
      }
    }
    return { body: onWall(body), screws: cuts };
  },
};

// ---------------------------------------------------------------- volume wedge

const wedge: HoldDef = {
  id: "volume-wedge",
  name: "Volume wedge",
  group: "Sculpted",
  blurb: "A small macro: a rounded prismoid whose top face is shrunk and shifted, giving one steep side and one slab.",
  params: [
    P("w", "Width", 40, 180, 110),
    P("h", "Height", 40, 150, 90),
    P("d", "Depth", 15, 80, 45),
    frac("topScale", "Top scale", 0.1, 0.9, 0.45, 0.05),
    P("shiftX", "Top shift (along wall)", -40, 40, 18),
    P("shiftY", "Top shift (up)", -40, 40, 10),
    P("rEdge", "Edge radius", 2, 15, 6),
  ],
  build(v, screws) {
    const r = max(1, min(v.rEdge, v.w / 2 - 2, v.h / 2 - 2, v.d - 3));
    const ts = max(0.1, min(v.topScale, 0.9));
    const sx = max(-v.w / 3, min(v.shiftX, v.w / 3)), sy = max(-v.h / 3, min(v.shiftY, v.h / 3));
    const a = v.w - 2 * r, b = v.h - 2 * r, H = v.d - r;
    // minkowski(convex prismoid, ball) = hull of balls on its corners
    const balls: Manifold[] = [];
    for (const ex of [-1, 1]) for (const ey of [-1, 1]) {
      balls.push(sphere(r, [(ex * a) / 2, (ey * b) / 2, -r], 32));
      balls.push(sphere(r, [sx + (ex * a * ts) / 2, sy + (ey * b * ts) / 2, H - r], 32));
    }
    const body = hull(balls);
    const cuts = screws ? [1, -1].map((s) => tScrew([sx + (s * (v.w * ts)) / 4, sy, v.d - 3])) : [];
    return { body: onWall(body), screws: cuts };
  },
};

// -------------------------------------------------- revolved, textured pinches

function sweptSlab(lp: number, wp: number, d: number, ctrl: [number, number][], tex: (p: V3) => number, texSize: number): Manifold {
  const prof: [number, number][] = [...bezier(ctrl, 32), [0, d], [0, 0]];
  const body = revolve(prof, 96).scale([lp / wp, 1, 1]);
  // fade the relief out at the wall so the back stays flat
  return displaceField(body, max(texSize / 5, 0.8), (p) => tex(p) * min(1, p[2] / 2));
}

const granite: HoldDef = {
  id: "granite-pinch",
  name: "Granite pinch",
  group: "Sculpted",
  blurb: "A low oval pinch with a lumpy, hilly granite skin baked into the shape.",
  params: [
    P("lp", "Length", 60, 160, 95),
    P("wp", "Width", 30, 70, 48),
    P("d", "Depth", 8, 19, 17),
    P("texSize", "Lump size", 4, 16, 8),
    frac("texDepth", "Lump depth", 0.4, 2, 1.2, 0.1),
  ],
  build(v, screws) {
    const wp = max(30, v.wp), lp = max(v.lp, wp), d = max(8, min(v.d, 20 - v.texDepth));
    const T = v.texSize;
    const hills = (p: V3) => v.texDepth * 0.5 * (1 + sin((2 * PI * p[0]) / T) * sin((2 * PI * p[1]) / T) * 0.6 + 0.4 * sin((2 * PI * (p[0] + p[1] * 0.6)) / (1.7 * T)));
    const body = sweptSlab(lp, wp, d, [[wp / 2, 0], [wp / 2, d * 0.55], [(wp / 2) * 0.75, d * 0.92], [(wp / 2) * 0.42, d]], hills, T);
    const cuts = screws ? [1, -1].map((s) => tScrew([s * lp * 0.28, 0, d * 0.8])) : [];
    return { body: onWall(body), screws: cuts };
  },
};

const ripple: HoldDef = {
  id: "ripple-slab",
  name: "Ripple slab",
  group: "Sculpted",
  blurb: "A flat oval slab with water-ripple ribs: a low foothold or a friction smear.",
  params: [
    P("ls", "Length", 60, 200, 120),
    P("ws", "Width", 50, 140, 85),
    P("d", "Depth", 8, 25, 14),
    P("rippleSize", "Ripple spacing", 5, 20, 10),
    frac("rippleDepth", "Ripple depth", 0.4, 2.5, 1.4, 0.1),
  ],
  build(v, screws) {
    const ws = v.ws, ls = max(v.ls, ws), d = max(8, min(v.d, 25));
    const T = v.rippleSize;
    const rib = (p: V3) => v.rippleDepth * 0.5 * (1 + cos((2 * PI * hypot(p[0] * (ws / ls), p[1])) / T));
    const body = sweptSlab(ls, ws, d, [[ws / 2, 0], [ws / 2, d * 0.7], [(ws / 2) * 0.8, d], [(ws / 2) * 0.4, d]], rib, T);
    const cuts = screws ? [1, -1].map((s) => tScrew([s * ls * 0.28, 0, d * 0.8])) : [];
    return { body: onWall(body), screws: cuts };
  },
};

// ------------------------------------------------------------------ alien face

function alien(v: Values, screws: boolean) {
  const { w, h, d } = v;
  const lift = v.lift, undercut = v.undercut, widest = 0.1;
  const cranium = v.cranium, jaw = v.jaw, face = 2.6;
  const eyeR = v.eyeR, eyeOut = v.eyeOut, eyeSep = v.eyeSep, eyeRise = 0.12;
  const noseLen = 20, noseW = 21, noseOut = v.noseOut, noseDrop = 0.3;
  const mouthW = 32, mouthH = 16, mouthOut = v.mouthOut, mouthSplit = 3.4, mouthDrop = 0.62;
  const screwScale = 0.7, screwDrop = 5.5, seatD = 15;

  const base = max(0.3, min((w / 2 - undercut) / (w / 2), 0.995));
  const az = max(d - lift, 1);
  const bz = max(lift, 0.01) / sqrt(max(1 - base * base, 1e-6));
  const hw = (t: number) =>
    t >= widest
      ? pow(max(0, 1 - pow((t - widest) / max(1 - widest, 1e-6), cranium)), 1 / cranium)
      : pow(max(0, 1 - pow((widest - t) / max(1 + widest, 1e-6), jaw)), 1 / jaw);
  const sAt = (z: number) => (z <= lift ? sqrt(max(0, 1 - ((lift - z) / bz) ** 2)) : pow(max(0, 1 - pow((z - lift) / az, face)), 1 / face));
  const zAt = (s: number) => lift + az * pow(max(0, 1 - pow(s, face)), 1 / face);
  const outr = (ang: number) => {
    let lo = -1, hi = 1;
    for (let n = 0; n < 32; n++) {
      const m = (lo + hi) / 2;
      if ((atan2((m * h) / 2, (hw(m) * w) / 2) * 180) / PI < ang) lo = m; else hi = m;
    }
    const t = (lo + hi) / 2;
    return hypot((hw(t) * w) / 2, (t * h) / 2);
  };
  const surfZ = (x: number, y: number) => {
    const r = hypot(x, y);
    return r < 0.001 ? d : zAt(min(r / outr((atan2(y, Math.abs(x)) * 180) / PI), 1));
  };

  // outline: right half bottom->top, then mirrored back down (CCW from +z)
  const nOut = 48;
  const half: [number, number][] = [];
  for (let i = 0; i <= nOut; i++) { const t = -1 + (2 * i) / nOut; half.push([(hw(t) * w) / 2, (t * h) / 2]); }
  const outline: [number, number][] = [...half];
  for (let i = nOut - 1; i >= 1; i--) outline.push([-half[i][0], half[i][1]]);
  const nb = 6, nf = 36;
  const levels: [number, number][] = [];
  for (let i = 0; i < nb; i++) { const z = (lift * i) / nb; levels.push([z, sAt(z)]); }
  for (let i = 0; i < nf; i++) { const s = 1 - i / nf; levels.push([zAt(s), s]); }
  const rings = levels.map(([z, s]) => outline.map(([x, y]) => [x * s, y * s, z] as V3));
  const head = loft(rings, [0, 0, d]);

  const feat = 32;
  const eyeX = (eyeSep * w) / 4, eyeY = (eyeRise * h) / 2;
  const eyes = [1, -1].map((m) => sphere(eyeR, [m * eyeX, eyeY, surfZ(eyeX, eyeY) - (eyeR - eyeOut)], 48));
  const yt = (-noseDrop * h) / 2 + noseLen / 2, ym = (-noseDrop * h) / 2, yb = (-noseDrop * h) / 2 - noseLen / 2;
  const rt = max(noseW * 0.17, 1.5), rm = max(noseW * 0.26, 1.8), rw = max(noseW * 0.32, 2);
  const nose = hull([
    sphere(rt, [0, yt, surfZ(0, yt) - (rt - noseOut * 0.35)], feat),
    sphere(rm, [0, ym, surfZ(0, ym) - (rm - noseOut * 0.8)], feat),
    ...[1, -1].map((m) => { const x = m * (noseW / 2 - rw); return sphere(rw, [x, yb, surfZ(x, yb) - (rw - noseOut)], feat); }),
  ]);
  const lipBalls = () => {
    const r = mouthH / 2, re = max(r * 0.5, 1.5), y = (-mouthDrop * h) / 2, sink = r - mouthOut;
    return [-1, 0, 1].map((i) => { const x = i * (mouthW / 2 - re); return sphere(i === 0 ? r : re, [x, y, surfZ(x, y) - sink], feat); });
  };
  const lips = hull(lipBalls());
  let body = union([head, ...eyes, nose, lips]);
  if (mouthSplit > 0) {
    const y = (-mouthDrop * h) / 2;
    body = body.subtract(hull(lipBalls()).translate([0, -y, 0]).scale([1.04, 0.18, 1]).translate([0, y, mouthOut + mouthH * 0.35 - mouthSplit]));
  }
  const cuts: Manifold[] = [];
  if (screws) {
    const sc = screwScale, scb = screwScale * 0.72;
    const top: V3 = [0, (h / 2) * sc, sc], bot: V3 = [0, (-h / 2) * scb, scb];
    const pair = (t: number): V3[] => [1, -1].map((m) => [(m * hw(t) * w) / 2 * sc, (t * h) / 2 * sc, sc] as V3);
    const pts = v.fasteners <= 2 ? [top, bot] : v.fasteners === 3 ? [top, ...pair(-0.55)] : [top, bot, ...pair(-0.2)];
    for (const p of pts) {
      const z = zAt(p[2]) - screwDrop;
      // the template also bores a wide seat pocket down to the screw
      cuts.push(K().Manifold.cylinder(d, seatD / 2, seatD / 2, 48).translate([p[0], p[1], z]).rotate([90, 0, 90]));
      cuts.push(tScrew([p[0], p[1], z]));
    }
  }
  return { body: onWall(body), screws: cuts };
}

const alienFace: HoldDef = {
  id: "alien-face",
  name: "Alien face",
  group: "Sculpted",
  blurb: "A superellipse skull with eyeballs, a nose and split lips: every feature is a hold.",
  params: [
    P("w", "Width", 70, 200, 124),
    P("h", "Height", 80, 240, 158),
    P("d", "Depth", 35, 110, 62),
    P("lift", "Wall lift", 2, 40, 13),
    P("undercut", "Undercut", 0, 18, 6),
    frac("cranium", "Cranium squareness", 1.6, 8, 4.2, 0.1),
    frac("jaw", "Jaw squareness", 1.1, 8, 3.2, 0.1),
    P("eyeR", "Eye radius", 6, 36, 19),
    P("eyeOut", "Eye bulge", 2, 28, 11),
    frac("eyeSep", "Eye spacing", 0.3, 1.3, 0.82),
    P("noseOut", "Nose bulge", 1, 16, 7),
    P("mouthOut", "Lip bulge", 1, 14, 6),
    P("fasteners", "Screws", 2, 4, 3, 1, ""),
  ],
  build: alien,
};

export const SCULPTED: HoldDef[] = [gradientEdge, systemEdge, rockFlake, tufaRib, waveRail, cobble, wedge, granite, ripple, alienFace];
