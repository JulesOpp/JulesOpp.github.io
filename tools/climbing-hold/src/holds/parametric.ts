// Equation-driven holds, ported from the desktop app's
// holdlib/parametric.py (the holds written for the app itself -- the
// Jeremy Till "spherical" edge, foot chip, star and crimp rail are not
// included). Every body has its flat back on x = 0.

import {
  cutWall, ellipsoid, hull, K, resize, revolve, roundedBox, sphere, sunkScrew, union,
  type Manifold,
} from "../kernel";
import { deg, frac, P, type HoldDef } from "../spec";

const { sqrt, sin, cos, atan2, exp, hypot, min, max, PI } = Math;
const rad = (d: number) => (d * PI) / 180;

// ---------------------------------------------------------------- jug (bucket)

/** Half cross-section of a unit bowl with a rim tangent to both surfaces. */
function bowlProfile(wall = 0.3, thick: number | null = null, lift = 0.1, lip = 0.13, steps = 40): [number, number][] {
  const R = 1;
  let Ri = R - wall;
  let K2 = (R - lip) ** 2 - (Ri + lip) ** 2;
  const v = (K2 + lift ** 2) / (2 * lift);
  if (thick !== null) {
    const fMin = 0.05;
    const cx0 = sqrt(max((R - lip) ** 2 - v * v, 0));
    let tLo = 0.06, tHi = min(R + lip - cx0 - 0.02, 0.6);
    const tF = ((R + lip) ** 2 + v * v - (R - lip) ** 2 - (v - fMin) ** 2) / (2 * (R + lip + v - fMin));
    if (v >= lift && tF + v >= fMin) tLo = max(tLo, tF);
    else if (v < lift && fMin - tF - v >= 0) tHi = min(tHi, tF);
    thick = min(max(thick, min(tLo, wall)), max(tHi, wall));
    Ri = R - thick;
    K2 = (R - lip) ** 2 - (Ri + lip) ** 2;
    const s = sqrt(max(v * v - K2, 0));
    lift = v >= lift ? v - s : v + s;
  }
  const cx = sqrt(max((R - lip) ** 2 - v ** 2, 1e-9));
  const cz = R + v;
  const T1 = [(cx * R) / (R - lip), R + (v * R) / (R - lip)];
  const T2 = [(cx * Ri) / (Ri + lip), R + lift + ((v - lift) * Ri) / (Ri + lip)];
  const pts: [number, number][] = [[0, 0]];
  const a1 = atan2(T1[0], R - T1[1]);
  for (let i = 1; i <= steps; i++) {
    const a = (a1 * i) / steps;
    pts.push([R * sin(a), R * (1 - cos(a))]);
  }
  const f1 = atan2(T1[1] - cz, T1[0] - cx);
  let f2 = atan2(T2[1] - cz, T2[0] - cx);
  if (f2 < f1) f2 += 2 * PI;
  for (let i = 1; i <= steps; i++) {
    const f = f1 + ((f2 - f1) * i) / steps;
    pts.push([cx + lip * cos(f), cz + lip * sin(f)]);
  }
  const b2 = atan2(T2[0], R + lift - T2[1]);
  for (let i = 1; i <= steps; i++) {
    const b = (b2 * (steps - i)) / steps;
    pts.push([Ri * sin(b), R + lift - Ri * cos(b)]);
  }
  return pts.map(([x, y]) => [max(x, 0), y]);
}

const jug: HoldDef = {
  id: "jug",
  name: "Jug (bucket)",
  group: "Parametric",
  blurb: "A deep revolved bowl, tilted back so the lip is incut. The rim is a true tangent fillet: no sharp edge anywhere on the grip.",
  params: [
    P("width", "Width", 80, 200, 135),
    P("depth", "Depth", 30, 90, 60),
    deg("tilt", "Lip incut", 0, 30, 15),
    frac("wall", "Rim height", 0.1, 1, 0.74),
    frac("thick", "Shell thickness", 0.08, 0.45, 0.3),
    frac("spread", "Screw spread", 0.1, 0.48, 0.38),
    frac("pairH", "Screw pair height", 0.1, 1.1, 0.7),
    frac("centreH", "Centre screw height", 0, 0.6, 0.16),
  ],
  build(v, screws) {
    const wall = 0.16 + 0.19 * v.wall;
    const t = rad(v.tilt), shift = 0.32;
    const R = v.depth / (1 - sin(t) + shift);
    const ky = v.width / (2 * R);
    const body = cutWall(
      revolve(bowlProfile(wall, v.thick)).scale(R).rotate([0, -v.tilt, 0]).translate([shift * R, 0, 0]).scale([1, ky, 1]),
    );
    const cuts: Manifold[] = [];
    if (screws) {
      const cfx = (shift - sin(t)) * R, cfz = cos(t) * R;
      for (const [fy, fz] of [[-v.spread, v.pairH], [v.spread, v.pairH], [0, v.centreH]]) {
        const y = fy * v.width, z = fz * R;
        const sx = cfx + R * sqrt(max(0.05, 1 - (y / (ky * R)) ** 2 - ((z - cfz) / R) ** 2));
        cuts.push(sunkScrew(min(sx, 38), y, z, 3.5));
      }
    }
    return { body, screws: cuts };
  },
};

// ------------------------------------------------------------- pinch blade

const pinch: HoldDef = {
  id: "pinch",
  name: "Pinch blade",
  group: "Parametric",
  blurb: "A tapered vertical blade to squeeze, thumb on one side and fingers on the other, with optional dishes scooped into each face.",
  params: [
    P("width", "Width", 30, 90, 53),
    P("depth", "Depth", 40, 100, 67),
    P("height", "Height", 60, 160, 106),
    deg("angle", "Side angle (90 = prism)", 65, 90, 78),
    P("dish", "Dish depth", 0, 10, 4),
    P("dishLength", "Dish length", 20, 80, 47),
    P("dishHeight", "Dish height", 20, 90, 58),
    P("radius", "Rounding", 6, 24, 14),
  ],
  build(v, screws) {
    const shrink = (0.7 * v.depth) / Math.tan(rad(min(v.angle, 89.9)));
    const tipW = max(v.width - 2 * shrink, 0.2 * v.width);
    const tipH = max(v.height - 2 * shrink, 0.2 * v.height);
    const root = roundedBox(v.depth * 0.6, v.width, v.height, v.radius);
    const tip = roundedBox(v.depth * 0.5, tipW, tipH, v.radius).translate([v.depth * 0.75, 0, 0]);
    let body = cutWall(hull([root, tip]));
    if (v.dish > 0) {
      const ry = max(2.5 * v.dish, 6);
      const xc = 0.55 * v.depth;
      const ySide = v.width / 2 - (max(0, xc - 0.3 * v.depth) / (0.7 * v.depth)) * max(0, (v.width - tipW) / 2);
      const cut = ellipsoid([v.dishLength / 2, ry, v.dishHeight / 2], [xc, ySide + ry - v.dish, 0]);
      body = body.subtract(cut).subtract(cut.mirror([0, 1, 0]));
    }
    const cuts = screws ? [sunkScrew(v.depth * 0.55, 0, v.height * 0.22, 0), sunkScrew(v.depth * 0.55, 0, -v.height * 0.22, 0)] : [];
    return { body, screws: cuts };
  },
};

// ------------------------------------------------------------- palm sloper

function bellProfile(sigma = 0.5, rim = 0.08, steps = 48): [number, number][] {
  const E = exp(-1 / sigma ** 2);
  const fr = rim * 0.9;
  const bell = (r: number) => rim + ((1 - rim) * (exp(-((r / sigma) ** 2)) - E)) / (1 - E);
  const r1 = 1 - fr;
  const pts: [number, number][] = [];
  for (let i = 0; i <= steps; i++) pts.push([(i / steps) * r1, bell((i / steps) * r1)]);
  const m = (bell(r1) - bell(r1 - 1e-4)) / 1e-4;
  const n = hypot(1, m);
  const cx = r1 + (fr * m) / n, cz = bell(r1) - fr / n;
  const fEnd = atan2(bell(r1) - cz, r1 - cx);
  const k = Math.floor(steps / 3);
  for (let i = 1; i <= k; i++) {
    const f = fEnd * (1 - i / k);
    pts.push([cx + fr * cos(f), cz + fr * sin(f)]);
  }
  pts.push([cx + fr, 0], [0, 0]);
  return pts;
}

const sloper: HoldDef = {
  id: "sloper",
  name: "Palm sloper",
  group: "Parametric",
  blurb: "A bell-shaped friction dome: a convex crown rolling into a concave skirt at the wall, like a resin sloper rather than half an egg.",
  params: [
    P("width", "Width", 90, 220, 160),
    P("depth", "Depth", 30, 90, 55),
    P("height", "Height", 70, 180, 125),
    frac("sigma", "Crown tightness", 0.3, 0.8, 0.5),
  ],
  build(v, screws) {
    const body = resize(revolve(bellProfile(v.sigma)).rotate([0, 90, 0]), [v.depth, v.width, v.height]);
    const cuts: Manifold[] = [];
    if (screws) {
      const E = exp(-1 / v.sigma ** 2);
      for (const [fy, fz] of [[0, 0.24], [-0.26, -0.18], [0.26, -0.18]]) {
        const rho = min(0.98, hypot(2 * fy, 2 * fz));
        const g = (exp(-((rho / v.sigma) ** 2)) - E) / (1 - E);
        cuts.push(sunkScrew(v.depth * (0.08 + 0.92 * g), fy * v.width, fz * v.height, 3));
      }
    }
    return { body: cutWall(body, 0), screws: cuts };
  },
};

// ------------------------------------------------------------------ pocket

const pocket: HoldDef = {
  id: "pocket",
  name: "Doughnut pocket",
  group: "Parametric",
  blurb: "A squashed torus sliced flat at the wall: slot one or two fingers through and pull on the rounded lower rim.",
  params: [
    P("width", "Width", 70, 170, 112),
    P("height", "Height", 50, 130, 88),
    P("depth", "Depth", 15, 55, 30),
    P("tube", "Ring thickness", 20, 50, 36),
    P("lift", "Wall lift", 0, 55, 2),
  ],
  build(v, screws) {
    const lift = min(max(v.lift, 1.5), v.tube / 2 - 4);
    const rc0 = (v.width - v.tube) / 2;
    const circle = K().CrossSection.circle(v.tube / 2, 48).translate([rc0, 0]);
    const ring = resize(K().Manifold.revolve(circle, 64).rotate([0, 90, 0]), [(v.depth - lift) * 2, v.width, v.height]).translate([lift, 0, 0]);
    const cuts: Manifold[] = [];
    if (screws) {
      const squish = v.height / v.width;
      for (const ang of [90, 210, 330]) {
        const a = rad(ang);
        cuts.push(sunkScrew(v.depth, rc0 * cos(a), rc0 * sin(a) * squish, 3.5));
      }
    }
    return { body: cutWall(ring), screws: cuts };
  },
};

// ------------------------------------------------------------------ sphere

const ball: HoldDef = {
  id: "ball",
  name: "Sphere",
  group: "Parametric",
  blurb: "A smooth ellipsoid dome lifted off the wall before the flat cut, so it bulges over its own base.",
  params: [
    P("width", "Width", 40, 200, 90),
    P("length", "Length (off wall)", 20, 100, 60),
    P("height", "Height", 40, 200, 90),
    P("lift", "Wall lift", 0, 40, 2),
  ],
  build(v, screws) {
    const lift = min(max(v.lift, 0), 0.42 * v.length);
    const ax = v.length - lift;
    const body = cutWall(ellipsoid([ax, v.width / 2, v.height / 2], [lift, 0, 0], 64));
    const cuts: Manifold[] = [];
    if (screws)
      for (const [fy, fz] of [[0, 0.26], [-0.26, -0.2], [0.26, -0.2]]) {
        const s = 1 - (2 * fy) ** 2 - (2 * fz) ** 2;
        cuts.push(sunkScrew(min(lift + ax * sqrt(max(0.05, s)), 38), fy * v.width, fz * v.height, 3));
      }
    return { body, screws: cuts };
  },
};

// -------------------------------------------------------------------- rung

const rung: HoldDef = {
  id: "rung",
  name: "Rung",
  group: "Parametric",
  blurb: "A rounded bar to wrap the whole hand over. Shifted off the wall before the flat cut, so the underside is incut.",
  params: [
    P("width", "Width (along wall)", 40, 220, 140),
    P("depth", "Depth (off wall)", 15, 60, 35),
    P("diameter", "Bar thickness", 25, 65, 45),
  ],
  build(v, screws) {
    const r = v.diameter / 2;
    const lift = max(min(v.depth - r, 0.8 * r), -0.45 * r);
    const Ly = max(v.width - 2 * r, 0.15 * v.width);
    const body = cutWall(hull([sphere(r, [lift, -Ly / 2, 0]), sphere(r, [lift, Ly / 2, 0])]));
    const cuts: Manifold[] = [];
    if (screws) {
      const ys = [Ly * 0.3, -Ly * 0.3];
      if (v.width > 200) ys.push(Ly * 0.46, -Ly * 0.46);
      for (const y of ys) cuts.push(sunkScrew(min(lift + r, 38), y, 0, 3));
    }
    return { body, screws: cuts };
  },
};

// -------------------------------------------------------------------- wave

const wave: HoldDef = {
  id: "wave",
  name: "Scalloped rail",
  group: "Parametric",
  blurb: "A crimp rail with a row of finger seats scooped into the top, so each finger settles into its own dish.",
  params: [
    P("width", "Width", 60, 220, 150),
    P("depth", "Depth", 10, 35, 26),
    P("height", "Height", 18, 50, 34),
    P("radius", "Lip radius", 3, 16, 10),
    P("scallops", "Finger seats", 1, 8, 4, 1, ""),
    P("scoop", "Seat depth", 0, 10, 6),
  ],
  build(v, screws) {
    const box = roundedBox(2 * v.depth, v.width, 2 * v.height, v.radius);
    let body = cutWall(box).trimByPlane([0, 0, 1], 0);
    const n = max(Math.round(v.scallops), 1);
    if (v.scoop > 0) {
      const usable = v.width * 0.82, spacing = usable / n;
      const sy = spacing * 0.42, sx = v.depth * 0.95, sz = max(v.scoop * 1.9, 7);
      const xc = v.depth * 0.72;
      const cutters: Manifold[] = [];
      for (let i = 0; i < n; i++) {
        const y = -usable / 2 + (i + 0.5) * spacing;
        cutters.push(ellipsoid([sx / 2, sy, sz], [xc, y, v.height + sz - v.scoop]));
      }
      body = body.subtract(union(cutters));
    }
    // the rail lives in z >= 0: centre it on the origin like the others
    body = body.translate([0, 0, -v.height / 2]);
    const cuts: Manifold[] = [];
    if (screws) {
      const z = v.height * 0.32 - v.height / 2;
      const ys = [v.width / 4, -v.width / 4];
      if (v.width > 200) ys.push(v.width / 2.6, -v.width / 2.6);
      for (const y of ys) cuts.push(sunkScrew(v.depth * 0.85, y, z, 2));
    }
    return { body, screws: cuts };
  },
};

// ------------------------------------------------------------------- hueco

const hueco: HoldDef = {
  id: "hueco",
  name: "Hueco",
  group: "Parametric",
  blurb: "A dome with a bowl scooped out through its face. The cave is wider inside than at the mouth, so it is incut all the way round.",
  params: [
    P("width", "Width", 80, 180, 120),
    P("height", "Height", 70, 170, 110),
    P("depth", "Depth", 30, 75, 45),
    P("rim", "Rim thickness", 6, 22, 12),
    frac("mouth", "Mouth size", 0.1, 1, 0.55),
  ],
  build(v, screws) {
    const lift = min(3, 0.35 * v.depth);
    const ax = v.depth - lift, ay = v.width / 2, az = v.height / 2;
    const outer = ellipsoid([ax, ay, az], [lift, 0, 0], 64);
    const rim = max(min(v.rim, 0.4 * ay, 0.4 * az), 3);
    const mouth = min(max(v.mouth, 0.1), 1);
    const back = 0.16 * ax;
    const cx = (0.55 + 0.32 * mouth) * ax, cy = ay - rim, cz = az - rim;
    const cav = ellipsoid([cx, cy, cz], [lift + back + cx, 0, 0], 64);
    const body = cutWall(outer.subtract(cav));
    const cuts: Manifold[] = [];
    if (screws)
      for (const ang of [90, 210, 330]) {
        const a = rad(ang);
        const y = 1.08 * cy * cos(a), z = 1.08 * cz * sin(a);
        const t = max(1 - (y / ay) ** 2 - (z / az) ** 2, 0.04);
        cuts.push(sunkScrew(min(lift + ax * sqrt(t), 40), y, z, 3.5));
      }
    return { body, screws: cuts };
  },
};

export const PARAMETRIC: HoldDef[] = [jug, pinch, sloper, pocket, ball, rung, wave, hueco];
