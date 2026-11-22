// The machinable edge: a hold designed to be cut from aluminium.
//
// The preview here is a manifold-3d mesh of the same construction the
// STEP builder (src/brep/machinable-edge.ts) makes in OpenCascade: hull
// of equal spheres, flat back, optional plinth, optional vice flats,
// three countersunk screws. The two small fillets (plinth crease, flat
// edges) exist only in the STEP; on the mesh they would change a few
// cubic millimetres and are left out.

import { cutWall, hull, K, screwNegative, sphere, type Manifold } from "../kernel";
import { coreSpheres, EDGE_DEFAULTS, screwAxes, validateEdge, type EdgeParams } from "../edge/params";
import { P, type HoldDef, type Values } from "../spec";

export function edgeParams(v: Values, screws: boolean): EdgeParams {
  return {
    width: v.width, depth: v.depth, lipHeight: v.lipHeight, lipLength: v.lipLength,
    lowerWidth: v.lowerWidth, drop: v.drop, edgeRadius: v.edgeRadius,
    sideFlat: v.sideFlat, flatFillet: v.flatFillet, plinth: v.plinth,
    plinthFillet: v.plinthFillet, screwSeat: v.screwSeat, screws,
  };
}

/** Mesh of the mounting-face section, extruded back by `h` (x in [-h, 0]). */
function plinthOf(body: Manifold, h: number): Manifold {
  const section = body.rotate([0, -90, 0]).slice(1e-4); // template x -> z
  return section.extrude(h).translate([0, 0, -h]).rotate([0, 90, 0]);
}

export function edgeMesh(p: EdgeParams): { body: Manifold; screws: Manifold[] } {
  validateEdge(p);
  const balls = coreSpheres(p).map((c) => sphere(p.edgeRadius, c, 64));
  let body = cutWall(hull(balls));
  if (p.plinth > 0) body = body.add(plinthOf(body, p.plinth)).translate([p.plinth, 0, 0]);
  if (p.sideFlat > 0) {
    const keep = p.width / 2 - p.sideFlat;
    body = body.intersect(K().Manifold.cube([1000, 2 * keep, 1000], true).translate([0, 0, 0]));
  }
  const screws = p.screws
    ? screwAxes(p).map(([x0, y, z]) => screwNegative(x0, y, z, p.depth + p.plinth + 10))
    : [];
  return { body, screws };
}

const d = EDGE_DEFAULTS;
export const MACHINABLE_EDGE: HoldDef = {
  id: "machinable-edge",
  name: "Machinable edge",
  group: "Machinable",
  step: true,
  blurb:
    "An incut edge built only from planes, cylinders and spheres, with vice flats on both ends and a straight base, so it can be milled from aluminium. Export it as STEP for a CNC shop.",
  params: [
    P("width", "Width before flats", 60, 200, d.width),
    P("sideFlat", "Vice flat (trim per side)", 0, 50, d.sideFlat, 1, " mm", "Trims each rounded end back to a flat a plain vice can grip. 0 = rounded ends."),
    P("depth", "Lip stand-off", 12, 50, d.depth),
    P("lipHeight", "Lip height", 0, 20, d.lipHeight, 0.5),
    P("lipLength", "Lip length", 30, 120, d.lipLength),
    P("lowerWidth", "Base width", 30, 90, d.lowerWidth),
    P("drop", "Drop", 15, 60, d.drop),
    P("edgeRadius", "Edge radius", 2, 8, d.edgeRadius, 0.1, " mm", "Every edge and corner is this ball: the smallest convex radius the shop has to cut."),
    P("plinth", "Plinth (straight base)", 0, 20, d.plinth, 1, " mm", "Extrudes the mounting face straight back: bigger flats and an easy base to machine."),
    P("plinthFillet", "Plinth crease fillet", 0, 12, d.plinthFillet, 0.5, " mm", "Rounds the one concave crease. A bigger radius lets the shop finish with a bigger ball cutter."),
    P("flatFillet", "Flat edge fillet", 0, 3, d.flatFillet, 0.5),
    P("screwSeat", "Screw seat (metal behind head)", 4, 30, d.screwSeat, 0.5),
  ],
  build(v, screws) {
    return edgeMesh(edgeParams(v, screws));
  },
};
