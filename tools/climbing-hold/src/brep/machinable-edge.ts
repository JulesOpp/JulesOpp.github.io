// The machinable edge as an analytic B-rep solid, for STEP export.
//
// Port of the construction in the desktop app's brep/system_edge.py:
// the convex hull of equal spheres is exactly the convex hull of their
// centres grown by the radius, and growing a convex polyhedron by a
// radius is a single OCCT offset (Arc join). The result is real planes,
// cylinders and spheres -- the clean surfaces a CAM package wants --
// where an STL would only hand the shop facets.
//
// Axes: +x stands off the wall, y runs along the wall, z is up. The
// mounting face is the x = 0 plane.
//
// Everything here needs OpenCascade loaded (replicad's setOC) and is
// only ever run inside the STEP worker or the node tests.

import {
  basicFaceExtrusion,
  getOC,
  makeBox,
  makePolygon,
  makeSolid,
  measureArea,
  measureVolume,
  revolution,
  Vector,
  type Edge,
  type Face,
  type Shape3D,
} from "replicad";
import { coreHullFaces } from "../edge/hull";
import {
  type EdgeParams,
  coreSpheres,
  screwAxes,
  SCREW,
  validateEdge,
} from "../edge/params";

const TOL = 1e-6;

function single(shape: Shape3D, what: string): Shape3D {
  const solids = (shape as unknown as { solids?: Shape3D[] }).solids;
  if (solids && solids.length !== 1)
    throw new Error(`${what}: expected one solid, got ${solids.length} -- check the parameters for a self-intersecting hull`);
  return solids && solids.length === 1 ? solids[0] : shape;
}

/** Grow a convex solid by `r`: Minkowski sum with a ball. */
function grow(solid: Shape3D, r: number): Shape3D {
  const oc = getOC();
  const mk = new oc.BRepOffsetAPI_MakeOffsetShape();
  const progress = new oc.Message_ProgressRange();
  mk.PerformByJoin(
    solid.wrapped, r, TOL,
    oc.BRepOffset_Mode.BRepOffset_Skin,
    false, false,
    oc.GeomAbs_JoinType.GeomAbs_Arc,
    true,
    progress,
  );
  if (!mk.IsDone()) throw new Error("OpenCascade could not offset the hull");
  // replicad's cast wraps a raw TopoDS_Shape in the right class
  const shape = castShape(mk.Shape());
  mk.delete();
  progress.delete();
  return shape;
}

// replicad exports `cast`; imported lazily to keep the type surface small
import { cast } from "replicad";
const castShape = (s: unknown) => cast(s as never) as unknown as Shape3D;

function isWallFace(f: Face): boolean {
  if (f.geomType !== "PLANE") return false;
  const n = f.normalAt();
  return Math.abs(f.center.x) < TOL && Math.abs(Math.abs(n.x) - 1) < TOL;
}

function mountingFace(shape: Shape3D): Face {
  const found = shape.faces.filter(isWallFace);
  if (found.length !== 1) throw new Error(`expected one face on the wall plane, found ${found.length}`);
  return found[0];
}

function facesOfEdge(shape: Shape3D, edge: Edge): Face[] {
  return shape.faces.filter((f) => f.edges.some((e) => e.isSame(edge)));
}

/** Edges lying in x = xAt whose two faces actually meet at a corner. */
function sharpEdgesInPlane(shape: Shape3D, xAt: number, minAngleDeg = 0.5): Edge[] {
  const out: Edge[] = [];
  for (const edge of shape.edges) {
    const onPlane = [0, 0.25, 0.5, 0.75, 1].every((t) => Math.abs(edge.pointAt(t).x - xAt) < 1e-5);
    if (!onPlane) continue;
    const faces = facesOfEdge(shape, edge);
    if (faces.length < 2) continue;
    const mid = edge.pointAt(0.5);
    const d = Math.abs(faces[0].normalAt(mid).dot(faces[1].normalAt(mid)));
    if ((Math.acos(Math.min(1, d)) * 180) / Math.PI >= minAngleDeg) out.push(edge);
  }
  return out;
}

function filletEdges(shape: Shape3D, edges: Edge[], r: number, what: string, hint: string): Shape3D {
  if (!edges.length || r <= 0) return shape;
  let rounded: Shape3D;
  try {
    rounded = shape.fillet((e: Edge) => (edges.some((x) => x.isSame(e)) ? r : null));
  } catch (err) {
    throw new Error(`OpenCascade could not round the ${what} at r = ${r} mm. ${hint} (${(err as Error).message ?? err})`);
  }
  return single(rounded, `${what} rounding`);
}

/** Extrude the mounting face straight back by `height`, round the crease. */
function wallPlinth(hold: Shape3D, height: number, fillet: number): Shape3D {
  const back = mountingFace(hold);
  const prism = basicFaceExtrusion(back, new Vector([-height, 0, 0]));
  let body = single(hold.fuse(prism), "wall plinth");
  if (fillet > 0) {
    body = filletEdges(body, sharpEdgesInPlane(body, 0), fillet, "plinth junction", "Try a smaller plinth fillet.");
  }
  return body.translate(height, 0, 0);
}

/** Trim the rounded ends to two parallel vice flats and break their edges. */
function sideFlats(hold: Shape3D, width: number, sideFlat: number, fillet: number, span: number): Shape3D {
  const keepY = width / 2 - sideFlat;
  if (keepY <= 0) throw new Error("side flat must leave material behind");
  const slab = makeBox([-span / 2, -keepY, -span / 2], [span / 2, keepY, span / 2]);
  const trimmed = single(hold.intersect(slab), "side flats");
  if (fillet <= 0) return trimmed;
  const flatEdges: Edge[] = [];
  let shortest = Infinity;
  for (const f of trimmed.faces) {
    if (f.geomType !== "PLANE") continue;
    if (Math.abs(Math.abs(f.normalAt().y) - 1) > TOL || Math.abs(Math.abs(f.center.y) - keepY) > TOL) continue;
    for (const e of f.edges) {
      flatEdges.push(e);
      shortest = Math.min(shortest, e.length);
    }
  }
  return filletEdges(trimmed, flatEdges, fillet, "side flats",
    `The shortest segment of the flat's outline is ${shortest.toFixed(2)} mm: try a flat fillet below that, or a bigger side flat.`);
}

/**
 * One screw's removed metal as a single revolve: counterbore, 90 deg
 * countersink, shank. The profile lives in the half-plane through the
 * screw axis (parallel to x), so every surface is an exact cylinder or cone.
 */
export function screwCutter(x0: number, y: number, z: number, reach: number): Shape3D {
  const { cbDia, shankDia, shankLen } = SCREW;
  const cone = (cbDia - shankDia) / 2; // 45 deg half angle
  const rc = cbDia / 2, rs = shankDia / 2;
  // (x, radial) profile, radial measured along +z from the axis
  const prof: [number, number][] = [
    [x0 - shankLen, 0],
    [x0 - shankLen, rs],
    [x0 - cone, rs],
    [x0, rc],
    [x0 + reach, rc],
    [x0 + reach, 0],
  ];
  const face = makePolygon(prof.map(([x, r]) => [x, y, z + r]));
  return revolution(face, [0, y, z], [1, 0, 0], 360) as Shape3D;
}

export interface EdgeSolidReport {
  solid: Shape3D;
  volume: number; // mm^3
  area: number; // mm^2
  faces: number;
  flatArea: number; // per side, mm^2 (0 without side flats)
  bbox: [number, number, number];
}

export function machinableEdgeSolid(p: EdgeParams): EdgeSolidReport {
  validateEdge(p);
  const r = p.edgeRadius;
  const faces = coreHullFaces(coreSpheres(p)).map((ring) => makePolygon(ring));
  const core = makeSolid(faces) as unknown as Shape3D;
  const grown = grow(core, r);

  const span = 4 * (Math.abs(p.drop) + Math.abs(p.lipHeight) + p.depth + p.width);
  const keep = makeBox([0, -span / 2, -span / 2], [span, span / 2, span / 2]);
  let hold = single(grown.intersect(keep), "flat back");

  if (p.plinth > 0) hold = wallPlinth(hold, p.plinth, p.plinthFillet);
  if (p.sideFlat > 0) hold = sideFlats(hold, p.width, p.sideFlat, p.flatFillet, span);
  if (p.screws) {
    const reach = p.depth + p.plinth + 10;
    for (const [x0, y, z] of screwAxes(p)) hold = hold.cut(screwCutter(x0, y, z, reach));
  }
  hold = single(hold, "finished hold");

  let flatArea = 0;
  if (p.sideFlat > 0) {
    const keepY = p.width / 2 - p.sideFlat;
    for (const f of hold.faces) {
      if (f.geomType === "PLANE" && f.center.y > 0 && Math.abs(f.center.y - keepY) < 1e-4 && Math.abs(Math.abs(f.normalAt().y) - 1) < 1e-6)
        flatArea += measureArea(f);
    }
  }
  const [lo, hi] = hold.boundingBox.bounds;
  return {
    solid: hold,
    volume: measureVolume(hold),
    area: measureArea(hold),
    faces: hold.faces.length,
    flatArea,
    bbox: [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]],
  };
}
