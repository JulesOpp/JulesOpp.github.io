// Parameters of the machinable edge, shared by the mesh preview
// (manifold-3d) and the STEP builder (OpenCascade), so the two can never
// disagree about where a sphere or a screw sits.

export interface EdgeParams {
  width: number; // overall length along the wall before flats, mm
  depth: number; // lip stand-off from the wall, mm
  lipHeight: number; // lip centre above the top corners, mm
  lipLength: number; // mm
  lowerWidth: number; // width of the bottom of the hold, mm
  drop: number; // how far the hold reaches below the top corners, mm
  edgeRadius: number; // the ball every edge and corner is rounded with, mm
  sideFlat: number; // trimmed off each end to leave vice flats (0 = none), mm
  flatFillet: number; // rounds the edge round each flat, mm
  plinth: number; // straight base extruded out of the wall face, mm
  plinthFillet: number; // rounds the crease where the base meets the body, mm
  screwSeat: number; // countersink bottom, forward of the wall face, mm
  screws: boolean;
}

/** The vice-friendly preset: flats, a 10 mm plinth, a 6 mm junction fillet. */
export const EDGE_DEFAULTS: EdgeParams = {
  width: 160,
  depth: 40,
  lipHeight: 10.5,
  lipLength: 70,
  lowerWidth: 50,
  drop: 45,
  edgeRadius: 5.3,
  sideFlat: 30,
  flatFillet: 2,
  plinth: 10,
  plinthFillet: 6,
  screwSeat: 18,
  screws: true,
};

/** #8 / M4 flat-head screw: counterbore, 90 deg countersink, shank. */
export const SCREW = { cbDia: 13, shankDia: 5, shankLen: 41.3 };

export function validateEdge(p: EdgeParams): void {
  if (p.edgeRadius <= 0) throw new Error("edge radius must be positive");
  if (p.depth <= p.edgeRadius) throw new Error("depth must exceed the edge radius");
  if (p.lowerWidth <= 2 * p.edgeRadius) throw new Error("lower width must exceed twice the edge radius");
  if (p.width <= 2 * p.edgeRadius) throw new Error("width must exceed twice the edge radius");
  if (p.sideFlat > 0 && p.width / 2 - p.sideFlat <= p.edgeRadius)
    throw new Error("side flats would cut away the whole hold");
  if (p.screws && p.screwSeat <= 2) throw new Error("screw seat must sit at least 2 mm forward of the wall face");
}

/** Sphere centres whose hull, grown by edgeRadius, is the hold. */
export function coreSpheres(p: EdgeParams): [number, number, number][] {
  const r = p.edgeRadius;
  const halfLip = Math.max(p.lipLength / 2 - r, 0);
  const pts: [number, number, number][] = [
    [0, p.width / 2 - r, 0],
    [0, -(p.width / 2 - r), 0],
    [0, p.lowerWidth / 2 - r, -(p.drop - r)],
    [0, -(p.lowerWidth / 2 - r), -(p.drop - r)],
    [p.depth - r, halfLip, p.lipHeight],
    [p.depth - r, -halfLip, p.lipHeight],
  ];
  const seen = new Set<string>();
  return pts.filter((q) => {
    const k = q.map((v) => v.toFixed(9)).join(",");
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
}

/**
 * Screw axes as (seat x, y, z). Three screws: a pair high up and one
 * low in the middle, spread across what is left between the flats. All
 * three seats sit at the same depth.
 */
export function screwAxes(p: EdgeParams): [number, number, number][] {
  const half = p.sideFlat > 0 ? p.width / 2 - p.sideFlat : p.width / 2;
  const pairY = Math.min(0.5 * half, p.width / 4);
  return [
    [p.screwSeat, pairY, -0.2 * p.drop],
    [p.screwSeat, -pairY, -0.2 * p.drop],
    [p.screwSeat, 0, -0.5 * p.drop],
  ];
}
