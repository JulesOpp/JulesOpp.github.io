import { beforeAll, describe, expect, it } from "vitest";
import { loadOC } from "./occ";
import { machinableEdgeSolid } from "../src/brep/machinable-edge";
import { EDGE_DEFAULTS, type EdgeParams } from "../src/edge/params";

// Reference volumes from the desktop app's brep/system_edge.py
// (cadquery), computed at the same core geometry with screws off.
const BASE: EdgeParams = {
  width: 100, depth: 20, lipHeight: 10.569620253164558, lipLength: 70, lowerWidth: 50,
  drop: 28, edgeRadius: 5.329113924050633, sideFlat: 0, flatFillet: 0, plinth: 0,
  plinthFillet: 0, screwSeat: 8.4, screws: false,
};
const FLATS: EdgeParams = {
  ...BASE, width: 160, sideFlat: 30, flatFillet: 2, depth: 40, drop: 45, plinth: 10,
  plinthFillet: 6, lipHeight: 10.5696, edgeRadius: 5.329,
};

beforeAll(loadOC, 60_000);

describe("machinable edge B-rep", () => {
  it("matches brep/ on the plain hull (36593.9 mm3, 19 faces)", () => {
    const r = machinableEdgeSolid(BASE);
    expect(r.volume).toBeCloseTo(36593.9, 0);
    expect(r.faces).toBe(19);
  });

  it("matches brep/ with plinth and side flats", () => {
    const r = machinableEdgeSolid(FLATS);
    // brep/: 150759.4 mm3, 44 faces; README: 628 mm2 per flat
    expect(Math.abs(r.volume - 150759.4) / 150759.4).toBeLessThan(1e-4);
    expect(r.faces).toBe(44);
    expect(r.flatArea).toBeCloseTo(628.3, 0);
    expect(r.bbox[1]).toBeCloseTo(100, 3); // across the jaws
  }, 60_000);

  it("builds the shipped default with screws as one solid", () => {
    const r = machinableEdgeSolid(EDGE_DEFAULTS);
    expect(r.volume).toBeGreaterThan(0);
    expect(r.flatArea).toBeGreaterThan(400);
  }, 60_000);
});
