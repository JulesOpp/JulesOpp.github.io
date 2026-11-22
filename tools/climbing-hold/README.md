# Climbing Hold Generator (web)

Source for `blog/climbing-hold-generator.html`. It is a browser port of the desktop
[Climbing Hold Generator](https://github.com/JulesOpp/Climbing-Hold-Generator), and it is
built into `assets/blog/climbing-hold/` (the built files are committed, as Pages has no
build step).

```
npm install
npm test          # vitest: every hold at default/min/max, textures, desktop parity, STEP parity
npm run typecheck
npm run build     # -> assets/blog/climbing-hold/{chg.js, chgw.js, chstepw.js, *.wasm, LICENSES.txt}
```

## Layout

- `src/kernel.ts`: manifold-3d setup and an `arena()` that frees all the WASM objects a
  build creates. It also has the primitives, the screw negative (Ø13 counterbore, 90°
  countersink, Ø5 shank) and loft/bezier/squircle.
- `src/holds/`: the holds.
  - `parametric.ts`: ported from `holdlib/parametric.py`.
  - `sculpted.ts`: rebuilds of the app's own `.scad` templates. Template frame →
    `rotate([90,0,90])`.
  - `machinable.ts`: the mesh preview of the CNC hold.
  - The Jeremy Till holds (spherical edge, foot chip, star, crimp rail, and the
    "Original geometry" templates) are deliberately **not** ported: there is no license
    for them. The one exception is the System edge, whose template credits Jeremy Till
    and Jules Oppenheim jointly; it is included at Jules's request.
- `src/edge/params.ts`: the machinable edge's parameters, sphere centres and screw axes.
  The mesh preview and the B-rep both use it.
- `src/brep/machinable-edge.ts`: the STEP solid in OpenCascade (replicad). It follows
  `brep/system_edge.py`: hull of the core points, an Arc-join offset, the flat back, the
  plinth with its crease fillet, the side flats with their fillet, and the screws.
- `src/texture.ts`: port of `geometry/texture.py`. Divots and grit use a seeded PRNG, so
  they do not match the Python patterns exactly.
- `src/workers/`:
  - `mesh.ts`: builds every preview and STL.
  - `step.ts`: OCCT. It is loaded only when someone clicks *Download STEP*, because the
    wasm is 22 MB (about 7 MB gzipped).
- `scripts/decimate_lafaye.py`: makes `assets/blog/climbing-hold/lafaye/`. It takes six of
  Jérémy LAFAYE's CC BY-SA 4.0 scans, decimates them to 14k triangles, re-orients them so
  the mounting face is x = 0, and writes the manifest and licence.
- `scripts/prepare_wupinstick.py`: makes `assets/blog/climbing-hold/wupinstick/`. It
  re-orients Wupinstick's *Big incut* (CC BY 4.0, thingiverse.com/thing:2398720) the same
  way, but does not decimate it, and writes the manifest and licence.
- `scripts/parity_reference.py`: renders the ported holds with the desktop app's own
  OpenSCAD pipeline. Its output is `test/parity-reference.json`.

## Checks behind the numbers

- `test/brep.test.ts` reproduces `brep/system_edge.py` exactly: 36593.9 mm³ and 19 faces
  for the plain hull, and 150759 mm³, 44 faces and 628 mm² flats for the flats+plinth
  preset.
- Re-importing the default STEP into cadquery gives one valid solid in millimetres with
  55 faces: 29 cylinders, 11 planes, 4 spheres, 3 cones, 2 tori and 6 B-spline fillet
  blends.
