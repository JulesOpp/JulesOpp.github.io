"""Volumes and extents of the ported holds at their defaults (screws off),
rendered by the desktop app's own OpenSCAD pipeline, as a JSON reference
for test/parity.test.ts.

    python3 scripts/parity_reference.py "<path to Climbing Hold Generator>" > test/parity-reference.json
"""
import json
import sys
import tempfile

sys.path.insert(0, sys.argv[1])
import trimesh  # noqa: E402
from holdlib.registry import HOLD_SPECS, build_model  # noqa: E402
from render.pipeline import render_stl  # noqa: E402

# web id -> desktop name
PAIRS = {
    "jug": "Jug (bucket)", "pinch": "Pinch", "sloper": "Sloper", "ball": "Sphere",
    "pocket": "Pocket", "rung": "Rung (bar)", "wave": "Wave (scalloped rail)", "hueco": "Hueco (cave)",
    "gradient-edge": "Gradient edge", "rock-flake": "Rock flake", "tufa-rib": "Tufa rib",
    "wave-rail": "Wave rail", "cobble": "Cobble sloper (BOSL2)", "volume-wedge": "Volume wedge (BOSL2)",
    "alien-face": "Alien face (BOSL2)", "system-edge": "System edge",
}

out = {}
with tempfile.TemporaryDirectory() as wd:
    for web, name in PAIRS.items():
        if name not in HOLD_SPECS:
            print(f"skip {name!r}: not in HOLD_SPECS", file=sys.stderr)
            continue
        spec = HOLD_SPECS[name]
        vals = {p.key: (p.map(p.default) if getattr(p, "map", None) else p.default) for p in spec.params}
        stl = render_stl(build_model(name, vals, screws=False), 48, wd, name=web)
        m = trimesh.load(stl)
        out[web] = {"volume": float(abs(m.volume)), "extents": [float(x) for x in m.extents]}
        print(web, round(out[web]["volume"]), [round(x, 1) for x in m.extents], file=sys.stderr)
print(json.dumps(out, indent=1))
