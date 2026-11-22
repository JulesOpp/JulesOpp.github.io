"""Decimate a handful of Jeremy LAFAYE's CC BY-SA 4.0 hold scans for the web page.

    python3 -I scripts/decimate_lafaye.py "<path to Climbing Hold Generator>/stl_holds"

Each pick is loaded, simplified to TARGET faces with fast_simplification
(quadric decimation), turned so its largest flat face is the mounting face
on x = 0 (the generator's house frame: +x off the wall, z up), centred on
the wall, checked watertight with manifold3d, and written as binary STL to
assets/blog/climbing-hold/lafaye/. A JSON manifest carries names and
attribution. The decimated meshes are adaptations of the originals and
stay under CC BY-SA 4.0.
"""

import json
import sys
from pathlib import Path

import fast_simplification
import manifold3d as m3d
import numpy as np
import trimesh

TARGET = 14000
HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent.parent / "assets" / "blog" / "climbing-hold" / "lafaye"
AUTHOR = "Jérémy LAFAYE"
AUTHOR_URL = "https://printables.com/@JeremyLAFAYE_2080610"
LICENSE = "CC BY-SA 4.0"
LICENSE_URL = "https://creativecommons.org/licenses/by-sa/4.0/"

PICKS = [
    # (source folder, file, display name, printables pack url)
    ("10x5-holds-sample-pack-model_files", "bh_10x5_01.stl", "10x5 hold #01", "https://printables.com/model/881409-10x5-holds-sample-pack"),
    ("10x5-holds-sample-pack-model_files", "bh_10x5_18.stl", "10x5 hold #18", "https://printables.com/model/881409-10x5-holds-sample-pack"),
    ("finger-pocket-hold-pack-model_files", "fh1_3_dec.stl", "Finger pocket macro #1", "https://printables.com/model/881406-finger-pocket-hold-pack"),
    ("jugs-sample-pack-model_files", "ch1_s.stl", "Jug #1 (S)", "https://printables.com/model/881407-jugs-sample-pack"),
    ("rock-sample-pack-model_files", "rockhold_t2_01.stl", "Rock hold T2-01", "https://printables.com/model/881408-rock-sample-pack"),
    ("rock-sample-pack-model_files", "rockhold_t1_04.stl", "Rock hold T1-04", "https://printables.com/model/881408-rock-sample-pack"),
]


def mounting_rotation(mesh):
    """Rotation taking the bbox face with the most coplanar area to x = 0,
    hold in +x."""
    lo, hi = mesh.bounds
    n, c, a = mesh.face_normals, mesh.triangles_center, mesh.area_faces
    best = None
    for axis in range(3):
        for sign, plane in ((-1, lo[axis]), (1, hi[axis])):
            on = (np.abs(c[:, axis] - plane) < 0.5) & (n[:, axis] * sign > 0.95)
            area = a[on].sum()
            if best is None or area > best[0]:
                best = (area, axis, sign)
    _, axis, sign = best
    outward = np.zeros(3)
    outward[axis] = sign
    # rotate the outward normal of the base onto -x
    return trimesh.geometry.align_vectors(outward, [-1.0, 0.0, 0.0])


def main(src):
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for folder, name, title, url in PICKS:
        mesh = trimesh.load(Path(src) / folder / name, force="mesh")
        ratio = max(0.0, 1.0 - TARGET / len(mesh.faces))
        v, f = fast_simplification.simplify(mesh.vertices, mesh.faces, ratio)
        mesh = trimesh.Trimesh(v, f, process=True)
        if mesh.volume < 0:          # some scans arrive wound inside-out
            mesh.invert()
        mesh.apply_transform(mounting_rotation(mesh))
        lo, hi = mesh.bounds
        mesh.apply_translation([-lo[0], -(lo[1] + hi[1]) / 2, -(lo[2] + hi[2]) / 2])
        gl = m3d.Mesh(vert_properties=np.ascontiguousarray(mesh.vertices, np.float32),
                      tri_verts=np.ascontiguousarray(mesh.faces, np.uint32))
        gl.merge()
        man = m3d.Manifold(gl)
        if man.status() != m3d.Error.NoError or man.volume() <= 0 or man.genus() < 0:
            raise SystemExit(f"{name}: not a closed solid after decimation ({man.status()})")
        out = man.to_mesh()
        clean = trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=False)
        slug = Path(name).stem.replace("_", "-")
        clean.export(OUT / f"{slug}.stl")
        size = (OUT / f"{slug}.stl").stat().st_size
        ext = clean.bounds[1] - clean.bounds[0]
        print(f"{slug}: {len(clean.faces)} faces, {size/1024:.0f} KB, {ext.round(1)} mm, {man.volume()/1000:.1f} cm3")
        manifest.append(dict(id=f"lafaye-{slug}", file=f"{slug}.stl", name=title, source=url,
                             author=AUTHOR, authorUrl=AUTHOR_URL, license=LICENSE, licenseUrl=LICENSE_URL,
                             original=name))
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    (OUT / "LICENSE.txt").write_text(
        f"The .stl files in this folder are decimated, re-oriented adaptations of climbing-hold\n"
        f"models by {AUTHOR} ({AUTHOR_URL}), published on Printables under the\n"
        f"Creative Commons Attribution-ShareAlike 4.0 International licence ({LICENSE_URL}).\n"
        f"Changes: mesh simplified to ~{TARGET} triangles, rotated so the mounting face is x = 0,\n"
        f"and centred. These adaptations are shared under the same licence. Sources:\n"
        + "".join(f"  {p[2]}: {p[3]} ({p[1]})\n" for p in PICKS))


if __name__ == "__main__":
    main(sys.argv[1])
