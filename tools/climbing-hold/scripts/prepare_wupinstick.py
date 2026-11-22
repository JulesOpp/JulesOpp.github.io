"""Re-orient Wupinstick's "Big incut" (CC BY 4.0) for the web page.

    python3 -I scripts/prepare_wupinstick.py "<path to Climbing Hold Generator>/stl_holds"

Same treatment as decimate_lafaye.py, minus the decimation (the file is
already small): largest flat face to x = 0, centred, checked closed with
manifold3d, written to assets/blog/climbing-hold/wupinstick/ with a
manifest and licence note.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import manifold3d as m3d  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402
from decimate_lafaye import mounting_rotation  # noqa: E402

OUT = HERE.parent.parent.parent / "assets" / "blog" / "climbing-hold" / "wupinstick"
AUTHOR = "Wupinstick"
AUTHOR_URL = "https://www.thingiverse.com/wupinstick/designs"
SOURCE = "https://www.thingiverse.com/thing:2398720"
LICENSE = "CC BY 4.0"
LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"
# (folder, file, name, default rotation on the wall in degrees)
PICKS = [("Wupinstick", "Big_incut.stl", "Big incut", 90)]


def main(src):
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for folder, name, title, rotation in PICKS:
        mesh = trimesh.load(Path(src) / folder / name, force="mesh")
        mesh.merge_vertices()
        if mesh.volume < 0:
            mesh.invert()
        mesh.apply_transform(mounting_rotation(mesh))
        lo, hi = mesh.bounds
        mesh.apply_translation([-lo[0], -(lo[1] + hi[1]) / 2, -(lo[2] + hi[2]) / 2])
        gl = m3d.Mesh(vert_properties=np.ascontiguousarray(mesh.vertices, np.float32),
                      tri_verts=np.ascontiguousarray(mesh.faces, np.uint32))
        gl.merge()
        man = m3d.Manifold(gl)
        if man.status() != m3d.Error.NoError or man.volume() <= 0 or man.genus() < 0:
            raise SystemExit(f"{name}: not a closed solid ({man.status()})")
        out = man.to_mesh()
        clean = trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=False)
        slug = Path(name).stem.lower().replace("_", "-")
        clean.export(OUT / f"{slug}.stl")
        ext = clean.bounds[1] - clean.bounds[0]
        print(f"{slug}: {len(clean.faces)} faces, {(OUT / f'{slug}.stl').stat().st_size/1024:.0f} KB, "
              f"{ext.round(1)} mm, {man.volume()/1000:.1f} cm3, genus {man.genus()}")
        manifest.append(dict(id=f"wupinstick-{slug}", file=f"{slug}.stl", name=title, source=SOURCE,
                             author=AUTHOR, authorUrl=AUTHOR_URL, license=LICENSE, licenseUrl=LICENSE_URL,
                             original=name, rotation=rotation))
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    (OUT / "LICENSE.txt").write_text(
        f"The .stl files in this folder are re-oriented copies of models from \"Rock Climbing Holds\"\n"
        f"by {AUTHOR} ({SOURCE}), licensed under the Creative Commons Attribution 4.0\n"
        f"International licence ({LICENSE_URL}).\n"
        f"Changes: rotated so the mounting face is x = 0, centred, and re-meshed as one closed solid.\n"
        + "".join(f"  {p[2]}: {p[1]}\n" for p in PICKS))


if __name__ == "__main__":
    main(sys.argv[1])
