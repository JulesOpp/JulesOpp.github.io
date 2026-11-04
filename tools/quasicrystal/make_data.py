"""Write the QC Sandbox post's precomputed data block.

The browser cannot run the phase-diagram code (it needs SciPy and
about a second per diagram), so the six preset diagrams, their label
positions and the colour tables are computed here by the Python core
and pasted into the page between /*QCDATA:BEGIN*/ and /*QCDATA:END*/.

Run with the Quasicrystal repo's interpreter (branch build/app):

    ~/GitHub/Quasicrystal/.venv/bin/python \
        tools/quasicrystal/make_data.py blog/quasicrystal-sandbox.html

QC_REPO overrides where the Quasicrystal checkout is.
"""
import base64
import json
import os
import re
import subprocess
import sys

REPO = os.path.expanduser(os.environ.get("QC_REPO",
                                         "~/GitHub/Quasicrystal"))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from qcsandbox import colormaps, phasediagram  # noqa: E402
from qcsandbox.params import Q_PRESETS  # noqa: E402

NX = NY = 160
BEGIN, END = "/*QCDATA:BEGIN*/", "/*QCDATA:END*/"


def rle(labels):
    # (value, run) byte pairs, runs capped at 255: the regions are a
    # handful of blobs, so this is a few KB instead of 25 KB a preset.
    flat = labels.astype(np.uint8).ravel()
    out = bytearray()
    i = 0
    while i < len(flat):
        j = i
        while j < len(flat) and flat[j] == flat[i] and j - i < 255:
            j += 1
        out += bytes((flat[i], j - i))
        i = j
    return bytes(out)


def unrle(data, size):
    out = np.empty(size, np.uint8)
    pos = 0
    for v, n in zip(data[::2], data[1::2], strict=True):
        out[pos:pos + n] = v
        pos += n
    assert pos == size
    return out


def build():
    diagrams = {}
    for name, q in Q_PRESETS.items():
        d = phasediagram.cached(name, q, NX, NY)
        packed = rle(d.labels)
        back = unrle(packed, NX * NY).reshape(NY, NX)
        assert np.array_equal(back, d.labels), name
        diagrams[name] = dict(
            q=q, nx=NX, ny=NY, names=list(d.names), note=d.note,
            colours=[colormaps.region_colour(n) for n in d.names],
            rle=base64.b64encode(packed).decode(),
            centres={k: [round(v, 5) for v in c] for k, c in
                     phasediagram.region_centres(d).items()})
    luts = {n: base64.b64encode(colormaps.lut(n).tobytes()).decode()
            for n in colormaps.NAMES}
    rev = subprocess.run(["git", "-C", REPO, "rev-parse", "--short",
                          "HEAD"], capture_output=True, text=True)
    return dict(source=f"Quasicrystal@{rev.stdout.strip()}",
                lut=luts, diagrams=diagrams)


def main(path):
    data = json.dumps(build(), separators=(",", ":"))
    html = open(path, encoding="utf-8").read()
    pat = re.compile(re.escape(BEGIN) + ".*?" + re.escape(END), re.S)
    if not pat.search(html):
        sys.exit(f"no {BEGIN} ... {END} block in {path}")
    html = pat.sub(lambda _: f"{BEGIN}{data}{END}", html, count=1)
    open(path, "w", encoding="utf-8").write(html)
    print(f"ok: {len(Q_PRESETS)} diagrams, {len(data)} bytes")


if __name__ == "__main__":
    main(sys.argv[1])
