"""Export reference fields from the Python QC Sandbox core.

tools/quasicrystal/parity.html runs the post's own JavaScript against
these, so the browser port is checked against the code it ports
rather than against itself.  Output goes to tools/quasicrystal/ref/
(git-ignored).

    ~/GitHub/Quasicrystal/.venv/bin/python \
        tools/quasicrystal/export_reference.py
"""
import dataclasses
import json
import os
import sys

REPO = os.path.expanduser(os.environ.get("QC_REPO",
                                         "~/GitHub/Quasicrystal"))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from qcsandbox import colormaps, seeds  # noqa: E402
from qcsandbox.params import Params  # noqa: E402
from qcsandbox.sim import (Simulation, _fourier_resample,  # noqa: E402
                           _wrap_window)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref")
manifest = {}


def save(name, a, **meta):
    a = np.ascontiguousarray(a)
    if np.iscomplexobj(a):
        a = np.stack([a.real, a.imag], -1)
    dtype = "uint8" if a.dtype == np.uint8 else "float32"
    a.astype(dtype).tofile(os.path.join(OUT, name + ".bin"))
    manifest[name] = dict(file=name + ".bin", dtype=dtype,
                          shape=list(a.shape), **meta)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(1)

    # FFT on its own: the full complex transform the shaders compute.
    x = rng.standard_normal((128, 128)).astype(np.float32)
    save("fft_in", x)
    save("fft_out", np.fft.fft2(x.astype(np.float64)))

    # 100 steps at T = 0 from the same start, both dynamics.
    for mode in ("overdamped", "inertial"):
        p = Params(mode=mode, friction=1.0)
        sim = Simulation(n=128, params=p, seed=0)
        rho0 = sim.field()
        for _ in range(100):
            sim.step(1)
        save(f"rho0_{mode}", rho0)
        save(f"rho100_{mode}", sim.field(), params=dataclasses.asdict(
            sim.params), q=sim.params.q, dx=sim.dx, steps=100)

    # The same from a quasicrystal seed, where the field is order one
    # and the cubic term matters.
    sim = Simulation(n=128, seed=0)
    sim.set_field(seeds.make("quasicrystal", sim))
    save("rho0_qc", sim.field())
    for _ in range(100):
        sim.step(1)
    save("rho100_qc", sim.field(), dx=sim.dx)

    # Phase-diagram labels, to check the post's embedded copy.
    from qcsandbox import phasediagram
    from qcsandbox.params import Q_PRESETS
    for name in Q_PRESETS:
        d = phasediagram.cached(name, Q_PRESETS[name], 160, 160)
        save(f"labels_{name}", d.labels.astype(np.uint8),
             names=list(d.names))

    # One noisy step from a flat field: the noise variance per cell.
    p = Params(eps=0.0, T=0.01)
    sim = Simulation(n=512, params=p, seed=0)
    sim.set_field(np.zeros((512, 512), np.float32))
    sim.step(1)
    manifest["noise_step"] = dict(var=float(np.var(sim.field())),
                                  params=dataclasses.asdict(p), n=512,
                                  dx=sim.dx)

    # Seeds, on both a small and the default grid.
    for n in (128, 512):
        for q_name in ("k12", "k5", "custom"):
            p = Params(q_name=q_name, q_custom=1.9)
            sim = Simulation(n=n, params=p, seed=0)
            for kind in ("hexagonal", "quasicrystal", "half"):
                if q_name == "custom" and kind == "hexagonal":
                    continue
                save(f"seed_{kind}_{q_name}_{n}", seeds.make(kind, sim),
                     kind=kind, q_name=q_name, q=p.q, n=n, dx=sim.dx,
                     eps=p.eps)

    # Resampling for grid and zoom changes.
    a = rng.standard_normal((128, 128))
    save("resample_in", a)
    save("resample_up", _fourier_resample(a, 256), m=256)
    save("resample_down", _fourier_resample(a, 64), m=64)
    save("resample_odd", _fourier_resample(_wrap_window(a, 85), 128),
         window=85, m=128)

    # Rendering: the colour lookup and the Fourier view.
    sim = Simulation(n=128, seed=0)
    for _ in range(300):
        sim.step(1)
    f = sim.field()
    save("render_field", f, dx=sim.dx)
    table = colormaps.lut("red-white-blue")
    save("render_rgb", colormaps.to_rgb(f, table, -0.3, 0.3),
         lo=-0.3, hi=0.3, cmap="red-white-blue")
    save("render_fourier", colormaps.fourier_rgb(f, table, sim.dx),
         cmap="red-white-blue")
    step = max(1, f.shape[0] // 256)
    manifest["render_p995"] = float(
        np.percentile(np.abs(f[::step, ::step]), 99.5))

    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    print(f"ok: {len(manifest)} entries in {OUT}")


if __name__ == "__main__":
    main()
