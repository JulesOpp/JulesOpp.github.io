#!/usr/bin/env python3
"""Independent brute-force check of ssg_descent.py.

Builds an explicit modulated structure obeying the superspace group (two general-position
atoms with random displacive modulation functions, symmetrised by applying every SSG
operation), cuts the 3D section at a given t0 in the N1 x N2 x N3 supercell, and lets spglib
find the space group of the resulting atom list.  This uses none of the coset / survival logic
of ssg_descent.py.

String transformation used (superspace basis, x_s4 = t + q.x):
  atom with average position xb and displacement u(y), y = internal coordinate of the
  average position, maps under g = (R, M, eps | tE, tI) to
      xb' = R xb + tE,   u'(y') = R u( eps^-1 (y' - M xb - tI) ).
  The 3D atoms of the section t0 are at  xb + L + u(t0 + Q (xb + L))  for lattice vectors L.
"""
import sys, itertools
import numpy as np
import spglib
from ssg_descent import parse_ssg, eval_q, Descent, run, fmt_t

rng = np.random.default_rng(1)

def build(ds, t0, cell_len, n_orbits=2, amp=0.02):
    d = ds.d
    Q = np.array([[float(x) for x in r] for r in ds.Q])
    ops = [(np.array([[float(x) for x in r] for r in o["R"]]),
            np.array([[float(x) for x in r] for r in o["M"]]) if d else None,
            np.array([[float(x) for x in r] for r in o["E"]]),
            np.array([float(x) for x in o["tE"]]),
            np.array([float(x) for x in o["tI"]])) for o in ds.ops]
    # harmonics: all integer vectors m in {-1,0,1}^d (plus a couple of 2nd order)
    harm = [np.array(m) for m in itertools.product([-1, 0, 1], repeat=d) if any(m)]
    harm += [2 * np.eye(d, dtype=int)[i] for i in range(d)]
    atoms = []  # (species, xb, function)
    for s in range(n_orbits):
        xb0 = rng.random(3)
        A = rng.normal(size=(len(harm), 3)) * amp
        B = rng.normal(size=(len(harm), 3)) * amp
        def u0(y, A=A, B=B):
            ph = 2 * np.pi * np.array([m @ y for m in harm])
            return np.cos(ph) @ A + np.sin(ph) @ B
        for R, M, E, tE, tI in ops:
            xb = R @ xb0 + tE
            Ei = np.linalg.inv(E)
            def u(y, R=R, M=M, Ei=Ei, tI=tI, xb0=xb0, u0=u0):
                return R @ u0(Ei @ (y - M @ xb0 - tI))
            atoms.append((s + 1, xb, u))
    N = ds.N
    pos, typ = [], []
    t0 = np.array([float(x) for x in t0])
    for sp, xb, u in atoms:
        for L in itertools.product(*[range(k) for k in N]):
            xa = xb + np.array(L)
            y = t0 + Q @ xa
            r = xa + u(y)
            pos.append(r / np.array(N)); typ.append(sp)
    pos = np.array(pos) % 1.0
    # remove exact duplicates (shouldn't happen for general positions)
    lat = np.diag([cell_len[i] * N[i] for i in range(3)])
    return lat, pos, np.array(typ)

def check(label, params, cell_len=(12.00, 12.21, 7.94), G=None, max_classes=99):
    ssg, Q, ds, Gu, out = run(label, params, G=G)
    print(f"## {label} {ssg['symbol']} params={params}")
    ok = True
    for r in out[:max_classes]:
        for t0 in r["t0s"][:2]:
            lat, pos, typ = build(ds, t0, cell_len)
            dset = spglib.get_symmetry_dataset((lat, pos, typ), symprec=1e-4)
            num = dset.number if hasattr(dset, "number") else dset["number"]
            sym = dset.international if hasattr(dset, "international") else dset["international"]
            flag = "OK" if num == r["sg_number"] else "MISMATCH"
            ok &= num == r["sg_number"]
            print(f"  t0={fmt_t(t0):>16}  predicted {r['sg_symbol']:>10} ({r['sg_number']:>3})  "
                  f"brute force {sym:>10} ({num:>3})  {flag}")
    return ok

if __name__ == "__main__":
    ok = True
    ok &= check("62.1.9.2", {"g": "1/3"})
    ok &= check("55.1.10.6", {"g": "1/3"})
    ok &= check("55.1.10.6", {"g": "1/4"})
    ok &= check("55.3.122.93", {"a1": "1/4", "b1": "1/4", "g2": "1/4"})
    ok &= check("55.3.122.96", {"a1": "1/4", "b1": "1/4", "g2": "1/4"})
    ok &= check("55.3.123.100", {"a1": "1/4", "b1": "1/4", "g2": "1/4"})
    ok &= check("14.3.16.6", {"a1": "1/2", "a2": "1/4", "b1": "1/4", "b2": "1/4", "g2": "1/4"},
                cell_len=(12.21, 7.94, 12.00))
    print("ALL OK" if ok else "SOME MISMATCHES")
