#!/usr/bin/env python3
"""Validate ssg_descent.py against the IUCr IT Superspace Group Finder (Orlov, Palatinus &
Chapuis), https://it.iucr.org/resources/finder/  -- table '3+1D to 3D group transformation
for commensurate modulation' (retrieved 2026-11-10 via the finder's gettransform.php).

Expected = set of 3D groups (symbols in the parental basis, as printed by the finder) over all
t0 columns (General, n/N, 1/4N+n/N, 1/2N+n/N, 3/4N+n/N).  We compare space-group TYPES
(number) and also the setting symbol produced by ssg_descent.hm_in_setting.
"""
import sys
from ssg_descent import run

# (BYU label, IT finder name, gamma, expected setting symbols from IT (spaces removed))
CASES = [
    ("62.1.9.2", "62.2 Pnma(00g)0s0", "1/3", {"P21", "P21/a", "P212121"}),
    ("62.1.9.2", "62.2 Pnma(00g)0s0", "2/5", {"Pn", "P21/n", "Pn21a"}),
    ("62.1.9.2", "62.2 Pnma(00g)0s0", "1/2", {"Pc", "P21/c", "P21ca"}),
    ("62.1.9.6", "62.6 Pmcn(00g)s00", "1/3", {"P21", "P21/n", "P212121"}),
    ("62.1.9.6", "62.6 Pmcn(00g)s00", "2/5", {"Pc", "P21/c", "P21cn"}),
    ("62.1.9.6", "62.6 Pmcn(00g)s00", "1/2", {"Pc", "P21/c", "Pc21n"}),
    ("55.1.9.2", "55.2 Pbam(00g)s00", "1/3", {"Pa", "P21/a", "P21am"}),
    ("55.1.9.2", "55.2 Pbam(00g)s00", "1/2", {"Pna21", "Pnam"}),
    ("55.1.9.3", "55.3 Pbam(00g)ss0", "1/3", {"P2", "P2/m", "P21212"}),
    ("55.1.9.3", "55.3 Pbam(00g)ss0", "1/2", {"Pnn2", "Pnnm"}),
    ("55.1.9.4", "55.4 Pcma(00g)", "2/5", {"Pcm21", "Pcma"}),
    ("55.1.9.4", "55.4 Pcma(00g)", "1/3", {"Pm", "P2/m", "P21ma"}),
    ("55.1.9.4", "55.4 Pcma(00g)", "1/2", {"Pm", "P2/m", "P21ma"}),
    ("55.1.9.5", "55.5 Pcma(00g)0s0", "1/3", {"P21", "P21/a", "P21221"}),
    ("55.1.9.5", "55.5 Pcma(00g)0s0", "1/2", {"Pc", "P2/c", "P21ca"}),
    ("55.1.9.5", "55.5 Pcma(00g)0s0", "2/5", {"Pc", "P21/c", "Pc2a"}),
    ("55.1.10.6", "55.6 Pcma(01/2g)", "1/4", {"Am", "A2/m", "A21ma"}),
    ("55.1.10.6", "55.6 Pcma(01/2g)", "2/5", {"Pcm21", "Pcma", "Pcmn"}),
    ("55.1.10.6", "55.6 Pcma(01/2g)", "1/3", {"Pnm21", "Pnmn", "Pnma"}),
]

def short(setting):
    """'P 2_1/n 2/m 2_1/n' -> 'P21/n2/m21/n' then reduce to IT short form."""
    cen, *ax = setting.split()
    ax = [a.replace("_", "") for a in ax]
    if sum(a != "1" for a in ax) == 1:            # monoclinic: keep only the unique axis
        ax = [a for a in ax if a != "1"]
    elif all("/" in a for a in ax):               # holohedral orthorhombic: short = planes
        ax = [a.split("/")[1] for a in ax]
    return cen + "".join(ax)

ok = True
for lab, itname, g, exp in CASES:
    ssg, Q, ds, G, out = run(lab, {"g": g})
    got = {short(r["setting"]) for r in out}
    flag = "OK " if got == exp else "MISMATCH"
    ok &= got == exp
    print(f"{flag} {itname:22s} g={g:4s} IT={sorted(exp)}  ours={sorted(got)}")
print("ALL OK" if ok else "SOME MISMATCHES")
sys.exit(0 if ok else 1)
