# Spike (2026-11-10): BYU ssg_tables.zip (supercentred ops) + lattice transform -> basic-setting ops.
# Compared with cached BYU HTML pages in analysis/ssg_descent/cache: 17/17 exact set matches.
# Usage: unzip ssg_tables.zip into $BYU_TABLES, then: python3 -I convert_check.py
import os
import re, sys, glob, os
from fractions import Fraction as Fr
sys.path.insert(0, "/Users/julesoppenheim/dials/In8IrSe9_Xtal4_FullRun/analysis/ssg_descent")
import ssg_descent as S
T = os.environ.get("BYU_TABLES", "byu/t") + "/"
def blocks(fn):
    return [b.strip().split("\n") for b in open(fn).read().split("\n \n") if b.strip()]
groups, lat = {}, {}
for d in (1,2,3):
    for b in blocks(T+f"ssg_group{d}_tables.txt"):
        lab, sym = b[0].split(None,1)
        rest = " ".join(b[2:])
        rest = re.split(r"[HKLMNP0-9-]+:", rest)[0]
        m = re.match(r"((?:\s*\([^)]*\),?)+)(.*)", rest)
        cen = re.findall(r"\(([^)]*)\)", m.group(1)); body = m.group(2)
        groups[lab] = dict(sym=sym, d=d, cen=cen, ops=[o.strip() for o in body.split(";") if o.strip()])
    for b in blocks(T+f"ssg_lattice{d}_tables.txt"):
        lab = b[0].split()[0]
        pl = " ".join(l for l in b[1:] if l.lstrip().startswith("A") or (not l.lstrip().startswith(("Q","where","(")) and "=" in l and "a" in l))
        lat[lab] = dict(P=[x.split("=")[1].strip() for x in re.split(r",\s*", pl.strip().rstrip(",")) if x.strip()], d=d)
def lin(e, names, n):
    c,k = S.parse_lin(e.replace(" ",""), 0) if False else (None,None)
def parse_vec(expr, names):
    e = expr.replace(" ","");  e = e if e[0] in "+-" else "+"+e
    coef=[Fr(0)]*len(names); const=Fr(0)
    for sg,body in re.findall(r"([+-])([^+-]+)", e):
        s=1 if sg=="+" else -1
        m=re.fullmatch(r"(\d+(?:/\d+)?)?\*?([A-Za-z]+\d*)", body)
        if m: coef[names.index(m.group(2))]+=s*(Fr(m.group(1)) if m.group(1) else 1)
        else: const+=s*Fr(body)
    return coef,const
import numpy as np
def inv(M):
    import sympy; return [[Fr(str(x)) for x in r] for r in sympy.Matrix(M).inv().tolist()]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def basic_ops(lab):
    g=groups[lab]; d=g["d"]; n=3+d
    cls=".".join(lab.split(".")[1:3]); L=lat[cls]
    an=[f"a{i+1}" for i in range(n)]
    P=[parse_vec(p,an)[0] for p in L["P"]]       # A_i = sum_j P_ij a_j
    PT=[list(r) for r in zip(*P)]; PTi=inv(PT)
    XN="XYZTUV"[:n]
    out=set()
    for c in g["cen"]:
        cv=[Fr(x) for x in c.split(",")]
        for o in g["ops"]:
            rows=[parse_vec(p,list(XN)) for p in o.split(",")]
            W=[r[0] for r in rows]; w=[r[1]+ci for r,ci in zip(rows,cv)]
            Wb=mm(mm(PT,W),PTi); wb=[x%1 for x in [r[0] for r in mm(PT,[[x] for x in w])]]
            out.add((tuple(map(tuple,Wb)),tuple(wb)))
    return out
def html_ops(lab):
    s=S.parse_ssg(lab); out=set()
    for c in s["centering"]:
        for W,w in s["ops"]:
            out.add((tuple(map(tuple,W)),tuple((a+b)%1 for a,b in zip(w,c))))
    return out
def closure_mod_lattice(ops):  # compare as sets of cosets incl. centring
    return ops
ok=0
for fn in sorted(glob.glob("/Users/julesoppenheim/dials/In8IrSe9_Xtal4_FullRun/analysis/ssg_descent/cache/*.html")):
    lab=os.path.basename(fn)[:-5]
    a,b=basic_ops(lab),html_ops(lab)
    r = a==b; ok+=r
    print(lab, "MATCH" if r else f"DIFF {len(a)} {len(b)}")
print(ok, len(groups))
