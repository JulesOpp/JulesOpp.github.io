#!/usr/bin/env python3
"""Dump identification fixtures for test/identify.test.ts (milestone M4).

For every cached BYU group (local/cache/*.html), every IT-finder case of
validate_it_finder.py and every brute-force (In8IrSe9) case of verify_bruteforce.py,
run the grid prototype (ssg_descent.run) and write, per class, the 3D operations
(`ops3`, conventional-supercell basis, centrings included as R = I operations)
together with spglib's number/symbol and the prototype's setting symbol.

Output: local/fixtures/identify/<label>__<params>.json  (BYU-derived => gitignored, D1).

Usage:  /opt/anaconda3/bin/python3 -I reference/make_identify_fixtures.py
"""
import ast, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ssg_descent import run  # noqa: E402

ROOT = os.path.join(HERE, "..")
CACHE = os.path.join(ROOT, "local", "cache")
OUT = os.path.join(ROOT, "local", "fixtures", "identify")

# Default parameters for cached groups that no case list mentions (the In8IrSe9 values).
DEFAULT_D3 = {"a1": "1/4", "b1": "1/4", "g2": "1/4"}
DEFAULT_D1 = {"g": "1/3"}


def it_cases():
    """CASES from validate_it_finder.py, read with ast (importing it would run the script)."""
    tree = ast.parse(open(os.path.join(HERE, "validate_it_finder.py")).read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "CASES" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("CASES not found")


def bruteforce_cases():
    """(label, params, cell_len) of every check(...) call in verify_bruteforce.py's main block."""
    tree = ast.parse(open(os.path.join(HERE, "verify_bruteforce.py")).read())
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "check":
            label = ast.literal_eval(node.args[0])
            params = ast.literal_eval(node.args[1])
            cell = (12.00, 12.21, 7.94)
            for kw in node.keywords:
                if kw.arg == "cell_len":
                    cell = ast.literal_eval(kw.value)
            out.append((label, params, cell))
    return out


def fr(x):
    return str(x)


def dump(label, params, cell_len, it=None):
    ssg, Q, ds, G, out = run(label, params)
    classes = []
    for r in out:
        classes.append(dict(
            t0=[fr(x) for x in r["t0"]], dim=int(r["dim"]), npts=int(r["npts"]),
            surv=[int(j) for j in r["surv"]],
            sg_number=int(r["sg_number"]), sg_symbol=r["sg_symbol"], setting=r["setting"],
            ops3=[dict(R=[[fr(x) for x in row] for row in R], t=[fr(x) for x in t]) for R, t in r["ops3"]],
        ))
    rec = dict(label=label, symbol=ssg["symbol"], d=ssg["d"], params=params,
               Q=[[fr(x) for x in row] for row in Q], N=list(ds.N), ncells=len(ds.cells),
               cell_len=list(cell_len), grid=G, classes=classes)
    if it is not None:
        rec["it"] = it
    tag = "_".join(f"{k}={v.replace('/', 'o')}" for k, v in sorted(params.items()))
    fn = os.path.join(OUT, f"{label}__{tag}.json")
    with open(fn, "w") as f:
        json.dump(rec, f, separators=(",", ":"))
    print(f"{fn}: {len(classes)} classes, N={ds.N}, index={ds.N[0]*ds.N[1]*ds.N[2]//len(ds.cells)}")
    return fn


def main():
    os.makedirs(OUT, exist_ok=True)
    done = set()
    for lab, itname, g, exp in it_cases():
        key = (lab, (("g", g),))
        if key in done:
            continue
        done.add(key)
        dump(lab, {"g": g}, (12.00, 12.21, 7.94), it=dict(name=itname, expected=sorted(exp)))
    for lab, params, cell in bruteforce_cases():
        key = (lab, tuple(sorted(params.items())))
        if key in done:
            continue
        done.add(key)
        dump(lab, params, cell)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".html"):
            continue
        lab = fn[:-5]
        if any(k[0] == lab for k in done):
            continue
        params = DEFAULT_D3 if lab.split(".")[1] == "3" else DEFAULT_D1
        done.add((lab, tuple(sorted(params.items()))))
        dump(lab, params, (12.00, 12.21, 7.94))


if __name__ == "__main__":
    main()
