"""Cross-checks of exact.py against the grid prototype (ssg_descent.py) and the independent
brute force (verify_bruteforce.build).  Used by test_exact.py and make_golden.py."""
import ast, itertools, math, os
from fractions import Fraction as Fr

import numpy as np

import exact
import intmat as im
import ssg_descent as proto

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- IT-finder cases
def it_cases():
    """CASES and short() from validate_it_finder.py, without running its module-level check."""
    src = open(os.path.join(HERE, "validate_it_finder.py")).read()
    tree = ast.parse(src)
    ns = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and node.targets[0].id == "CASES":
            ns["CASES"] = ast.literal_eval(node.value)
        if isinstance(node, ast.FunctionDef) and node.name == "short":
            exec(compile(ast.Module([node], []), "validate_it_finder.py", "exec"), ns)
    return ns["CASES"], ns["short"]


IN8_CASES = [(f"55.3.122.{i}", {"a1": "1/4", "b1": "1/4", "g2": "1/4"}) for i in range(93, 99)] + \
            [(f"55.3.123.{i}", {"a1": "1/4", "b1": "1/4", "g2": "1/4"}) for i in range(99, 102)] + \
            [("14.3.16.6", {"a1": "1/2", "a2": "1/4", "b1": "1/4", "b2": "1/4", "g2": "1/4"})]


# ---------------------------------------------------------------- prototype bridge
def proto_descent(ex):
    """Prototype Descent for the same group and Q.  Uses the BYU page when it is cached
    (independent parse), otherwise the canonical operator list of the record."""
    label = ex.rec.get("label")
    if label and os.path.exists(os.path.join(proto.CACHE, label + ".html")):
        ssg = proto.parse_ssg(label)
    else:
        nv = 3 + ex.d
        ssg = dict(d=ex.d, qsym=ex.rec["q"], centering=[tuple(Fr(0) for _ in range(nv))],
                   ops=[([[Fr(x) for x in r] for r in W], [Fr(x) for x in w]) for W, w in ex.ops])
    ds = proto.Descent(ssg, ex.Q)
    # memoise Descent.offsets per op (pure function of the op; survive() recomputes it
    # with a D^3 loop on every call, which dominates the run time for cubic groups)
    cache, raw = {}, ds.offsets
    ds.offsets = lambda op: cache[id(op)] if id(op) in cache else cache.setdefault(id(op), raw(op))
    return ds


def op_map(ex, ds):
    """prototype op index -> exact op index"""
    idx = {(tuple(tuple(r) for r in W), tuple(w)): k for k, (W, w) in enumerate(ex.ops)}
    out = []
    for o in ds.ops:
        W = [list(map(int, r)) + [0] * ex.d for r in o["R"]] + \
            [list(map(int, m)) + list(map(int, e)) for m, e in zip(o["M"], o["E"])]
        w = [exact.frac1(x) for x in list(o["tE"]) + list(o["tI"])]
        out.append(idx[(tuple(tuple(r) for r in W), tuple(w))])
    return out


def proto_point(ex, ds, omap, t):
    """(surviving [(exact idx, n)], sg number) computed with the prototype at exact t."""
    surv = []
    for j, o in enumerate(ds.ops):
        n = ds.survive(o, list(t))
        if n is not None:
            surv.append((omap[j], tuple(n)))
    surv = sorted(set(surv))
    o3 = ds.group3d(tuple(t))
    o3 = [([[int(x) for x in r] for r in R], tt) for R, tt in o3]
    return surv, exact.identify(o3, ex.lattice()).number, proto.hm_in_setting(o3)


def default_grid(ex):
    D = exact.den_lcm([x for r in ex.Q for x in r])
    return 8 * D if ex.d == 3 else 16 * D


def grid_compare(ex, G=None):
    """T4: every grid class maps to exactly one exact class with the same surviving set,
    n and type; every exact class is hit by the grid or checked at its exact point.
    Returns a dict with lists of problems (empty = agreement)."""
    ds = proto_descent(ex)
    omap = op_map(ex, ds)
    d = ex.d
    G = G or default_grid(ex)
    grid = np.array(list(itertools.product(range(G), repeat=d)), dtype=np.int64)
    # ---- prototype keys per point (same arithmetic as Descent.scan)
    offs = [ds.offsets(o) for o in ds.ops]
    K = G * exact.den_lcm([x for o in offs for v in o for x in v] + [Fr(1, G)])
    gk = grid * (K // G)
    keys = np.full((len(grid), len(ds.ops)), -1, dtype=np.int64)
    nlist = []
    for j, o in enumerate(ds.ops):
        E = np.array([[int(x) for x in r] for r in o["E"]], dtype=np.int64)
        v = (-(gk @ E.T - gk)) % K
        code = np.zeros(len(grid), dtype=np.int64)
        for i in range(d):
            code = code * K + v[:, i]
        al = list(offs[j].items())
        nlist.append([n for _, n in al])
        allowed = np.array([sum(int(x * K) * K ** (d - 1 - i) for i, x in enumerate(vv)) for vv, _ in al],
                           dtype=np.int64)
        order = np.argsort(allowed)
        pos = np.clip(np.searchsorted(allowed[order], code), 0, len(allowed) - 1)
        hit = allowed[order][pos] == code
        keys[hit, j] = order[pos[hit]]
    # ---- exact classification per point (vectorised flat membership)
    cls_list = ex.classes
    member = np.zeros((len(grid), len(cls_list)), dtype=bool)
    for ci, c in enumerate(cls_list):
        F = c["flat"]
        if not F.C:
            member[:, ci] = True
            continue
        T = im.matmul([list(r) for r in F.C], ex.Binv)
        L = exact.den_lcm([Fr(x, G) for r in T for x in r] + list(F.f))
        Ti = np.array([[int(Fr(x, G) * L) for x in r] for r in T], dtype=np.int64)
        fi = np.array([int(x * L) for x in F.f], dtype=np.int64)
        member[:, ci] = np.all(((grid @ Ti.T) - fi) % L == 0, axis=1)
    dims = np.array([c["dim"] for c in cls_list])
    big = np.where(member, dims[None, :], 99)
    cls_of = np.argmin(big, axis=1)
    assert np.all((big == big[np.arange(len(grid)), cls_of][:, None]).sum(axis=1) == 1)
    # ---- compare
    problems, hit, shared_keys = [], set(), []
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    ngrid_classes = len(uniq)
    for u in range(len(uniq)):
        pts = np.nonzero(inv == u)[0]
        ecs = set(cls_of[pts].tolist())
        hit |= ecs
        psurv = sorted((omap[j], tuple(nlist[j][uniq[u][j]])) for j in range(len(ds.ops)) if uniq[u][j] >= 0)
        psurv = sorted(set(psurv))
        # A grid key (surviving ops + n at that t) can be shared by several exact classes:
        # n depends on the Λ-translate of t, so different s-flats can produce identical keys
        # (e.g. 55.1.9.4, g=2/5: t=0 and t=1/2).  They must then agree in ops and type.
        orbits = {cls_list[e]["orbit"] for e in ecs}
        opsets = {tuple(k for k, _ in cls_list[e]["surviving"]) for e in ecs}
        types = {cls_list[e]["sg_number"] for e in ecs}
        if len(opsets) != 1 or len(types) != 1:
            problems.append(f"grid class {psurv} spans inconsistent exact classes {sorted(ecs)}")
            continue
        if len(ecs) > 1:
            shared_keys.append((sorted(ecs), len(orbits) > 1))
        t0 = tuple(Fr(int(x), G) for x in grid[pts[0]])
        c = cls_list[cls_of[pts[0]]]
        if psurv != ex.surviving_at_t(t0):
            problems.append(f"n mismatch at {t0}: grid {psurv} exact {ex.surviving_at_t(t0)}")
        if [k for k, _ in psurv] != [k for k, _ in c["surviving"]]:
            problems.append(f"surviving mismatch: grid {psurv} exact {c['surviving']}")
        _, num, _ = proto_point(ex, ds, omap, t0)
        if num != c["sg_number"]:
            problems.append(f"type mismatch at {t0}: grid {num} exact {c['sg_number']}")
    missed = [ci for ci in range(len(cls_list)) if ci not in hit]
    for ci in missed:
        c = cls_list[ci]
        surv, num, _ = proto_point(ex, ds, omap, c["example_t"])
        if surv != ex.surviving_at_t(c["example_t"]) or \
                [k for k, _ in surv] != [k for k, _ in c["surviving"]] or num != c["sg_number"]:
            problems.append(f"missed class {ci} disagrees at its example point")
    return dict(problems=problems, grid_classes=ngrid_classes, missed=missed, G=G,
                shared_keys=shared_keys)


# ---------------------------------------------------------------- brute force (T3)
def bruteforce(ex, n_extra=2, seed=1):
    """For every class: cut an explicit symmetrised modulated structure at the generic
    point and n_extra other points of the same flat, identify with spglib."""
    import spglib
    import verify_bruteforce as vb
    vb.rng = np.random.default_rng(seed)
    rng = np.random.default_rng(seed + 1)
    ds = proto_descent(ex)
    # generic invariant metric of the basic point group
    G0 = rng.normal(size=(3, 3)); G0 = G0 @ G0.T + 3 * np.eye(3)
    Rs = [np.array(o["R"], dtype=float) for o in ex.opinfo]
    Gb = sum(R.T @ G0 @ R for R in Rs) / len(Rs)
    Nd = np.diag(ex.N).astype(float)
    lat = np.linalg.cholesky(Nd @ Gb @ Nd)
    problems = []
    for ci, c in enumerate(ex.classes):
        F = c["flat"]
        pts = [c["generic_t"]]
        a = 0
        while len(pts) < 1 + n_extra and a < 50:
            a += 1
            s = F.point([Fr(int(k), 997) for k in rng.integers(1, 997, size=F.dim)])
            if ex.in_open_part(F, s):
                pts.append(ex.t_of_s(s))
            if F.dim == 0:
                break
        for t in pts:
            # spglib occasionally fails (returns None) on an unlucky random structure with
            # near-coincident atoms; draw a fresh structure then (a wrong type is never retried)
            for _attempt in range(3):
                _, pos, typ = vb.build(ds, t, (1, 1, 1))
                ds_ = spglib.get_symmetry_dataset((lat, pos, typ), symprec=1e-4)
                if ds_ is not None:
                    break
            if ds_ is None or ds_.number != c["sg_number"]:
                problems.append(f"class {ci} at t={[str(x) for x in t]}: brute force "
                                f"{None if ds_ is None else ds_.number}, exact {c['sg_number']}")
    return problems


# ---------------------------------------------------------------- Λ-shift invariance (T5)
def shift_invariance(ex, seed=2):
    rng = np.random.default_rng(seed)
    ds = proto_descent(ex)
    omap = op_map(ex, ds)
    problems = []
    for ci, c in enumerate(ex.classes):
        n = [int(x) for x in rng.integers(-5, 6, size=3)]
        t2 = [exact.frac1(a + b) for a, b in zip(c["generic_t"], im.matvec(ex.Q, n))]
        # n rule: t -> t + Q m  shifts the matched n by (R - 1) m (mod L_super)
        want = []
        for k, n0 in c["surviving"]:
            R = ex.opinfo[k]["R"]
            nn = [n0[i] + sum(R[i][j] * n[j] for j in range(3)) - n[i] for i in range(3)]
            want.append((k, ex.nlookup[tuple(exact.frac1(x) for x in im.matvec(ex.Q, nn))]))
        if ex.surviving_at_t(t2) != want:
            problems.append(f"class {ci}: n at t+Qn does not follow n + (R-1)m")
        if ex.classify_s(ex.s_of_t(t2)) is not c["flat"]:
            problems.append(f"class {ci}: shifted point left the flat")
        surv, num, setting = proto_point(ex, ds, omap, t2)
        if [k for k, _ in surv] != [k for k, _ in c["surviving"]] or num != c["sg_number"] \
                or setting != c["setting"]:
            problems.append(f"class {ci}: shift by Qn={n} changed the group")
    return problems
