"""pytest suite for the exact reference algorithm (PLAN.md §7: T1-T5).

BYU-dependent tests are skipped when tools/ssg-finder/local/ (cache, data) is absent.
Run:  /opt/anaconda3/bin/python3 -I -m pytest tools/ssg-finder/reference
"""
import functools, glob, json, os, random, sys
from fractions import Fraction as Fr

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import intmat as im  # noqa: E402

LOCAL = os.path.join(HERE, "..", "local")
CACHE = os.path.join(LOCAL, "cache")
GOLDEN = os.path.join(LOCAL, "fixtures", "golden")
byu = pytest.mark.skipif(not os.path.isdir(CACHE), reason="local/ BYU cache absent")


# ================================================================ T5: integer normal forms
def _rand_mat(rng, m, n):
    return [[rng.randint(-6, 6) * rng.choice([0, 1, 1]) for _ in range(n)] for _ in range(m)]


def test_snf_hnf_properties():
    rng = random.Random(0)
    for _ in range(2000):
        m, n = rng.randint(1, 5), rng.randint(1, 5)
        A = _rand_mat(rng, m, n)
        D, U, V, Vi = im.snf(A)
        assert im.matmul(im.matmul(U, A), V) == D
        assert im.matmul(V, Vi) == im.eye(n) and im.is_unimodular(U)
        assert all(D[i][j] == 0 for i in range(m) for j in range(n) if i != j)
        ds = [D[i][i] for i in range(min(m, n)) if D[i][i]]
        assert all(x > 0 for x in ds) and all(b % a == 0 for a, b in zip(ds, ds[1:]))
        H, W = im.hnf_with_transform(A)
        assert im.matmul(W, A) == H and im.is_unimodular(W)
        Hn = im.hnf(A)
        piv = [next(j for j in range(n) if r[j]) for r in Hn]
        assert piv == sorted(set(piv))
        for i, p in enumerate(piv):
            assert Hn[i][p] > 0 and all(0 <= Hn[k][p] < Hn[i][p] for k in range(i))
        # HNF is canonical: invariant under a unimodular change of generators
        assert im.hnf(im.matmul(W, A)) == Hn
        S = im.saturate(A)
        assert len(S) == im.rank(A)
        if S:
            D2, *_ = im.snf(S)
            assert all(D2[i][i] == 1 for i in range(len(S)))


def test_flat_solve_and_intersections():
    import exact
    rng = random.Random(3)
    for _ in range(300):
        d = rng.randint(1, 3)
        G = _rand_mat(rng, rng.randint(1, 3), d)
        h = [Fr(rng.randint(0, 7), 8) for _ in G]
        comps = exact.solve(G, h, d)
        for F in comps:
            p = F.point([Fr(rng.randint(0, 96), 97) for _ in range(F.dim)])
            assert all(exact.frac1(sum(a * b for a, b in zip(r, p)) - x) == 0 for r, x in zip(G, h))
            assert F.contains(p)
        assert len({F.key for F in comps}) == len(comps)
        if len(comps) >= 2:  # components are pairwise disjoint; intersection is commutative/idempotent
            A, B = comps[0], comps[1]
            assert exact.intersect(A, B) == []
            assert [F.key for F in exact.intersect(A, A)] == [A.key]


# ================================================================ BYU-dependent
@functools.lru_cache(maxsize=None)
def _ex(label, params):
    import exact
    return exact.run(label, dict(params))


def _key(p):
    return tuple(sorted(p.items()))


def _cases():
    import compare
    cases, short = compare.it_cases()
    return cases, short, compare.IN8_CASES


if os.path.isdir(CACHE):
    IT_CASES, SHORT, IN8 = _cases()
else:  # pragma: no cover
    IT_CASES, SHORT, IN8 = [], None, []
ALL = [(lab, {"g": g}) for lab, _, g, _ in IT_CASES] + list(IN8)
IDS = [f"{lab}:{','.join(f'{k}={v}' for k, v in p.items())}" for lab, p in ALL]


@byu
@pytest.mark.parametrize("lab,itname,g,expected", IT_CASES, ids=[f"{c[0]}:g={c[2]}" for c in IT_CASES])
def test_it_finder_rows(lab, itname, g, expected):
    ex = _ex(lab, _key({"g": g}))
    assert {SHORT(c["setting"]) for c in ex.classes} == expected


@byu
@pytest.mark.parametrize("lab,params", ALL, ids=IDS)
def test_grid_prototype(lab, params):
    import compare
    ex = _ex(lab, _key(params))
    r = compare.grid_compare(ex)
    assert r["problems"] == []
    # every exact class of dimension >= 1 is hit by the grid
    assert all(ex.classes[ci]["dim"] == 0 for ci in r["missed"])


@byu
@pytest.mark.parametrize("lab,params", ALL, ids=IDS)
def test_bruteforce(lab, params):
    import compare
    assert compare.bruteforce(_ex(lab, _key(params))) == []


@byu
@pytest.mark.parametrize("lab,params", ALL, ids=IDS)
def test_lambda_shift_invariance(lab, params):
    import compare
    assert compare.shift_invariance(_ex(lab, _key(params))) == []


@byu
@pytest.mark.parametrize("lab,params", ALL, ids=IDS)
def test_structure(lab, params):
    """Closure/orbit invariants and the 3D group closing at every generic point."""
    import exact
    ex = _ex(lab, _key(params))
    assert ex.ops[0][1] == tuple(Fr(0) for _ in range(3 + ex.d))
    assert len(ex.ops) == ex.rec["nops"]
    keys = [c["flat"].key for c in ex.classes]
    assert len(set(keys)) == len(keys)
    for c in ex.classes:
        # generic point lies in exactly its own flat's open part
        assert ex.classify_s(c["generic_s"]) is c["flat"]
        assert ex.classify_s(ex.s_of_t(c["example_t"])) is c["flat"]
        # surviving set closes under multiplication (it is the 3D group's coset set)
        o3 = {(tuple(map(tuple, R)), t) for R, t in c["ops3"]}
        for R1, t1 in c["ops3"]:
            for R2, t2 in c["ops3"][:8]:
                R = im.matmul(R1, R2)
                t = tuple(exact.frac1(sum(R1[i][k] * t2[k] for k in range(3)) + t1[i]) for i in range(3))
                assert (tuple(map(tuple, R)), t) in o3
        # SSG images of a class have the same type
        for o in ex.opinfo:
            img = ex.by_key[ex.image(c["flat"], o).key]
            assert img["orbit"] == c["orbit"] and img["sg_number"] == c["sg_number"]
    # intersections of closure flats stay in the closure
    fl = list(ex.flats.values())
    for A in fl[:12]:
        for B in fl[:12]:
            AB = {F.key for F in exact.intersect(A, B)}
            assert AB <= set(ex.flats) and AB == {F.key for F in exact.intersect(B, A)}


# ---------------------------------------------------------------- T2: In8IrSe9 README §4
AXES = {"55.3.122": "zxy", "55.3.123": "zxy", "14.3.16": "yzx"}
#           plane t3≡0  plane t3≡1/8  line          point t3≡0     point t3≡1/8
README = {
    93: ("A m 1 1", "A m 1 1", "A 2 1 1", "A 2/m 1 1", "A 2/m 1 1"),
    94: ("A m 1 1", "A m 1 1", "A 2_1 1 1", "A 2_1/m 1 1", "A 2_1/m 1 1"),
    95: ("A m 1 1", "A m 1 1", "A 2 1 1", "A 2/m 1 1", "A 2/m 1 1"),
    96: ("A e 1 1", "A e 1 1", "A 2 1 1", "A 2/e 1 1", "A 2/e 1 1"),
    97: ("A e 1 1", "A e 1 1", "A 2_1 1 1", "A 2_1/e 1 1", "A 2_1/e 1 1"),
    98: ("A e 1 1", "A e 1 1", "A 2 1 1", "A 2/e 1 1", "A 2/e 1 1"),
    99: ("A m 1 1", "A e 1 1", "A 2 1 1", "A 2/m 1 1", "A 2/e 1 1"),
    100: ("A m 1 1", "A e 1 1", "A 2_1 1 1", "A 2_1/m 1 1", "A 2_1/e 1 1"),
    101: ("A m 1 1", "A e 1 1", "A 2 1 1", "A 2/m 1 1", "A 2/e 1 1"),
}


def _our_setting(ops3, axes):
    import ssg_descent as proto
    P = [[1 if axes[i] == "xyz"[j] else 0 for j in range(3)] for i in range(3)]
    Pi = im.transpose(P)
    o = [(im.matmul(im.matmul(P, R), Pi), tuple(im.matvec(P, t))) for R, t in ops3]
    return proto.hm_in_setting(o)


@byu
@pytest.mark.parametrize("lab,params", IN8, ids=[c[0] for c in IN8])
def test_in8irse9_readme(lab, params):
    ex = _ex(lab, _key(params))
    axes = AXES[lab.rsplit(".", 1)[0]]
    dims = {c["dim"] for c in ex.classes}
    gen = [c for c in ex.classes if c["dim"] == 3]
    assert len(gen) == 1 and gen[0]["sg_number"] == 1 and _our_setting(gen[0]["ops3"], axes) == "A 1 1 1"
    if lab == "14.3.16.6":
        assert dims == {3, 0}
        assert all(c["sg_number"] == 2 for c in ex.classes if c["dim"] == 0)
        return
    p0, p8, line, pt0, pt8 = README[int(lab.rsplit(".", 1)[1])]
    assert dims == {3, 2, 1, 0}
    for c in ex.classes:
        s = _our_setting(c["ops3"], axes)
        t1, t2, t3 = c["generic_t"]
        if c["dim"] == 2:
            assert c["locus"] in ("t3 ≡ 0 (mod 1/4)", "t3 ≡ 1/8 (mod 1/4)")
            assert s == (p0 if c["locus"].endswith("0 (mod 1/4)") else p8)
        elif c["dim"] == 1:
            assert s == line and "t3" not in c["locus"]
            assert (8 * t1).denominator == 1 and (8 * t2).denominator == 1 and (4 * (t1 + t2)).denominator == 1
        elif c["dim"] == 0:
            assert s == (pt0 if (4 * t3).denominator == 1 else pt8)
            assert (8 * t3).denominator == 1


# ---------------------------------------------------------------- golden fixtures are current
@byu
@pytest.mark.skipif(not os.path.isdir(GOLDEN), reason="no golden fixtures written yet")
def test_golden_fixtures_reproduce():
    import exact, make_golden
    files = sorted(glob.glob(os.path.join(GOLDEN, "*__*.json")))
    rng = random.Random(5)
    sample = [f for f in files if not os.path.basename(f)[0].isdigit() or True]
    sample = rng.sample(sample, min(25, len(sample)))
    recs = make_golden.load_data() if os.path.isdir(make_golden.DATA) else {}
    for fn in sample:
        j = json.load(open(fn))
        lab = j["label"]
        if os.path.exists(os.path.join(CACHE, lab + ".html")):
            rec = lab
        elif lab in recs:
            rec = recs[lab]
        else:
            continue
        ex = exact.run(rec, j["params"])
        assert json.loads(json.dumps(ex.to_json(lab), ensure_ascii=False)) == j, fn
