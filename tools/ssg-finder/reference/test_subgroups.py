"""pytest suite for the (3+1)D equivalence oracle and the k-index-1 subgroup enumeration.

Needs the BYU tables in tools/ssg-finder/local/data (licence rule D1); skipped without them.
The default run covers a representative subset (every Bravais class, plus the groups named
below).  Set SSG_FULL=1 for the sweep over all 775 (3+1)D records:

  SSG_FULL=1 /opt/anaconda3/bin/python3 -I -m pytest tools/ssg-finder/reference/test_subgroups.py
"""
import collections, functools, json, os, sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ssg_equiv import load_records  # noqa: E402

DATA = os.path.join(HERE, "..", "local", "data")
FULL = os.environ.get("SSG_FULL") == "1"
byu = pytest.mark.skipif(not os.path.isdir(DATA), reason="local/ BYU data absent")

# groups that must always be covered, whatever the subset
NAMED = ["14.1.6.3", "15.1.7.3", "165.1.24.1", "194.1.24.1"]
NAMED_SG = {14, 15, 62, 165, 194}


def rec(label):
    return load_records(1)[label]


def sg_of(label):
    return int(label.split(".")[0])


@functools.lru_cache(maxsize=None)
def sample():
    """Labels under test: all of them with SSG_FULL=1, else the first record of every
    Bravais class, all records of NAMED_SG and NAMED."""
    recs = load_records(1)
    if FULL:
        return list(recs)
    seen, out = set(), []
    for label, r in recs.items():
        if r["cls"] not in seen or sg_of(label) in NAMED_SG or label in NAMED:
            seen.add(r["cls"])
            out.append(label)
    return out


@byu
def test_sample_covers_named_groups():
    s = sample()
    for label in NAMED:
        assert label in s
    assert any(sg_of(l) == 62 for l in s)
    assert len({rec(l)["cls"] for l in s}) == len({r["cls"] for r in load_records(1).values()})


@byu
def test_self_equivalent():
    from ssg_equiv import equivalent, q_generic
    bad = [l for l in sample() if equivalent(rec(l), q_generic(rec(l)), rec(l)) is None]
    assert bad == []


@byu
def test_distinct_not_equivalent():
    from ssg_equiv import equivalent, q_generic
    by_sg = collections.defaultdict(list)
    for l in sample():
        by_sg[sg_of(l)].append(l)
    hits = []
    for labels in by_sg.values():
        for a in labels:
            ra = rec(a)
            qa = q_generic(ra)
            for b in labels:
                if a != b and equivalent(ra, qa, rec(b)) is not None:
                    hits.append((a, b))
    assert hits == []


@byu
def test_equivalence_map_is_valid():
    """The returned map really conjugates A onto B (checked independently of the search)."""
    from ssg_equiv import equivalent, q_generic, check_map
    for l in ["165.1.24.1", "14.1.6.3", "15.1.7.3"]:
        r = rec(l)
        res = equivalent(r, q_generic(r), r)
        assert res is not None and check_map(r, q_generic(r), r, res)
    # subgroup maps, incl. P onto a centred type (P4bm ⊃ Cmm2, Bm: det P = 1/2)
    import exact
    from subgroups import t_subgroups, subgroup_record
    for l in ["100.1.19.1", "165.1.24.1"]:
        r = rec(l)
        for c in t_subgroups(r):
            sub = subgroup_record(r, exact.build_ops(r), c["pg_ops"])
            assert check_map(sub, q_generic(r), rec(c["byu"]), c["map"]), (l, c["pg_ops"])


def _rec(sym, gen, q="0,0,g", cen=()):
    import exact
    r = dict(sym=sym, cls="", d=1, q=[q], params=exact.params_of([q]), cen=list(cen), gen=list(gen))
    r["nops"] = len(exact.build_ops(r))
    return r


def _transform(r, P, p, eta, m, p4):
    """The record r rewritten in the setting x' = P x + p, x4' = m·x + η x4 + p4."""
    from fractions import Fraction as Fr
    import exact, intmat as im
    from ssg_equiv import _S, _Sinv, mm, mv, vm, mat
    P, m = mat(P), [Fr(x) for x in m]
    s = [Fr(x) for x in p] + [Fr(p4)]
    S, Si = _S(P, m, eta), _Sinv(im.inv_rat(P), m, eta)
    gen = []
    for g in r["gen"]:
        W, w = exact.parse_op(g, 4)
        Wp = mm(mm(S, mat(W)), Si)
        wp = [a + b - c for a, b, c in zip(mv(S, w), s, mv(Wp, s))]
        gen.append(exact.format_op([[int(x) for x in row] for row in Wp], [exact.frac1(x) for x in wp]))
    cen = [",".join(exact.fstr(exact.frac1(x)) for x in mv(S, [Fr(c) for c in v.split(",")])) for v in r["cen"]]
    zero = {k: 0 for k in r["params"]}
    b0 = exact.eval_q(r, zero)[0]
    Pi = im.inv_rat(P)
    c0 = vm([eta * x + y for x, y in zip(b0, m)], Pi)
    comps = [exact.fstr(x) for x in c0]
    for k in r["params"]:
        bd = [x - y for x, y in zip(exact.eval_q(r, dict(zero, **{k: 1}))[0], b0)]
        for j, x in enumerate(vm([eta * v for v in bd], Pi)):
            if x:
                comps[j] = ("" if comps[j] == "0" else comps[j] + "+") + f"{exact.fstr(x)}*{k}"
    return dict(r, gen=gen, cen=cen, q=[",".join(c.replace("+-", "-") for c in comps)])


def test_screw_phase_removed_by_q_shift():
    """P2_1(0,0,g)s ~ P2_1(0,0,g)0 via q -> q + c* (needs the integer freedom in m)."""
    from ssg_equiv import equivalent, q_generic, check_map
    a = _rec("P2_1(0,0,g)s", ["-x,-y,z+1/2,t+1/2"])
    b = _rec("P2_1(0,0,g)0", ["-x,-y,z+1/2,t"])
    c = _rec("P2(0,0,g)s", ["-x,-y,z,t+1/2"])
    res = equivalent(a, q_generic(a), b)
    assert res is not None and check_map(a, q_generic(a), b, res)
    assert equivalent(c, q_generic(c), _rec("P2(0,0,g)0", ["-x,-y,z,t"])) is None


@byu
def test_transformed_copy_equivalent():
    """A record rewritten in another superspace setting is found equivalent, and the map checks."""
    from ssg_equiv import equivalent, q_generic, check_map
    settings = [([[0, 1, 0], [0, 0, 1], [1, 0, 0]], ["1/4", "0", "1/8"], -1, [1, 0, 0], "1/3"),
                ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], ["0", "1/2", "0"], 1, [0, 0, 1], "1/2"),
                ([[1, 1, 0], [0, 1, 0], [0, 0, 1]], ["1/3", "0", "0"], -1, [0, 0, 0], "0")]
    bad = []
    for l in NAMED + ["62.1.9.2", "2.1.1.1"]:
        r = rec(l)
        for st in settings:
            t = _transform(r, *st)
            res = equivalent(r, q_generic(r), t)
            if res is None or not check_map(r, q_generic(r), t, res):
                bad.append((l, st))
    assert bad == []


@byu
def test_p3c1_classes():
    from subgroups import t_subgroups
    cl = t_subgroups(rec("165.1.24.1"))
    assert len(cl) == 10
    assert sorted(c["sg"] for c in cl) == [1, 2, 5, 9, 15, 143, 147, 150, 158, 165]
    assert sorted(c["index"] for c in cl) == [1, 2, 2, 2, 3, 4, 6, 6, 6, 12]
    assert {c["sg"]: c["class_size"] for c in cl} == {165: 1, 158: 1, 150: 1, 147: 1, 15: 3, 143: 1,
                                                     2: 1, 5: 3, 9: 3, 1: 1}
    c150 = [c for c in cl if c["sg"] == 150]
    assert len(c150) == 1 and c150[0]["byu"] == "150.1.24.1"
    c165 = [c for c in cl if c["sg"] == 165]
    assert c165[0]["byu"] == "165.1.24.1" and c165[0]["index"] == 1


@byu
def test_every_subgroup_identified():
    from subgroups import t_subgroups
    bad = []
    for l in sample():
        for c in t_subgroups(rec(l)):
            if len(c["matches"]) != 1:
                bad.append((l, c["pg_ops"], c["sg"], c["matches"]))
    assert bad == []


# Public hygiene checks use synthetic records and fixtures, never local BYU data.
def test_cap_jobs():
    from make_subgroup_fixtures import cap_jobs
    assert [cap_jobs(n) for n in (16, 0, 3)] == [4, 1, 3]


def test_worker_env():
    import subprocess
    env = dict(os.environ, SSG_JOBS="16", **{k: "8" for k in (
        "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
    code = ("import sys; sys.path.insert(0, " + repr(HERE) + "); "
            "import make_subgroup_fixtures as m; assert m.MAX_JOBS == 4; "
            "assert all(m.os.environ[k] == '1' for k in "
            "('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')); "
            "m.os.environ['OMP_NUM_THREADS'] = '8'; m.worker_init(); "
            "assert m.os.environ['OMP_NUM_THREADS'] == '1'")
    subprocess.run([sys.executable, "-I", "-c", code], env=env, check=True)


def test_d_guard_raises():
    import ssg_equiv as se
    from subgroups import t_subgroups
    for fn in (t_subgroups, se.q_generic, se.Prep):
        with pytest.raises(ValueError, match="d = 1"):
            fn({"d": 4})


def test_load_records(tmp_path):
    from ssg_equiv import load_records, load_d1
    (tmp_path / "sg001.json").write_text(json.dumps({"groups": {
        "1.1.1.1": {"d": 1}, "1.2.1.1": {"d": 2}}}))
    records = load_records(1, tmp_path)
    assert records == {"1.1.1.1": {"d": 1, "label": "1.1.1.1"}}
    assert load_records(1, tmp_path) is records
    assert load_d1(tmp_path) is records
    assert list(load_records(2, tmp_path)) == ["1.2.1.1"]


def test_summary_seconds_after_resume(tmp_path, monkeypatch):
    import make_subgroup_fixtures as m

    class SerialPool:
        def __init__(self, jobs, initializer):
            assert jobs == 4
            initializer()
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def imap_unordered(self, fn, labels):
            return map(fn, labels)

    monkeypatch.setattr(m.multiprocessing, "Pool", SerialPool)
    calls = []
    def fixture(label):
        calls.append(label)
        return dict(label=label, sym="synthetic", order=1, classes=[])
    monkeypatch.setattr(m, "fixture", fixture)
    ticks = iter([10, 11, 13, 15, 20, 23])
    monkeypatch.setattr(m.time, "time", lambda: next(ticks))
    m.main(["1.1.1.1", "--jobs", "16"], out=tmp_path)
    prior = json.loads((tmp_path / "_summary.json").read_text())
    m.main(["1.1.1.1", "--jobs", "16", "--resume"], out=tmp_path)
    summary = json.loads((tmp_path / "_summary.json").read_text())
    assert prior["seconds"] == prior["wall_seconds_this_run"] == 5
    assert summary["wall_seconds_this_run"] == 3
    assert summary["seconds"] == prior["seconds"] + summary["wall_seconds_this_run"] == 8
    assert summary["jobs"] == 4
    assert calls == ["1.1.1.1"]
