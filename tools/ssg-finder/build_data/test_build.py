"""Tests for the BYU data pipeline (byu_tables.py, verify_html.py).

The parser edge-case tests use invented blocks in BYU's text format (no BYU content).
Everything that reads real BYU data skips when tools/ssg-finder/local/ is absent
(decision D1: BYU data is never committed).

    /opt/anaconda3/bin/python3 -I -m pytest tools/ssg-finder/build_data
"""
import json, os, random, sys
from fractions import Fraction as Fr

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import byu_tables as B  # noqa: E402

LOCAL = B.LOCAL
ZIP = B.ZIP_DEFAULT
DATA = os.path.join(LOCAL, "data")
CACHE = os.path.join(LOCAL, "cache")
CACHED_17 = ["14.3.16.6", "55.1.9.2", "55.1.9.3", "55.1.9.4", "55.1.9.5", "55.1.10.6",
             "55.3.122.93", "55.3.122.94", "55.3.122.95", "55.3.122.96", "55.3.122.97",
             "55.3.122.98", "55.3.123.99", "55.3.123.100", "55.3.123.101", "62.1.9.2", "62.1.9.6"]

need_zip = pytest.mark.skipif(not os.path.exists(ZIP), reason="local/raw/ssg_tables.zip absent (BYU data is local-only)")
need_data = pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "index.json")),
                               reason="local/data absent: run build_data/byu_tables.py")


# ------------------------------------------------------------ synthetic parser edge cases
def test_lin_parsing():
    assert B.lin("-2a1+a4", None) == ({"a1": Fr(-2), "a4": Fr(1)}, Fr(0))
    assert B.lin("g+1/2", {"g"}) == ({"g": Fr(1)}, Fr(1, 2))
    assert B.lin("2*a3", None) == ({"a3": Fr(2)}, Fr(0))
    assert B.lin("-1/2", set()) == ({}, Fr(-1, 2))
    assert B.lin_vec("-X+1/2", list("XYZT")) == ([-1, 0, 0, 0], Fr(1, 2))
    with pytest.raises(B.BuildError):
        B.lin("a1+q", {"a1"})
    with pytest.raises(B.BuildError):
        B.lin("a1*+2", None)


def test_group_block_wrapped_with_condition_prefixes():
    # invented group: wrapped centring list, operators wrapped mid-list and mid-operator,
    # reflection conditions with digit/minus prefixes, wrapped too
    lines = """999.3.1.7 Pfake(a1,b1,0)000(-a1,b1,0)00s(0,0,g2)000
(A1,0,0); (0,B1,0);
(0,0,G2)
(0,0,0,0,0,0), (0,0,0,1/2,1/2,0),
(1/2,1/2,0,0,0,0)
X,Y,Z,T,U,V; -X,-Y,Z,-T,-U+1/2,V;
X+1/2,-Y,-Z,T,-U,
-V+1/2; -X+1/2,Y,-Z,-T,U+1/2,-V+1/2
HKLMNP:M+N=2n; H-H00N0:N=2n; 2KKLM-MP:P=2n;
0KL2NNP:P=2n""".split("\n")
    g = B.parse_group_block(lines, 3)
    assert g["label"] == "999.3.1.7"
    assert g["sym"] == "Pfake(a1,b1,0)000(-a1,b1,0)00s(0,0,g2)000"
    assert g["qs"] == [["A1", "0", "0"], ["0", "B1", "0"], ["0", "0", "G2"]]
    assert len(g["cen"]) == 3 and g["cen"][1] == [0, 0, 0, Fr(1, 2), Fr(1, 2), 0]
    assert g["ops"] == ["X,Y,Z,T,U,V", "-X,-Y,Z,-T,-U+1/2,V",
                        "X+1/2,-Y,-Z,T,-U,-V+1/2", "-X+1/2,Y,-Z,-T,U+1/2,-V+1/2"]
    assert g["refl"].startswith("HKLMNP:") and "H-H00N0:" in g["refl"]
    assert "2KKLM-MP:" in g["refl"] and g["refl"].endswith("0KL2NNP:P=2n")


def test_group_block_without_conditions_and_d1():
    g = B.parse_group_block(["998.1.1.1 Pfake(a,b,g)0", "(A,B,G)", "(0,0,0,0)", "X,Y,Z,T; -X,-Y,-Z,-T"], 1)
    assert g["ops"] == ["X,Y,Z,T", "-X,-Y,-Z,-T"] and g["refl"] == ""
    assert B.label_q(g["sym"], 1) == (["a,b,g"], ["a", "b", "g"])


def test_group_block_condition_starting_with_digit():
    g = B.parse_group_block(["997.1.1.1 Pfake(0,0,g)s", "(0,0,G)", "(0,0,0,0)",
                             "X,Y,Z,T; -X,-Y,Z,T+1/2", "00LM:M=2n"], 1)
    assert g["ops"] == ["X,Y,Z,T", "-X,-Y,Z,T+1/2"] and g["refl"] == "00LM:M=2n"


def test_group_block_rejects_garbage():
    with pytest.raises(B.BuildError):
        B.parse_group_block(["996.1.1.1 Pfake(a,b,g)0", "(A,B,G)", "(0,0,0,0)", "X,Y,Z; -X,-Y,-Z,-T"], 1)
    with pytest.raises(B.BuildError):
        B.parse_group_block(["notalabel Pfake(a,b,g)0", "(A,B,G)", "(0,0,0,0)", "X,Y,Z,T"], 1)


def test_lattice_block_wrapped():
    lines = """3.999 Pfake(a1,b1,1/2)(-a1,b1,1/2)(0,0,g2)
A1=a1, A2=a2, A3=2a3+a4+a5,
A4=a4-a5, A5=a4+a5, A6=a6
Q1=(A1,0,0), Q2=(0,B1,0),
Q3=(0,0,G2),
where A1=a1, B1=b1,
G2=2g2
(0,0,1/2,0,1/2,0), (0,0,0,1/2,1/2,0),
(0,0,1/2,1/2,0,0)""".split("\n")
    L = B.parse_lattice_block(lines, 3)
    assert L["P"][2] == [0, 0, 2, 1, 1, 0] and L["P"][3] == [0, 0, 0, 1, -1, 0]
    assert L["Q"] == [["A1", "0", "0"], ["0", "B1", "0"], ["0", "0", "G2"]]
    assert L["where"] == {"A1": "a1", "B1": "b1", "G2": "2g2"}
    assert len(L["cen"]) == 3


def test_label_q_params_order():
    qs, ps = B.label_q("Pfake(1/2,b,g)(a,1/2,g+1/2)(a+1/2,b+1/2,1/2)", 3)
    assert qs == ["1/2,b,g", "a,1/2,g+1/2", "a+1/2,b+1/2,1/2"] and ps == ["b", "g", "a"]
    with pytest.raises(B.BuildError):
        B.label_q("Pfake(a,b,g)0", 2)


def test_fmt_op_roundtrip():
    W = [[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, -1]]
    s = B.fmt_op(W, [Fr(1, 2), 0, Fr(3, 4), 0], 4)
    assert s == "-x+1/2,y,z+3/4,z-t"
    assert B.parse_op(s, 4) == (W, [Fr(1, 2), 0, Fr(3, 4), 0])
    assert B.parse_op("x-1/4,y,z,t", 4)[1][0] == Fr(3, 4)  # reduced to [0,1)


def test_inverse():
    A = [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]
    assert B.matmul(A, B.inverse(A)) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


# ------------------------------------------------------------ real tables (local only)
@pytest.fixture(scope="module")
def tables():
    return B.load_zip(ZIP)


@need_zip
def test_counts(tables):
    groups, lats = tables
    assert len(groups) == 16697 and len(lats) == 322
    assert [sum(g["d"] == d for g in groups) for d in (1, 2, 3)] == [775, 3338, 12584]
    assert [sum(L["d"] == d for L in lats.values()) for d in (1, 2, 3)] == [24, 83, 215]


@need_zip
def test_real_condition_prefixes_split_cleanly(tables):
    groups, _ = tables
    for prefix in ("H-H00N0:", "2KKLM-MP:", "0KL2NNP:"):
        hits = [g for g in groups if prefix in g["refl"]]
        assert hits, prefix
        for g in hits:
            assert all(":" not in o and len(o.split(",")) == 3 + g["d"] for o in g["ops"])


@need_zip
def test_real_wrapped_blocks(tables):
    groups, _ = tables
    # the largest groups wrap their operator lists over many lines
    big = max(groups, key=lambda g: len(g["ops"]))
    assert len(big["ops"]) >= 48


@need_zip
def test_convert_sample_invariants(tables):
    """Full invariant check on a deterministic sample (the build checks all groups)."""
    groups, lats = tables
    rng = random.Random(7)
    sample = groups[::97] + [g for g in groups if g["label"] in CACHED_17]
    for g in sample:
        rec, info = B.convert(g, lats[".".join(g["label"].split(".")[1:3])], rng)
        assert rec["nops"] == len(B.canonical_ops(rec))
        assert info["sg"] == int(g["label"].split(".")[0])


@need_zip
def test_invariant_violation_is_detected(tables):
    groups, lats = tables
    g = dict(next(x for x in groups if x["label"] == "55.3.122.96"))
    L = lats["3.122"]
    # drop an operator -> not closed
    bad = dict(g, ops=g["ops"][:-1])
    with pytest.raises(B.BuildError):
        B.convert(bad, L, random.Random(1))
    # corrupt the label q -> M + eps Q != Q R
    bad = dict(g, sym=g["sym"].replace("(-a1,b1,0)", "(a1,b1,0)"))
    with pytest.raises(B.BuildError):
        B.convert(bad, L, random.Random(1))


@need_data
def test_written_data():
    idx = json.load(open(os.path.join(DATA, "index.json")))
    assert idx["count"] == 16697 == len(idx["labels"]) == len(set(idx["labels"]))
    d55 = json.load(open(os.path.join(DATA, "sg055.json")))
    assert d55["sg"] == 55
    r = d55["groups"]["55.3.122.96"]
    assert set(r) == {"sym", "cls", "d", "q", "params", "cen", "gen", "nops"}
    assert r["cls"] == "3.122" and r["d"] == 3 and r["params"] == ["a1", "b1", "g2"]
    assert r["q"] == ["a1,b1,0", "-a1,b1,0", "0,0,g2"]
    # basic setting is primitive here: the supercentring (0,0,0,1/2,1/2,0) maps to an integer vector
    assert r["cen"] == [] and r["nops"] == 8
    ops = B.canonical_ops(r)
    assert B.fmt_op(*ops[0], 6) == "x,y,z,t,u,v" and len(ops) == 8
    # Fmmm: the F-centring survives into the basic setting
    f = json.load(open(os.path.join(DATA, "sg069.json")))["groups"]
    for c in f.values():
        assert {"0,1/2,1/2", "1/2,0,1/2", "1/2,1/2,0"} <= {",".join(x.split(",")[:3]) for x in c["cen"]}
        assert c["nops"] == len(B.canonical_ops(c)) and c["nops"] % (1 + len(c["cen"])) == 0


@need_data
@pytest.mark.parametrize("label", CACHED_17)
def test_cached_html_matches(label):
    if not os.path.exists(os.path.join(CACHE, label + ".html")):
        pytest.skip("page not cached")
    import verify_html as V
    sg = int(label.split(".")[0])
    rec = json.load(open(os.path.join(DATA, f"sg{sg:03d}.json")))["groups"][label]
    assert V.compare(label, rec) == []
