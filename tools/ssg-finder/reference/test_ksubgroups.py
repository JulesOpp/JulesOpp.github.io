"""Public hand groups and optional local d=1 k-index oracle checks."""
import os
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import subgroups as ts
import ksubgroups as ks


def hand(name):
    """Synthetic primitive groups; no catalogue data."""
    diagonal = {
        'P1': [], 'P-1': [(-1, -1, -1, -1)],
        'P2/m': [(-1, -1, 1, 1), (1, 1, -1, -1)],
    }
    gen = [ts.exact.format_op([[int(i == j)*a[i] for j in range(4)]
                              for i in range(4)], [0]*4) for a in diagonal[name]]
    return dict(sym=name, d=1, q=['0,0,g'] if name == 'P2/m' else ['a,b,g'],
                params=['g'] if name == 'P2/m' else ['a', 'b', 'g'],
                cen=[], gen=gen)


@pytest.mark.parametrize('n,k,count', [(4,2,15), (4,3,40), (6,2,63)])
def test_sublattices(n, k, count):
    assert len(ks.sublattices(n, k)) == count


@pytest.mark.parametrize('k,count', [(2,15), (3,40)])
def test_p1(k, count):
    assert len(ks.k_subgroups(hand('P1'), k, {})) == count


def test_inversion():
    rows = ks.k_subgroups(hand('P-1'), 2, {})
    # The plan's 30 counts t_index=1 (total index 2), the IT A1 analogue.
    assert len([r for r in rows if r['t_index'] == 1]) == 30
    assert len([r for r in rows if r['t_index'] == 2]) == 15
    assert len(rows) == 45


def test_brute():
    rec = hand('P2/m')
    assert ks.brute_subgroups(rec, 2) == {ks.row_key(rec, r) for r in ks.k_subgroups(rec, 2, {})}


def test_no_conjugate_classes():
    rec = hand('P2/m')
    rows = ks.k_subgroups(rec, 2, {})
    keys = {ks.row_key(rec, r) for r in rows}
    assert len(keys) == len(rows)
    for row in rows:
        orbit = ks.full_orbit(rec, row, 2)
        assert len(orbit) == row['class_size']
        assert len(orbit & keys) == 1


@pytest.mark.parametrize('name', ['P1', 'P-1', 'P2/m'])
def test_k1(name):
    rec = hand(name)
    old = ts.t_subgroups(rec, {})
    new = ks.k_subgroups(rec, 1, {})
    assert [{k:r[k] for k in old[0]} for r in new] == old


@pytest.mark.parametrize('param', ['a', 'a1', 'g2', 'A1'])
def test_rebase_implicit_q_multiplication(param):
    rec = dict(sym='synthetic', d=1, q=[f'-2{param},{param},1/3'],
               params=[param], cen=[], gen=[])
    rebased, S = ks.rebase(rec, ks.exact.build_ops(rec), ks.im.eye(4))
    assert S == ks.im.eye(4)
    # Identity rebasing preserves every component exactly, with explicit '*'.
    assert rebased['q'] == [f'-2*{param},{param},1/3']
    from fractions import Fraction
    from ssg_descent import eval_q
    expected = [[Fraction(-4, 5), Fraction(2, 5), Fraction(1, 3)]]
    assert eval_q(rec['q'], {param: Fraction(2, 5)}) == expected
    assert ks.exact.eval_q(rebased, {param: Fraction(2, 5)}) == expected


HAVE = (Path(__file__).resolve().parents[1]/'local/data').is_dir()


@pytest.mark.skipif(not HAVE, reason='local corpus absent')
def test_local_all_k2():
    db = ts.se.load_records(1)
    for rec in db.values():
        assert all(len(r['matches']) == 1 for r in ks.k_subgroups(rec, 2, db))


@pytest.mark.skipif(not HAVE, reason='local corpus absent')
def test_local_brute():
    db = ts.se.load_records(1)
    eligible = [r for r in db.values() if len({W for W,w in ts.exact.build_ops(r)})*16 <= 2048]
    for rec in random.Random(10).sample(eligible, 20):
        assert ks.brute_subgroups(rec, 2) == {ks.row_key(rec,r) for r in ks.k_subgroups(rec,2,{})}
