"""Exact matrix oracle regressions; local corpus checks skip when unavailable."""
import os
import sys
from itertools import combinations
sys.path.insert(0, os.path.dirname(__file__))
import pytest
import ssg_equiv as se

byu = pytest.mark.skipif(not os.path.isdir(se.DATA), reason='local BYU records unavailable')


def hand():
    return dict(d=2, params=['a1','b1','g1','a2','b2','g2'],
                q=['a1,b1,g1','a2,b2,g2'], gen=[], cen=[], sym='P1', cls='2.1')


def test_parallel_q_mixing():
    a = hand()
    b = dict(a, q=['a1+a2,b1+b2,g1+g2','a2,b2,g2'])
    m = se.equivalent(a, se.q_generic(a), b)
    assert m is not None and abs(se.im.det(m['E'])) == 1
    assert se.check_map(a, se.q_generic(a), b, m)


def test_generic_cycle_and_symbolic_rank():
    r = hand()
    q = se.q_generic(r)
    assert len(q) == 2 and q[1][0] == se.Fr(4010, 10007)
    assert se.rank_Mstar(r, q) == 5
    bad = dict(r, q=['a1,b1,g1','2a1,2b1,2g1'])
    assert se.rank_Mstar(bad, se.q_generic(bad)) == 4


def test_forced_q_mixing():
    a = dict(hand(),params=['a','b'],q=['a,0,0','0,b,0'])
    b = dict(a,q=['a,b,0','0,b,0'])
    q = se.q_generic(a)
    m = se.equivalent(a,q,b)
    assert m is not None and se.check_map(a,q,b,m)
    assert any(m['E'][i][j] for i in range(2) for j in range(2) if i != j)
    known = dict(P=se.im.eye(3),p=[0]*3,E=[[1,1],[0,1]],M=[[0]*3 for _ in range(2)],p4=[0,0])
    assert se.check_map(a,q,b,known)


def test_affine_transformed_d23():
    # Synthetic inversion representations with phase/origin and q-basis shear.
    for d in (2,3):
        a = dict(d=d,params=['a','b','c'][:d],q=['a,0,0','0,b,0','0,0,c'][:d],
                 gen=['-x,-y,-z,'+','.join('-'+t for t in 'tuv'[:d])],cen=[],sym='P-1',cls='')
        E = [[1,1],[0,1]] if d==2 else [[1,1,0],[0,1,1],[0,0,1]]
        M = [[se.Fr(i==j) for j in range(3)] for i in range(d)]
        p,p4 = [se.Fr(1,8),se.Fr(1,6),se.Fr(1,4)],[se.Fr(i+1,12) for i in range(d)]
        S=se._S(se.im.eye(3),M,E); Si=se._Sinv(se.im.eye(3),M,E)
        gen=[]
        for W,w in se.exact.build_ops(a)[1:]:
            Wp=se.mm(se.mm(S,W),Si); shift=p+p4
            wp=[se.exact.frac1(x+y-z) for x,y,z in zip(se.mv(S,w),shift,se.mv(Wp,shift))]
            gen.append(se.exact.format_op(Wp,wp))
        # Q' = E Q + M, expressions constructed exactly.
        qs=[]
        for i in range(d):
            qs.append(','.join(str(int(M[i][j]))+''.join('+'+str(E[i][k])+'*'+a['params'][k] for k in range(d) if E[i][k] and k==j) for j in range(3)))
        b=dict(a,gen=gen,q=qs)
        q=se.q_generic(a)
        assert se.check_map(a,q,b,dict(P=se.im.eye(3),p=p,E=E,M=M,p4=p4))
        m=se.equivalent(a,q,b)
        assert m is not None and se.check_map(a,q,b,m)


def test_affine_family_torsion_obstruction():
    a=dict(hand(),params=['a','b'],q=['a,0,0','b,0,0'])
    b=dict(a,q=['a,1/2,0','b,0,0'])
    assert se._family_torsion(se.prep(a)) == (2,1)
    assert se._family_torsion(se.prep(b)) == (2,2)
    assert se.equivalent(a,se.q_generic(a),b) is None
    # A narrower subgroup family may sit inside a larger canonical family.
    c=dict(hand(),params=['a','b','c'],q=['a,c,0','b,0,0'])
    m=se.equivalent(b,se.q_generic(b),c)
    assert m is not None and se.check_map(b,se.q_generic(b),c,m)


def test_d1_unchanged():
    r=dict(d=1,params=['g'],q=['0,0,g'],gen=['-x,-y,z,t'],cen=[],sym='P2',cls='')
    q=se.q_generic(r)
    m=se.equivalent(r,q,r)
    assert len(q)==1 and m['E']==[[1]] and m['M']==[[0,0,0]] and m['p4']==[0]
    assert 'eta' not in m and 'm' not in m
    assert se.check_map(r,q,r,m)


def test_zero_parameter_family():
    assert se.solve_lin([], [se.Fr(1)], 0) is None
    r=dict(d=1,params=[],q=['0,0,1/2'],gen=[],cen=[],sym='P1',cls='')
    with pytest.raises(ValueError,match='outside'):
        se.rank_Mstar(r,[[0,0,0]])
    known=dict(P=se.im.eye(3),p=[0,0,0],E=[[1]],M=[[0,0,0]],p4=[0])
    assert not se.check_map(r,[[0,0,0]],r,known)


def test_external_orientation():
    r=hand()
    reversing=dict(P=[[-1,0,0],[0,1,0],[0,0,1]],p=[0]*3,E=se.im.eye(2),M=[[0]*3 for _ in range(2)],p4=[0,0])
    assert not se.check_map(r,se.q_generic(r),r,reversing)


@byu
def test_chiral_pair_distinct():
    db=records(2)
    a,b=db['98.2.71.7'],db['98.2.71.8']
    assert se.equivalent(a,se.q_generic(a),b) is None


def records(d):
    db = se.load_records(d)
    if not db:
        pytest.skip('local BYU records unavailable')
    return db


def sample(db):
    """First record of every class and every basic space group, stable label order."""
    chosen = {}
    for label, r in db.items():
        chosen.setdefault(('cls', r['cls']), label)
        chosen.setdefault(('sg', label.split('.')[0]), label)
    return {label: db[label] for label in sorted(set(chosen.values()))}


@byu
@pytest.mark.parametrize('d', [2, 3])
def test_rank_Mstar(d):
    for label, r in records(d).items():
        assert se.rank_Mstar(r, se.q_generic(r)) == 3+d, label


@byu
@pytest.mark.parametrize('d', [2, 3])
def test_self_equivalent(d):
    db = records(d)
    if d == 3 and os.environ.get('SSG_FULL') != 'd3':
        db = sample(db)
    for label, r in db.items():
        q = se.q_generic(r)
        m = se.equivalent(r, q, r)
        assert m is not None and se.check_map(r, q, r, m), label


@byu
@pytest.mark.parametrize('d', [2, 3])
def test_distinct_not_equivalent(d):
    db = records(d)
    if os.environ.get('SSG_FULL') != f'd{d}':
        db = sample(db)
    buckets = {}
    for label, r in db.items():
        buckets.setdefault((label.split('.')[0],se.invariant(r, se.q_generic(r))), []).append(label)
    pairs = [pair for labels in buckets.values() for pair in combinations(labels, 2)]
    if os.environ.get('SSG_FULL') != f'd{d}' and len(pairs) > 100:
        # Keep default regressions quick across the complete ordered pair list.
        # Acceptance tests every pair via check_equiv_d23.py's resumable harness.
        pairs = [pairs[i * len(pairs) // 100] for i in range(100)]
    for a,b in pairs:
        assert se.equivalent(db[a], se.q_generic(db[a]), db[b]) is None, (a,b)


@byu
def test_star_permutation():
    r = records(3)['221.3.213.10']
    q = se.q_generic(r)
    # The [111] threefold is a group operator, hence an explicit self-map.
    P = [[0,0,1],[1,0,0],[0,1,0]]
    for W,w in se.exact.build_ops(r):
        if [list(row[:3]) for row in W[:3]] == P:
            m = dict(P=P,p=list(w[:3]), E=[list(row[3:]) for row in W[3:]],
                     M=[list(row[:3]) for row in W[3:]],p4=list(w[3:]))
            # Three chosen arms of a four-arm tetrahedral star: the omitted
            # arm is -sum(q_i), hence E need not be a permutation matrix.
            assert abs(se.im.det(m['E'])) == 1
            assert se.mm(se.mm(m['E'],m['E']),m['E']) == se.im.eye(3)
            assert m['E'] != se.im.eye(3)
            assert se.check_map(r,q,r,m)
            return
    pytest.fail('threefold operator missing')


def test_resume_binds_arithmetic_and_tolerates_interruption(tmp_path,monkeypatch):
    import json
    import check_equiv_d23 as runner
    root=tmp_path/'reference'
    root.mkdir()
    monkeypatch.setattr(se,'HERE',str(root))
    monkeypatch.setattr(se,'load_records',lambda d: {'synthetic':hand()})
    arithmetic=tmp_path/'intmat.py'
    arithmetic.write_text('arithmetic version one')
    monkeypatch.setattr(se.im,'__file__',str(arithmetic))
    calls=[]
    check=runner._check
    def record(case):
        calls.append(case)
        return check(case)
    monkeypatch.setattr(runner,'_check',record)
    runner.run(2,'self')
    checkpoint=tmp_path/'local/fixtures/subgroups/_oracle_d2_self_full.jsonl'
    assert len(calls)==1
    with checkpoint.open('a') as f:
        f.write('{interrupted\n')
    runner.run(2,'self')
    assert len(calls)==1
    arithmetic.write_text('arithmetic version two')
    runner.run(2,'self')
    assert len(calls)==2
    rows=[json.loads(line) for line in checkpoint.read_text().splitlines() if line.startswith('{"case"')]
    assert len(rows)==2 and rows[0]['signature'] != rows[1]['signature']


def test_operator_family_constant_is_exact():
    r=dict(d=1,params=['a'],q=['1/2,0,a'],gen=['-x,-y,z,t'],cen=[],sym='P2',cls='')
    with pytest.raises(AssertionError,match='q-family constant'):
        se.Prep(r)
    valid=dict(r,gen=['-x,-y,z,t-x'])
    assert se.Prep(valid).d==1


def test_check_map_validates_source_q_family_and_shape():
    a=dict(d=1,params=['a'],q=['a,0,0'],gen=[],cen=[],sym='P1',cls='')
    b=dict(a,params=['a','b'],q=['a,b,0'])
    identity=dict(P=se.im.eye(3),p=[0,0,0],E=[[1]],M=[[0,0,0]],p4=[0])
    assert not se.check_map(a,[[0,1,0]],b,identity)
    assert se.check_map(a,[[1,0,0]],b,identity)
    for Q in ([], [[1,0]], [[1,0,0],[0,0,0]]):
        assert not se.check_map(a,Q,b,identity)


def _parallel_d3(two_offsets=False,centred=False):
    return dict(d=3,params=['a','b','c'],
                q=['a,1/2,0','b,0,1/2' if two_offsets else 'b,0,0','c,0,0'],
                gen=[],cen=['1/2,0,0,0,0,0'] if centred else [],sym='P1',cls='')


def test_physical_offset_quotient_rejects_before_search(monkeypatch):
    a,b=_parallel_d3(),_parallel_d3(two_offsets=True)
    assert se._family_torsion(se.prep(a)) == se._family_torsion(se.prep(b)) == (3,2)
    assert se._physical_offset_quotient(se.prep(a)) == (1,(2,))
    assert se._physical_offset_quotient(se.prep(b)) == (1,(2,2))
    def forbidden(*args):
        pytest.fail('different physical offset quotients must reject before normalizer/E enumeration')
    monkeypatch.setattr(se,'normaliser_candidates',forbidden)
    assert se.equivalent(a,se.q_generic(a),b) is None


def test_physical_offset_quotient_centred_unimodular_invariance():
    a=_parallel_d3(two_offsets=True,centred=True)
    P=[[1,0,0],[0,1,1],[0,0,1]]
    E=[[1,1,0],[0,1,1],[0,0,1]]
    M=[[2,0,1],[0,1,0],[2,1,1]]
    Pi=se.im.inv_rat(P)
    b0,dirs=se.family(a)
    c=se.mm([[x+y for x,y in zip(r,t)] for r,t in zip(se.mm(E,b0),M)],Pi)
    bd=[se.mm(se.mm(E,b),Pi) for b in dirs]
    q=[]
    for i in range(3):
        components=[]
        for j in range(3):
            expression=se.exact.fstr(c[i][j])
            for parameter,b in zip(a['params'],bd):
                if b[i][j]: expression+='+'+se.exact.fstr(b[i][j])+'*'+parameter
            components.append(expression.replace('+-','-'))
        q.append(','.join(components))
    b=dict(a,q=q)
    known=dict(P=P,p=[0,0,0],E=E,M=M,p4=[0,0,0])
    assert se._embedded_offset_compatible(se.prep(a),se.prep(b),P,Pi)
    assert se.check_map(a,se.q_generic(a),b,known)
    assert se._physical_offset_quotient(se.prep(a)) == (1,(2,2))
    assert se._physical_offset_quotient(se.prep(b)) == (1,(2,2))


def test_physical_offset_quotient_narrower_family_and_skip_scope():
    from types import SimpleNamespace
    a=_parallel_d3()
    b=dict(a,params=['a','b','c','h'],q=['a,h,0','b,0,0','c,0,0'])
    assert se._physical_offset_quotient(se.prep(a)) == (1,(2,))
    assert se._physical_offset_quotient(se.prep(b)) == (2,())
    m=se.equivalent(a,se.q_generic(a),b)
    assert m is not None and se.check_map(a,se.q_generic(a),b,m)
    mixed=SimpleNamespace(rec=dict(a,cen=['0,0,0,1/2,0,0']))
    assert se._physical_offset_quotient(mixed) is None


def test_rational_generator_lattice_has_no_implicit_ambient_units():
    B=se._rational_generator_lattice([[2,0],[0,3]],2)
    assert abs(se.im.det(B))==6
    B=se._rational_generator_lattice([[se.Fr(1,2),0],[0,se.Fr(3,2)]],2)
    assert abs(se.im.det(B))==se.Fr(3,4)


def test_physical_offset_quotient_narrower_same_physical_rank():
    a=dict(d=2,params=['a'],q=['a,0,1/2','0,a,0'],gen=[],cen=[],sym='P1',cls='')
    b=dict(a,params=['a','b'],q=['a,0,1/2','0,b,0'])
    assert se._physical_offset_quotient(se.prep(a)) == (2,(2,))
    assert se._physical_offset_quotient(se.prep(b)) == (2,(2,))
    known=dict(P=se.im.eye(3),p=[0,0,0],E=se.im.eye(2),M=[[0,0,0],[0,0,0]],p4=[0,0])
    q=se.q_generic(a)
    assert se._embedded_offset_compatible(se.prep(a),se.prep(b),se.im.eye(3),se.im.eye(3))
    assert se.check_map(a,q,b,known)
    result=se.equivalent(a,q,b)
    assert result is not None and se.check_map(a,q,b,result)


def _parallel_offset_records_p222():
    a=dict(d=3,params=['a','b','c'],q=['a,1/2,0','b,0,0','c,0,0'],
           gen=['x,-y,-z,t-y,u,v','-x,y,-z,-t+y,-u,-v'],cen=[],sym='P222',cls='')
    b=dict(a,q=['a,1/2,1/2','b,0,0','c,0,0'],
           gen=['x,-y,-z,t-y-z,u,v','-x,y,-z,-t+y,-u,-v'])
    return a,b


def test_embedded_offsets_reject_p222_before_e_search(monkeypatch):
    a,b=_parallel_offset_records_p222()
    assert se._physical_offset_quotient(se.prep(a)) == se._physical_offset_quotient(se.prep(b)) == (1,(2,))
    assert not se._embedded_offset_compatible(se.prep(a),se.prep(b),se.im.eye(3),se.im.eye(3))
    def forbidden(*args):
        pytest.fail('embedded offset mismatch must reject every external candidate before E enumeration')
    monkeypatch.setattr(se,'_e_candidates',forbidden)
    assert se.equivalent(a,se.q_generic(a),b) is None


def test_embedded_offsets_allow_external_shear_in_p1():
    a,b=_parallel_offset_records_p222()
    a,b=dict(a,gen=[],sym='P1'),dict(b,gen=[],sym='P1')
    P=[[1,0,0],[0,1,-1],[0,0,1]]
    Pi=se.im.inv_rat(P)
    assert se._embedded_offset_compatible(se.prep(a),se.prep(b),P,Pi)
    known=dict(P=P,p=[0,0,0],E=se.im.eye(3),M=[[0]*3 for _ in range(3)],p4=[0,0,0])
    q=se.q_generic(a)
    assert se.check_map(a,q,b,known)
    result=se.equivalent(a,q,b)
    assert result is not None and se.check_map(a,q,b,result)


def test_embedded_offsets_require_primitive_external_unimodularity():
    a=_parallel_d3(centred=True)
    b=_parallel_d3()
    A,B=se.prep(a),se.prep(b)
    assert not se._embedded_offset_compatible(A,B,se.im.eye(3),se.im.eye(3))
    P=[[2,0,0],[0,1,0],[0,0,1]]
    assert se._embedded_offset_compatible(A,B,P,se.im.inv_rat(P))
    P=[[4,0,0],[0,1,0],[0,0,1]]
    assert not se._embedded_offset_compatible(A,B,P,se.im.inv_rat(P))
    from types import SimpleNamespace
    mixed=SimpleNamespace(rec=dict(a,cen=['0,0,0,1/2,0,0']))
    assert se._embedded_offset_compatible(mixed,B,se.im.eye(3),se.im.eye(3))
