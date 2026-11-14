"""Public exact controls for the optional simultaneous affine E sieve."""
import os
import sys
from fractions import Fraction as F

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ssg_equiv as se
import affine_search as af

# spglib announces its error-handling migration on valid synthetic calls.
pytestmark = pytest.mark.filterwarnings("ignore:Set OLD_ERROR_HANDLING.*:DeprecationWarning")


def parallel(screw=False, phase=False):
    return dict(sym="public parallel P2", cls="", d=3,
                q=["0,0,a", "0,0,b", "0,0,c"], params=["a", "b", "c"],
                gen=["-x,-y," + ("z+1/2" if screw else "z") + "," +
                     ("t+1/2" if phase else "t") + ",u,v"], cen=[], nops=2)


def sieve(a, b, P=None):
    P = se.mat(se.im.eye(3)) if P is None else P
    return af.feasible_e_lattice(se.prep(a), se.prep(b), P, se.im.inv_rat(P))


def test_affine_projection_retains_congruences_without_saturation():
    # 2*x - 4*y = 2 gives x=1 mod 2, rather than all integers.
    lattice = af.integer_projection([[F(2), F(-4)]], [F(2)], 1, 1)
    assert lattice.contains([[1]]) and lattice.contains([[-3]])
    assert not lattice.contains([[0]]) and not lattice.contains([[2]])
    assert lattice.has_unimodular_residue(2)


def test_affine_projection_detects_inconsistent_diophantine_system():
    lattice = af.integer_projection([[F(2)]], [F(1)], 1, 1)
    assert lattice.empty
    assert not lattice.contains([[1]])


@pytest.mark.parametrize("prime", [2, 3])
def test_affine_projection_detects_singular_residue_coset(prime):
    lattice = af.AffineELattice([0, 0, 0, 1], [[prime, 0, 0, 0]], 2)
    assert not lattice.has_unimodular_residue(prime)
    assert lattice.contains([[prime, 0], [0, 1]])


def test_simultaneous_sieve_rejects_phase_without_solution():
    lattice = sieve(parallel(), parallel(phase=True))
    assert lattice.empty


def test_simultaneous_sieve_rejects_only_singular_internal_maps():
    lattice = sieve(parallel(phase=True), parallel())
    assert not lattice.empty
    assert not lattice.has_unimodular_residue(2)
    assert not lattice.contains(se.im.eye(3))


def test_simultaneous_sieve_keeps_screw_phase_compensation():
    a, b = parallel(screw=True), parallel(screw=True, phase=True)
    assert sieve(a,b).contains(se.im.eye(3))
    mp = dict(P=se.mat(se.im.eye(3)),p=[F(0)]*3,E=se.im.eye(3),
              M=[[F(0),F(0),F(1)],[F(0)]*3,[F(0)]*3],p4=[F(0)]*3)
    assert se.check_map(a,se.q_generic(a),b,mp)


@pytest.mark.parametrize("centred", ["none", "source", "target"])
def test_simultaneous_sieve_preserves_shear_origins_narrow_family_and_rational_P(centred):
    a = dict(sym="public narrow inversion",cls="",d=3,
             q=["a,0,0","0,a,0","0,0,a"],params=["a"],
             gen=["-x,-y,-z,-t,-u,-v"],cen=[],nops=2)
    b = dict(a,sym="public wide inversion",q=["a,b,c","d,e,f","g,h,i"],
             params=list("abcdefghi"),gen=["-x,-y,-z,-t+1/2,-u+1/3,-v+1/4"])
    if centred != "none":
        r = a if centred == "source" else b
        r.update(cen=["1/2,0,0,0,0,0"],nops=4)
    P=se.mat(se.im.eye(3))
    P[0][0]=F(2) if centred == "source" else F(1,2) if centred == "target" else F(1)
    E=[[1,1,0],[0,1,1],[0,0,1]]
    M=[[F(2),F(0),F(0)],[F(0)]*3,[F(0)]*3]
    mp=dict(P=P,p=[F(0)]*3,E=E,M=M,p4=[F(1,4),F(1,6),F(1,8)])
    assert se.check_map(a,se.q_generic(a),b,mp)
    lattice=sieve(a,b,P)
    assert lattice.contains(E)
    assert lattice.has_unimodular_residue(2) and lattice.has_unimodular_residue(3)


def test_simultaneous_sieve_skips_coupled_lattice():
    a=parallel()
    a.update(cen=["1/2,0,0,1/2,0,0"],nops=4)
    assert sieve(a,a) is None


def test_d1_sieved_stream_is_original_bounded_stream_subsequence(monkeypatch):
    a=dict(sym='synthetic screw',cls='',d=1,q=['0,0,a'],params=['a'],
           gen=['-x,-y,z+1/2,t+1/2'],cen=[])
    b=dict(a,gen=['-x,-y,z+1/2,t'])
    A,B=se.prep(a),se.prep(b)
    P=se.mat(se.im.eye(3))
    candidates=[[[1]],[[-1]]]
    monkeypatch.setattr(se,"_e_candidates",lambda *args: iter(candidates))
    lattice=sieve(a,b)
    assert list(af.e_candidates(A,B,P,P)) == [E for E in candidates if lattice.contains(E)]


def test_equivalent_optional_callback_and_default(monkeypatch):
    import numpy as np
    a,b=parallel(screw=True),parallel(screw=True,phase=True)
    monkeypatch.setattr(se,"normaliser_candidates",lambda *args: [np.eye(3,dtype=int)])
    calls=[]
    def callback(*args):
        calls.append(True)
        return af.e_candidates(*args)
    plain=se.equivalent(a,se.q_generic(a),b)
    filtered=se.equivalent(a,se.q_generic(a),b,e_candidates=callback)
    assert calls and plain is not None and filtered is not None
    assert se.check_map(a,se.q_generic(a),b,filtered)


def test_fixture_elapsed_time_uses_monotonic_clock(tmp_path, monkeypatch):
    import make_subgroup_fixtures as m
    monkeypatch.setattr(m,"fixture",lambda label: {"label":label,"classes":[]})
    monkeypatch.setattr(m,"source_signature",lambda d: "public-clock-test")
    ticks=iter([100.0,103.5])
    monkeypatch.setattr(m.time,"monotonic",lambda: next(ticks))
    monkeypatch.setattr(m.time,"time",lambda: (_ for _ in ()).throw(AssertionError("wall clock used")))
    result=m.work("1.1.1.1",out=tmp_path)
    assert result[3] == 3.5


@pytest.mark.parametrize('d',[2,3])
def test_public_mirror_affine_coset(d, monkeypatch):
    # Independently constructed mirror: thirds in q and half-period phase.
    # Integrality forces E[0,0]=2 (mod 3), E[1,0]=0 (mod 3);
    # the mirror phase forces that first column to be even,odd,(even).
    tails=['t','u','v'][:d]
    a=dict(sym='synthetic mirror',cls='',d=d,params=list('abc')[:d],
           q=['1/3,0,a','0,0,b','0,0,c'][:d],cen=[],
           gen=['x,-y,z,'+','.join(x+('+1/2' if i==0 else '') for i,x in enumerate(tails))])
    b=dict(a,q=['2/3,0,a','0,0,b','0,0,c'][:d],
           gen=['x,-y,z,'+','.join(x+('+1/2' if i==1 else '') for i,x in enumerate(tails))])
    E=se.im.eye(d); E[0][:2]=[2,1];E[1][:2]=[3,1]
    M=[[0]*3 for _ in range(d)]; M[1][0]=-1
    known=dict(P=se.im.eye(3),p=[0]*3,E=E,M=M,p4=[0]*d)
    q=se.q_generic(a)
    assert se.check_map(a,q,b,known)
    A,B=se.prep(a),se.prep(b)
    P=se.mat(se.im.eye(3))
    lattice=af.feasible_e_lattice(A,B,P,P)
    assert lattice.contains(E)
    # Isolate the exact internal search at the known external identity.
    import numpy as np
    N=se.mm(B.T,se.im.inv_rat(A.T))
    monkeypatch.setattr(se,'normaliser_candidates',lambda *args:[np.array(N,dtype=int)])
    result=se.equivalent(a,q,b,e_candidates=af.e_candidates)
    assert result is not None
    assert se.check_map(a,q,b,result)
