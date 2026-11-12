#!/usr/bin/env python3
"""Exact bounded-search oracle for affine (3+d)D equivalence, d=1,2,3.

Maps have x'=Px+p, t'=Mx+Et+p4. Q is a d×3 exact matrix.
E comes from the saturated integer lattice of intertwining and symbolic
family constraints, with lattice coefficients bounded by two. P uses the
existing bounded external normaliser. Every accepted map is checked exactly;
float spglib motif standardisation is only used to propose external settings.
"""
import functools, glob, itertools, json, math, os, sys

for _thread_var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_thread_var] = "1"
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # also under python -I
import exact  # noqa: E402
import intmat as im  # noqa: E402

NBOUND = 2
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "local", "data")


# ================================================================ small exact linear algebra
def mat(A):
    return [[Fr(x) for x in r] for r in A]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def mv(A, v):
    return [sum(a * x for a, x in zip(r, v)) for r in A]


def vm(v, A):
    return [sum(v[k] * A[k][j] for k in range(len(v))) for j in range(len(A[0]))]


def tr(A):
    return [list(r) for r in zip(*A)]


def is_int(x):
    return Fr(x).denominator == 1


def den(xs):
    return exact.lcm(1, *[Fr(x).denominator for x in xs])


def rref(A):
    """Reduced row echelon form of a rational matrix; returns (R, pivot columns)."""
    M = [[Fr(x) for x in r] for r in A]
    piv, r = [], 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    return M[:r], piv


def nullspace(A, ncol):
    """Rational basis (list of vectors) of {x : A x = 0}."""
    if not A:
        return [[Fr(int(i == j)) for j in range(ncol)] for i in range(ncol)]
    R, piv = rref(A)
    out = []
    for f in (c for c in range(ncol) if c not in piv):
        x = [Fr(0)] * ncol
        x[f] = Fr(1)
        for row, p in zip(R, piv):
            x[p] = -row[f]
        out.append(x)
    return out


def solve_lin(A, b, ncol):
    """One rational solution of A x = b (free variables 0), or None."""
    if not A:
        return None if any(b) else [Fr(0)] * ncol
    R, piv = rref([list(r) + [Fr(bi)] for r, bi in zip(A, b)])
    if ncol in piv:
        return None
    x = [Fr(0)] * ncol
    for row, p in zip(R, piv):
        x[p] = row[ncol]
    return x


def solve_int(M, r):
    """One integer solution u of M u = r (M, r rational), or None."""
    n = len(M[0])
    D = den([x for row in M for x in row])
    Mi = [[int(x * D) for x in row] for row in M]
    ri = [Fr(x) * D for x in r]
    if not all(is_int(x) for x in ri):
        return None
    Dm, U, V, _ = im.snf(Mi)
    Ur = [sum(u * x for u, x in zip(Ur_, ri)) for Ur_ in U]
    w = [Fr(0)] * n
    for i in range(len(Mi)):
        d = Dm[i][i] if i < n else 0
        if d:
            if Ur[i] % d:
                return None
            w[i] = Ur[i] / d
        elif Ur[i] != 0:
            return None
    return [int(x) for x in mv(V, w)]


def mixed_solve(G, H, h, a, b):
    """Solve G s + H t ≡ h (mod Z^k) with s real (len a) and t integer (len b).
    G: k×a, H: k×b (rational; rows may be empty lists), h: k.  Returns (s, t) or None."""
    k = len(h)
    if k == 0:
        return [Fr(0)] * a, [0] * b
    # left null space of G: rows n with n G = 0; then v ∈ colspace(G) + Z^k  <=>  N v ∈ N Z^k
    Nn = nullspace(tr(G), k) if a else [[Fr(int(i == j)) for j in range(k)] for i in range(k)]
    if Nn:
        NH = mm(Nn, H) if b else [[] for _ in Nn]
        Mx = [list(n) + list(nh) for n, nh in zip(Nn, NH)]
        u = solve_int(Mx, mv(Nn, h))
        if u is None:
            return None
        z, t = u[:k], u[k:]
    else:
        z, t = [0] * k, [0] * b
    rhs = [h[i] - z[i] - sum(H[i][j] * t[j] for j in range(b)) for i in range(k)]
    s = solve_lin(G, rhs, a) if a else []
    assert s is not None, "mixed congruence solution must satisfy the real linear system"
    return s, t


def int_lattice_of_real(G):
    """Basis (columns, as list of vectors) of {y : G y ∈ Z^k}, G rational k×a of full column rank."""
    a = len(G[0])
    D = den([x for r in G for x in r])
    Gi = [[int(x * D) for x in r] for r in G]
    Dg, U, V, _ = im.snf(Gi)
    assert all(Dg[i][i] for i in range(a)), "G not of full column rank"
    return [[Fr(V[r][i] * D, Dg[i][i]) for r in range(a)] for i in range(a)]


# ================================================================ group data
def q_generic(rec):
    """Exact d×3 generic Q, cycling seeds with +1009 per three parameters."""
    if rec["d"] not in (1, 2, 3):
        raise ValueError("supported dimensions are d = 1, 2, 3")
    return exact.eval_q(rec, {p: Fr(exact.GENERIC_SEEDS[i % 3] + 1009*(i//3),
                                   exact.P_GENERIC) for i,p in enumerate(rec["params"])})


def family(rec):
    return _family(tuple(rec["q"]), tuple(rec["params"]))


@functools.lru_cache(maxsize=2048)
def _family(q, params):
    rec = dict(q=q, params=params)
    zero = dict.fromkeys(rec["params"], 0)
    b0 = exact.eval_q(rec, zero)
    return b0, [[[x-y for x,y in zip(r,s)] for r,s in zip(
        exact.eval_q(rec, dict(zero, **{p:1})), b0)] for p in rec["params"]]


def rank_Mstar(rec, Q):
    """Q-rank of the *symbolic family* module generated by e_j and q_i.

    Rational evaluations always have rank <=3 and cannot express incommensurate
    rank 3+d. Stack each independent formal parameter's coefficients instead;
    verify the supplied exact evaluation has the advertised shape/family.
    """
    d = rec["d"]
    if len(Q) != d or any(len(r) != 3 for r in Q):
        raise ValueError("Q must have d rows and three columns")
    b0, dirs = family(rec)
    if solve_lin(tr([sum(b,[]) for b in dirs]),
                 [x-y for x,y in zip(sum(mat(Q),[]),sum(b0,[]))],len(dirs)) is None:
        raise ValueError("Q is outside the record's affine family")
    rows = [[Fr(int(i==j)) for j in range(3)] + [Fr(0)]*(3*len(dirs)) for i in range(3)]
    rows += [list(b0[i])+[x for b in dirs for x in b[i]] for i in range(d)]
    return len(rref(rows)[1])


def invariant(rec, Q):
    """Cheap exact conjugacy invariants, including dimension and family rank."""
    return _invariant(rec_key(rec), rec["d"], rank_Mstar(rec,Q))


@functools.lru_cache(maxsize=20000)
def _invariant(key, d, rank):
    gens, cen, q = key
    rec = dict(gen=gens, cen=cen, q=q, d=d)
    ops = exact.build_ops(rec)
    chars = []
    seen = set()
    for W,_ in ops:
        if W in seen: continue
        seen.add(W)
        R = [list(r[:3]) for r in W[:3]]
        power = im.eye(3)
        for order in range(1, 49):
            power = mm(power,R)
            if power == im.eye(3): break
        else: raise ValueError("noncrystallographic point operation")
        chars.append((order, im.det(R),sum(R[i][i] for i in range(3)),
                      sum(W[i][i] for i in range(3,3+rec['d']))))
    # Raw nops changes under primitive/centred basis changes; |K| does not.
    return d, len(seen), rank, tuple(sorted(chars))


def lattice_basis(vecs, n):
    """Columns basis (n×n rational) of Z^n + Σ Z v."""
    vs = [[Fr(int(i == j)) for j in range(n)] for i in range(n)] + [[Fr(x) for x in v] for v in vecs]
    D = den([x for v in vs for x in v])
    H = im.hnf([[int(x * D) for x in v] for v in vs])
    assert len(H) == n, "lattice basis must have full rank"
    return tr([[Fr(x, D) for x in r] for r in H])


def mul4(A, B):
    n = len(A)
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)) for i in range(n))


class Prep:
    """Everything about one record that the search needs, computed once."""

    def __init__(self, rec):
        if rec["d"] not in (1,2,3):
            raise ValueError("supported dimensions are d = 1, 2, 3")
        self.d, self.n = rec["d"], 3+rec["d"]
        self.rec = rec
        self.ops = exact.build_ops(rec)
        self.lin = {}
        for W, w in self.ops:
            self.lin.setdefault(W, w)
        self.K = {tuple(tuple(r[:3]) for r in W[:3]): tuple(tuple(r[3:]) for r in W[3:]) for W in self.lin}   # R -> ε
        self.gens = [exact.parse_op(g, self.n) for g in rec["gen"]]
        cen = [[Fr(x) for x in c.split(",")] for c in rec["cen"]]
        self.Lam = lattice_basis(cen, self.n)
        self.Laminv = im.inv_rat(self.Lam)
        self.detLam = im.det(self.Lam)
        self.L3 = lattice_basis([c[:3] for c in cen], 3)
        self.b0, self.Bd = family(rec)
        for W in self.lin:
            R, E = [r[:3] for r in W[:3]], [r[3:] for r in W[3:]]
            constant = [[x-y for x,y in zip(r,t)] for r,t in zip(mm(self.b0,R),mm(E,self.b0))]
            assert constant == [list(r[:3]) for r in W[3:]], "q-family constant does not match operator M"
            for b in self.Bd:
                assert mm(b,R) == mm(E,b), rec.get("label")
        assert len(self.K) == len(self.lin), "point action is not faithful over the external group"
        self._std()

    def _std(self):
        """spglib standardisation: T (x_std = T x), space-group number, standard K and lattice."""
        import spglib
        Rs = list(self.K)
        G0 = np.array([[1.0, 0.13, 0.07], [0.13, 1.31, 0.11], [0.07, 0.11, 1.73]])
        G = sum(np.array(R, float).T @ G0 @ np.array(R, float) for R in Rs)
        L = np.linalg.cholesky(G)                                  # rows a, b, c: L Lᵀ = G
        motif = [np.array([0.1137, 0.2371, 0.3719]), np.array([0.4219, 0.0613, 0.2903]),
                 np.array([0.3391, 0.4787, 0.1531])]          # three species: no accidental symmetry
        pos, num = [], []
        for sp, x0 in enumerate(motif):
            orb = []
            for W, w in self.ops:
                x = (np.array([r[:3] for r in W[:3]], float) @ x0 + np.array([float(t) for t in w[:3]])) % 1.0
                # Float motif deduplication feeds spglib's basic-group standardisation.
                # Candidate conjugation and lattice membership are checked exactly in
                # equivalent(), _solve_m(), and _solve_shift(); check_map() verifies maps.
                if not any(np.allclose((x - y + 0.5) % 1.0 - 0.5, 0, atol=1e-8) for y in orb):
                    orb.append(x)
            pos += orb
            num += [sp + 1] * len(orb)
        ds = spglib.get_symmetry_dataset((L, np.array(pos), num), symprec=1e-5)
        assert ds is not None, "spglib failed"
        self.sg = int(ds.number)
        T = [[Fr(float(x)).limit_denominator(24) for x in r] for r in ds.transformation_matrix]
        Ti = im.inv_rat(T)
        Kstd = set()
        for R in Rs:
            Rp = mm(mm(T, mat(R)), Ti)
            assert all(is_int(x) for r in Rp for x in r), "bad spglib transformation"
            Kstd.add(tuple(tuple(int(x) for x in r) for r in Rp))
        assert len(self.ops) == len(ds.rotations), (len(self.ops), len(ds.rotations))
        self.T, self.Ti = T, Ti
        self.Kstd = frozenset(Kstd)
        Ls = mm(T, self.L3)
        self.Lstd = tuple(tuple(r) for r in lattice_basis([c for c in tr(Ls)], 3))


def rec_key(rec):
    return (tuple(rec["gen"]), tuple(rec["cen"]), tuple(rec["q"]))


_PREP = {}


def prep(rec):
    k = rec_key(rec)
    if k not in _PREP:
        _PREP[k] = Prep(rec)
    return _PREP[k]


# ================================================================ normaliser candidates
@functools.lru_cache(maxsize=1)
def _unimodular():
    """All integer 3×3 matrices with entries in [-NBOUND, NBOUND] and det ±1 (n×3×3 int64),
    ordered by (Σ|N − I|, entries) so the identity comes first."""
    vals = np.arange(-NBOUND, NBOUND + 1)
    rows = np.array(np.meshgrid(vals, vals, vals, indexing="ij")).reshape(3, -1).T      # all rows
    nr = len(rows)
    I0, I1 = np.meshgrid(np.arange(nr), np.arange(nr), indexing="ij")
    I0, I1 = I0.ravel(), I1.ravel()
    cr = np.cross(rows[I0], rows[I1])                                                   # r0 × r1
    out = []
    for r2 in rows:                                     # det = (r0 × r1) · r2
        dd = cr @ r2
        sel = np.abs(dd) == 1
        out.append(np.stack([rows[I0[sel]], rows[I1[sel]], np.broadcast_to(r2, (sel.sum(), 3))], axis=1))
    N = np.concatenate(out)
    cost = np.abs(N - np.eye(3, dtype=int)).reshape(len(N), -1).sum(1)
    order = np.lexsort(tuple(N.reshape(len(N), -1).T[::-1]) + (cost,))
    return N[order]


@functools.lru_cache(maxsize=None)
def normaliser_candidates(Kstd, Lstd):
    """Integer N (bounded, det ±1) with N K N^-1 = K and N L = L, one per coset N·K (m×3×3 array)."""
    N = _unimodular()
    # Float determinants only recover the sign of these integer unimodular matrices;
    # _unimodular() selected det ±1 with integer arithmetic. equivalent() checks det(P)
    # and the superspace lattice volume exactly before accepting a candidate.
    det = np.round(np.linalg.det(N)).astype(np.int64)
    adj = np.stack([np.cross(N[:, 1], N[:, 2]), np.cross(N[:, 2], N[:, 0]), np.cross(N[:, 0], N[:, 1])], axis=2)
    Ninv = adj * det[:, None, None]                    # inverse of a det ±1 matrix
    code = lambda M: ((M + 1) * np.array([3 ** i for i in range(9)]).reshape(3, 3)).sum((-1, -2))
    kcodes = np.array(sorted(code(np.array(R)) for R in Kstd))
    ok = np.ones(len(N), bool)
    for R in Kstd:
        C = N @ np.array(R) @ Ninv
        ok &= (np.abs(C).max((1, 2)) <= 1)
        ok &= np.isin(code(np.clip(C, -1, 1)), kcodes)
    L = np.array([[float(x) for x in r] for r in Lstd])
    Li = np.linalg.inv(L)
    X = Li @ N[ok] @ L
    # Float prefilter rejects N that do not preserve the standard centred lattice.
    # _solve_m() checks superspace lattice integrality exactly, and equivalent()
    # checks its determinant/volume; check_map() independently verifies both.
    good = np.all(np.abs(X - np.round(X)) < 1e-9, axis=(1, 2))
    # One N per coset N·K: if S is an equivalence with linear part P R (R ∈ K_A), then S∘g^-1
    # (g ∈ A with linear part R) is one with linear part P.  Keep the first (simplest) of each.
    cand = N[ok][good]
    # Superspace type equivalence preserves external orientation (2013 eq.7).
    # Filter before coset representatives: an improper K element changes sign.
    cand = cand[np.linalg.det(cand) > 0]
    if len(Kstd) == 1:
        return cand
    prods = np.stack([cand @ np.array(R) for R in Kstd], axis=1)          # m×|K|×3×3
    codes = (prods + 8).reshape(len(cand), len(Kstd), 9) @ (17 ** np.arange(9))     # |entries| <= 6
    codes.sort(axis=1)
    _, first = np.unique(codes, axis=0, return_index=True)
    return cand[np.sort(first)]


# ================================================================ the equivalence test
def _S(P, M, E):
    if isinstance(E, int): M, E = [M], [[E]]  # legacy d1 private helper callers
    d = len(E)
    return [list(r)+[Fr(0)]*d for r in P]+[list(m)+list(e) for m,e in zip(M,E)]


def _Sinv(Pi, M, E):
    if isinstance(E, int): M, E = [M], [[E]]
    Ei = im.inv_rat(E)
    return _S(Pi, [[-x for x in r] for r in mm(mm(Ei,M),Pi)],Ei)


def _det_polynomial_content(basis, d):
    """GCD of exact determinant-polynomial coefficients; zero means singular span.

    A content >1 excludes unimodular points without enumerating the lattice.
    This necessary obstruction is exact and never rejects det ±1.
    """
    coefficients = {}
    for perm in itertools.permutations(range(d)):
        sign = (-1)**sum(perm[i]>perm[j] for i in range(d) for j in range(i+1,d))
        for inds in itertools.product(range(len(basis)),repeat=d):
            value = sign * math.prod(basis[k][i*d+perm[i]] for i,k in enumerate(inds))
            if value:
                key=tuple(sorted(inds))
                coefficients[key]=coefficients.get(key,0)+value
    return math.gcd(*coefficients.values()) if coefficients else 0


def _e_candidates(A,B,P,Pi):
    """Saturated integer lattice of intertwining, family-preserving E.

    Enumerate lattice coefficients in [-2,2]; this is a bounded oracle search,
    like the external normaliser bound, not a decision procedure for arbitrary
    settings. Identity and signed basis vectors are visited before larger sums.
    """
    d = A.d
    constraints = []
    units = [[[Fr(int(i*d+j == k)) for j in range(d)] for i in range(d)] for k in range(d*d)]
    for R, EA in A.K.items():
        Rp = mm(mm(P,mat(R)),Pi)
        key = tuple(tuple(int(x) for x in r) for r in Rp) if all(is_int(x) for r in Rp for x in r) else None
        if key not in B.K: return
        EB = B.K[key]
        values = [[x-y for x,y in zip(sum(mm(U,EA),[]),sum(mm(EB,U),[]))] for U in units]
        constraints += tr(values)
    orth = nullspace([sum(b,[]) for b in B.Bd],d*3)
    for b in A.Bd:
        values = [mv(orth,sum(mm(mm(U,b),Pi),[])) for U in units]
        constraints += tr(values)
    basis = nullspace(constraints,d*d)
    if not basis: return
    # Clearing denominators alone misses integral points: saturate first.
    integer = im.saturate([[int(x*den(v)) for x in v] for v in basis])
    integer = im.hnf(integer)
    if _det_polynomial_content(integer,d) != 1: return
    eye = im.eye(d)
    def allowed(E):
        v=sum(E,[])
        return all(sum(a*b for a,b in zip(row,v)) == 0 for row in constraints)
    yielded = set()
    for E in (eye,[[-x for x in r] for r in eye]):
        if allowed(E):
            yielded.add(tuple(sum(E,[]))); yield E
    # Iterate by support so useful sparse maps precede the large free-family cube.
    k=len(integer)
    for support in range(1,k+1):
        for inds in itertools.combinations(range(k),support):
            for coeffs in itertools.product((1,-1,2,-2),repeat=support):
                v=[sum(c*integer[i][j] for i,c in zip(inds,coeffs)) for j in range(d*d)]
                key=tuple(v)
                if key in yielded: continue
                E=[v[i*d:(i+1)*d] for i in range(d)]
                if abs(im.det(E)) == 1:
                    yielded.add(key); yield E


def _solve_m(A,B,P,Q,E):
    d,n=A.d,A.n
    c=[[x-y for x,y in zip(r,s)] for r,s in zip(mm(B.b0,P),mm(E,Q))]
    X0=mm(mm(B.Laminv,_S(P,c,E)),A.Lam)
    BdP=[mm(b,P) for b in B.Bd]
    k=len(BdP)
    Xl=[mm(mm(B.Laminv,[[Fr(0)]*n for _ in range(3)]+[list(r)+[Fr(0)]*d for r in b]),A.Lam) for b in BdP]
    G=[[Xl[l][i][j] for l in range(k)] for i in range(n) for j in range(n)]
    h=[-x for r in X0 for x in r]
    res=mixed_solve(G,[[] for _ in h],h,k,0)
    if res is None: return None
    y0,_=res
    M0=[[c[i][j]+sum(y0[l]*BdP[l][i][j] for l in range(k)) for j in range(3)] for i in range(d)]
    Y=int_lattice_of_real(G) if k else []
    dirs=[[[sum(y[l]*BdP[l][i][j] for l in range(k)) for j in range(3)] for i in range(d)] for y in Y]
    return M0,dirs


def _solve_shift(A,B,P,Pi,M0,dirs,E):
    n,d=A.n,A.d
    S,Si=_S(P,M0,E),_Sinv(Pi,M0,E)
    G,H,h=[],[],[]
    for W,w in A.gens:
        Wp=mm(mm(S,mat(W)),Si)
        key=tuple(tuple(int(x) for x in r) for r in Wp) if all(is_int(x) for r in Wp for x in r) else None
        if key not in B.lin: return None
        # Family directions intertwine, so changing M along dirs leaves Wp fixed.
        IW=[[Fr(int(i==j))-Wp[i][j] for j in range(n)] for i in range(n)]
        G+=mm(B.Laminv,IW)
        ext=[[Fr(0)]*3+mv(D,w[:3]) for D in dirs]
        H+=tr([mv(B.Laminv,v) for v in ext]) if dirs else [[] for _ in range(n)]
        h+=mv(B.Laminv,[x-y for x,y in zip(B.lin[key],mv(S,w))])
    res=mixed_solve(G,H,h,n,len(dirs))
    if res is None: return None
    shift,t=res
    M=[[M0[i][j]+sum(ti*D[i][j] for ti,D in zip(t,dirs)) for j in range(3)] for i in range(d)]
    return shift,M


def _family_torsion(A):
    """Order of the affine family in matrices/(direction span + integer lattice).

    Only used for product lattices L3 ⊕ Z^d and equally dimensional families.
    Then E and L3B^-1 P L3A are unimodular integer matrices, and M L3A is
    integral, so the induced automorphism preserves this torsion order.
    Narrow parent families in subgroup matching deliberately skip comparison.
    """
    if hasattr(A, '_torsion'): return A._torsion
    if any(any(Fr(x) for x in c.split(',')[3:]) for c in A.rec['cen']):
        A._torsion = None
        return None
    directions=[sum(mm(b,A.L3),[]) for b in A.Bd]
    dimension=len(rref(directions)[1])
    N=nullspace(directions,A.d*3)
    if not N:
        A._torsion=(dimension,1)
        return A._torsion
    scale=den(sum(N,[]))
    Ni=[[int(x*scale) for x in row] for row in N]
    D,U,_,_=im.snf(Ni)
    value=mv(U,mv(Ni,sum(mm(A.b0,A.L3),[])))
    order=den([value[i]/D[i][i] for i in range(len(N))])
    A._torsion=(dimension,order)
    return A._torsion



def _rational_generator_lattice(generators, dimension):
    """Column basis of exactly the supplied full-rank rational lattice.

    Unlike lattice_basis(), no ambient unit vectors are added here.
    """
    scale=den([x for row in generators for x in row])
    H=im.hnf([[int(Fr(x)*scale) for x in row] for row in generators])
    assert len(H)==dimension, "rational generator lattice must have full rank"
    return tr([[Fr(x,scale) for x in row] for row in H])


def _physical_offset_quotient(A):
    """Physical direction rank and finite quotient of projected constant rows.

    For product lattice L3 x Z^d, use primitive Q=Q L3. Let V be the span of
    all parameter-coefficient rows, N a basis of its annihilator, L0=N Z^3,
    and L1=L0+sum_i Z N Q0[i]^T. Return the Smith factors >1 of L1/L0.

    Under Q'=(E Q+integer) U^-1, E,U unimodular, equal physical row ranks
    turn direction inclusion into equality. Annihilators obey N'=T N U^T,
    so both lattices transform by T: integer shifts contribute only L0 and
    unimodular E preserves the generated constant-row subgroup. Unequal
    physical ranks (narrower parent families) and nonproduct lattices skip
    comparison in equivalent(). The quotient is trivial when V spans Q^3.
    """
    if hasattr(A,'_physical_quotient'): return A._physical_quotient
    if any(any(Fr(x) for x in c.split(',')[3:]) for c in A.rec['cen']):
        A._physical_quotient=None
        return None
    rows=[row for b in A.Bd for row in mm(b,A.L3)]
    physical_rank=len(rref(rows)[1])
    N=nullspace(rows,3)
    if not N:
        A._physical_quotient=(physical_rank,())
        return A._physical_quotient
    base=tr(N)
    projected=[mv(N,row) for row in mm(A.b0,A.L3)]
    B0=_rational_generator_lattice(base,len(N))
    B1=_rational_generator_lattice(base+projected,len(N))
    inclusion=mm(im.inv_rat(B1),B0)
    assert all(is_int(x) for row in inclusion for x in row)
    D,_,_,_=im.snf([[int(x) for x in row] for row in inclusion])
    factors=tuple(abs(D[i][i]) for i in range(len(N)) if abs(D[i][i])>1)
    A._physical_quotient=(physical_rank,factors)
    return A._physical_quotient



def _embedded_offset_data(B):
    """Cache B's annihilator and embedded projected-offset superlattice."""
    if hasattr(B,'_embedded_offsets'): return B._embedded_offsets
    if any(any(Fr(x) for x in c.split(',')[3:]) for c in B.rec['cen']):
        B._embedded_offsets=None
        return None
    rows=[row for b in B.Bd for row in mm(b,B.L3)]
    N=nullspace(rows,3)
    if not N:
        B._embedded_offsets=([],[])
    else:
        projected=[mv(N,row) for row in mm(B.b0,B.L3)]
        B._embedded_offsets=(N,_rational_generator_lattice(tr(N)+projected,len(N)))
    return B._embedded_offsets


def _embedded_offset_compatible(A,B,P,Pi):
    """Necessary product-lattice conditions for a fixed external P, before E.

    U=L3B^-1 P L3A must be unimodular integral. Target annihilator N kills
    every transformed source direction (invertible E cannot change this).
    The exact embedded lattices N Z^3 + sum_i Z N Q0[i]^T must agree:
    B-family terms vanish under N, primitive M is integral, and unimodular E
    preserves the projected-row subgroup. No equal-family/rank gate is needed.
    """
    if any(any(Fr(x) for x in c.split(',')[3:]) for c in A.rec['cen']): return True
    data=_embedded_offset_data(B)
    if data is None: return True
    U=mm(mm(im.inv_rat(B.L3),P),A.L3)
    if not all(is_int(x) for row in U for x in row) or abs(im.det(U))!=1: return False
    N,target=data
    if not N: return True
    right=mm(Pi,B.L3)
    for b in A.Bd:
        if any(any(mv(N,row)) for row in mm(b,right)): return False
    projected=[mv(N,row) for row in mm(A.b0,right)]
    return _rational_generator_lattice(tr(N)+projected,len(N)) == target


def equivalent(rec_a,Q_a,rec_b):
    """Map {P,p,E,M,p4} from A at exact Q_a to B, or None (bounded search)."""
    if rec_a['d'] != rec_b['d']: return None
    Q=mat(Q_a)
    if invariant(rec_a,Q) != invariant(rec_b,q_generic(rec_b)): return None
    A,B=prep(rec_a),prep(rec_b)
    if A.sg != B.sg or A.Kstd != B.Kstd or A.Lstd != B.Lstd: return None
    ta,tb=_family_torsion(A),_family_torsion(B)
    if ta is not None and tb is not None and ta[0] == tb[0] and ta[1] != tb[1]: return None
    qa,qb=_physical_offset_quotient(A),_physical_offset_quotient(B)
    if qa is not None and qb is not None and qa[0] == qb[0] and qa[1] != qb[1]: return None
    # Exact self-map avoids constructing the normaliser for corpus validation.
    if rec_key(rec_a) == rec_key(rec_b):
        return dict(P=mat(im.eye(3)),p=[Fr(0)]*3,E=im.eye(A.d),M=[[Fr(0)]*3 for _ in range(A.d)],p4=[Fr(0)]*A.d)
    for N in normaliser_candidates(B.Kstd,B.Lstd):
        P=mm(mm(B.Ti,mat(N.tolist())),A.T)
        if im.det(P) <= 0 or im.det(P)*abs(A.detLam) != abs(B.detLam): continue
        Pi=im.inv_rat(P)
        if not _embedded_offset_compatible(A,B,P,Pi): continue
        for E in _e_candidates(A,B,P,Pi):
            sm=_solve_m(A,B,P,Q,E)
            if sm is None: continue
            sh=_solve_shift(A,B,P,Pi,*sm,E)
            if sh is None: continue
            shift,M=sh
            result=dict(P=P,p=shift[:3],E=E,M=M,p4=shift[3:])
            if check_map(rec_a,Q,rec_b,result): return result
    return None


def check_map(rec_a,Q_a,rec_b,res):
    """Independently verify lattice, affine family and every operation exactly."""
    if rec_a['d'] != rec_b['d']: return False
    A,B=prep(rec_a),prep(rec_b)
    if len(Q_a) != A.d or any(len(row) != 3 for row in Q_a): return False
    Q_a=mat(Q_a)
    source_rhs=[x-y for x,y in zip(sum(Q_a,[]),sum(A.b0,[]))]
    if solve_lin(tr([sum(b,[]) for b in A.Bd]),source_rhs,len(A.Bd)) is None: return False
    # Read legacy fixture maps during the d1 fixture migration.
    E=res.get('E',[[res.get('eta',1)]])
    if not all(is_int(x) for row in E for x in row): return False
    M=mat(res['M'] if 'M' in res else [res['m']])
    p4=res['p4'] if isinstance(res['p4'],list) else [res['p4']]
    P=mat(res['P'])
    if im.det(P) <= 0: return False
    Pi=im.inv_rat(P)
    s=[Fr(x) for x in res['p']+p4]
    S,Si=_S(P,M,E),_Sinv(Pi,M,E)
    X=mm(mm(B.Laminv,S),A.Lam)
    if not all(is_int(x) for r in X for x in r) or abs(im.det(X)) != 1 or abs(im.det(E)) != 1: return False
    qp=mm([[x+y for x,y in zip(r,t)] for r,t in zip(mm(E,mat(Q_a)),M)],Pi)
    rhs=[x-y for x,y in zip(sum(qp,[]),sum(B.b0,[]))]
    if solve_lin(tr([sum(b,[]) for b in B.Bd]),rhs,len(B.Bd)) is None: return False
    # Family directions must transform into B's symbolic parameter span too.
    for b in A.Bd:
        if solve_lin(tr([sum(v,[]) for v in B.Bd]),sum(mm(mm(E,b),Pi),[]),len(B.Bd)) is None: return False
    for W,w in A.ops:
        Wp=mm(mm(S,mat(W)),Si)
        wp=[a+b-c for a,b,c in zip(mv(S,w),s,mv(Wp,s))]
        if not any(mat(WB)==Wp and all(is_int(x) for x in mv(B.Laminv,[a-b for a,b in zip(wp,wB)])) for WB,wB in B.ops): return False
    return len(A.ops)*abs(A.detLam)==len(B.ops)*abs(B.detLam)


# ================================================================ data access / CLI
@functools.lru_cache(maxsize=None)
def load_records(d: int, data_dir=DATA) -> dict[str, dict]:
    """Cached {label: record} for dimension d, in file / label order."""
    out = {}
    for f in sorted(glob.glob(os.path.join(data_dir, "sg*.json"))):
        with open(f) as source:
            groups = json.load(source)["groups"]
        for label, r in groups.items():
            if r["d"] == d:
                out[label] = dict(r, label=label)
    return out


def load_d1(data_dir=DATA): return load_records(1, data_dir)


if __name__ == "__main__":
    db = load_d1()
    a, b = db[sys.argv[1]], db[sys.argv[2]]
    r = equivalent(a, q_generic(a), b)
    print(None if r is None else {k: (exact.fstr(v) if isinstance(v, Fr) else v) for k, v in r.items()})
