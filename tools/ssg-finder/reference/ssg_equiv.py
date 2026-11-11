#!/usr/bin/env python3
"""ssg_equiv.py -- exact equivalence test between (3+1)D superspace groups (offline oracle).

Two groups A, B (records as in CONTRACTS §2, d = 1) are equivalent iff a superspace affine map

    x' = P x + p,    x4' = m·x + η x4 + p4        (S = [[P, 0], [m, η]], s = (p, p4))

sends the superspace lattice Σ_A onto Σ_B, the operations of A onto those of B (mod Σ_B), and
the given q of A into the q family of B *exactly*:  q' = (η q + m) P^-1 = q_B(y) for some y.
(Records satisfy M = qR − εq exactly at generic q, so "mod L*" is carried by m and the lattice.)

Search.  Both basic groups are standardised with spglib (x_std = T x).  P runs over
T_B^-1 N T_A, with N an integer matrix (entries in -NBOUND..NBOUND, det ±1) that normalises the
standard point group and preserves the standard centred lattice.  For each P and η:
  * m:  m = (b0 + y Bd) P − η q, y real, with Λ_B^-1 S Λ_A integral (an affine congruence in y);
        the integer freedom in y is kept as unknowns t ∈ Z^k for the next step.
  * (p, p4, t):  S w + (1 − W') s ≡ w_B(W') (mod Σ_B) for the generators of A.
Both steps are mixed real/integer congruences, solved exactly with the Smith normal form.

Usage:  python3 ssg_equiv.py <label_a> <label_b>
"""
import functools, glob, json, math, os, sys
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
        return [Fr(0)] * ncol
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
    assert s is not None
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
    """q of a d = 1 record at the generic parameter values of CONTRACTS §3.1."""
    assert rec["d"] == 1, "d > 1 not supported"
    params = {p: Fr(exact.GENERIC_SEEDS[i], exact.P_GENERIC) for i, p in enumerate(rec["params"])}
    return exact.eval_q(rec, params)[0]


def lattice_basis(vecs, n):
    """Columns basis (n×n rational) of Z^n + Σ Z v."""
    vs = [[Fr(int(i == j)) for j in range(n)] for i in range(n)] + [[Fr(x) for x in v] for v in vecs]
    D = den([x for v in vs for x in v])
    H = im.hnf([[int(x * D) for x in v] for v in vs])
    assert len(H) == n
    return tr([[Fr(x, D) for x in r] for r in H])


def mul4(A, B):
    n = len(A)
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)) for i in range(n))


class Prep:
    """Everything about one record that the search needs, computed once."""

    def __init__(self, rec):
        assert rec["d"] == 1, "only (3+1)D groups (d = 1) are supported"
        self.rec = rec
        self.ops = exact.build_ops(rec)
        self.lin = {}
        for W, w in self.ops:
            self.lin.setdefault(W, w)
        self.K = {tuple(tuple(r[:3]) for r in W[:3]): W[3][3] for W in self.lin}   # R -> ε
        self.gens = [exact.parse_op(g, 4) for g in rec["gen"]]
        cen = [[Fr(x) for x in c.split(",")] for c in rec["cen"]]
        self.Lam = lattice_basis(cen, 4)
        self.Laminv = im.inv_rat(self.Lam)
        self.detLam = im.det(self.Lam)
        self.L3 = lattice_basis([c[:3] for c in cen], 3)
        # q family: q(y) = b0 + y Bd
        zero = {p: 0 for p in rec["params"]}
        b0 = exact.eval_q(rec, zero)[0]
        self.b0 = b0
        self.Bd = [[x - y for x, y in zip(exact.eval_q(rec, dict(zero, **{p: 1}))[0], b0)] for p in rec["params"]]
        for W in self.lin:   # family directions must satisfy b (R − ε) = 0 for every op
            R = [r[:3] for r in W[:3]]
            for b in self.Bd:
                assert all(x - W[3][3] * y == 0 for x, y in zip(vm(b, R), b)), rec.get("label")
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
    good = np.all(np.abs(X - np.round(X)) < 1e-9, axis=(1, 2))
    # One N per coset N·K: if S is an equivalence with linear part P R (R ∈ K_A), then S∘g^-1
    # (g ∈ A with linear part R) is one with linear part P.  Keep the first (simplest) of each.
    cand = N[ok][good]
    if len(Kstd) == 1:
        return cand
    prods = np.stack([cand @ np.array(R) for R in Kstd], axis=1)          # m×|K|×3×3
    codes = (prods + 8).reshape(len(cand), len(Kstd), 9) @ (17 ** np.arange(9))     # |entries| <= 6
    codes.sort(axis=1)
    _, first = np.unique(codes, axis=0, return_index=True)
    return cand[np.sort(first)]


# ================================================================ the equivalence test
def _S(P, m, eta):
    return [list(P[0]) + [Fr(0)], list(P[1]) + [Fr(0)], list(P[2]) + [Fr(0)], list(m) + [Fr(eta)]]


def _Sinv(Pi, m, eta):
    mi = [-eta * x for x in vm(m, Pi)]
    return _S(Pi, mi, eta)


def _solve_m(A, B, P, q, eta):
    """m = (b0 + y Bd) P − η q with Λ_B^-1 S Λ_A integral.  Returns (m0, dirs) with dirs the
    m-shifts of a basis of the integer freedom in y, or None."""
    c = [x - eta * y for x, y in zip(vm(B.b0, P), q)]
    X0 = mm(mm(B.Laminv, _S(P, c, eta)), A.Lam)
    BdP = [vm(b, P) for b in B.Bd]
    k = len(BdP)
    Xl = []
    for row in BdP:
        Z = [[Fr(0)] * 4 for _ in range(3)] + [list(row) + [Fr(0)]]
        Xl.append(mm(mm(B.Laminv, Z), A.Lam))
    G = [[Xl[l][i][j] for l in range(k)] for i in range(4) for j in range(4)]
    h = [-X0[i][j] for i in range(4) for j in range(4)]
    res = mixed_solve(G, [[] for _ in h], h, k, 0)
    if res is None:
        return None
    y0, _ = res
    m0 = [ci + sum(y0[l] * BdP[l][j] for l in range(k)) for j, ci in enumerate(c)]
    Y = int_lattice_of_real(G)
    dirs = [[sum(yv[l] * BdP[l][j] for l in range(k)) for j in range(3)] for yv in Y]
    return m0, dirs


def _solve_shift(A, B, P, Pi, m0, dirs, eta):
    """(p, p4) and integers t with S w_A + (1 − W') s ≡ w_B(W') (mod Σ_B) for A's generators."""
    S, Si = _S(P, m0, eta), _Sinv(Pi, m0, eta)
    G, H, h = [], [], []
    for W, w in A.gens:
        Wp = mm(mm(S, mat(W)), Si)
        key = tuple(tuple(int(x) for x in r) for r in Wp) if all(is_int(x) for r in Wp for x in r) else None
        if key not in B.lin:
            return None
        wB = B.lin[key]
        Sw = mv(S, w)
        IW = [[Fr(int(i == j)) - Wp[i][j] for j in range(4)] for i in range(4)]
        G += mm(B.Laminv, IW)
        ext = [sum(d[j] * w[j] for j in range(3)) for d in dirs]        # 4th-component shift per t
        H += [[B.Laminv[i][3] * e for e in ext] for i in range(4)]
        h += mv(B.Laminv, [x - y for x, y in zip(wB, Sw)])
    res = mixed_solve(G, H, h, 4, len(dirs))
    if res is None:
        return None
    s, t = res
    m = [m0[j] + sum(ti * d[j] for ti, d in zip(t, dirs)) for j in range(3)]
    return s, m


def equivalent(rec_a, q_a, rec_b):
    """Map {P, p, eta, m, p4} taking A (at q = q_a) onto B, or None."""
    A, B = prep(rec_a), prep(rec_b)
    if A.sg != B.sg or len(A.K) != len(B.K) or A.Kstd != B.Kstd or A.Lstd != B.Lstd:
        return None
    q = [Fr(x) for x in q_a]
    volA, volB = A.detLam, B.detLam
    Kgens = list(A.K.items())
    for N in normaliser_candidates(B.Kstd, B.Lstd):
        P = mm(mm(B.Ti, mat(N.tolist())), A.T)
        if abs(im.det(P)) * volA != abs(volB):
            continue
        Pi = im.inv_rat(P)
        ok = True
        for R, e in Kgens:
            Rp = mm(mm(P, mat(R)), Pi)
            key = tuple(tuple(int(x) for x in r) for r in Rp) if all(is_int(x) for r in Rp for x in r) else None
            if B.K.get(key) != e:
                ok = False
                break
        if not ok:
            continue
        for eta in (1, -1):
            sm = _solve_m(A, B, P, q, eta)
            if sm is None:
                continue
            sh = _solve_shift(A, B, P, Pi, sm[0], sm[1], eta)
            if sh is None:
                continue
            s, m = sh
            return dict(P=P, p=s[:3], eta=eta, m=m, p4=s[3])
    return None


def check_map(rec_a, q_a, rec_b, res):
    """Independent check that `res` maps A (at q_a) onto B; True/False."""
    A, B = prep(rec_a), prep(rec_b)
    P, eta, m = mat(res["P"]), res["eta"], [Fr(x) for x in res["m"]]
    s = [Fr(x) for x in res["p"]] + [Fr(res["p4"])]
    S, Si = _S(P, m, eta), _Sinv(im.inv_rat(P), m, eta)
    X = mm(mm(B.Laminv, S), A.Lam)
    if not all(is_int(x) for r in X for x in r) or abs(im.det(X)) != 1:
        return False
    qp = vm([eta * Fr(x) + y for x, y in zip(q_a, m)], im.inv_rat(P))
    if solve_lin(tr(B.Bd), [x - y for x, y in zip(qp, B.b0)], len(B.Bd)) is None if B.Bd else qp != B.b0:
        return False
    for W, w in A.ops:
        Wp = mm(mm(S, mat(W)), Si)
        wp = [a + b - c for a, b, c in zip(mv(S, w), s, mv(Wp, s))]
        if not any(mat(WB) == Wp and all(is_int(x) for x in mv(B.Laminv, [a - b for a, b in zip(wp, wB)]))
                   for WB, wB in B.ops):
            return False
    return len(A.ops) * A.detLam == len(B.ops) * B.detLam      # |ops| = |K| / det Λ


# ================================================================ data access / CLI
@functools.lru_cache(maxsize=None)
def load_d1(data_dir=DATA):
    """{label: record} for every d = 1 record in local/data."""
    out = {}
    for f in sorted(glob.glob(os.path.join(data_dir, "sg*.json"))):
        for label, r in json.load(open(f))["groups"].items():
            if r["d"] == 1:
                out[label] = dict(r, label=label)
    return out


if __name__ == "__main__":
    db = load_d1()
    a, b = db[sys.argv[1]], db[sys.argv[2]]
    r = equivalent(a, q_generic(a), b)
    print(None if r is None else {k: (exact.fstr(v) if isinstance(v, Fr) else v) for k, v in r.items()})
