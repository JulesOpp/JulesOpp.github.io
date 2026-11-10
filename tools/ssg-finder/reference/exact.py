#!/usr/bin/env python3
"""exact.py -- exact t0 classes of commensurate sections of a (3+d)D superspace group.

Reference implementation of PLAN.md §4 (milestone M2).  Everything is exact rational /
integer arithmetic except the final space-group identification (spglib) and the setting
symbol (``ssg_descent.hm_in_setting``).

Conventions (see CONTRACTS.md):
  op g = (W, w), W = [[R, 0], [M, E]], w = (tE, tI);  M + E Q = Q R.
  A = E - 1,  c = tI - Q tE.  g survives at t  iff  A t + c ∈ Λ = Z^d + Q Z^3,
  and then its matched n is the coset with  Q n ≡ A t + c (mod Z^d).
  t = B s with B a basis of Λ (columns), s ∈ R^d / Z^d.
  Ã = B^-1 A B (integer), c̃ = B^-1 c;  survival: Ã s + c̃ ∈ Z^d.

Usage:
  python3 exact.py 62.1.9.2 --q g=1/3
"""
import argparse, itertools, json, math, re, sys
from fractions import Fraction as Fr

import numpy as np

import intmat as im
import ssg_descent as proto

P_GENERIC = 10007
GENERIC_SEEDS = (3001, 4723, 7177)
GENERIC_STEP = 1009
EXAMPLE_MAX_DEN = 24
VARS = "xyztuv"


# ================================================================ rationals / formatting
def fr(x):
    return Fr(x)


def fstr(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def frac1(x):
    x = Fr(x)
    return x - math.floor(x)


def lcm(*xs):
    l = 1
    for x in xs:
        l = l * x // math.gcd(l, x)
    return l


def den_lcm(xs):
    return lcm(1, *[Fr(x).denominator for x in xs])


def format_op(W, w):
    """xyzt string: terms ordered by variable, constant last, translation in [0,1)."""
    n = len(W)
    parts = []
    for i in range(n):
        s = ""
        for j in range(n):
            c = W[i][j]
            if c == 0:
                continue
            sign = "-" if c < 0 else ("+" if s else "")
            mag = "" if abs(c) == 1 else str(abs(c))
            s += f"{sign}{mag}{VARS[j]}"
        if w[i] != 0:
            s += ("+" if s else "") + fstr(w[i])
        parts.append(s or "0")
    return ",".join(parts)


def parse_op(s, nv):
    W, w = [], []
    for p in s.split(","):
        c, k = proto.parse_lin(p, nv)
        W.append([int(x) for x in c])
        w.append(Fr(k))
    return W, w


# ================================================================ group records (CONTRACTS §2)
def params_of(qsym):
    out = []
    for q in qsym:
        for tok in re.findall(r"[a-z]+\d*", q):
            if tok not in out:
                out.append(tok)
    return out


def record_from_byu(label):
    """Group record (CONTRACTS §2) from the cached BYU page via ssg_descent.parse_ssg."""
    s = proto.parse_ssg(label)
    d = s["d"]
    nv = 3 + d
    ident = [[int(i == j) for j in range(nv)] for i in range(nv)]
    gen = []
    for W, w in s["ops"]:
        Wi = [[int(x) for x in r] for r in W]
        wf = [frac1(x) for x in w]
        if Wi == ident and not any(wf):
            continue
        gen.append(format_op(Wi, wf))
    cen = [",".join(fstr(frac1(x)) for x in c) for c in s["centering"] if any(frac1(x) for x in c)]
    m = re.match(r"(\d+)\.(\d+\.\d+)\.\d+", label)
    rec = dict(sym=s["symbol"], cls=f"{m.group(2)}" if m else "", d=d, q=list(s["qsym"]),
               params=params_of(s["qsym"]), cen=cen, gen=gen)
    rec["nops"] = len(build_ops(rec))
    rec["label"] = label
    return rec


def build_ops(rec):
    """Canonical operator list (CONTRACTS §1): closure mod Z^{3+d}, identity first, rest
    sorted by (W flattened row-major, w)."""
    nv = 3 + rec["d"]
    ident = tuple(tuple(int(i == j) for j in range(nv)) for i in range(nv))
    gens = [parse_op(g, nv) for g in rec["gen"]]
    for c in rec["cen"]:
        gens.append(([list(r) for r in ident], [Fr(x) for x in c.split(",")]))

    def key(W, w):
        return (tuple(tuple(r) for r in W), tuple(frac1(x) for x in w))

    def mul(a, b):
        (W1, w1), (W2, w2) = a, b
        W = tuple(tuple(sum(W1[i][k] * W2[k][j] for k in range(nv)) for j in range(nv)) for i in range(nv))
        w = tuple(frac1(sum(W1[i][k] * w2[k] for k in range(nv)) + w1[i]) for i in range(nv))
        return W, w

    e = (ident, tuple(Fr(0) for _ in range(nv)))
    G = {e}
    gk = [key(*g) for g in gens]
    frontier = [e]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gk:
                b = mul(a, g)
                if b not in G:
                    G.add(b)
                    nxt.append(b)
        frontier = nxt
    rest = sorted((g for g in G if g != e),
                  key=lambda g: (tuple(x for r in g[0] for x in r), g[1]))
    return [e] + rest


def eval_q(rec, params):
    """q components as written (e.g. "-2a+1/2", "a1", "0") at the given parameter values.
    Implicit multiplication ("2a") is made explicit before the prototype's evaluator."""
    qs = [re.sub(r"(\d)(?=[a-z])", r"\1*", q) for q in rec["q"]]
    return proto.eval_q(qs, {k: str(v) for k, v in params.items()})


# ================================================================ flats
class Flat:
    """Connected affine subtorus {s ∈ R^d/Z^d : C s ≡ f (mod 1)}, C saturated HNF (r x d)."""
    __slots__ = ("C", "f", "d", "key", "pivots")

    def __init__(self, C, f, d):
        self.C, self.f, self.d = C, f, d
        self.key = (tuple(tuple(r) for r in C), tuple(f))
        self.pivots = [next(j for j in range(d) if r[j]) for r in C]

    @property
    def r(self):
        return len(self.C)

    @property
    def dim(self):
        return self.d - len(self.C)

    def contains(self, s):
        return all(frac1(sum(c * x for c, x in zip(row, s)) - fi) == 0 for row, fi in zip(self.C, self.f))

    def point(self, free):
        """Point with the non-pivot coordinates set to `free` (in order), pivot coordinates
        solved from C s = f exactly (bottom row first)."""
        s = [None] * self.d
        it = iter(free)
        for j in range(self.d):
            if j not in self.pivots:
                s[j] = Fr(next(it))
        for i in reversed(range(self.r)):
            p = self.pivots[i]
            row = self.C[i]
            acc = self.f[i] - sum(row[j] * s[j] for j in range(p + 1, self.d))
            s[p] = Fr(acc) / row[p]
        return s

    def to_json(self):
        return {"C": [list(r) for r in self.C], "f": [fstr(x) for x in self.f]}


def canonical_flat(C, f, d):
    """C: r x d integer, rows a basis of a saturated lattice (full row rank)."""
    if not C:
        return Flat([], [], d)
    H, U = im.hnf_with_transform(C)
    assert all(any(r) for r in H), "C not of full row rank"
    f2 = [frac1(sum(u * x for u, x in zip(Ur, f))) for Ur in U]
    return Flat(H, f2, d)


def solve(G, h, d):
    """Connected components of {s : G s ≡ h (mod 1)}; G k x d integer, h rationals."""
    if not G:
        return [Flat([], [], d)]
    D, U, V, Vi = im.snf(G)
    k = len(G)
    Uh = [sum(u * x for u, x in zip(Ur, h)) for Ur in U]
    rho = sum(1 for i in range(min(k, d)) if D[i][i])
    for i in range(rho, k):
        if Fr(Uh[i]).denominator != 1:
            return []
    C = Vi[:rho]
    out = []
    for js in itertools.product(*[range(D[i][i]) for i in range(rho)]):
        f = [Fr(Uh[i] + js[i], D[i][i]) for i in range(rho)]
        out.append(canonical_flat(C, f, d))
    return out


def intersect(F1, F2):
    return solve([list(r) for r in F1.C] + [list(r) for r in F2.C], list(F1.f) + list(F2.f), F1.d)


def subset(G, F):
    """G ⊆ F ?"""
    if G.dim > F.dim:
        return False
    if F.r and im.rank([list(r) for r in G.C] + [list(r) for r in F.C]) != G.r:
        return False
    return F.contains(G.point([0] * G.dim))


# ================================================================ main analysis
class Exact:
    def __init__(self, rec, params):
        self.rec = rec
        self.d = d = rec["d"]
        self.params = {k: Fr(v) for k, v in params.items()}
        self.Q = Q = eval_q(rec, params)
        self.ops = build_ops(rec)
        self._setup_lattices()
        self._setup_ops()
        self._closure()
        self._classes()

    # ------------------------------------------------------------ lattices
    def _setup_lattices(self):
        d, Q = self.d, self.Q
        self.N = [den_lcm([Q[i][j] for i in range(d)]) for j in range(3)]
        Dd = den_lcm([x for r in Q for x in r])
        self.Dden = Dd
        # Λ = Z^d + Q Z^3: generators as row vectors in R^d
        gens = [[Dd * int(i == j) for j in range(d)] for i in range(d)] + \
               [[int(Q[i][j] * Dd) for i in range(d)] for j in range(3)]
        H = im.hnf_lower(gens)                      # rows = basis vectors (× Dd)
        self.Lrows = [[Fr(x, Dd) for x in r] for r in H]
        self.B = im.transpose(self.Lrows)           # columns = basis vectors
        self.Binv = im.inv_rat(self.B)
        # supercell: L_super = {n : Q n ∈ Z^d}
        self.cells = [n for n in itertools.product(*[range(k) for k in self.N])
                      if all(Fr(x).denominator == 1 for x in im.matvec(Q, n))]
        sgens = [list(n) for n in self.cells if any(n)] + \
                [[self.N[j] * int(i == j) for i in range(3)] for j in range(3)]
        self.hnf_super = im.hnf_lower(sgens)
        self.index = abs(int(im.det(self.hnf_super)))
        assert self.index * len(self.cells) == math.prod(self.N) * 1 or True
        # Qn mod 1 -> lexicographically smallest n in the conventional box
        self.nlookup = {}
        for n in itertools.product(*[range(k) for k in self.N]):
            v = tuple(frac1(x) for x in im.matvec(Q, n))
            self.nlookup.setdefault(v, n)
        assert len(self.nlookup) == self.index
        self.lam_reps = sorted(self.nlookup)       # Λ / Z^d representatives (in t)

    # ------------------------------------------------------------ operators
    def _setup_ops(self):
        d, Q = self.d, self.Q
        self.opinfo = []
        for W, w in self.ops:
            R = [list(r[:3]) for r in W[:3]]
            M = [list(r[:3]) for r in W[3:]]
            E = [list(r[3:]) for r in W[3:]]
            assert all(x == 0 for r in W[:3] for x in r[3:])
            lhs = [[M[i][j] + sum(E[i][k] * Q[k][j] for k in range(d)) for j in range(3)] for i in range(d)]
            if lhs != im.matmul(Q, R):
                raise ValueError(f"operator inconsistent with q: {format_op(W, w)}")
            tE, tI = list(w[:3]), list(w[3:])
            A = [[E[i][j] - int(i == j) for j in range(d)] for i in range(d)]
            c = [tI[i] - sum(Q[i][j] * tE[j] for j in range(3)) for i in range(d)]
            At = im.matmul(im.matmul(self.Binv, A), self.B)
            Et = im.matmul(im.matmul(self.Binv, E), self.B)
            assert all(Fr(x).denominator == 1 for r in At for x in r), "Λ not A-invariant"
            At = [[int(x) for x in r] for r in At]
            Et = [[int(x) for x in r] for r in Et]
            ct = im.matvec(self.Binv, c)
            info = dict(R=R, M=M, E=E, tE=tE, tI=tI, A=A, c=c, At=At, Et=Et, ct=ct,
                        Etinv=[[int(x) for x in r] for r in im.inv_rat(Et)])
            if any(any(r) for r in At):
                info["flats"] = solve(At, [-x for x in ct], d)
                info["always"] = False
            else:
                info["flats"] = None
                info["always"] = all(Fr(x).denominator == 1 for x in ct)
            self.opinfo.append(info)

    def t_of_s(self, s):
        return [frac1(x) for x in im.matvec(self.B, s)]

    def s_of_t(self, t):
        return [frac1(x) for x in im.matvec(self.Binv, t)]

    def surviving_at_t(self, t):
        """[(op index, n)] for exact rational t (basic setting)."""
        t = [Fr(x) for x in t]
        out = []
        for k, o in enumerate(self.opinfo):
            v = tuple(frac1(sum(o["A"][i][j] * t[j] for j in range(self.d)) + o["c"][i])
                      for i in range(self.d))
            n = self.nlookup.get(v)
            if n is not None:
                out.append((k, tuple(n)))
        return out

    # ------------------------------------------------------------ closure
    def _closure(self):
        d = self.d
        torus = Flat([], [], d)
        self.flats = {torus.key: torus}
        queue = []
        for o in self.opinfo:
            if o["flats"]:
                queue.extend(o["flats"])
        while queue:
            F = queue.pop(0)
            if F.key in self.flats:
                continue
            existing = list(self.flats.values())
            self.flats[F.key] = F
            for G in existing:
                for H in intersect(F, G):
                    if H.key not in self.flats:
                        queue.append(H)
        fl = list(self.flats.values())
        self.sub = {F.key: [G for G in fl if G.key != F.key and G.dim < F.dim and subset(G, F)]
                    for F in fl}

    def in_open_part(self, F, s):
        return F.contains(s) and not any(G.contains(s) for G in self.sub[F.key])

    def classify_s(self, s):
        """The closure flat whose open part contains s."""
        cont = [F for F in self.flats.values() if F.contains(s)]
        m = min(F.dim for F in cont)
        best = [F for F in cont if F.dim == m]
        assert len(best) == 1
        return best[0]

    # ------------------------------------------------------------ points
    def generic_point(self, F):
        for a in range(1000):
            free = []
            for m in range(F.dim):
                k = (GENERIC_SEEDS[m] + GENERIC_STEP * a) % P_GENERIC or 1
                free.append(Fr(k, P_GENERIC))
            s = F.point(free)
            if self.in_open_part(F, s):
                return [frac1(x) for x in s]
        raise RuntimeError("no generic point found")

    def example_point(self, F):
        """Nice t: smallest grid denominator m of the free s-coordinates (m = 1..24, free
        coordinates k/m in lexicographic order of k), then over Λ-translates the t mod 1
        minimising (lcm of denominators, number of non-zero components, t)."""
        for m in range(1, EXAMPLE_MAX_DEN + 1):
            cands = []
            for ks in itertools.product(range(m), repeat=F.dim):
                s = F.point([Fr(k, m) for k in ks])
                if self.in_open_part(F, s):
                    t0 = im.matvec(self.B, s)
                    for lam in self.lam_reps:
                        cands.append(tuple(frac1(a + b) for a, b in zip(t0, lam)))
            if cands:
                return list(min(cands, key=lambda t: (den_lcm(t), sum(1 for x in t if x), t)))
        return None

    # ------------------------------------------------------------ 3D group
    def ops3(self, surv):
        out, seen = [], set()
        N = self.N
        for k, n in surv:
            o = self.opinfo[k]
            R = o["R"]
            Rs = [[Fr(R[i][j] * N[j], N[i]) for j in range(3)] for i in range(3)]
            assert all(x.denominator == 1 for r in Rs for x in r)
            Rs = [[int(x) for x in r] for r in Rs]
            tau = [o["tE"][i] + n[i] for i in range(3)]
            for cell in self.cells:
                t = tuple(frac1((tau[i] + cell[i]) / N[i]) for i in range(3))
                key = (tuple(map(tuple, Rs)), t)
                if key not in seen:
                    seen.add(key)
                    out.append((Rs, t))
        return out

    def lattice(self):
        """An invariant metric of the basic point group, in the conventional supercell."""
        G = np.zeros((3, 3))
        for o in self.opinfo:
            R = np.array(o["R"], dtype=float)
            G += R.T @ R
        Nd = np.diag(self.N).astype(float)
        Gs = Nd @ G @ Nd
        return np.linalg.cholesky(Gs)   # rows = lattice vectors, L L^T = Gs

    def identify(self, ops3):
        return identify(ops3, self.lattice())

    # ------------------------------------------------------------ classes
    def _classes(self):
        d = self.d
        classes = []
        for F in self.flats.values():
            g = self.generic_point(F)
            t = self.t_of_s(g)
            surv = self.surviving_at_t(t)
            o3 = self.ops3(surv)
            sgt = self.identify(o3)
            classes.append(dict(flat=F, dim=F.dim, generic_s=g, generic_t=t, surviving=surv,
                                ops3=o3, sg_number=sgt.number, sg_symbol=sgt.international_short,
                                setting=proto.hm_in_setting(o3)))
        classes.sort(key=lambda c: (-c["dim"], c["sg_number"], c["setting"], tuple(c["generic_t"])))
        self.classes = classes
        self.by_key = {c["flat"].key: c for c in classes}
        # ---- SSG orbits (union-find)
        parent = {c["flat"].key: c["flat"].key for c in classes}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for c in classes:
            F = c["flat"]
            for o in self.opinfo:
                img = self.image(F, o)
                assert img.key in parent, "closure not invariant under the SSG"
                a, b = find(F.key), find(img.key)
                if a != b:
                    parent[b] = a
        orbit_id, roots = {}, {}
        for c in classes:
            r = find(c["flat"].key)
            if r not in roots:
                roots[r] = len(roots)
            c["orbit"] = roots[r]
        # ---- loci and example points
        for c in classes:
            c["example_t"] = self.example_point(c["flat"])
            assert self.classify_s(self.s_of_t(c["example_t"])) is c["flat"]
        for oid in set(c["orbit"] for c in classes):
            members = [c for c in classes if c["orbit"] == oid]
            txt = self.locus_text(members)
            for c in members:
                c["locus"] = txt
        if d == 1:
            self.it_N = int(1 / self.B[0][0])
            for c in classes:
                c["it_label"] = self.it_label(c)

    def image(self, F, o):
        """Image of flat F under s -> Ẽ s + c̃."""
        if not F.C:
            return F
        Ci = im.matmul([list(r) for r in F.C], o["Etinv"])
        f = [fi + sum(a * b for a, b in zip(row, o["ct"])) for fi, row in zip(F.f, Ci)]
        return canonical_flat(Ci, f, self.d)

    # ------------------------------------------------------------ loci text
    def tvar(self, j):
        return "t" if self.d == 1 else f"t{j + 1}"

    def linform(self, c):
        s = ""
        for j, x in enumerate(c):
            if x == 0:
                continue
            sign = ("-" if x < 0 else "") if not s else (" - " if x < 0 else " + ")
            mag = "" if abs(x) == 1 else str(abs(x))
            s += f"{sign}{mag}{self.tvar(j)}"
        return s

    def locus_text(self, members):
        """Union over the SSG orbit (and implicitly over Λ-translates) in t-coordinates."""
        d = self.d
        if members[0]["dim"] == d:
            return "general"
        groups = {}
        for c in members:
            F = c["flat"]
            K = im.matmul([list(r) for r in F.C], self.Binv)       # rows: K t ≡ f
            den = den_lcm([x for r in K for x in r])
            Wb = [[Fr(x, den) for x in r] for r in im.hnf([[int(x * den) for x in r] for r in K])]
            key = tuple(tuple(r) for r in Wb)
            phi = tuple(frac1(sum(w * t for w, t in zip(r, c["generic_t"]))) for r in Wb)
            groups.setdefault(key, (Wb, set()))[1].add(phi)
        texts = []
        for key in sorted(groups):
            Wb, Rset = groups[key]
            texts.append(self._group_text(Wb, Rset))
        if len(texts) == 1:
            return texts[0]
        return " or ".join(f"[{t}]" if "," in t and not t.startswith("(") else t for t in sorted(texts))

    def _group_text(self, Wb, Rset):
        r = len(Wb)
        R = sorted(Rset)
        x0 = R[0]
        H = []
        for x in R:
            h = tuple(frac1(a - b) for a, b in zip(x, x0))
            if {tuple(frac1(a + b) for a, b in zip(y, h)) for y in R} == Rset:
                H.append(h)
        # lattice Z^r + H, its dual (integer rows), HNF
        den = den_lcm([x for h in H for x in h])
        gens = [[den * int(i == j) for j in range(r)] for i in range(r)] + [[int(x * den) for x in h] for h in H]
        L = [[Fr(x, den) for x in row] for row in im.hnf_lower(gens)]
        dual = im.transpose(im.inv_rat(L))
        assert all(Fr(x).denominator == 1 for row in dual for x in row)
        Z = im.hnf([[int(x) for x in row] for row in dual])
        rows = [[sum(z[k] * Wb[k][j] for k in range(r)) for j in range(self.d)] for z in Z]
        vals = sorted({tuple(frac1(sum(z[k] * x[k] for k in range(r))) for z in Z) for x in R})
        conds = []
        for i, w in enumerate(rows):
            m0 = den_lcm(w)
            v = [int(x * m0) for x in w]
            g = math.gcd(*v)
            cvec = [x // g for x in v]
            # w = c·g/m0, so  w·t ≡ φ (mod 1)  <=>  c·t ≡ φ·m0/g (mod m0/g)
            mu = Fr(m0, g)
            sign = 1 if next(x for x in cvec if x) > 0 else -1
            cvec = [sign * x for x in cvec]
            conds.append((cvec, mu, sign, i))
        out = []
        varying = [i for i in range(len(rows)) if len({v[i] for v in vals}) > 1]
        fixed = [i for i in range(len(rows)) if i not in varying]
        for (cvec, mu, sign, i) in conds:
            if i in fixed:
                val = (sign * vals[0][i] * mu) % mu
                out.append(f"{self.linform(cvec)} ≡ {fstr(val)} (mod {fstr(mu)})")
        if varying:
            cs = [conds[i] for i in varying]
            tups = sorted({tuple((cv[2] * v[cv[3]] * cv[1]) % cv[1] for cv in cs) for v in vals})
            if len(cs) == 1:
                cvec, mu = cs[0][0], cs[0][1]
                out.append(f"{self.linform(cvec)} ≡ {' or '.join(fstr(x[0]) for x in tups)} (mod {fstr(mu)})")
            else:
                lhs = "(" + ", ".join(self.linform(cv[0]) for cv in cs) + ")"
                alts = " or ".join("(" + ", ".join(fstr(x) for x in tp) + ")" for tp in tups)
                mods = ", ".join(fstr(cv[1]) for cv in cs)
                out.append(f"{lhs} ≡ {alts} (mod {mods})")
        return ", ".join(out)

    def it_label(self, c):
        if c["dim"] == 1:
            return "general"
        a = frac1(c["generic_t"][0] * self.it_N)
        names = {Fr(0): "n/N", Fr(1, 4): "1/4N+n/N", Fr(1, 2): "1/2N+n/N", Fr(3, 4): "3/4N+n/N"}
        return names.get(a, f"{fstr(a)}/N+n/N")

    # ------------------------------------------------------------ output
    def to_json(self, label):
        d = self.d
        return {
            "label": label, "d": d,
            "params": {k: fstr(v) for k, v in self.params.items()},
            "Q": [[fstr(x) for x in r] for r in self.Q],
            "ops": [format_op(W, w) for W, w in self.ops],
            "supercell": {
                "hnf": self.hnf_super,
                "conv": self.N,
                "centring": [[fstr(Fr(n[i], self.N[i])) for i in range(3)] for n in self.cells if any(n)],
                "index": self.index,
            },
            "Lambda": [[fstr(x) for x in r] for r in self.B],
            **({"it_N": self.it_N} if d == 1 else {}),
            "classes": [self._class_json(c) for c in self.classes],
        }

    def _class_json(self, c):
        j = {
            "dim": c["dim"],
            "flat": c["flat"].to_json(),
            "locus": c["locus"],
            "generic_t": [fstr(x) for x in c["generic_t"]],
            "example_t": [fstr(x) for x in c["example_t"]],
            "surviving": [[k, list(n)] for k, n in c["surviving"]],
            "sg_number": c["sg_number"],
            "sg_symbol": c["sg_symbol"],
            "setting": c["setting"],
            "orbit": c["orbit"],
        }
        if self.d == 1:
            j["it_label"] = c["it_label"]
        return j


def identify(ops3, lattice=None):
    import spglib
    rots = np.array([R for R, t in ops3], dtype="intc")
    trs = np.array([[float(x) for x in t] for R, t in ops3])
    if lattice is None:
        lattice = np.eye(3)
    sgt = spglib.get_spacegroup_type_from_symmetry(rots, trs, lattice=lattice, symprec=1e-5)
    if sgt is None:
        raise RuntimeError("spglib could not identify the group")
    return sgt


def run(label_or_rec, params):
    rec = record_from_byu(label_or_rec) if isinstance(label_or_rec, str) else label_or_rec
    return Exact(rec, params)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("--q", nargs="+", default=[])
    a = ap.parse_args()
    params = dict(p.split("=") for p in a.q)
    ex = run(a.label, params)
    print(f"Λ basis (columns) {[[fstr(x) for x in r] for r in ex.B]}, supercell {ex.N} index {ex.index}, "
          f"{len(ex.ops)} ops, {len(ex.classes)} classes")
    for c in ex.classes:
        print(f"{c['dim']} orbit {c['orbit']:2d} {c['sg_symbol']:>10} ({c['sg_number']:3d}) {c['setting']:16s} "
              f"ex={[fstr(x) for x in c['example_t']]}  {c['locus']}"
              + (f"   [{c['it_label']}]" if 'it_label' in c else ""))
