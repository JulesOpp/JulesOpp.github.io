"""Necessary simultaneous affine feasibility, before the bounded internal search.

Only product lattices L3 x Z^d are handled. Other presentations use the original
oracle unchanged. The affine integer E lattice is never saturated: its finite
congruences are essential. This module proposes no maps; check_map remains the
independent authority for every accepted result.
"""
import itertools
import math

import ssg_equiv as se


def _det(v, d):
    if d == 1:
        return v[0]
    if d == 2:
        return v[0]*v[3]-v[1]*v[2]
    return (v[0]*(v[4]*v[8]-v[5]*v[7])
            -v[1]*(v[3]*v[8]-v[5]*v[6])
            +v[2]*(v[3]*v[7]-v[4]*v[6]))


class AffineELattice:
    """An exact coset offset + row_Z(basis), or the empty set."""

    def __init__(self, offset, basis, d):
        self.d = d
        self.empty = offset is None
        self.offset = offset
        self.basis = se.im.hnf(basis) if basis else []
        self.pivots = [(next(i for i,x in enumerate(row) if x), row)
                       for row in self.basis]

    def contains(self, E):
        if self.empty:
            return False
        v = [x-y for x,y in zip(sum(E, []), self.offset)]
        for pivot, row in self.pivots:
            q, r = divmod(v[pivot], row[pivot])
            if r:
                return False
            if q:
                v = [x-q*y for x,y in zip(v, row)]
        return not any(v)

    def has_unimodular_residue(self, prime):
        """Whether this coset has any nonsingular matrix over F_prime.

        Nonzero determinant is necessary for det +/-1. Enumeration is finite:
        at most prime^(d*d) <= 19683 points for the primes 2 and 3 used here.
        """
        if self.empty:
            return False
        offset = [x % prime for x in self.offset]
        if _det(offset, self.d) % prime:
            return True
        basis = []
        pivots = []
        for row in self.basis:
            v = [x % prime for x in row]
            for pivot, b in zip(pivots, basis):
                if v[pivot]:
                    c = v[pivot]
                    v = [(x-c*y) % prime for x,y in zip(v,b)]
            if any(v):
                pivot = next(i for i,x in enumerate(v) if x)
                inverse = pow(v[pivot], -1, prime)
                v = [(x*inverse) % prime for x in v]
                pivots.append(pivot)
                basis.append(v)
        # HNF pivots need not stay ordered after reduction modulo prime.
        # Elimination by the insertion-order pivots still gives an independent
        # spanning basis: each new row vanishes at every earlier pivot.
        def visit(i, v):
            if i == len(basis):
                return bool(_det(v, self.d) % prime)
            for _ in range(prime):
                if visit(i+1, v):
                    return True
                v = [(x+y) % prime for x,y in zip(v,basis[i])]
            return False
        return visit(0, offset)


def integer_projection(H, h, nproject, d):
    """Project all integer solutions H*x=h onto the first nproject entries.

    Rational row operations and denominator clearing preserve the equations.
    Smith free columns generate the full integer kernel; projecting and taking
    HNF preserves the resulting (possibly nonsaturated) congruence lattice.
    """
    if not H:
        return AffineELattice([0]*nproject, se.im.eye(nproject), d)
    rows, rhs = [], []
    seen = set()
    for row, value in zip(H,h):
        scale = se.den(list(row)+[value])
        v = [int(x*scale) for x in row]+[int(value*scale)]
        common = math.gcd(*v)
        if not common:
            continue
        v = [x//common for x in v]
        if not any(v[:-1]):
            return AffineELattice(None, [], d)
        key = tuple(v)
        if key not in seen:
            seen.add(key)
            rows.append(v[:-1])
            rhs.append(v[-1])
    if not rows:
        return AffineELattice([0]*nproject, se.im.eye(nproject), d)
    D,U,V,_ = se.im.snf(rows)
    ur = se.im.matvec(U, rhs)
    nvars = len(rows[0])
    rank = sum(bool(D[i][i]) for i in range(min(len(rows),nvars)))
    w = [0]*nvars
    for i in range(len(rows)):
        diagonal = D[i][i] if i < rank else 0
        if diagonal:
            q,r = divmod(ur[i], diagonal)
            if r:
                return AffineELattice(None, [], d)
            w[i] = q
        elif ur[i]:
            return AffineELattice(None, [], d)
    offset = se.im.matvec(V[:nproject], w)
    basis = [[V[i][j] for i in range(nproject)] for j in range(rank,nvars)]
    return AffineELattice(offset, basis, d)


def _point_generators(A):
    """Smallest generating set of the finite faithful external point group."""
    if hasattr(A, '_affine_point_generators'):
        return A._affine_point_generators
    rotations = list(A.K)
    index = {R:i for i,R in enumerate(rotations)}
    identity = index[tuple(map(tuple,se.im.eye(3)))]
    table = [[index[se.mul4(R,T)] for T in rotations] for R in rotations]
    def closure(gens):
        seen = {identity}
        pending = [identity]
        while pending:
            a = pending.pop()
            for g in gens:
                b = table[a][g]
                if b not in seen:
                    seen.add(b)
                    pending.append(b)
        return seen
    greedy, generated = [], {identity}
    for i in range(len(rotations)):
        if i not in generated:
            greedy.append(i)
            generated = closure(greedy)
    chosen = greedy
    pool = [i for i in range(len(rotations)) if i != identity]
    for count in range(1,len(greedy)):
        smaller = next((list(gs) for gs in itertools.combinations(pool,count)
                        if len(closure(gs)) == len(rotations)), None)
        if smaller is not None:
            chosen = smaller
            break
    lifts = {tuple(tuple(r[:3]) for r in W[:3]): (W,w) for W,w in A.lin.items()}
    A._affine_point_generators = [lifts[rotations[i]] for i in chosen]
    return A._affine_point_generators


def feasible_e_lattice(A, B, P, Pi):
    """All E satisfying the simultaneous necessary product-lattice equations.

    Return None to skip this sieve for nonproduct lattices. Integer variables
    are vec(E), vec(J=M L3A), and one target-lattice lift per point generator.
    Real variables are the target family parameters and the common origin.
    Source b0 is used because transformed source directions can be absorbed in
    target parameters; their homogeneous inclusion constraints are also imposed.
    """
    for obj in (A,B):
        if any(any(se.Fr(x) for x in c.split(',')[3:]) for c in obj.rec['cen']):
            return None
    d,n = A.d,A.n
    U = se.mm(se.mm(se.im.inv_rat(B.L3),P),A.L3)
    if not all(se.is_int(x) for row in U for x in row) or abs(se.im.det(U)) != 1:
        return AffineELattice(None, [], d)
    generators = _point_generators(A)
    target = {tuple(tuple(r[:3]) for r in W[:3]): (W,w) for W,w in B.lin.items()}
    pairs = []
    for WA,wa in generators:
        Rp = se.mm(se.mm(P,[r[:3] for r in WA[:3]]),Pi)
        if not all(se.is_int(x) for row in Rp for x in row):
            return AffineELattice(None, [], d)
        key = tuple(tuple(int(x) for x in row) for row in Rp)
        if key not in target:
            return AffineELattice(None, [], d)
        WB,wb = target[key]
        pairs.append((WA,wa,WB,wb))
    ne,nj,ny = d*d,d*3,len(B.Bd)
    ni,nr = ne+nj+n*len(pairs),ny+n
    Li = se.im.inv_rat(A.L3)
    BP = se.mm(B.b0,P)
    directions = [se.mm(b,P) for b in B.Bd]
    equations = []
    def add(h, g, rhs=0):
        if any(h) or any(g) or rhs:
            equations.append(g+h+[rhs])
    # (B0 + sum y_k Bk) P = E A0 + J L3A^-1.
    for i in range(d):
        for j in range(3):
            h,g = [0]*ni,[0]*nr
            for k in range(d):
                h[i*d+k] = -A.b0[k][j]
            for k in range(3):
                h[ne+i*3+k] = -Li[k][j]
            for k in range(ny):
                g[k] = directions[k][i][j]
            add(h,g,-BP[i][j])
    # E A_k P^-1 belongs to the target direction span, for every source k.
    orth = se.nullspace([sum(b,[]) for b in B.Bd],d*3)
    for b in A.Bd:
        bPi = se.mm(b,Pi)
        for row in orth:
            h = [0]*ni
            for i in range(d):
                for k in range(d):
                    h[i*d+k] = sum(row[i*3+j]*bPi[k][j] for j in range(3))
            add(h,[0]*nr)
    for number,(WA,wa,WB,wb) in enumerate(pairs):
        RA,EA,MA = [r[:3] for r in WA[:3]],[r[3:] for r in WA[3:]],[r[:3] for r in WA[3:]]
        EB,MB = [r[3:] for r in WB[3:]],[r[:3] for r in WB[3:]]
        LiRA,MBP = se.mm(Li,RA),se.mm(MB,P)
        # E EA = EB E and M RA + E MA = MB P + EB M.
        for i in range(d):
            for j in range(d):
                h = [0]*ni
                for k in range(d):
                    h[i*d+k] += EA[k][j]
                    h[k*d+j] -= EB[i][k]
                add(h,[0]*nr)
            for j in range(3):
                h = [0]*ni
                for k in range(d):
                    h[i*d+k] += MA[k][j]
                for k in range(3):
                    h[ne+i*3+k] += LiRA[k][j]
                    for ell in range(d):
                        h[ne+ell*3+k] -= EB[i][ell]*Li[k][j]
                add(h,[0]*nr,MBP[i][j])
        # L_B^-1 (S wA + (I-WB) s - wB) is an integer vector.
        IW = [[int(i==j)-WB[i][j] for j in range(n)] for i in range(n)]
        shift = se.mm(B.Laminv,IW)
        primitive_w = se.mv(Li,wa[:3])
        const = se.mv(P,wa[:3])+[0]*d
        rhs = se.mv(B.Laminv,[x-y for x,y in zip(wb,const)])
        for row in range(n):
            h,g = [0]*ni,[0]*nr
            for i in range(d):
                for k in range(d):
                    h[i*d+k] += B.Laminv[row][3+i]*wa[3+k]
                for k in range(3):
                    h[ne+i*3+k] += B.Laminv[row][3+i]*primitive_w[k]
            h[ne+nj+number*n+row] = -1
            g[ny:] = shift[row]
            add(h,g,rhs[row])
    # RREF puts real-variable pivots first. Rows with zero real block are
    # exactly the residual compatibility equations on the integer variables.
    reduced,_ = se.rref(equations)
    residual = [row[nr:] for row in reduced if not any(row[:nr])]
    return integer_projection([r[:-1] for r in residual],[r[-1] for r in residual],ne,d)


def e_candidates(A, B, P, Pi):
    """Filter the original bounded E stream, preserving its order and bounds."""
    lattice = feasible_e_lattice(A,B,P,Pi)
    if lattice is None:
        yield from se._e_candidates(A,B,P,Pi)
        return
    if lattice.empty or not all(lattice.has_unimodular_residue(p) for p in (2,3)):
        return
    for E in se._e_candidates(A,B,P,Pi):
        if lattice.contains(E):
            yield E
