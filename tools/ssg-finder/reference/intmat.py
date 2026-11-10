"""Exact integer matrix normal forms (pure Python ints, no numpy).

Matrices are lists of rows of Python ints.  All functions return new lists.

* ``hnf(A)``          row-style Hermite normal form: ``H = U A`` with U unimodular,
                      H in row-echelon form (zero rows dropped), pivots positive and
                      entries above each pivot reduced into ``[0, pivot)``.
* ``hnf_lower(A)``    the same lattice, but lower-triangular for a full-rank square
                      result (pivots on the diagonal, entries *below* each pivot
                      reduced into ``[0, pivot)``).  Rows = basis vectors.
* ``snf(A)``          Smith normal form ``U A V = D`` with U, V unimodular and
                      ``d_1 | d_2 | ...``.  Also returns ``V^-1``.
* ``saturate(A)``     basis (in HNF) of ``span_Q(rows A) ∩ Z^n``.
"""
from fractions import Fraction


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def matmul(A, B):
    if not A:
        return []
    m = len(B[0]) if B else 0
    return [[sum(a * B[k][j] for k, a in enumerate(row)) for j in range(m)] for row in A]


def matvec(A, v):
    return [sum(a * b for a, b in zip(row, v)) for row in A]


def transpose(A, ncols=None):
    if not A:
        return [[] for _ in range(ncols or 0)]
    return [list(r) for r in zip(*A)]


def det(A):
    """Exact determinant (Fractions or ints) by fraction-free Gaussian elimination."""
    n = len(A)
    M = [[Fraction(x) for x in r] for r in A]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d


def inv_rat(A):
    """Exact inverse of a square rational matrix (list of rows), as Fractions."""
    n = len(A)
    M = [[Fraction(x) for x in r] + [Fraction(int(i == j)) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [r[n:] for r in M]


def rank(A):
    if not A:
        return 0
    M = [[Fraction(x) for x in r] for r in A]
    rk, ncol = 0, len(M[0])
    for c in range(ncol):
        p = next((r for r in range(rk, len(M)) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / M[rk][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[rk])]
        rk += 1
    return rk


# ---------------------------------------------------------------- HNF
def hnf_with_transform(A, ncols=None):
    """Return (H, U) with H = U A (U square unimodular, len(A) x len(A)).

    H keeps all len(A) rows; the non-zero rows come first, in echelon form.
    """
    m = len(A)
    n = len(A[0]) if A else (ncols or 0)
    H = [list(map(int, r)) for r in A]
    U = eye(m)
    row = 0
    for col in range(n):
        if row >= m:
            break
        # Euclid on column `col`, rows row..m-1
        while True:
            nz = [i for i in range(row, m) if H[i][col] != 0]
            if not nz:
                break
            piv = min(nz, key=lambda i: abs(H[i][col]))
            if piv != row:
                H[row], H[piv] = H[piv], H[row]
                U[row], U[piv] = U[piv], U[row]
            done = True
            for i in range(row + 1, m):
                if H[i][col] != 0:
                    q = H[i][col] // H[row][col]
                    H[i] = [a - q * b for a, b in zip(H[i], H[row])]
                    U[i] = [a - q * b for a, b in zip(U[i], U[row])]
                    if H[i][col] != 0:
                        done = False
            if done:
                break
        if row < m and H[row][col] != 0:
            if H[row][col] < 0:
                H[row] = [-a for a in H[row]]
                U[row] = [-a for a in U[row]]
            p = H[row][col]
            for i in range(row):
                q = H[i][col] // p
                if q:
                    H[i] = [a - q * b for a, b in zip(H[i], H[row])]
                    U[i] = [a - q * b for a, b in zip(U[i], U[row])]
            row += 1
    return H, U


def hnf(A, ncols=None):
    """Row HNF, zero rows dropped."""
    H, _ = hnf_with_transform(A, ncols)
    return [r for r in H if any(r)]


def hnf_lower(A):
    """Lower-triangular HNF basis of the row lattice of A (full column rank expected).

    Pivots are on the diagonal; entries below a pivot (same column) lie in [0, pivot).
    Obtained by reversing the column order, taking the row HNF, and reversing back.
    """
    n = len(A[0])
    R = [list(reversed(r)) for r in A]
    H = hnf(R)
    H = [list(reversed(r)) for r in H]
    H.reverse()
    return H


# ---------------------------------------------------------------- SNF
def snf(A):
    """Smith normal form.  Returns (D, U, V, Vinv) with U A V = D (m x n)."""
    m = len(A)
    n = len(A[0]) if A else 0
    D = [list(map(int, r)) for r in A]
    U, V, Vi = eye(m), eye(n), eye(n)

    def swap_rows(i, j):
        D[i], D[j] = D[j], D[i]
        U[i], U[j] = U[j], U[i]

    def swap_cols(i, j):
        for r in D:
            r[i], r[j] = r[j], r[i]
        for r in V:
            r[i], r[j] = r[j], r[i]
        Vi[i], Vi[j] = Vi[j], Vi[i]

    def add_row(dst, src, k):  # row dst += k row src
        D[dst] = [a + k * b for a, b in zip(D[dst], D[src])]
        U[dst] = [a + k * b for a, b in zip(U[dst], U[src])]

    def add_col(dst, src, k):  # col dst += k col src ; Vinv: row src -= k row dst
        for r in D:
            r[dst] += k * r[src]
        for r in V:
            r[dst] += k * r[src]
        Vi[src] = [a - k * b for a, b in zip(Vi[src], Vi[dst])]

    t = 0
    while t < min(m, n):
        nz = [(abs(D[i][j]), i, j) for i in range(t, m) for j in range(t, n) if D[i][j] != 0]
        if not nz:
            break
        _, i, j = min(nz)
        if i != t:
            swap_rows(i, t)
        if j != t:
            swap_cols(j, t)
        while True:
            changed = False
            for i in range(t + 1, m):
                if D[i][t]:
                    q = D[i][t] // D[t][t]
                    add_row(i, t, -q)
                    if D[i][t]:
                        changed = True
            for j in range(t + 1, n):
                if D[t][j]:
                    q = D[t][j] // D[t][t]
                    add_col(j, t, -q)
                    if D[t][j]:
                        changed = True
            if changed:
                cand = [(abs(D[i][t]), i, t) for i in range(t + 1, m) if D[i][t]] + \
                       [(abs(D[t][j]), t, j) for j in range(t + 1, n) if D[t][j]]
                _, i, j = min(cand)
                if i != t:
                    swap_rows(i, t)
                if j != t:
                    swap_cols(j, t)
                continue
            # divisibility
            bad = next(((i, j) for i in range(t + 1, m) for j in range(t + 1, n)
                        if D[i][j] % D[t][t]), None)
            if bad is None:
                break
            add_row(t, bad[0], 1)
        if D[t][t] < 0:
            D[t] = [-a for a in D[t]]
            U[t] = [-a for a in U[t]]
        t += 1
    return D, U, V, Vi


def saturate(A, ncols=None):
    """HNF basis of the saturation of the row lattice of A."""
    if not A or not any(any(r) for r in A):
        return []
    D, U, V, Vi = snf(A)
    rho = sum(1 for i in range(min(len(D), len(D[0]))) if D[i][i])
    return hnf(Vi[:rho])


def is_unimodular(A):
    return abs(det(A)) == 1
