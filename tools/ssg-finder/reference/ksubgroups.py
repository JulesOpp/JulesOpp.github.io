#!/usr/bin/env python3
"""Exact finite-index subgroup oracle, with lattices represented by row lower HNF.

All affine computations use column coordinates. A sigma row basis H therefore
has column basis H.T in the parent's full translation-lattice coordinates.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
import itertools
import sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import subgroups as ts

se, exact, im = ts.se, ts.exact, ts.se.im


def frozen(A):
    return tuple(tuple(r) for r in A)


def sublattices(n, k):
    """All row lower HNFs of determinant k, without a prime-index assumption."""
    if n < 1 or k < 1:
        raise ValueError('positive dimension and index required')
    def diagonals(m, value):
        if m == 1:
            yield (value,)
        else:
            for a in range(1, value+1):
                if value % a == 0:
                    for rest in diagonals(m-1, value//a):
                        yield (a,)+rest
    out = []
    for diag in diagonals(n, k):
        places = [(i,j) for i in range(n) for j in range(i)]
        for values in itertools.product(*(range(diag[j]) for i,j in places)):
            H = [[diag[i] if i == j else 0 for j in range(n)] for i in range(n)]
            for (i,j),v in zip(places,values):
                H[i][j] = v
            out.append(H)
    return sorted(out)


def reduce(v, H):
    """Canonical affine vector modulo the row lattice H, in ambient coordinates."""
    v = list(map(F,v))
    for i in reversed(range(len(H))):
        q = v[i]//H[i][i]
        v = [x-q*h for x,h in zip(v,H[i])]
    return tuple(v)


def lattice(rec):
    return se.lattice_basis([[F(x) for x in c.split(',')] for c in rec['cen']], 3+rec['d'])


def context(rec):
    L = lattice(rec)
    Li = im.inv_rat(L)
    ops = exact.build_ops(rec)
    first = {}
    for i,(W,w) in enumerate(ops):
        first.setdefault(W,i)
    idx = list(first.values())
    lifts = [(frozen(se.mm(se.mm(Li,ops[i][0]),L)),tuple(se.mv(Li,ops[i][1]))) for i in idx]
    pos = {W:i for i,(W,w) in enumerate(lifts)}
    table = [[pos[frozen(se.mm(a[0],b[0]))] for b in lifts] for a in lifts]
    return L, ops, idx, lifts, table


def point_subgroups(table):
    found = {}
    for r in range(4):
        for gens in itertools.combinations(range(1,len(table)),r):
            found.setdefault(ts._closure(gens,table),gens)
    return sorted(found.items(),key=lambda item:(len(item[0]),tuple(sorted(item[0]))))


def multiply(a,b,H):
    W,w = a
    V,v = b
    return frozen(se.mm(W,V)), reduce([x+y for x,y in zip(w,se.mv(W,v))],H)


def close(gens,H,limit=None):
    n = len(H)
    identity = frozen(im.eye(n)), (F(0),)*n
    seen, frontier = {identity}, [identity]
    while frontier:
        a = frontier.pop()
        for g in gens:
            b = multiply(a,g,H)
            if b not in seen:
                seen.add(b)
                if limit is not None and len(seen) > limit:
                    return None
                frontier.append(b)
    return frozenset(seen)


def pair_key(H,ops):
    return frozen(H), tuple(sorted({(W,reduce(w,H)) for W,w in ops}))


def conjugate(key,g):
    H,ops = key
    V,v = g
    Vi = im.inv_rat(V)
    B = se.mm(V,se.tr(H))
    J = im.hnf_lower([[int(x) for x in row] for row in se.tr(B)])
    out = []
    for W,w in ops:
        R = frozen(se.mm(se.mm(V,W),Vi))
        u = [a+b-c for a,b,c in zip(v,se.mv(V,w),se.mv(R,v))]
        out.append((R,u))
    return pair_key(J,out)


def orbit(key,gens):
    seen, todo = {key}, [key]
    while todo:
        a = todo.pop()
        for g in gens:
            b = conjugate(a,g)
            if b not in seen:
                seen.add(b)
                todo.append(b)
    return seen


def kind(rec,H):
    """Compare exact internal intersections and external projections."""
    L = lattice(rec)
    B = se.mm(L,se.tr(H))
    def block_basis(M):
        scale = se.den([x for row in M for x in row])
        return [[F(x,scale) for x in row] for row in im.hnf([[int(x*scale) for x in row] for row in se.tr(M)])]
    parent, child = block_basis(L), block_basis(B)
    if parent[3:] == child[3:]:
        return '3d'
    def projection(rows):
        scale = se.den([x for row in rows for x in row[:3]])
        return tuple(tuple(F(x,scale) for x in r) for r in im.hnf([[int(x*scale) for x in row[:3]] for row in rows]))
    return 'internal' if projection(parent) == projection(child) else 'mixed'


def rebase(rec, ops, sigma):
    """Return a primitive subgroup record and S=B^-1, with B block triangular.

    ops are actual affine operators in parent coordinates; sigma is a row lower
    HNF in the full parent lattice basis. q is transformed as D^-1(q A-C).
    """
    import sympy as sp
    n = 3+rec['d']
    B0 = se.mm(lattice(rec),se.tr(sigma))
    scale = se.den([x for r in B0 for x in r])
    rows = im.hnf([[int(x*scale) for x in row] for row in se.tr(B0)])
    B = se.tr([[F(x,scale) for x in row] for row in rows])
    # Positive orientation in the physical space is mandatory.
    A = [r[:3] for r in B[:3]]
    if sp.Matrix(A).det() < 0:
        for row in B:
            row[0] = -row[0]
    S = im.inv_rat(B)
    gen = [exact.format_op(se.mm(se.mm(S,W),B),se.mv(S,w)) for W,w in ops]
    A = sp.Matrix([r[:3] for r in B[:3]])
    C = sp.Matrix([r[:3] for r in B[3:]])
    D = sp.Matrix([r[3:] for r in B[3:]])
    symbols = {p:sp.Symbol(p) for p in rec['params']}
    Q = sp.Matrix([[sp.sympify(x,locals=symbols) for x in row.split(',')] for row in rec['q']])
    Q2 = D.inv()*(Q*A-C)
    q = [','.join(str(sp.expand(Q2[i,j])) for j in range(3)) for i in range(rec['d'])]
    rec2 = dict(sym='<rebased subgroup>',cls='',d=rec['d'],q=q,params=list(rec['params']),cen=[],gen=gen)
    rec2['nops'] = len(exact.build_ops(rec2))
    return rec2,S


def row_key(rec,row):
    L,ops,idx,lifts,table = context(rec)
    pos = {i:j for j,i in enumerate(idx)}
    Li = im.inv_rat(L)
    return pair_key(row['sigma'],[(lifts[pos[i]][0],se.mv(Li,[F(x) for x in w])) for i,w in zip(row['pg_ops'],row['tw'])])


def full_orbit(rec,row,k):
    """Independent exhaustive conjugation over every element of G/kΣ."""
    L,ops,idx,lifts,table = context(rec)
    key = row_key(rec,row)
    return {conjugate(key,(W,[x+y for x,y in zip(w,t)]))
            for W,w in lifts for t in itertools.product(range(k),repeat=len(L))}


def k_subgroups(rec,k,db=None,*,stats=None):
    """All G-conjugacy classes at lattice index k, identified after rebasing.

    db={} suppresses catalogue identification for public synthetic tests. stats
    receives invariant_lattices (sum over P) and complement_candidates counts.
    """
    if rec['d'] != 1:
        raise ValueError('this implementation is gated to d=1')
    if k < 1:
        raise ValueError('positive index required')
    L,ops,idx,lifts,table = context(rec)
    n = len(L)
    if k == 1:
        out = ts.t_subgroups(rec,db)
        for r in out:
            r.update(sigma=im.eye(n),kind='3d',tw=[list(ops[i][1]) for i in r['pg_ops']],t_index=r['index'],k_index=1)
        return out
    db = se.load_records(1) if db is None else db
    gens = lifts+[(frozen(im.eye(n)),tuple(int(i==j) for j in range(n))) for i in range(n)]
    seen, out = set(), []
    counters = dict(invariant_lattices=0,complement_candidates=0)
    pointsubs = point_subgroups(table)
    for H in sublattices(n,k):
        reps = list(itertools.product(*(range(H[i][i]) for i in range(n))))
        for P,pgens in pointsubs:
            if any(any(reduce(se.mv(lifts[g][0],h),H)) for g in pgens for h in H):
                continue
            counters['invariant_lattices'] += 1
            for choices in itertools.product(reps,repeat=len(pgens)):
                counters['complement_candidates'] += 1
                gg = [(lifts[g][0],[a+b for a,b in zip(lifts[g][1],t)]) for g,t in zip(pgens,choices)]
                sub = close(gg,H,len(P))
                if sub is None:
                    continue
                key = pair_key(H,sub)
                if key in seen:
                    continue
                orb = orbit(key,gens)
                seen.update(orb)
                J,aff = min(orb)
                by_W = {W:i for i,(W,w) in zip(idx,lifts)}
                ordered = sorted(aff,key=lambda a:by_W[a[0]])
                physical = [(frozen(se.mm(se.mm(L,W),im.inv_rat(L))),tuple(se.mv(L,w))) for W,w in ordered]
                row = dict(pg_ops=[by_W[W] for W,w in ordered],sigma=[list(r) for r in J],
                           kind=kind(rec,J),tw=[list(w) for W,w in physical],t_index=len(lifts)//len(P),
                           k_index=k,index=k*len(lifts)//len(P),class_size=len(orb),matches=[])
                if db:
                    rebased,S = rebase(rec,physical,J)
                    q = se.q_generic(rebased)
                    sg = se.prep(rebased).sg
                    row['sg'] = sg
                    maps = []
                    for label,cand in ts._candidates(db,sg,ts._type_key(rebased,q)):
                        if not ts._phase_compatible(rebased,cand):
                            continue
                        result = se.equivalent(rebased,q,cand,e_candidates=ts.affine_search.e_candidates)
                        if result is not None:
                            row['matches'].append(label)
                            maps.append(result)
                    if len(maps) == 1:
                        row.update(byu=row['matches'][0],byu_sym=db[row['matches'][0]]['sym'],map=maps[0])
                out.append(row)
    if stats is not None:
        stats.update(counters)
    return sorted(out,key=lambda r:(r['index'],tuple(int(x) for x in r.get('byu','0').split('.')),r['sigma'],r['pg_ops'],r['tw']))


def brute_subgroups(rec,k):
    """Independent subgroup joining in G/kΣ, filtered by translation index.

    Enumerates finite subgroups by adjoining cyclic generators, without using
    point-group generators, invariant lattices, or complement enumeration.
    """
    L,ops,idx,lifts,table = context(rec)
    n = len(L)
    H = [[k*int(i==j) for j in range(n)] for i in range(n)]
    elements = sorted({(W,reduce([a+b for a,b in zip(w,t)],H)) for W,w in lifts
                       for t in itertools.product(range(k),repeat=n)})
    ident = frozen(im.eye(n))
    zero = ident,(F(0),)*n
    pos = {e:i for i,e in enumerate(elements)}
    identity = pos[zero]
    translations = {i for i,(W,w) in enumerate(elements) if W == ident}
    allowed = k**(n-1)
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def product(a,b):
        return pos[multiply(elements[a],elements[b],H)]
    def join(gens):
        found, todo, nt = {identity}, [identity], 1
        while todo:
            a = todo.pop()
            for g in gens:
                b = product(a,g)
                if b not in found:
                    found.add(b)
                    nt += b in translations
                    if nt > allowed:
                        return None
                    todo.append(b)
        return frozenset(found)
    initial = frozenset([identity])
    seen, todo, qualifying = {initial}, [(initial,())], []
    while todo:
        group, gens = todo.pop()
        if len(group & translations) == allowed:
            qualifying.append(group)
        for e in range(len(elements)):
            if e in group:
                continue
            joined = join(gens+(e,))
            if joined is not None and joined not in seen:
                seen.add(joined)
                todo.append((joined,gens+(e,)))
    # Quotient-group conjugation is independent of the lattice/orbit algorithm.
    inverses = {g:next(h for h in range(len(elements)) if product(g,h) == identity)
                for g in range(len(elements))}
    conjugate_sets, result = set(), set()
    for group in qualifying:
        if group in conjugate_sets:
            continue
        conjugates = {frozenset(product(product(g,h),inverses[g]) for h in group)
                      for g in range(len(elements))}
        conjugate_sets.update(conjugates)
        keys = []
        for sub in conjugates:
            trans = [elements[i][1] for i in sub & translations]
            J = im.hnf_lower([[int(x) for x in w] for w in trans]+H)
            keys.append(pair_key(J,[elements[i] for i in sub]))
        result.add(min(keys))
    return result
