#!/usr/bin/env python3
"""subgroups.py -- k-index-1 (translationengleiche) subgroups of a (3+d)D superspace group (d=1,2,3).

For incommensurate q the point group K_s of G is isomorphic to the basic point group K, and a
subgroup with the same superspace lattice Σ is the full preimage of a subgroup of K.  So the
k-index-1 subgroups up to conjugacy in G are the subgroups of K up to K-conjugacy.  Each class
is identified against the BYU records with the same basic space-group number using
``ssg_equiv.equivalent`` at the parent's generic q (CONTRACTS §5).

Usage:  python3 subgroups.py 165.1.24.1
"""
import functools, itertools, os, sys

# Direct CLI invocation must enforce the resource budget before exact/numpy.
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # also under python -I
import exact  # noqa: E402
import ssg_equiv as se  # noqa: E402
import affine_search  # noqa: E402


def _closure(gens, table):
    G = {0}
    frontier = [0]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = table[a][g]
                if b not in G:
                    G.add(b)
                    nxt.append(b)
        frontier = nxt
    return frozenset(G)


def subgroup_record(rec, ops, idx):
    """Group record (CONTRACTS §2) of the preimage of the point operations at op indices `idx`."""
    gen = [exact.format_op([list(r) for r in ops[i][0]], list(ops[i][1])) for i in idx if i != 0]
    r = dict(sym=f"<subgroup of {rec['sym']}>", cls="", d=rec["d"], q=list(rec["q"]),
             params=list(rec["params"]), cen=list(rec["cen"]), gen=gen)
    r["nops"] = len(exact.build_ops(r))
    return r


# A process handles at most three dimensions. Retain the database object as well
# as its id so a caller's short-lived synthetic database cannot alias a cache.
_candidate_cache = []


# This invariant uses the full superspace translation lattice, including coupled
# external/internal centring. It makes no product-lattice or q-family assumption.
@functools.lru_cache(maxsize=8192)
def _torus_smith(W, L):
    n = len(W)
    Li = se.im.inv_rat(L)
    A = se.mm(se.mm(Li, [[int(i == j)-W[i][j] for j in range(n)] for i in range(n)]), L)
    if any(x.denominator != 1 for row in A for x in row):
        raise ValueError('point action does not preserve the full translation lattice')
    D, U, V, Vi = se.im.snf([[int(x) for x in row] for row in A])
    factors = tuple(abs(D[i][i]) for i in range(n) if D[i][i])
    zero_rows = tuple(i for i in range(n) if not D[i][i])
    return factors, tuple(tuple(row) for row in U), zero_rows


@functools.lru_cache(maxsize=20000)
def _fixed_torus_profile(gens, centrings, d):
    from fractions import Fraction
    n = 3+d
    L = tuple(tuple(row) for row in se.lattice_basis(
        [[Fraction(x) for x in c.split(',')] for c in centrings], n))
    Li = se.im.inv_rat(L)
    rec = dict(d=d, gen=gens, cen=centrings)
    first = {}
    for W, w in exact.build_ops(rec):
        first.setdefault(W, w)
    profile = []
    for W, w in first.items():
        factors, U, zero_rows = _torus_smith(W, L)
        Ub = se.mv(U, se.mv(Li, w))
        nonempty = all(Ub[i].denominator == 1 for i in zero_rows)
        profile.append((factors, nonempty))
    return tuple(sorted(profile))


def fixed_torus_profile(rec):
    """Multiset of Smith factors and affine fixed-set nonemptiness on R^n/Σ.

    Nonempty entries have dimension n-len(factors) and product(factors)
    connected components (empty product one). Representatives with equal W
    differ by Σ, so the result does not depend on the chosen affine lift.
    """
    return _fixed_torus_profile(tuple(rec['gen']), tuple(rec['cen']), rec['d'])


@functools.lru_cache(maxsize=20000)
def _phase_character_profile(gens, centrings, d):
    """Finite phase quotient for scalar internal action and product lattice.

    None means ineligible, never an obstructing value. The only allowed phase
    comparison is between two eligible profiles; see CONTRACTS proof.
    """
    from fractions import Fraction as F
    if any(F(x).denominator != 1 for c in centrings for x in c.split(',')[3:]):
        return None
    first = {}
    for W, w in exact.build_ops(dict(d=d, gen=gens, cen=centrings)):
        first.setdefault(W, w)
    if len(first) > 24:
        return None  # conservative complexity guard, never reject on ineligibility
    L3 = se.lattice_basis([[F(x) for x in c.split(',')[:3]] for c in centrings], 3)
    L3i = se.im.inv_rat(L3)
    rotations, epsilons, translations, phases = [], [], [], []
    for W, w in first.items():
        eps = W[3][3]
        if eps not in (-1, 1):
            return None
        if (any(W[i][j] for i in range(3,3+d) for j in range(3)) or
                any(W[i][j] for i in range(3) for j in range(3,3+d)) or
                any(W[i][3+j] != eps*int(i==3+j) for i in range(3,3+d) for j in range(d))):
            return None
        R = se.mm(se.mm(L3i, [list(row[:3]) for row in W[:3]]), L3)
        if any(x.denominator != 1 for row in R for x in row):
            return None
        rotations.append(R)
        epsilons.append(eps)
        translations.append(se.mv(L3i, w[:3]))
        phases.append(w[3:])
    m = len(rotations)
    constraints = [[R[i][j]-eps*int(i==j) for i in range(3)]
                   for R,eps in zip(rotations,epsilons) for j in range(3)]
    basis = se.nullspace(constraints, 3)
    scale = se.den([x for row in basis for x in row])
    covariant = se.im.saturate([[int(x*scale) for x in row] for row in basis])
    N = se.nullspace([[F(1-eps) for eps in epsilons]], m)
    if not N:
        return ()
    generators = [list(col) for col in se.tr(N)]
    for u in covariant:
        h = [sum(x*y for x,y in zip(u, w)) for w in translations]
        generators.append(se.mv(N,h))
    B0 = se._rational_generator_lattice(generators, len(N))
    B1 = se._rational_generator_lattice(
        generators+[se.mv(N,[w[i] for w in phases]) for i in range(d)], len(N))
    inclusion = se.mm(se.im.inv_rat(B1), B0)
    if any(x.denominator != 1 for row in inclusion for x in row):
        raise ValueError('phase quotient lattices do not nest integrally')
    D, _, _, _ = se.im.snf([[int(x) for x in row] for row in inclusion])
    return tuple(abs(D[i][i]) for i in range(len(D)) if abs(D[i][i]) > 1)


def phase_character_profile(rec):
    return _phase_character_profile(tuple(rec['gen']), tuple(rec['cen']), rec['d'])


def _phase_compatible(a, b):
    pa, pb = phase_character_profile(a), phase_character_profile(b)
    return pa is None or pb is None or pa == pb


def _type_key(rec, q):
    d, order, rank, characters = se.invariant(rec, q)
    # Family dimension/module rank is deliberately absent: a subgroup of a
    # narrower parent family can have a canonical type with more free parameters.
    return d, order, characters, fixed_torus_profile(rec)


def _candidates(db, sg, type_key):
    """Cheap SG catalogue first; build expensive profiles only for requested SGs."""
    for source, by_sg, indexed in _candidate_cache:
        if source is db:
            break
    else:
        by_sg, indexed = {}, {}
        for label, cand in db.items():
            by_sg.setdefault(int(label.split('.')[0]), []).append((label,cand))
        _candidate_cache.append((db, by_sg, indexed))
        del _candidate_cache[:-3]
    if sg not in indexed:
        buckets = {}
        for label, cand in by_sg.get(sg, []):
            key = _type_key(cand, se.q_generic(cand))
            buckets.setdefault(key, []).append((label,cand))
        indexed[sg] = buckets
    return indexed[sg].get(type_key, [])


def t_subgroups(rec, db=None):
    """K-conjugacy classes of subgroups of K, each with its BYU identification.

    Returns a list of dicts sorted by (index, BYU label as a tuple of ints, pg_ops):
      pg_ops      op indices (canonical list of rec) of one op per point operation, sorted
      index       |K| / |P|
      class_size  number of K-conjugates
      sg          basic space-group number of the subgroup
      matches     BYU labels equivalent to the subgroup at the parent's generic q
      byu, byu_sym, map   the unique match and the map from the subgroup onto it (if unique)
    """
    if rec["d"] not in (1, 2, 3):
        raise ValueError("supported dimensions are d = 1, 2, 3")
    db = db if db is not None else se.load_records(rec["d"])
    ops = exact.build_ops(rec)
    first = {}
    for i, (W, w) in enumerate(ops):
        first.setdefault(W, i)
    K = list(first)                                   # point operations, in canonical-list order
    n = len(K)
    pos = {W: i for i, W in enumerate(K)}
    table = [[pos[se.mul4(a, b)] for b in K] for a in K]
    inv = [next(j for j in range(n) if table[i][j] == 0) for i in range(n)]
    assert K[0] == ops[0][0], "canonical point-operation order must begin with the identity"

    subs = set()
    for r in range(4):                                # every point group is generated by <= 3 ops
        for gens in itertools.combinations(range(1, n), r):
            subs.add(_closure(gens, table))

    def conj(S, g):
        return frozenset(table[table[g][h]][inv[g]] for h in S)

    def op_idx(S):
        return tuple(sorted(first[K[h]] for h in S))

    classes, seen = [], set()
    for S in sorted(subs, key=lambda S: (len(S), op_idx(S))):
        if S in seen:
            continue
        orbit = {conj(S, g) for g in range(n)}
        seen |= orbit
        rep = min(orbit, key=op_idx)
        classes.append((rep, len(orbit)))

    q = se.q_generic(rec)
    out = []
    for S, size in classes:
        idx = op_idx(S)
        sub = subgroup_record(rec, ops, idx)
        sg = se.prep(sub).sg
        matches, maps = [], []
        for label, cand in _candidates(db, sg, _type_key(sub, q)):
            if not _phase_compatible(sub, cand):
                continue
            res = se.equivalent(sub, q, cand, e_candidates=affine_search.e_candidates)
            if res is not None:
                matches.append(label)
                maps.append(res)
        c = dict(pg_ops=list(idx), index=n // len(S), class_size=size, sg=sg, matches=matches)
        if len(matches) == 1:
            c.update(byu=matches[0], byu_sym=db[matches[0]]["sym"], map=maps[0])
        out.append(c)
    out.sort(key=lambda c: (c["index"], tuple(int(x) for x in c.get("byu", "0").split(".")), c["pg_ops"]))
    return out


if __name__ == "__main__":
    db = se.load_records(int(sys.argv[1].split(".")[1]))
    rec = db[sys.argv[1]]
    print(rec["sym"])
    for c in t_subgroups(rec, db):
        print(f"  [{c['index']:2d}] x{c['class_size']} sg {c['sg']:3d} {c['matches']} "
              f"{c.get('byu_sym', '')}  ops {c['pg_ops']}")
