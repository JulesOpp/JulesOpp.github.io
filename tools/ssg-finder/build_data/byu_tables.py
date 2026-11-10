#!/usr/bin/env python3
"""
byu_tables.py -- parse BYU ISO(3+d)D `ssg_tables.zip`, convert every superspace group
from the supercentred to the basic setting, check invariants, and write
`local/data/sgNNN.json` + `local/data/index.json` (CONTRACTS.md §2).

Everything here is exact (fractions.Fraction / ints); no sympy.

Conversion (PLAN §2.2): basis A = P a  =>  coordinates x_basic = Pᵀ X_super, so
    W_b = Pᵀ W_s P⁻ᵀ,   w_b = Pᵀ w_s (mod 1),   c_b = Pᵀ c_s (mod 1).

Licence (decision D1): the input and every output are BYU-derived and must stay
under tools/ssg-finder/local/ (gitignored).

Usage:
    /opt/anaconda3/bin/python3 -I tools/ssg-finder/build_data/byu_tables.py [--zip Z] [--out DIR] [--no-write]
"""
import argparse, hashlib, json, math, os, random, re, sys, time, zipfile
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
LOCAL = os.path.normpath(os.path.join(HERE, "..", "local"))
ZIP_DEFAULT = os.path.join(LOCAL, "raw", "ssg_tables.zip")
OUT_DEFAULT = os.path.join(LOCAL, "data")

SUPER_VARS = "XYZTUV"
BASIC_VARS = "xyztuv"


class BuildError(Exception):
    pass


# ============================================================ small exact helpers
def lin(expr, names):
    """Parse a linear expression such as '-2a1+a4', 'g+1/2', 'X+1/2', '2*a3', '0'.

    Returns (coef dict name->Fraction, constant Fraction).  `names` is a set (or
    None to accept any identifier [A-Za-z]+\\d*).
    """
    e = expr.replace(" ", "")
    if not e:
        raise BuildError(f"empty expression")
    if e[0] not in "+-":
        e = "+" + e
    coef, const = {}, Fr(0)
    pos = 0
    for m in re.finditer(r"([+-])([^+-]+)", e):
        if m.start() != pos:
            raise BuildError(f"cannot parse {expr!r}")
        pos = m.end()
        s = 1 if m.group(1) == "+" else -1
        body = m.group(2)
        t = re.fullmatch(r"(\d+(?:/\d+)?)?\*?([A-Za-z]+\d*)", body)
        if t:
            v = t.group(2)
            if names is not None and v not in names:
                raise BuildError(f"unknown variable {v!r} in {expr!r}")
            coef[v] = coef.get(v, Fr(0)) + s * (Fr(t.group(1)) if t.group(1) else 1)
        elif re.fullmatch(r"\d+(?:/\d+)?", body):
            const += s * Fr(body)
        else:
            raise BuildError(f"cannot parse term {body!r} in {expr!r}")
    if pos != len(e):
        raise BuildError(f"cannot parse {expr!r}")
    return coef, const


def lin_vec(expr, names):
    """Linear expression over an ordered variable list -> (coef list, const)."""
    c, k = lin(expr, set(names))
    return [c.get(v, Fr(0)) for v in names], k


def lin_eval(expr, values):
    c, k = lin(expr, set(values))
    return k + sum(a * values[v] for v, a in c.items())


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def matvec(A, v):
    return [sum(a * b for a, b in zip(r, v)) for r in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def inverse(A):
    """Gauss-Jordan inverse over Fractions."""
    n = len(A)
    M = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            raise BuildError("singular matrix")
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def mod1(x):
    return x - math.floor(x)


def fmt_rat(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_op(W, w, n):
    """CONTRACTS §1 xyzt string: terms ordered by variable, constant last, w in [0,1)."""
    rows = []
    for i in range(n):
        s = ""
        for j in range(n):
            c = W[i][j]
            if c:
                a = abs(c)
                s += ("+" if c > 0 else "-") + ("" if a == 1 else fmt_rat(a)) + BASIC_VARS[j]
        if w[i]:
            s += "+" + fmt_rat(w[i])
        s = s[1:] if s.startswith("+") else s
        rows.append(s or "0")
    return ",".join(rows)


def parse_op(s, n, names=BASIC_VARS):
    """'-x+1/2,y,...' -> (W list of int lists, w list of Fractions mod 1)."""
    parts = s.split(",")
    if len(parts) != n:
        raise BuildError(f"operator {s!r}: {len(parts)} components, expected {n}")
    W, w = [], []
    for p in parts:
        c, k = lin_vec(p, list(names[:n]))
        W.append(c)
        w.append(mod1(k))
    return W, w


# ============================================================ table parsing
def split_blocks(text):
    """Blocks are separated by a line containing only whitespace."""
    out, cur = [], []
    for line in text.splitlines():
        if line.strip():
            cur.append(line.rstrip())
        elif cur:
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


# reflection-condition prefixes: HKLM:, H-H00N0:, 2KKLM-MP:, 0KL2NNP: ...
REFL_RE = re.compile(r"(?:(?<=\s)|^)[HKLMNP0-9-]+:")
LABEL_RE = re.compile(r"(\d+)\.(\d)\.(\d+)\.(\d+)$")


def parse_group_block(lines, d):
    n = 3 + d
    head = lines[0].split(None, 1)
    if len(head) != 2:
        raise BuildError(f"bad header {lines[0]!r}")
    label, sym = head[0], head[1].strip()
    if not LABEL_RE.fullmatch(label):
        raise BuildError(f"bad label {label!r}")
    rest = " ".join(l.strip() for l in lines[1:])
    # leading run of parenthesised tuples: d supercentred q's (3 comps), then centrings (n comps)
    qs, cen = [], []
    pos = 0
    while True:
        m = re.compile(r"\s*[;,]?\s*\(([^()]*)\)").match(rest, pos)
        if not m:
            break
        comps = [c.strip() for c in m.group(1).split(",")]
        if len(comps) == 3 and len(qs) < d and not cen:
            qs.append(comps)
        elif len(comps) == n and len(qs) == d:
            cen.append([Fr(c) for c in comps])
        else:
            raise BuildError(f"{label}: unexpected tuple ({m.group(1)})")
        pos = m.end()
    if len(qs) != d or not cen:
        raise BuildError(f"{label}: {len(qs)} q's / {len(cen)} centrings")
    body = rest[pos:]
    r = REFL_RE.search(body)
    ops_txt, refl = (body[:r.start()], body[r.start():].strip()) if r else (body, "")
    ops = [o for o in re.sub(r"\s+", "", ops_txt).split(";") if o]
    for o in ops:
        if len(o.split(",")) != n or not re.fullmatch(r"[-+XYZTUV0-9/,]+", o):
            raise BuildError(f"{label}: bad operator {o!r}")
    return dict(label=label, sym=sym, d=d, qs=qs, cen=cen, ops=ops, refl=refl)


def parse_lattice_block(lines, d):
    n = 3 + d
    head = lines[0].split(None, 1)
    cls, sym = head[0], head[1].strip()
    if not re.fullmatch(rf"{d}\.\d+", cls):
        raise BuildError(f"bad class {cls!r}")
    rest = " ".join(l.strip() for l in lines[1:])
    iq, iw = rest.find("Q1="), rest.find("where")
    if iq < 0:
        raise BuildError(f"{cls}: no Q1=")
    a_txt = rest[:iq]
    A = dict(re.findall(r"A(\d)=([^,]+)", a_txt))
    if sorted(A) != [str(i + 1) for i in range(n)]:
        raise BuildError(f"{cls}: A lines {sorted(A)}")
    names = [f"a{j + 1}" for j in range(n)]
    P = [lin_vec(A[str(i + 1)].strip(), names) for i in range(n)]
    if any(k != 0 for _, k in P):
        raise BuildError(f"{cls}: constant in A_i")
    P = [c for c, _ in P]
    q_end = iw if iw >= 0 else len(rest)
    Q = re.findall(r"Q(\d)=\(([^()]*)\)", rest[iq:q_end])
    if [int(i) for i, _ in Q] != list(range(1, d + 1)):
        raise BuildError(f"{cls}: Q's {Q}")
    Q = [[c.strip() for c in q.split(",")] for _, q in Q]
    where, cen = {}, []
    if iw >= 0:
        wt = rest[iw + len("where"):]
        ic = wt.find("(")
        w_txt, c_txt = (wt[:ic], wt[ic:]) if ic >= 0 else (wt, "")
        for item in w_txt.split(","):
            if item.strip():
                k, v = item.split("=")
                where[k.strip()] = v.strip()
        for c in re.findall(r"\(([^()]*)\)", c_txt):
            comps = c.split(",")
            if len(comps) != n:
                raise BuildError(f"{cls}: centring {c}")
            cen.append([Fr(x) for x in comps])
        if re.sub(r"[\s,]|\([^()]*\)", "", c_txt):
            raise BuildError(f"{cls}: trailing text {c_txt!r}")
    return dict(cls=cls, sym=sym, d=d, P=P, Q=Q, where=where, cen=cen)


def label_q(sym, d):
    """Basic-setting q's from the label symbol, components as written."""
    qs = re.findall(r"\(([^()]*)\)", sym)
    if len(qs) != d or any(len(q.split(",")) != 3 for q in qs):
        raise BuildError(f"cannot read {d} q's from {sym!r}")
    qs = [q.replace(" ", "") for q in qs]
    params = []
    for q in qs:
        for v in re.findall(r"[a-z]+\d*", q):
            if v not in params:
                params.append(v)
    return qs, params


def load_zip(path):
    with zipfile.ZipFile(path) as z:
        files = {os.path.basename(i.filename): z.read(i).decode("ascii")
                 for i in z.infolist() if not i.is_dir()}
    groups, lats = [], {}
    for d in (1, 2, 3):
        for b in split_blocks(files[f"ssg_lattice{d}_tables.txt"]):
            L = parse_lattice_block(b, d)
            if L["cls"] in lats:
                raise BuildError(f"duplicate class {L['cls']}")
            lats[L["cls"]] = L
        for b in split_blocks(files[f"ssg_group{d}_tables.txt"]):
            groups.append(parse_group_block(b, d))
    return groups, lats


# ============================================================ integer group machinery
_MM, _ROWS = {}, {}


def matmul_flat(n, A, B):
    """Product of flat n x n integer matrices, memoised (point groups are tiny)."""
    key = (A, B)
    C = _MM.get(key)
    if C is None:
        if len(_MM) > 200000:
            _MM.clear()
            _ROWS.clear()
        C = _MM[key] = tuple(sum(A[i * n + k] * B[k * n + j] for k in range(n))
                             for i in range(n) for j in range(n))
    return C


class IntOps:
    """Operators as (W flat int tuple, w int tuple mod D) for fast closure."""

    def __init__(self, n, D):
        self.n, self.D = n, D

    def enc(self, W, w):
        flat = []
        for row in W:
            for x in row:
                if Fr(x).denominator != 1:
                    raise BuildError(f"non-integer W entry {x}")
                flat.append(int(x))
        ww = []
        for x in w:
            y = mod1(Fr(x)) * self.D
            if y.denominator != 1:
                raise BuildError(f"translation {x} not in (1/{self.D})Z")
            ww.append(int(y))
        return tuple(flat), tuple(ww)

    def dec(self, op):
        W, w = op
        n = self.n
        return [list(W[i * n:(i + 1) * n]) for i in range(n)], [Fr(x, self.D) for x in w]

    def mul(self, a, b):
        n, D = self.n, self.D
        A, ta = a
        B, tb = b
        C = matmul_flat(n, A, B)
        rows = _ROWS.get(A)
        if rows is None:
            rows = _ROWS[A] = [A[i * n:(i + 1) * n] for i in range(n)]
        t = tuple((sum(x * y for x, y in zip(r, tb)) + c) % D for r, c in zip(rows, ta))
        return C, t

    def identity(self):
        n = self.n
        return tuple(int(i == j) for i in range(n) for j in range(n)), (0,) * n

    def closure(self, gens, limit=5000):
        e = self.identity()
        seen = {e}
        frontier = [e]
        while frontier:
            nxt = []
            for a in frontier:
                for g in gens:
                    c = self.mul(a, g)
                    if c not in seen:
                        seen.add(c)
                        nxt.append(c)
                        if len(seen) > limit:
                            raise BuildError("closure exceeds limit (not a finite group mod lattice?)")
            frontier = nxt
        return seen


def point_closure(n, gens):
    """Closure of integer matrices (flat tuples) under multiplication."""
    e = tuple(int(i == j) for i in range(n) for j in range(n))
    seen = {e}
    frontier = [e]
    while frontier:
        nxt = []
        for A in frontier:
            for B in gens:
                C = matmul_flat(n, A, B)
                if C not in seen:
                    seen.add(C)
                    nxt.append(C)
        frontier = nxt
    return seen


def lcm(*xs):
    l = 1
    for x in xs:
        l = l * x // math.gcd(l, x)
    return l


def rat_key(op):
    """CONTRACTS §1 sort key: (W flattened row-major, w), numeric comparison."""
    W, w = op
    return tuple(W), tuple(w)


# ============================================================ conversion
def convert(g, L, rng):
    """Convert one group to the basic setting and check all invariants.

    Returns (record dict for JSON, info dict).  Raises BuildError on any violation.
    """
    label, d = g["label"], g["d"]
    n = 3 + d
    sg, ld, cnum, gnum = (int(x) for x in LABEL_RE.fullmatch(label).groups())
    # --- label invariants
    if ld != d:
        raise BuildError(f"{label}: label d={ld} but in group{d} table")
    m = re.match(r"[A-Z][^(]*", g["sym"])
    if not m:
        raise BuildError(f"{label}: symbol {g['sym']!r}")
    cls = f"{ld}.{cnum}"
    if L["cls"] != cls:
        raise BuildError(f"{label}: class mismatch")

    P = L["P"]
    PT = transpose(P)
    PTi = inverse(PT)
    if any(x.denominator != 1 for r in P for x in r):
        raise BuildError(f"{label}: non-integer P")

    # --- supercentred operators and centrings
    ops_s = [parse_op(o, n, SUPER_VARS) for o in g["ops"]]
    cen_s = g["cen"]
    zero = [Fr(0)] * n
    if not any(c == zero for c in cen_s):
        raise BuildError(f"{label}: centring list lacks 0")
    # group centrings vs lattice-table centrings (supercentred, modulo 1)
    cset_g = {tuple(mod1(x) for x in c) for c in cen_s} - {tuple(zero)}
    cset_l = {tuple(mod1(x) for x in c) for c in L["cen"]}
    if cset_g != cset_l:
        raise BuildError(f"{label}: group centrings {cset_g} != lattice centrings {cset_l}")
    # supercentred q's of group table == lattice table Q's
    if g["qs"] != L["Q"]:
        raise BuildError(f"{label}: supercentred q's {g['qs']} != lattice {L['Q']}")

    # --- convert
    cen_b_all = [[mod1(x) for x in matvec(PT, c)] for c in cen_s]
    cen_b = sorted({tuple(c) for c in cen_b_all} - {tuple(zero)})
    conv = []
    for W, w in ops_s:
        Wb = matmul(matmul(PT, W), PTi)
        wb = [mod1(x) for x in matvec(PT, w)]
        conv.append((Wb, wb))

    # --- upper-right zero, integral, unimodular-ish
    for Wb, _ in conv:
        for i in range(n):
            for j in range(n):
                if Wb[i][j].denominator != 1:
                    raise BuildError(f"{label}: non-integer W_b {Wb}")
        if any(Wb[i][j] != 0 for i in range(3) for j in range(3, n)):
            raise BuildError(f"{label}: upper-right block non-zero {Wb}")

    D = lcm(*(x.denominator for _, w in conv for x in w), *(x.denominator for c in cen_b for x in c))
    I = IntOps(n, D)
    cen_i = [I.enc([[int(i == j) for j in range(n)] for i in range(n)], list(c)) for c in cen_b]
    full = set()
    for Wb, wb in conv:
        for c in [zero] + [list(c) for c in cen_b]:
            full.add(I.enc(Wb, [a + b for a, b in zip(wb, c)]))
    npg = len({W for W, _ in full})
    # BYU quirk: 5,748 groups list some coset representatives more than once, either
    # verbatim or differing by a centring vector (e.g. 24.1.12.2 lists X+1/2,-Y+1/2,-Z,-T
    # and X,-Y,-Z+1/2,-T, which differ by (1/2,1/2,1/2,0)).  These are harmless iff
    # every linear part's listed translations lie in a single centring coset, which the
    # |set| == |PG| x |cen| check below enforces.
    nredundant = len(ops_s) - npg
    nops = len(full)
    if nops != npg * (len(cen_b) + 1):
        raise BuildError(f"{label}: |set| {nops} != |PG| {npg} x |cen| {len(cen_b) + 1}")
    # (1) closure of all table ops + centrings has order nops
    tab_i = [I.enc(Wb, wb) for Wb, wb in conv]
    cl = I.closure(tab_i + cen_i, limit=4 * nops + 10)
    if cl != full:
        raise BuildError(f"{label}: closure order {len(cl)} != nops {nops}")

    # --- generators: greedy on the point group, largest subgroup first, ties by table order
    pg_order = []
    for op in tab_i:
        if op[0] not in pg_order:
            pg_order.append(op[0])
    e = I.identity()[0]
    chosen, H = [], {e}
    while len(H) < npg:
        best = None
        for W in pg_order:
            if W in H:
                continue
            Hn = point_closure(n, chosen + [W])
            if best is None or len(Hn) > len(best[1]):
                best = (W, Hn)
        chosen.append(best[0])
        H = best[1]
    # minimise: drop redundant gens (later gens can make earlier ones redundant)
    k = 0
    while k < len(chosen):
        trial = chosen[:k] + chosen[k + 1:]
        if len(point_closure(n, trial)) == npg:
            chosen = trial
        else:
            k += 1
    gens = []
    for W in chosen:
        reps = sorted((op for op in full if op[0] == W), key=lambda op: op[1])
        gens.append(reps[0])
    # (5) generators + centrings regenerate the full set
    if I.closure(gens + cen_i, limit=4 * nops + 10) != full:
        raise BuildError(f"{label}: generators do not regenerate the operator set")

    # --- label q's / params; (3) M + eps Q = Q R for generic params
    qs, params = label_q(g["sym"], d)
    for _trial in range(2):
        vals = {p: Fr(rng.choice(PRIMES), rng.choice(PRIMES[:40])) * rng.choice((1, -1)) for p in params}
        Q = [[lin_eval(c, vals) for c in q.split(",")] for q in qs]
        for W, _ in conv:
            R = [r[:3] for r in W[:3]]
            M = [r[:3] for r in W[3:]]
            E = [r[3:] for r in W[3:]]
            lhs = [[M[i][j] + sum(E[i][k] * Q[k][j] for k in range(d)) for j in range(3)] for i in range(d)]
            if lhs != matmul(Q, R):
                raise BuildError(f"{label}: M + eps Q != Q R for W={W}, params {vals}")
        # (extra) basic Q transforms into the lattice-table supercentred Q:
        #   Q_s = (P_IIᵀ - Q P_IEᵀ)^-1 (Q P_EEᵀ - P_EIᵀ)
        P_EE = [r[:3] for r in P[:3]]
        P_EI = [r[3:] for r in P[:3]]
        P_IE = [r[:3] for r in P[3:]]
        P_II = [r[3:] for r in P[3:]]
        left = [[a - b for a, b in zip(r1, r2)] for r1, r2 in zip(transpose(P_II), matmul(Q, transpose(P_IE)))]
        right = [[a - b for a, b in zip(r1, r2)] for r1, r2 in zip(matmul(Q, transpose(P_EE)), transpose(P_EI))]
        Qs = matmul(inverse(left), right)
        wvals = {}
        for k, v in L["where"].items():
            wvals[k] = lin_eval(v, vals)
        Qs_tab = [[lin_eval(c, wvals) if c not in ("0",) else Fr(0) for c in q] for q in L["Q"]]
        if Qs != Qs_tab:
            raise BuildError(f"{label}: Pᵀ-transformed Q {Qs} != table supercentred Q {Qs_tab}")

    rec = dict(
        sym=g["sym"], cls=cls, d=d, q=qs, params=params,
        cen=[",".join(fmt_rat(x) for x in c) for c in cen_b],
        gen=[fmt_op(*I.dec(op), n) for op in gens],
        nops=nops,
    )
    # round-trip check of the written strings
    gen_rt = [I.enc(*parse_op(s, n)) for s in rec["gen"]]
    cen_rt = [I.enc([[int(i == j) for j in range(n)] for i in range(n)], [Fr(x) for x in c.split(",")]) for c in rec["cen"]]
    if I.closure(gen_rt + cen_rt, limit=4 * nops + 10) != full:
        raise BuildError(f"{label}: written record does not round-trip")
    info = dict(sg=sg, npg=npg, nredundant=nredundant, ngen=len(gens), ncen=len(cen_b) + 1, full=full, I=I)
    return rec, info


def canonical_ops(rec):
    """CONTRACTS §1 canonical list from a record: closure of gen + cen, sorted, identity first."""
    n = 3 + rec["d"]
    gens = [parse_op(s, n) for s in rec["gen"]]
    cens = [[Fr(x) for x in c.split(",")] for c in rec["cen"]]
    D = lcm(*(x.denominator for _, w in gens for x in w), *(x.denominator for c in cens for x in c))
    I = IntOps(n, D)
    eye = [[int(i == j) for j in range(n)] for i in range(n)]
    cl = I.closure([I.enc(W, w) for W, w in gens] + [I.enc(eye, c) for c in cens])
    ops = [I.dec(op) for op in cl]
    ops.sort(key=lambda o: (tuple(x for r in o[0] for x in r), tuple(o[1])))
    e = (eye, [Fr(0)] * n)
    ops.remove(e)
    return [e] + ops


PRIMES = [p for p in range(1009, 20000) if all(p % q for q in range(2, int(p ** 0.5) + 1))][:400]


# ============================================================ driver
def build(zip_path=ZIP_DEFAULT, out_dir=OUT_DEFAULT, write=True, seed=20261110, verbose=True):
    t0 = time.time()
    sha = hashlib.sha256(open(zip_path, "rb").read()).hexdigest()
    groups, lats = load_zip(zip_path)
    t_parse = time.time() - t0
    rng = random.Random(seed)
    by_sg, errors = {}, []
    stats = dict(groups=0, classes=len(lats), by_d={1: 0, 2: 0, 3: 0}, max_nops=0, max_gen=0,
                 gen_hist={}, cen_hist={}, redundant_groups=0, redundant_ops=0)
    seen = set()
    for g in groups:
        if g["label"] in seen:
            errors.append(f"duplicate label {g['label']}")
            continue
        seen.add(g["label"])
        cls = ".".join(g["label"].split(".")[1:3])
        try:
            if cls not in lats:
                raise BuildError(f"{g['label']}: no lattice class {cls}")
            rec, info = convert(g, lats[cls], rng)
        except BuildError as ex:
            errors.append(str(ex))
            continue
        by_sg.setdefault(info["sg"], {})[g["label"]] = rec
        stats["groups"] += 1
        stats["by_d"][g["d"]] += 1
        stats["max_nops"] = max(stats["max_nops"], rec["nops"])
        stats["max_gen"] = max(stats["max_gen"], info["ngen"])
        stats["redundant_groups"] += info["nredundant"] > 0
        stats["redundant_ops"] += info["nredundant"]
        stats["gen_hist"][info["ngen"]] = stats["gen_hist"].get(info["ngen"], 0) + 1
        stats["cen_hist"][info["ncen"]] = stats["cen_hist"].get(info["ncen"], 0) + 1
    # SG coverage
    missing = [s for s in range(1, 231) if s not in by_sg]
    if missing:
        errors.append(f"space groups with no superspace group: {missing}")
    stats["time_parse_s"] = round(t_parse, 1)
    stats["time_total_s"] = round(time.time() - t0, 1)
    if write and not errors:
        if os.path.commonpath([os.path.abspath(out_dir), LOCAL]) != LOCAL:
            raise BuildError(f"refusing to write BYU-derived data outside {LOCAL}")
        os.makedirs(out_dir, exist_ok=True)
        for sg in sorted(by_sg):
            with open(os.path.join(out_dir, f"sg{sg:03d}.json"), "w") as f:
                json.dump({"sg": sg, "groups": by_sg[sg]}, f, separators=(",", ":"))
                f.write("\n")
        index = {"source": "BYU ISO(3+d)D ssg_tables.zip", "sha256": sha,
                 "count": stats["groups"],
                 "labels": [l for sg in sorted(by_sg) for l in by_sg[sg]]}
        with open(os.path.join(out_dir, "index.json"), "w") as f:
            json.dump(index, f, separators=(",", ":"))
            f.write("\n")
    if verbose:
        print(f"zip sha256 {sha}")
        print(f"classes {stats['classes']}  groups {stats['groups']} {stats['by_d']}  "
              f"max nops {stats['max_nops']}  max #gen {stats['max_gen']}")
        print(f"#gen histogram {dict(sorted(stats['gen_hist'].items()))}  "
              f"#centring (incl. 0) histogram {dict(sorted(stats['cen_hist'].items()))}")
        print(f"table quirk: {stats['redundant_groups']} groups list {stats['redundant_ops']} redundant "
              f"coset representatives (verbatim or differing by a centring); all consistent")
        print(f"invariants checked: closure order == nops, upper-right blocks zero, integral W_b, "
              f"M+eQ=QR (2 generic param draws), label SG/d/class, generators regenerate set, "
              f"Pᵀ-transformed Q == table Q_s, group/lattice centrings agree, JSON round-trip")
        print(f"time: parse {stats['time_parse_s']} s, total {stats['time_total_s']} s")
        print(f"VIOLATIONS: {len(errors)}")
        for e in errors[:50]:
            print("  ", e)
    return by_sg, stats, errors


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", default=ZIP_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args()
    _, stats, errors = build(a.zip, a.out, write=not a.no_write)
    if errors:
        print(f"FAILED: {len(errors)} invariant violations; nothing written", file=sys.stderr)
        sys.exit(1)
    if stats["groups"] != 16697 or stats["classes"] != 322:
        print(f"FAILED: expected 16697 groups / 322 classes", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
