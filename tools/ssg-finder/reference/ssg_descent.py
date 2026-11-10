#!/usr/bin/env python3
"""
ssg_descent.py -- 3D space groups of commensurate sections of a (3+d)D superspace group.

Given a superspace group from the BYU ISO(3+d)D tables (operators in the BASIC
space-group setting), rational values of its modulation-vector parameters, and an
optional axis permutation, enumerate the global phase t0 = (t1..td) on a fine
rational grid, and for every distinct set of surviving superspace operations
report the 3D space group of the commensurate supercell.

Theory (conventions of the BYU/Jana tables, superspace basis):
  superspace op  g = (W, w),  W = [[R, 0], [M, eps]],  w = (tau_E, tau_I)
     x'   = R x + tau_E
     x_I' = M x + eps x_I + tau_I
  3D crystal at phase t0 = section  x_I = Q x + t0   (rows of Q = q_j, x_s4 = t + q.x).
  Consistency (exact, as polynomials in the q parameters):  M + eps Q = Q R.
  g maps the t0-section onto itself, combined with the basic translation n and
  internal lattice translation m, iff
        (eps - 1) t0 + tau_I - Q (tau_E + n)  in  Z^d          (*)
  and then (R | tau_E + n) is a symmetry operation of the 3D structure.
  Pure translations: (1 | n) is a 3D symmetry iff Q n in Z^d  -> supercell lattice.
  t0 is defined modulo Z^d + Q Z^3 (an origin shift by a basic lattice vector n
  changes t0 by -Q n).

Usage examples
  python3 ssg_descent.py 62.1.9.2 --q g=1/3
  python3 ssg_descent.py 55.3.122.96 --q a1=1/4 b1=1/4 g2=1/4 --axes zxy
"""
import argparse, itertools, math, os, re, sys, html, subprocess, json
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get("SSGF_CACHE", os.path.join(HERE, "..", "local", "cache"))

# ---------------------------------------------------------------- fetch/parse
def fetch(label):
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, label + ".html")
    if not os.path.exists(fn) or os.path.getsize(fn) < 500:
        subprocess.run(["curl", "-s", "-A", "Mozilla/5.0", "-X", "POST", "-d",
                        f"label={label}&notation=xyzt&supercenter=true",
                        "https://iso.byu.edu/ssginfoform.php", "-o", fn], check=True)
    t = open(fn).read()
    t = html.unescape(re.sub(r"<[^>]+>", "\n", t))
    return re.sub(r"\s+", " ", t)

VARS = "xyztuv"

def parse_lin(expr, nvar):
    """'z-u+1/2' -> (coeff vector over x,y,z,t,u,v[:nvar], constant)"""
    e = expr.replace(" ", "")
    if e[0] not in "+-":
        e = "+" + e
    coef = [Fr(0)] * nvar
    const = Fr(0)
    for sgn, body in re.findall(r"([+-])([^+-]+)", e):
        s = 1 if sgn == "+" else -1
        m = re.fullmatch(r"(\d+(?:/\d+)?)?\*?([xyztuv])", body)
        if m:
            c = Fr(m.group(1)) if m.group(1) else Fr(1)
            coef[VARS.index(m.group(2))] += s * c
        else:
            const += s * Fr(body)
    return coef, const

def parse_ssg(label):
    t = fetch(label)
    m = re.search(r"Superspace group: (\S+) +(\S+)", t)
    symbol = m.group(2) if m else "?"
    # restrict to the basic-space-group setting block
    i0 = t.find("Independent modulation vectors:")
    i1 = t.find("SUPERCENTERED SETTING")
    if i1 < 0:
        i1 = t.find("Reflection conditions:")
    blk = t[i0:i1]
    qs = re.findall(r"q\d+=\(([^)]*)\)", blk)
    d = len(qs)
    nv = 3 + d
    cen = re.search(r"Centering: (.*?) Non-lattice generators", blk).group(1)
    cen = [tuple(Fr(x) for x in c.split(",")) for c in re.findall(r"\(([^)]*)\)", cen)]
    opstr = re.search(r"Non-lattice operators: (.*)", blk).group(1)
    ops = []
    for o in re.findall(r"\(([^)]*)\)", opstr):
        parts = o.split(",")
        assert len(parts) == nv, (o, nv)
        W, w = [], []
        for p in parts:
            c, k = parse_lin(p, nv)
            W.append(c); w.append(k)
        ops.append((W, w))
    return dict(label=label, symbol=symbol, d=d, qsym=qs, centering=cen, ops=ops)

def eval_q(qsym, params):
    Q = []
    for q in qsym:
        row = []
        for comp in q.split(","):
            comp = comp.strip()
            expr = comp
            for k in sorted(params, key=len, reverse=True):
                expr = re.sub(rf"(?<![a-z0-9]){k}(?![0-9])", f"({params[k]})", expr)
            row.append(Fr(_fr_eval(expr)))
        Q.append(row)
    return Q

def _fr_eval(expr):
    # evaluate with Fractions: wrap every number
    expr2 = re.sub(r"(\d+)", r"Fr(\1)", expr)
    return eval(expr2, {"Fr": Fr})

# ---------------------------------------------------------------- linear algebra helpers
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]

def matvec(A, v):
    return [sum(a * b for a, b in zip(r, v)) for r in A]

def frac1(v):
    return tuple(x - math.floor(x) for x in v)

def lcm_den(xs):
    l = 1
    for x in xs:
        l = l * Fr(x).denominator // math.gcd(l, Fr(x).denominator)
    return l

# ---------------------------------------------------------------- core
class Descent:
    def __init__(self, ssg, Q):
        self.ssg, self.Q, self.d = ssg, Q, ssg["d"]
        d = self.d
        # split ops; include centering translations as extra ops
        self.ops = []
        for W, w in ssg["ops"]:
            for c in ssg["centering"]:
                ww = [a + b for a, b in zip(w, c)]
                R = [row[:3] for row in W[:3]]
                M = [row[:3] for row in W[3:]]
                E = [row[3:] for row in W[3:]]
                assert all(row[3:] == [0] * d for row in W[:3])
                # consistency  M + E Q == Q R
                lhs = [[M[i][j] + sum(E[i][k] * Q[k][j] for k in range(d)) for j in range(3)]
                       for i in range(d)]
                rhs = matmul(Q, R)
                if lhs != rhs:
                    raise ValueError(f"inconsistent op/q: {W} {lhs} {rhs}")
                self.ops.append(dict(R=R, M=M, E=E, tE=ww[:3], tI=ww[3:]))
        self.D = lcm_den([x for r in Q for x in r])
        # supercell lattice: n in [0,D)^3 with Q n integer
        # conventional supercell N1 x N2 x N3 (setting axes) and its centring vectors
        self.N = [lcm_den([Q[i][j] for i in range(self.d)]) for j in range(3)]
        self.cells = [n for n in itertools.product(*[range(k) for k in self.N])
                      if all(x.denominator == 1 for x in matvec(Q, n))]

    def offsets(self, op):
        """set of  tau_I - Q(tau_E + n) mod 1  over basic translations n"""
        out = {}
        for n in itertools.product(range(self.D), repeat=3):
            v = frac1([a - b for a, b in zip(op["tI"],
                       matvec(self.Q, [x + y for x, y in zip(op["tE"], n)]))])
            out.setdefault(v, n)
        return out

    def survive(self, op, t0, offs=None):
        offs = offs or self.offsets(op)
        v = frac1([-(sum(op["E"][i][k] * t0[k] for k in range(self.d)) - t0[i])
                   for i in range(self.d)])
        return offs.get(v)  # n or None

    def scan(self, G):
        """grid t0 in (1/G)Z^d mod 1; vectorised exact integer arithmetic"""
        d = self.d
        offs = [self.offsets(op) for op in self.ops]
        K = G * lcm_den([x for o in offs for v in o for x in v] + [Fr(1, G)])
        grid = np.array(list(itertools.product(range(G), repeat=d)), dtype=np.int64) * (K // G)
        # keys[:, j] = index of the matched offset (i.e. of the coset n mod supercell) or -1
        keys = np.full((len(grid), len(self.ops)), -1, dtype=np.int32)
        for j, op in enumerate(self.ops):
            E = np.array([[int(x) for x in r] for r in op["E"]], dtype=np.int64)
            v = (-(grid @ E.T - grid)) % K
            code = np.zeros(len(grid), dtype=np.int64)
            for i in range(d):
                code = code * K + v[:, i]
            allowed = np.array([sum(int(x * K) * K ** (d - 1 - i) for i, x in enumerate(o))
                                for o in offs[j]], dtype=np.int64)
            order = np.argsort(allowed)
            pos = np.clip(np.searchsorted(allowed[order], code), 0, len(allowed) - 1)
            hit = allowed[order][pos] == code
            keys[hit, j] = order[pos[hit]]
        uniq, inv = np.unique(keys, axis=0, return_inverse=True)
        inv = np.asarray(inv).ravel()
        res = []
        for u in range(len(uniq)):
            idxs = np.nonzero(inv == u)[0]
            mask = uniq[u] >= 0
            # 'nice' representative: fewest non-zero components, then smallest
            best = min(idxs, key=lambda i: (int((grid[i] > 0).sum()), tuple(grid[i])))
            t0 = tuple(Fr(int(x), K) for x in grid[best])
            surv = [j for j in range(len(self.ops)) if mask[j]]
            rank = np.linalg.matrix_rank(np.vstack(
                [np.array([[float(x) for x in r] for r in self.ops[j]["E"]]) - np.eye(d)
                 for j in surv])) if surv else 0
            res.append(dict(t0=t0, surv=surv, npts=len(idxs), dim=d - rank, mask=mask))
        return res

    def group3d(self, t0, axes=None):
        """3D ops (R, tau) in the D x D x D supercell basis, all coset reps incl. centring"""
        ops3 = []
        for op in self.ops:
            n = self.survive(op, t0)
            if n is None:
                continue
            tau = [Fr(a + b, 1) for a, b in zip(op["tE"], n)]
            for c in self.cells:
                ops3.append((op["R"], tuple(((x + y) / k) % 1 for x, y, k in zip(tau, c, self.N))))
        # dedupe
        seen, out = set(), []
        for R, t in ops3:
            k = (tuple(map(tuple, R)), t)
            if k not in seen:
                seen.add(k); out.append((R, t))
        # rotation part in the supercell basis: diag(N)^-1 R diag(N)
        out = [([[R[i][j] * Fr(self.N[j], self.N[i]) for j in range(3)] for i in range(3)], t)
               for R, t in out]
        if axes:  # permutation: our coord i = setting coord axes[i]
            P = [[1 if axes[i] == "xyz"[j] else 0 for j in range(3)] for i in range(3)]
            Pi = [list(r) for r in zip(*P)]
            out = [(matmul(matmul(P, R), Pi), tuple(matvec(P, t))) for R, t in out]
        return out

# ---------------------------------------------------------------- identification
def identify(ops3, cellpar=(1, 1, 1)):
    import spglib
    rots = np.array([[[int(x) for x in r] for r in R] for R, t in ops3], dtype="intc")
    trs = np.array([[float(x) for x in t] for R, t in ops3])
    lat = np.diag(cellpar).astype(float)
    sgt = spglib.get_spacegroup_type_from_symmetry(rots, trs, lattice=lat, symprec=1e-5)
    return sgt

def hm_in_setting(ops3):
    """Hermann-Mauguin-like symbol in the GIVEN axes of the supercell (orthogonal axes only)."""
    E = np.eye(3, dtype=int)
    trans = [np.array([float(x) for x in t]) for R, t in ops3 if np.array_equal(np.array(R, dtype=int), E)]
    cen = ""
    tset = {tuple(np.round(t % 1, 6)) for t in trans}
    lat = []
    for name, v in (("A", (0, .5, .5)), ("B", (.5, 0, .5)), ("C", (.5, .5, 0)), ("I", (.5, .5, .5))):
        if v in tset:
            lat.append(name)
    cen = "F" if {"A", "B", "C"} <= set(lat) else (lat[0] if len(lat) == 1 else ("".join(lat) or "P"))
    if len(tset) > {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4}.get(cen, 0):
        cen += f"[{len(tset)}pts]"
    pts = {}
    for R, t in ops3:
        Ri = np.array([[int(x) for x in r] for r in R])
        pts.setdefault(tuple(Ri.flatten()), []).append(np.array([float(x) for x in t]))
    sym = []
    for ax in range(3):
        rot2 = -np.eye(3, dtype=int); rot2[ax, ax] = 1
        mir = np.eye(3, dtype=int); mir[ax, ax] = -1
        r = ""
        if tuple(rot2.flatten()) in pts:
            scr = [round(t[ax] % 1, 6) for t in pts[tuple(rot2.flatten())]]
            r = "2" if 0.0 in scr else "2_1" if 0.5 in scr else "2(?)"
        m = ""
        if tuple(mir.flatten()) in pts:
            gl = set()
            for t in pts[tuple(mir.flatten())]:
                g = tuple(round(x % 1, 6) if i != ax else 0 for i, x in enumerate(t))
                gl.add(g)
            names = []
            for g in gl:
                nz = [i for i in range(3) if g[i] != 0]
                if not nz: names.append("m")
                elif all(g[i] == .5 for i in nz):
                    names.append("abc"[nz[0]] if len(nz) == 1 else "n")
                elif all(g[i] == .25 or g[i] == .75 for i in nz): names.append("d")
                else: names.append("g?")
            order = "mebcandg?"
            nm = set(names)
            if {"a", "b", "c"} & nm and len({"a", "b", "c"} & nm) >= 2:
                m = "e"
            else:
                m = min(names, key=lambda s: order.index(s[0]) if s[0] in order else 99)
        sym.append(f"{r}/{m}" if r and m else (r or m or "1"))
    return cen + " " + " ".join(sym)

# ---------------------------------------------------------------- driver
def run(label, params, axes=None, G=None, cellpar=(1, 1, 1), verbose=False):
    ssg = parse_ssg(label)
    Q = eval_q(ssg["qsym"], params)
    ds = Descent(ssg, Q)
    if G is None:
        G = 8 * ds.D if ssg["d"] == 3 else 16 * ds.D
    res = ds.scan(G)
    agg = {}
    for r in res:
        ops3 = ds.group3d(r["t0"], axes)
        sgt = identify(ops3, cellpar)
        key = (tuple(r["surv"]), sgt.number, hm_in_setting(ops3))
        a = agg.setdefault(key, dict(t0s=[], npts=0, dim=r["dim"], surv=r["surv"],
                                     sg_number=sgt.number, sg_symbol=sgt.international_short,
                                     setting=key[2], nops3=len(ops3), ops3=ops3))
        a["t0s"].append(r["t0"]); a["npts"] += r["npts"]
    out = []
    for a in agg.values():
        a["t0s"].sort(key=lambda t: (sum(1 for x in t if x), t))
        a["t0"] = a["t0s"][0]
        out.append(a)
    out.sort(key=lambda r: (-r["dim"], -r["npts"], r["sg_number"]))
    return ssg, Q, ds, G, out

def fmt_t(t0):
    return "(" + ",".join(str(x) for x in t0) + ")"

def report(label, params, axes=None, G=None, cellpar=(1, 1, 1)):
    ssg, Q, ds, G, out = run(label, params, axes, G, cellpar)
    qtxt = ", ".join("(" + ",".join(str(x) for x in r) + ")" for r in Q)
    print(f"## {label}  {ssg['symbol']}   q = {qtxt}  (setting axes)  "
          f"supercell {'x'.join(map(str, ds.N))}/{len(ds.cells)} (index {math.prod(ds.N) // len(ds.cells)}), grid 1/{G}")
    print(f"{'t0 (sample)':>24} {'dim':>3} {'#grid':>7} {'#ops':>4}  {'type':>12}  setting(axes {axes or 'as SSG'})")
    for r in out:
        print(f"{fmt_t(r['t0']):>24} {r['dim']:>3} {r['npts']:>7} {len(r['surv']):>4}  "
              f"{r['sg_symbol']+' ('+str(r['sg_number'])+')':>12}  {r['setting']}")
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("--q", nargs="+", default=[], help="param=value, e.g. g=1/3")
    ap.add_argument("--axes", default=None, help="output-axis permutation, e.g. zxy")
    ap.add_argument("--grid", type=int, default=None)
    a = ap.parse_args()
    params = dict(p.split("=") for p in a.q)
    report(a.label, params, a.axes, a.grid)
