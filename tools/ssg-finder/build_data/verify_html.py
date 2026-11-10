#!/usr/bin/env python3
"""
verify_html.py -- check the converted BYU data (local/data/sg*.json, written by
byu_tables.py) against BYU's own basic-setting pages from ssginfoform.php.

One group per Bravais class (the lowest-numbered group of each of the 322 classes),
plus every page already in local/cache/.  Pages are fetched politely: a generic
User-Agent, >= 3 s between requests, and each page at most once (a page already in
the cache is never fetched again, even if it looks broken).

The comparison is exact set equality of the basic-setting operators, including all
centring combinations, modulo Z^(3+d); the q strings, symbol and centring vectors are
compared too.  HTML parsing is reference/ssg_descent.py::parse_ssg (imported as is).

Licence (decision D1): cached pages stay in local/cache/ (gitignored).

Usage:
    python3 -I build_data/verify_html.py --fetch     # fetch missing pages (~20 min), then compare
    python3 -I build_data/verify_html.py             # compare cached pages only
"""
import argparse, glob, json, os, sys, time, urllib.parse, urllib.request
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
LOCAL = os.path.join(ROOT, "local")
DATA = os.path.join(LOCAL, "data")
CACHE = os.path.join(LOCAL, "cache")
URL = "https://iso.byu.edu/ssginfoform.php"
USER_AGENT = "Mozilla/5.0"
DELAY_S = 3.5
MIN_PAGE = 500  # parse_ssg re-fetches pages smaller than this; never let it

sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "reference"))
os.environ.setdefault("SSGF_CACHE", CACHE)
import byu_tables as B  # noqa: E402


def load_records():
    recs = {}
    for fn in sorted(glob.glob(os.path.join(DATA, "sg*.json"))):
        recs.update(json.load(open(fn))["groups"])
    if not recs:
        sys.exit(f"no data in {DATA}; run byu_tables.py first")
    return recs


def class_representatives(recs):
    """Lowest-numbered group label of every Bravais class."""
    best = {}
    for lab, r in recs.items():
        key = tuple(int(x) for x in lab.split("."))
        if r["cls"] not in best or key < best[r["cls"]][0]:
            best[r["cls"]] = (key, lab)
    return [lab for _, lab in sorted(best.values())]


def cache_path(label):
    return os.path.join(CACHE, label + ".html")


def fetch_missing(labels, log=print):
    os.makedirs(CACHE, exist_ok=True)
    todo = [l for l in labels if not os.path.exists(cache_path(l))]
    log(f"fetch: {len(labels) - len(todo)} cached, {len(todo)} to fetch "
        f"(~{len(todo) * (DELAY_S + 1) / 60:.0f} min)")
    last = 0.0
    for i, lab in enumerate(todo):
        wait = DELAY_S - (time.time() - last)
        if wait > 0:
            time.sleep(wait)
        last = time.time()
        body = urllib.parse.urlencode({"label": lab, "notation": "xyzt", "supercenter": "true"}).encode()
        req = urllib.request.Request(URL, data=body, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                page = resp.read()
        except Exception as ex:  # record and move on: never hammer the server with retries
            log(f"  [{i + 1}/{len(todo)}] {lab}: FETCH ERROR {ex}")
            continue
        tmp = cache_path(lab) + ".part"
        with open(tmp, "wb") as f:
            f.write(page)
        os.replace(tmp, cache_path(lab))
        log(f"  [{i + 1}/{len(todo)}] {lab}: {len(page)} bytes")


def html_ops(ssg):
    out = set()
    for c in ssg["centering"]:
        for W, w in ssg["ops"]:
            out.add((tuple(tuple(Fr(x) for x in r) for r in W),
                     tuple(B.mod1(Fr(a) + b) for a, b in zip(w, c))))
    return out


def record_ops(rec):
    return {(tuple(tuple(Fr(x) for x in r) for r in W), tuple(w)) for W, w in B.canonical_ops(rec)}


def compare(label, rec):
    """Return a list of differences (empty = exact match)."""
    import ssg_descent as S
    fn = cache_path(label)
    if not os.path.exists(fn):
        return ["not cached"]
    if os.path.getsize(fn) < MIN_PAGE:
        return [f"cached page too small ({os.path.getsize(fn)} bytes)"]
    try:
        h = S.parse_ssg(label)
    except Exception as ex:
        return [f"HTML parse error: {ex!r}"]
    diffs = []
    if h["symbol"] != rec["sym"]:
        diffs.append(f"symbol {h['symbol']!r} != {rec['sym']!r}")
    hq = [q.replace(" ", "") for q in h["qsym"]]
    if hq != rec["q"]:
        diffs.append(f"q {hq} != {rec['q']}")
    hc = {tuple(B.mod1(x) for x in c) for c in h["centering"]} - {(Fr(0),) * (3 + rec["d"])}
    rc = {tuple(Fr(x) for x in c.split(",")) for c in rec["cen"]}
    if hc != rc:
        diffs.append(f"centring {sorted(hc)} != {sorted(rc)}")
    a, b = html_ops(h), record_ops(rec)
    if a != b:
        n = 3 + rec["d"]
        only_h = sorted(B.fmt_op(W, w, n) for W, w in a - b)[:3]
        only_r = sorted(B.fmt_op(W, w, n) for W, w in b - a)[:3]
        diffs.append(f"ops: html {len(a)} vs data {len(b)}; only html {only_h}; only data {only_r}")
    if len(b) != rec["nops"]:
        diffs.append(f"nops {rec['nops']} != closure {len(b)}")
    return diffs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true", help="fetch missing class representatives first")
    a = ap.parse_args()
    t0 = time.time()
    recs = load_records()
    reps = class_representatives(recs)
    if a.fetch:
        fetch_missing(reps, log=lambda s: print(s, flush=True))
    cached = sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(CACHE, "*.html")))
    labels = sorted(set(reps) | set(cached), key=lambda l: tuple(int(x) for x in l.split(".")))
    ok, bad, missing, classes_ok = 0, [], [], set()
    for lab in labels:
        if lab not in recs:
            bad.append((lab, ["label not in converted data"]))
            continue
        if not os.path.exists(cache_path(lab)):
            missing.append(lab)
            continue
        diffs = compare(lab, recs[lab])
        if diffs:
            bad.append((lab, diffs))
        else:
            ok += 1
            classes_ok.add(recs[lab]["cls"])
    print(f"classes: {len(reps)}; class representatives cached {len(reps) - len(missing)}; "
          f"classes with a matching page {len(classes_ok)}")
    print(f"pages compared {ok + len(bad)}: MATCH {ok}, MISMATCH {len(bad)}; not cached {len(missing)}")
    for lab, diffs in bad[:40]:
        print(f"  {lab}: " + "; ".join(diffs))
    print(f"time {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
