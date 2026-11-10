#!/usr/bin/env python3
"""Write golden fixtures (CONTRACTS §3) to local/fixtures/golden/.

Cases:
  * the 19 IT-finder rows of validate_it_finder.py;
  * the In8IrSe9 (3+3)D cases (compare.IN8_CASES);
  * if local/data/ exists: ~300 random (group, q) pairs stratified over Bravais classes,
    each q component k/N with one N in 2..8 per case.  For these the grid prototype is run
    too (T4) and any disagreement aborts.

Usage:  python3 -I make_golden.py [--random 300] [--seed 1] [--bruteforce] [--no-random]
"""
import argparse, glob, json, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import compare
import exact

LOCAL = os.path.join(HERE, "..", "local")
OUT = os.path.join(LOCAL, "fixtures", "golden")
DATA = os.path.join(LOCAL, "data")


def fixture_name(label, params):
    p = ",".join(f"{k}={str(v).replace('/', '_')}" for k, v in params.items())
    return f"{label}__{p}.json"


def write(ex, label, out=OUT):
    os.makedirs(out, exist_ok=True)
    fn = os.path.join(out, fixture_name(label, {k: exact.fstr(v) for k, v in ex.params.items()}))
    with open(fn, "w") as f:
        json.dump(ex.to_json(label), f, indent=1, ensure_ascii=False)
        f.write("\n")
    return fn


def load_data():
    recs = {}
    for fn in sorted(glob.glob(os.path.join(DATA, "sg*.json"))):
        j = json.load(open(fn))
        for label, rec in j["groups"].items():
            rec = dict(rec, label=label)
            recs[label] = rec
    return recs


def random_cases(recs, n, seed):
    rng = random.Random(seed)
    by_cls = {}
    for label, rec in recs.items():
        by_cls.setdefault(rec.get("cls", "?"), []).append(label)
    classes = sorted(by_cls)
    rng.shuffle(classes)
    out = []
    while len(out) < n and classes:
        for c in classes:
            if len(out) >= n:
                break
            label = rng.choice(sorted(by_cls[c]))
            rec = recs[label]
            N = rng.randint(2, 8)
            params = {p: exact.fstr(exact.Fr(rng.randint(1, N - 1), N)) for p in rec["params"]}
            out.append((label, params))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--random", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bruteforce", action="store_true")
    ap.add_argument("--no-random", action="store_true")
    a = ap.parse_args()
    summary = []
    cases, _ = compare.it_cases()
    fixed = [(lab, {"g": g}) for lab, _, g, _ in cases] + compare.IN8_CASES
    jobs = [(lab, p, None) for lab, p in fixed]
    if os.path.isdir(DATA) and not a.no_random:
        recs = load_data()
        jobs += [(lab, p, recs[lab]) for lab, p in random_cases(recs, a.random, a.seed)]
    elif not a.no_random:
        print("local/data/ not found: skipping the random cases")
    failures = 0
    for lab, p, rec in jobs:
        t = time.time()
        try:
            ex = exact.run(rec if rec is not None else lab, p)
        except Exception as e:  # report and continue
            print(f"ERROR {lab} {p}: {e!r}")
            failures += 1
            continue
        fn = write(ex, lab)
        row = dict(label=lab, params=p, d=ex.d, classes=len(ex.classes),
                   orbits=len({c["orbit"] for c in ex.classes}), file=os.path.basename(fn))
        if rec is not None or a.bruteforce:
            r = compare.grid_compare(ex)
            row.update(grid_classes=r["grid_classes"], missed=r["missed"],
                       shared_keys=len(r["shared_keys"]), problems=r["problems"])
            if r["problems"]:
                failures += 1
        if a.bruteforce:
            row["bruteforce"] = compare.bruteforce(ex)
            failures += bool(row["bruteforce"])
        row["seconds"] = round(time.time() - t, 2)
        summary.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
    with open(os.path.join(OUT, "_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
    print(f"{len(summary)} fixtures, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
