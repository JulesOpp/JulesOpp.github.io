#!/usr/bin/env python3
"""make_subgroup_fixtures.py -- write local/fixtures/subgroups/<label>.json (CONTRACTS §5).

One file per (3+1)D BYU record: its k-index-1 subgroup classes, each identified with the BYU
type and the map onto it.  BYU-derived, so local only (licence rule D1).

Usage:  python3 make_subgroup_fixtures.py [labels ...] [--jobs N] [--resume]

Uses at most 4 worker processes by default (override with --jobs or SSG_JOBS).  --resume keeps
every fixture that is already complete and only computes the missing ones.
"""
import argparse, json, multiprocessing, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # also under python -I
import exact  # noqa: E402
import ssg_equiv as se  # noqa: E402
import subgroups  # noqa: E402

MAX_JOBS = int(os.environ.get("SSG_JOBS", "4"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "local", "fixtures", "subgroups")


def fixture(label):
    db = se.load_d1()
    rec = db[label]
    classes = []
    for c in subgroups.t_subgroups(rec, db):
        j = dict(pg_ops=c["pg_ops"], index=c["index"], class_size=c["class_size"])
        if len(c["matches"]) == 1:
            mp = c["map"]
            assert all(x.denominator == 1 for x in mp["m"]), (label, c["pg_ops"], mp["m"])
            j.update(byu=c["byu"], byu_sym=c["byu_sym"],
                     P=[[exact.fstr(x) for x in r] for r in mp["P"]],
                     p=[exact.fstr(x) for x in mp["p"]], eta=mp["eta"],
                     m=[int(x) for x in mp["m"]], p4=exact.fstr(mp["p4"]))
        else:                                  # not expected; kept so the run log shows it
            j.update(byu=None, matches=c["matches"], sg=c["sg"])
        classes.append(j)
    order = len({W for W, w in exact.build_ops(rec)})
    return dict(label=label, sym=rec["sym"], order=order, classes=classes)


def work(label):
    t = time.time()
    fx = fixture(label)
    head = json.dumps({k: v for k, v in fx.items() if k != "classes"})[:-1]
    body = ",\n  ".join(json.dumps(c) for c in fx["classes"])
    path = os.path.join(OUT, f"{label}.json")
    with open(path + ".tmp", "w") as f:
        f.write(f'{head}, "classes": [\n  {body}\n]}}\n')
    os.replace(path + ".tmp", path)                    # never leave a truncated fixture behind
    bad = [c for c in fx["classes"] if c["byu"] is None]
    return label, len(fx["classes"]), bad, time.time() - t


def existing(label):
    """The fixture already on disk for `label` if it is complete and readable, else None."""
    try:
        with open(os.path.join(OUT, f"{label}.json")) as f:
            fx = json.load(f)
        return fx if fx.get("label") == label and isinstance(fx.get("classes"), list) else None
    except (OSError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("labels", nargs="*")
    ap.add_argument("--jobs", type=int, default=MAX_JOBS)
    ap.add_argument("--resume", action="store_true", help="keep complete fixtures, compute the rest")
    a = ap.parse_args()
    labels = a.labels or list(se.load_d1())
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    log = []
    if a.resume:
        todo = []
        for label in labels:
            fx = existing(label)
            if fx is None:
                todo.append(label)
            else:
                bad = [c for c in fx["classes"] if c["byu"] is None]
                log.append(dict(label=label, classes=len(fx["classes"]), unidentified=bad, seconds=None))
        print(f"resume: {len(log)} complete fixtures kept, {len(todo)} to compute", flush=True)
        labels = todo
    with multiprocessing.Pool(max(1, a.jobs)) as pool:
        for label, n, bad, dt in pool.imap_unordered(work, labels):
            log.append(dict(label=label, classes=n, unidentified=bad, seconds=round(dt, 2)))
            print(f"{label:14s} {n:3d} classes {dt:6.1f}s" + (f"  UNIDENTIFIED {bad}" if bad else ""), flush=True)
    log.sort(key=lambda r: tuple(int(x) for x in r["label"].split(".")))
    summary = dict(groups=len(log), classes=sum(r["classes"] for r in log),
                   unidentified=[r for r in log if r["unidentified"]], seconds=round(time.time() - t0, 1),
                   jobs=a.jobs, runs=log)
    with open(os.path.join(OUT, "_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"{summary['groups']} groups, {summary['classes']} classes, "
          f"{len(summary['unidentified'])} with unidentified classes, {summary['seconds']}s")


if __name__ == "__main__":
    main()
