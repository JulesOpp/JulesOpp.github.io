#!/usr/bin/env python3
"""make_subgroup_fixtures.py -- write local/fixtures/subgroups/<label>.json (CONTRACTS §5).

One file per (3+1)D BYU record: its k-index-1 subgroup classes, each identified with the BYU
type and the map onto it.  BYU-derived, so local only (licence rule D1).

Usage:  python3 make_subgroup_fixtures.py [labels ...] [--jobs N] [--resume]

Uses at most 4 worker processes (--jobs or SSG_JOBS can select fewer).  --resume keeps
every fixture that is already complete and only computes the missing ones.
"""
import argparse, functools, json, multiprocessing, os, sys, time


def worker_init():
    """Enforce the one-thread BLAS budget before numpy loads, including spawned workers."""
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"  # override inherited settings: the four-core budget is a hard cap


worker_init()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # also under python -I
import exact  # noqa: E402
import ssg_equiv as se  # noqa: E402
import subgroups  # noqa: E402


def cap_jobs(n: int) -> int:
    return max(1, min(n, 4))


MAX_JOBS = cap_jobs(int(os.environ.get("SSG_JOBS", "4")))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "local", "fixtures", "subgroups")


def fixture(label):
    db = se.load_records(1)
    rec = db[label]
    classes = []
    for c in subgroups.t_subgroups(rec, db):
        j = dict(pg_ops=c["pg_ops"], index=c["index"], class_size=c["class_size"])
        if len(c["matches"]) == 1:
            mp = c["map"]
            assert all(x.denominator == 1 for x in mp["M"][0]), (label, c["pg_ops"], mp["M"][0])
            j.update(byu=c["byu"], byu_sym=c["byu_sym"],
                     P=[[exact.fstr(x) for x in r] for r in mp["P"]],
                     p=[exact.fstr(x) for x in mp["p"]], eta=mp["E"][0][0],
                     m=[int(x) for x in mp["M"][0]], p4=exact.fstr(mp["p4"][0]))
        else:                                  # not expected; kept so the run log shows it
            j.update(byu=None, matches=c["matches"], sg=c["sg"])
        classes.append(j)
    order = len({W for W, w in exact.build_ops(rec)})
    return dict(label=label, sym=rec["sym"], order=order, classes=classes)


def work(label, out=OUT):
    t = time.time()
    fx = fixture(label)
    head = json.dumps({k: v for k, v in fx.items() if k != "classes"})[:-1]
    body = ",\n  ".join(json.dumps(c) for c in fx["classes"])
    path = os.path.join(out, f"{label}.json")
    with open(path + ".tmp", "w") as f:
        f.write(f'{head}, "classes": [\n  {body}\n]}}\n')
    os.replace(path + ".tmp", path)                    # never leave a truncated fixture behind
    bad = [c for c in fx["classes"] if c["byu"] is None]
    return label, len(fx["classes"]), bad, time.time() - t


def existing(label, out=OUT):
    """The fixture already on disk for `label` if it is complete and readable, else None."""
    try:
        with open(os.path.join(out, f"{label}.json")) as f:
            fx = json.load(f)
        return fx if fx.get("label") == label and isinstance(fx.get("classes"), list) else None
    except (OSError, ValueError):
        return None


def main(argv=None, out=OUT):
    ap = argparse.ArgumentParser()
    ap.add_argument("labels", nargs="*")
    ap.add_argument("--jobs", type=int, default=MAX_JOBS)
    ap.add_argument("--resume", action="store_true", help="keep complete fixtures, compute the rest")
    a = ap.parse_args(argv)
    a.jobs = cap_jobs(a.jobs)
    labels = a.labels or list(se.load_records(1))
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    log = []
    prior_seconds = 0
    if a.resume:
        try:
            with open(os.path.join(out, "_summary.json")) as f:
                prior_seconds = json.load(f).get("seconds", 0)
        except (OSError, ValueError):
            pass
        todo = []
        for label in labels:
            fx = existing(label, out)
            if fx is None:
                todo.append(label)
            else:
                bad = [c for c in fx["classes"] if c["byu"] is None]
                log.append(dict(label=label, classes=len(fx["classes"]), unidentified=bad, seconds=None))
        print(f"resume: {len(log)} complete fixtures kept, {len(todo)} to compute", flush=True)
        labels = todo
    with multiprocessing.Pool(a.jobs, initializer=worker_init) as pool:
        for label, n, bad, dt in pool.imap_unordered(functools.partial(work, out=out), labels):
            log.append(dict(label=label, classes=n, unidentified=bad, seconds=round(dt, 2)))
            print(f"{label:14s} {n:3d} classes {dt:6.1f}s" + (f"  UNIDENTIFIED {bad}" if bad else ""), flush=True)
    log.sort(key=lambda r: tuple(int(x) for x in r["label"].split(".")))
    wall_seconds = round(time.time() - t0, 1)
    summary = dict(groups=len(log), classes=sum(r["classes"] for r in log),
                   unidentified=[r for r in log if r["unidentified"]], seconds=round(prior_seconds + wall_seconds, 1),
                   wall_seconds_this_run=wall_seconds,
                   jobs=a.jobs, runs=log)
    with open(os.path.join(out, "_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(f"{summary['groups']} groups, {summary['classes']} classes, "
          f"{len(summary['unidentified'])} with unidentified classes, {summary['seconds']}s")


if __name__ == "__main__":
    main()
