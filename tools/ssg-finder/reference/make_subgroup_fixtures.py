#!/usr/bin/env python3
"""Generate source-bound, resumable local subgroup fixtures for d=1,2,3."""
import argparse, functools, hashlib, json, multiprocessing, os, sys, time
from collections import defaultdict
from pathlib import Path


def worker_init():
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"


worker_init()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact
import ssg_equiv as se
import subgroups


def cap_jobs(n):
    return max(1, min(n, 4))


MAX_JOBS = cap_jobs(int(os.environ.get("SSG_JOBS", "4")))
OUT = os.path.join(se.HERE, "..", "local", "fixtures", "subgroups")


def dimension(label):
    return int(label.split(".")[1])


def fixture(label):
    d = dimension(label)
    db = se.load_records(d)
    rec = db[label]
    ops = exact.build_ops(rec)
    q = se.q_generic(rec)
    classes = []
    for c in subgroups.t_subgroups(rec, db):
        j = dict(pg_ops=c["pg_ops"], index=c["index"], class_size=c["class_size"])
        if len(c["matches"]) == 1:
            mp = c["map"]
            sub = subgroups.subgroup_record(rec, ops, c["pg_ops"])
            if not se.check_map(sub, q, db[c["byu"]], mp):
                raise ValueError((label, c["pg_ops"], "invalid oracle map"))
            j.update(byu=c["byu"], byu_sym=c["byu_sym"],
                     P=[[exact.fstr(x) for x in row] for row in mp["P"]],
                     p=[exact.fstr(x) for x in mp["p"]],
                     E=[[int(x) for x in row] for row in mp["E"]],
                     M=[[exact.fstr(x) for x in row] for row in mp["M"]],
                     p4=[exact.fstr(x) for x in mp["p4"]])
        else:
            j.update(byu=None, matches=c["matches"], sg=c["sg"])
        classes.append(j)
    return dict(label=label, d=d, sym=rec["sym"], order=len({W for W,w in ops}), classes=classes)


@functools.lru_cache(maxsize=3)
def source_signature(d):
    digest = hashlib.sha256()
    for name in (__file__, subgroups.__file__, subgroups.affine_search.__file__, se.__file__, exact.__file__, se.im.__file__):
        digest.update(Path(name).read_bytes())
        digest.update(b"\0")
    digest.update(json.dumps(se.load_records(d), sort_keys=True).encode())
    return digest.hexdigest()


def valid_seconds(value):
    """Accept finite nonnegative JSON numeric timings, excluding booleans."""
    import math
    try:
        return type(value) in (int, float) and math.isfinite(value) and value >= 0
    except OverflowError:
        return False


def content_hash(fx):
    return hashlib.sha256(json.dumps({k:v for k,v in fx.items() if k != "_content_sha256"},
                                     sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=1) + "\n")
    os.replace(temp, path)


def work(label, out=OUT, signature=None):
    start = time.monotonic()
    fx = fixture(label)
    fx["_signature"] = signature or source_signature(dimension(label))
    fx["_seconds"] = time.monotonic() - start
    fx["_content_sha256"] = content_hash(fx)
    atomic_json(Path(out) / (label + ".json"), fx)
    bad = [c for c in fx["classes"] if c["byu"] is None]
    return label, len(fx["classes"]), bad, fx["_seconds"]


def valid_fixture(fx, label):
    """Check the on-disk schema in addition to its source/content bindings.

    Completeness includes both the trivial and full point subgroup, which every
    finite group has. The content checksum attests to the initial exact map
    checks; malformed freshly rewritten artifacts still fail this schema check.
    """
    from fractions import Fraction
    try:
        d = dimension(label)
        order = fx['order']
        classes = fx['classes']
        if (fx['label'] != label or fx['d'] != d or
                not isinstance(fx['sym'], str) or
                type(order) is not int or order < 1 or
                not isinstance(classes, list) or not classes):
            return False

        def vector(value, n):
            return (isinstance(value, list) and len(value) == n and
                    all(isinstance(x, str) and exact.fstr(Fraction(x)) == x for x in value))

        def matrix(value, n, m):
            return isinstance(value, list) and len(value) == n and all(vector(row, m) for row in value)

        trivial = full = False
        seen = set()
        for c in classes:
            idx, size, pg = c['index'], c['class_size'], c['pg_ops']
            if (type(idx) is not int or idx < 1 or order % idx or
                    type(size) is not int or size < 1 or order % size or
                    not isinstance(pg, list) or not pg or
                    any(type(i) is not int or i < 0 for i in pg) or
                    pg != sorted(set(pg)) or pg[0] != 0 or
                    len(pg) != order // idx or tuple(pg) in seen):
                return False
            seen.add(tuple(pg))
            if pg == [0] and idx == order and size == 1:
                trivial = True
            if idx == 1 and size == 1 and len(pg) == order:
                full = True
            if (not isinstance(c['byu'], str) or dimension(c['byu']) != d or
                    not isinstance(c['byu_sym'], str) or
                    not matrix(c['P'], 3, 3) or not vector(c['p'], 3) or
                    not matrix(c['M'], d, 3) or not vector(c['p4'], d)):
                return False
            E = c['E']
            if (not isinstance(E, list) or len(E) != d or
                    any(not isinstance(row, list) or len(row) != d or
                        any(type(x) is not int for x in row) for row in E) or
                    abs(se.im.det(E)) != 1 or
                    se.im.det([[Fraction(x) for x in row] for row in c['P']]) <= 0):
                return False
        return trivial and full
    except (KeyError, TypeError, ValueError, IndexError, ZeroDivisionError):
        return False


def existing(label, out=OUT, signature=None):
    """Reuse only a complete, untampered artifact bound to current source and corpus.

    Every map was independently check_map-validated before its atomic write;
    the checksum retains that validation on resume without repeating the search.
    """
    try:
        fx = json.loads((Path(out) / (label + ".json")).read_text())
        if not isinstance(fx, dict):
            return None
        expected = signature or source_signature(dimension(label))
        if (fx.get("label") != label or fx.get("_signature") != expected or
                fx.get("_content_sha256") != content_hash(fx) or
                not valid_fixture(fx, label) or
                not valid_seconds(fx.get("_seconds")) or
                any(not c.get("byu") for c in fx["classes"])):
            return None
        return fx
    except (OSError, ValueError, TypeError, KeyError):
        return None


def sample_labels(db):
    # Shared stratification: every Bravais class and every basic SG.
    from check_equiv_d23 import sample
    return list(sample(db))


def estimate_runtime(db, runs, jobs):
    """Descriptive stratified projection, not a statistical confidence bound.

    Each unmeasured label uses its SG/class intersection, then class, SG, or
    overall sample mean. Record the selected method and supporting labels.
    """
    measured = {r["label"]: r["seconds"] for r in runs if r.get("seconds") is not None}
    buckets = {"both": defaultdict(list), "cls": defaultdict(list), "sg": defaultdict(list)}
    for label in measured:
        sg, cls = label.split('.')[0], db[label]['cls']
        buckets['both'][(sg, cls)].append(label)
        buckets['cls'][cls].append(label)
        buckets['sg'][sg].append(label)
    means = {}
    for kind, groups in buckets.items():
        for key, support in groups.items():
            means[(kind, key)] = sum(measured[l] for l in support) / len(support)
    sample_mean = sum(measured.values()) / len(measured) if measured else None
    per_label = []
    for label, rec in db.items():
        if label in measured:
            method, support = "measured", [label]
        else:
            sg, cls = label.split('.')[0], rec['cls']
            choices = [("basic_sg_and_class_mean", "both", (sg,cls)),
                       ("class_mean", "cls", cls), ("basic_sg_mean", "sg", sg)]
            chosen = next(((m,k,v) for m,k,v in choices if v in buckets[k]), None)
            if chosen:
                method, kind, key = chosen
                support = buckets[kind][key]
            else:
                method, support = ("sample_mean", list(measured)) if measured else ("unmeasured", [])
        seconds = (measured[label] if method == 'measured' else
                   means[(kind,key)] if method not in ('sample_mean','unmeasured') else sample_mean)
        per_label.append(dict(label=label, seconds=seconds, projection_method=method, supporting_labels=support))
    complete = all(r['seconds'] is not None for r in per_label)
    total = sum(r['seconds'] for r in per_label) if complete else None
    return dict(method="SG/class stratified arithmetic means; wall estimate assumes balanced workers and excludes startup/stragglers",
                measured_labels=len(measured), total_labels=len(db), per_label=per_label,
                estimated_full_worker_seconds=total,
                estimated_full_wall_seconds=total / jobs if total is not None else None,
                jobs=jobs, full_d3_deferred=True)


def main(argv=None, out=OUT):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('labels', nargs='*')
    ap.add_argument('--d', type=int, choices=(1,2,3), default=1)
    ap.add_argument('--sample', action='store_true')
    ap.add_argument('--jobs', type=int, default=MAX_JOBS)
    ap.add_argument('--resume', action='store_true')
    a = ap.parse_args(argv)
    a.jobs = cap_jobs(a.jobs)
    if a.d == 3 and not a.sample and not a.labels:
        ap.error('full d3 sweeps are deferred; use --sample or explicit labels')
    db = se.load_records(a.d)
    labels = a.labels or (sample_labels(db) if a.sample else list(db))
    if any(dimension(l) != a.d or l not in db for l in labels):
        ap.error('all labels must belong to --d and the local corpus')
    # Explicit output is useful for synthetic tests; production always gets dN.
    out = Path(out) / f'd{a.d}' if str(out) == str(OUT) else Path(out)
    out.mkdir(parents=True, exist_ok=True)
    signature = source_signature(a.d)
    mono0 = time.monotonic()
    prior_seconds = 0
    if a.resume:
        try:
            prior = json.loads((out / '_summary.json').read_text())
            if (isinstance(prior, dict) and prior.get('signature') == signature
                    and valid_seconds(prior.get('seconds'))):
                prior_seconds = prior.get('seconds', 0)
        except (OSError, ValueError):
            pass
    log, todo = [], []
    for label in labels:
        fx = existing(label, out, signature) if a.resume else None
        if fx is None:
            todo.append(label)
        else:
            log.append(dict(label=label, classes=len(fx['classes']), unidentified=[], seconds=fx['_seconds']))
    print(json.dumps(dict(d=a.d, sample=a.sample, total=len(labels), remaining=len(todo), jobs=a.jobs)), flush=True)

    def save(wall, final=False):
        rows = sorted(log, key=lambda r: tuple(map(int, r['label'].split('.'))))
        summary = dict(d=a.d, signature=signature, sample=a.sample, expected_groups=len(labels),
                       groups=len(rows), complete=len(rows)==len(labels),
                       classes=sum(r['classes'] for r in rows),
                       unidentified=[r for r in rows if r['unidentified']],
                       seconds=prior_seconds+wall, wall_seconds_this_run=wall,
                       jobs=a.jobs, runs=rows,
                       estimate=(dict(status="complete", **estimate_runtime(db, rows, a.jobs))
                                 if final and len(rows)==len(labels)
                                 else dict(status="pending", reason="detailed projection is written only after all requested labels finish")))
        atomic_json(out / '_summary.json', summary)
        return summary

    with multiprocessing.Pool(a.jobs, initializer=worker_init) as pool:
        for label, n, bad, dt in pool.imap_unordered(functools.partial(work, out=out, signature=signature), todo):
            log.append(dict(label=label, classes=n, unidentified=bad, seconds=dt))
            save(time.monotonic()-mono0)
            print(f'{label}: {n} classes, {dt:.2f}s' + (' UNIDENTIFIED' if bad else ''), flush=True)
    summary = save(time.monotonic()-mono0, final=True)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('runs','estimate')}), flush=True)
    if summary['unidentified']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
