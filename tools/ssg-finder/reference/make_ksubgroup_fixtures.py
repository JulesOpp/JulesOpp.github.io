#!/usr/bin/env python3
"""Resumable d=1 k-index sample and fixture generation (local output only)."""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent))
import ksubgroups as ks

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT/'local'


def encoded(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def atomic(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp = path.with_suffix(path.suffix+'.partial')
    temp.write_bytes(encoded(value))
    temp.replace(path)


def signature(db):
    h = hashlib.sha256(encoded(db))
    for name in ('ksubgroups.py','make_ksubgroup_fixtures.py','subgroups.py','ssg_equiv.py','intmat.py','exact.py','affine_search.py'):
        h.update((ROOT/'reference'/name).read_bytes())
    return h.hexdigest()


def crystal(label):
    sg = int(label.split('.')[0])
    return next(i for i,end in enumerate((2,15,74,142,167,194,230)) if sg <= end)


def sample_labels(db):
    """Thirty parents, round-robin over present crystal systems, largest |K| first."""
    order = {label:len({W for W,w in ks.exact.build_ops(rec)}) for label,rec in db.items()}
    buckets = {}
    for label in db:
        buckets.setdefault(crystal(label),[]).append(label)
    for labels in buckets.values():
        labels.sort(key=lambda x:(-order[x],tuple(map(int,x.split('.')))))
    result = []
    while len(result) < min(30,len(db)):
        for system in sorted(buckets):
            if buckets[system] and len(result) < 30:
                result.append(buckets[system].pop(0))
    assert max(order[x] for x in result) == max(order.values())
    return result


_DB = None
_SIGNATURE = None


def init(db,sig):
    global _DB,_SIGNATURE
    _DB,_SIGNATURE = db,sig


def worker(item):
    label,k,path = item
    rec = _DB[label]
    start = time.monotonic()
    stats = {}
    rows = ks.k_subgroups(rec,k,_DB,stats=stats)
    failures = [dict(label=label,k=k,kind=r['kind'],HNF=r['sigma'],matches=r['matches'])
                for r in rows if len(r['matches']) != 1]
    counts = {name:sum((len(r['matches']) == 1 if name == 'exactly_one' else
                       len(r['matches']) == 0 if name == 'none' else len(r['matches']) > 1)
                      for r in rows) for name in ('exactly_one','none','several')}
    classes = []
    for r in rows:
        c = dict(r)
        c.update(c.pop('map',{}))
        c.setdefault('byu',None)
        classes.append(c)
    value = dict(label=label,d=1,k_index=k,sym=rec['sym'],
                 order=len({W for W,w in ks.exact.build_ops(rec)}),classes=classes,
                 _signature=_SIGNATURE,_seconds=time.monotonic()-start,
                 _stats=stats,_identification=counts,_failures=failures)
    value['_content_sha256'] = digest(value)
    atomic(Path(path),value)
    return dict(label=label,k=k,seconds=value['_seconds'],classes=len(rows),
                stats=stats,identification=counts,failures=failures)


def cached(path,sig):
    try:
        data = json.loads(path.read_text())
        checksum = data.pop('_content_sha256')
        if data['_signature'] != sig or digest(data) != checksum:
            return None
        return dict(label=data['label'],k=data['k_index'],seconds=data['_seconds'],
                    classes=len(data['classes']),stats=data['_stats'],
                    identification=data['_identification'],failures=data['_failures'])
    except (OSError,ValueError,KeyError):
        return None


def notes(results,db,jobs,complete):
    lines = ['# Task 10 d=1 measured cost checkpoint','',
             'Sample: 30 parents, every crystal system present in d=1, including maximum |K|.',
             'Projection: crystal-system mean worker seconds × population / workers.',
             'This excludes startup, scheduling imbalance and sample bias; it is not a bound.','',
             '| d | k | parents measured | classes | measured worker seconds | projected wall seconds |',
             '|---|---|---|---|---|---|']
    projections = {}
    for k in (2,3,4):
        measured = [r for r in results if r['k'] == k]
        if not measured:
            continue
        groups = {}
        for r in measured:
            groups.setdefault(crystal(r['label']),[]).append(r['seconds'])
        fallback = sum(r['seconds'] for r in measured)/len(measured)
        seconds = sum(sum(groups[crystal(label)])/len(groups[crystal(label)])
                      if crystal(label) in groups else fallback for label in db)/jobs
        projections[k] = seconds
        lines.append(f'| 1 | {k} | {len(measured)} | {sum(r["classes"] for r in measured)} | {sum(r["seconds"] for r in measured):.3f} | {seconds:.3f} |')
    lines += ['',f'Sample complete: {complete}. Effective workers: {jobs}.',
              'Metrics and exact identification failures are in sample checkpoint files.']
    failures = [f for r in results for f in r['failures']]
    if failures:
        lines += ['', 'STOP: unidentified or multiply identified classes (first five):']
        lines += ['- '+json.dumps(f,sort_keys=True) for f in failures[:5]]
    if any(s > 7200 for s in projections.values()):
        lines += ['', 'STOP: a projected d=1 k sweep exceeds two hours.']
    (LOCAL/'task10-logs/NOTES.md').write_text('\n'.join(lines)+'\n')
    return projections


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',action='store_true')
    p.add_argument('--resume',action='store_true',required=True)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--k',type=int,choices=(2,3,4),nargs='+',default=[2,3,4])
    args = p.parse_args()
    jobs = max(1,min(args.jobs,4))
    db = ks.se.load_records(1)
    sig = signature(db)
    labels = sample_labels(db) if args.sample else sorted(db,key=lambda x:tuple(map(int,x.split('.'))))
    log = LOCAL/'task10-logs'
    log.mkdir(parents=True,exist_ok=True)
    results,pending = [],[]
    for k in args.k:
        for label in labels:
            base = log/'sample'/f'k{k}' if args.sample else LOCAL/'fixtures/ksubgroups/d1'/f'k{k}'
            path = base/f'{label}.json'
            old = cached(path,sig)
            if old and not old['failures']:
                results.append(old)
            else:
                pending.append((label,k,str(path)))
    # A full sweep is allowed only after all 90 matching sample checkpoints pass.
    if not args.sample:
        sample = [cached(log/'sample'/f'k{k}'/f'{label}.json',sig)
                  for k in (2,3,4) for label in sample_labels(db)]
        if any(r is None or r['failures'] for r in sample):
            p.error('complete successful sample required before full generation')
        if any(s > 7200 for s in notes(sample,db,jobs,True).values()):
            p.error('projected wall time exceeds two-hour gate')
    start = time.monotonic()
    # Bound in-flight work, so a failure stops at most four active parents.
    with concurrent.futures.ProcessPoolExecutor(jobs,initializer=init,initargs=(db,sig)) as pool:
        work = iter(pending)
        active = {}
        for _ in range(jobs):
            item = next(work,None)
            if item is not None:
                active[pool.submit(worker,item)] = item
        stopped = False
        while active:
            done,_ = concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                item = active.pop(future)
                result = future.result()
                results.append(result)
                print(json.dumps(result,sort_keys=True),flush=True)
                if result['failures']:
                    stopped = True
            projections = notes(results,db,jobs,not active and len(results) == len(labels)*len(args.k)) if args.sample else {}
            atomic(log/('sample-summary.json' if args.sample else 'fixture-summary.json'),
                   dict(results=results,jobs=jobs,seconds=time.monotonic()-start,
                        complete=len(results) == len(labels)*len(args.k),signature=sig))
            if not stopped:
                for _ in range(len(done)):
                    item = next(work,None)
                    if item is not None:
                        active[pool.submit(worker,item)] = item
        if stopped:
            print('STOP: identification gate failed; see NOTES.md',flush=True)
            return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
