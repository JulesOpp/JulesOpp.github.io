#!/usr/bin/env python3
"""Resumable local corpus validation. Each successful case is flushed to JSONL.

Run a stratified --limit before long sweeps. d3 defaults to the sample; request
--full explicitly only after estimating its cost. Worker pools are capped at four processes; each worker uses one BLAS thread.
"""
import argparse
from collections import defaultdict
from itertools import combinations
import hashlib
import json
import multiprocessing
import os
import sys
import time
sys.path.insert(0,os.path.dirname(__file__))
import ssg_equiv as se


def sample(db):
    chosen={}
    for label,r in db.items():
        chosen.setdefault(('cls',r['cls']),label)
        chosen.setdefault(('sg',label.split('.')[0]),label)
    return {label:db[label] for label in sorted(set(chosen.values()))}


def pairs(db):
    buckets=defaultdict(list)
    for label,r in db.items():
        buckets[(label.split('.')[0],se.invariant(r,se.q_generic(r)))].append(label)
    return [pair for labels in buckets.values() for pair in combinations(labels,2)]


_DB = None
_KIND = None


def _worker_init(d,kind,full):
    global _DB,_KIND
    _DB=se.load_records(d)
    if d==3 and not full: _DB=sample(_DB)
    _KIND=kind


def _check(case):
    start=time.monotonic();a=_DB[case[0]];q=se.q_generic(a)
    b=a if _KIND=='self' else _DB[case[1]]
    result=se.equivalent(a,q,b)
    ok=result is not None and se.check_map(a,q,b,result) if _KIND=='self' else result is None
    return dict(case=case,ok=ok,seconds=time.monotonic()-start)


def run(d,kind,full=False,limit=None,jobs=1):
    db=se.load_records(d)
    if d==3 and not full: db=sample(db)
    cases=[(l,) for l in db] if kind=='self' else pairs(db)
    if limit is not None:
        # Stratify over the whole ordered case list, rather than easy first SGs.
        cases=[cases[i*len(cases)//min(limit,len(cases))] for i in range(min(limit,len(cases)))]
    path=os.path.join(se.HERE,'..','local','fixtures','subgroups',f'_oracle_d{d}_{kind}_{"full" if full or d==2 else "sample"}.jsonl')
    # Bind checkpoints to arithmetic, map checking and the case-selection logic.
    digest=hashlib.sha256()
    for source in (se.__file__,__file__,se.exact.__file__,se.im.__file__):
        with open(source,'rb') as f: digest.update(f.read())
        digest.update(b'\0')
    digest.update(json.dumps(db,sort_keys=True).encode())
    signature=digest.hexdigest()
    done=set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    row=json.loads(line)
                except json.JSONDecodeError:
                    continue  # interrupted final write: rerun the unverified case
                if row['ok'] and row.get('signature')==signature: done.add(tuple(row['case']))
    print(json.dumps(dict(d=d,kind=kind,total=len(cases),remaining=sum(c not in done for c in cases),path=path)),flush=True)
    os.makedirs(os.path.dirname(path),exist_ok=True)
    pending=[case for case in cases if case not in done]
    jobs=max(1,min(int(jobs),4))
    def save(rows):
        with open(path,'a') as f:
            for i,row in enumerate(rows,1):
                row['signature']=signature
                f.write(json.dumps(row)+'\n');f.flush()
                if i%100==0 or row['seconds']>5 or not row['ok']:
                    print(json.dumps({k:v for k,v in row.items() if k!='signature'}),flush=True)
                if not row['ok']: raise AssertionError(row['case'])
    start=time.monotonic()
    if jobs==1:
        _worker_init(d,kind,full)
        save(map(_check,pending))
    else:
        with multiprocessing.get_context('spawn').Pool(jobs,initializer=_worker_init,initargs=(d,kind,full)) as pool:
            save(pool.imap_unordered(_check,pending,chunksize=1))
    print(json.dumps(dict(completed=len(pending),seconds=time.monotonic()-start,jobs=jobs)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--d',type=int,choices=[2,3],required=True)
    parser.add_argument('--kind',choices=['self','pairs'],default='pairs')
    parser.add_argument('--full',action='store_true')
    parser.add_argument('--limit',type=int)
    parser.add_argument('--jobs',type=int,default=1,help='1–4 worker processes (clamped)')
    args=parser.parse_args()
    if args.limit is not None and args.limit<1: parser.error('--limit must be positive')
    run(args.d,args.kind,args.full,args.limit,args.jobs)
