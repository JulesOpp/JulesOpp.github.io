#!/usr/bin/env python3
"""Count exact invariant buckets; corpus output is confined to local/."""
import collections
import json
import os
import statistics
import sys
sys.path.insert(0, os.path.dirname(__file__))
import ssg_equiv as se


def counts():
    out = {}
    for d in (2,3):
        records = se.load_records(d)
        buckets = collections.Counter(se.invariant(r,se.q_generic(r)) for r in records.values())
        same_sg=collections.Counter((label.split('.')[0],se.invariant(r,se.q_generic(r))) for label,r in records.items())
        sizes=list(buckets.values())
        out[str(d)] = dict(records=len(records),bravais=dict(sorted(collections.Counter(r['cls'] for r in records.values()).items())),
                           buckets={repr(k):v for k,v in sorted(buckets.items())},
                           pairs=sum(n*(n-1)//2 for n in sizes),same_sg_pairs=sum(n*(n-1)//2 for n in same_sg.values()),median=statistics.median(sizes) if sizes else 0,
                           maximum=max(sizes,default=0))
    return out


if __name__ == '__main__':
    out=counts()
    destination=os.path.join(se.HERE,'..','local','fixtures','subgroups','_counts.json')
    os.makedirs(os.path.dirname(destination),exist_ok=True)
    with open(destination,'w') as f: json.dump(out,f,indent=2); f.write('\n')
    print(json.dumps(out,indent=2))
