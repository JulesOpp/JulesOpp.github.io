import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "convert_check.py")).read().split("ok=0")[0])
import json, gzip, collections, time
t=time.time()
def fmt(W,w,n):
    V="xyztuv"[:n]; rows=[]
    for i in range(n):
        s=""
        for j in range(n):
            c=W[i][j]
            if c: s+=("+" if c>0 else "-")+("" if abs(c)==1 else str(abs(c)))+V[j]
        if w[i]: s+="+"+str(w[i])
        rows.append(s.lstrip("+") or "0")
    return ",".join(rows)
by=collections.defaultdict(list); nops=[]
for lab,g in groups.items():
    n=3+g["d"]; ops=basic_ops(lab); nops.append(len(ops))
    by[lab.split(".")[0]].append([lab,g["sym"],[fmt(W,w,n) for W,w in sorted(ops)]])
tot=json.dumps(by,separators=(",",":")).encode()
print("groups",len(groups),"max ops(incl centring)",max(nops),"time",round(time.time()-t,1))
print("all raw",len(tot)//1024,"KB gz",len(gzip.compress(tot))//1024,"KB")
sz=sorted((len(json.dumps(v,separators=(",",":"))),k) for k,v in by.items())
print("per-SG raw KB: min",sz[0][0]//1024,"median",sz[len(sz)//2][0]//1024,"max",sz[-1][0]//1024,sz[-1][1], "gz max",len(gzip.compress(json.dumps(by[sz[-1][1]]).encode()))//1024)
