import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
// Exact independent map checks and deterministic setting changes; no oracle data.
import { matInv, matMul, Q, type QMat } from '../src/rational';
import { imatDet } from '../src/intmat';
import { buildOps, paramsOf, parseLinear, type GroupRecord, type SOp } from '../src/ssg';
import { formatQ, transformOps, type Standardised } from '../src/standardise';
import { crystalSystemOf } from '../src/symbol';
import { ssgSymbolInAxes } from '../src/ssgsymbol';
export const I = (n:number) => Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>+(i===j)));
export const rng=(seed:number)=>()=>{seed=(seed+0x6d2b79f5)>>>0;let t=Math.imul(seed^(seed>>>15),seed|1);t^=t+Math.imul(t^(t>>>7),t|61);return ((t^(t>>>14))>>>0)/4294967296;};
const qm=(A:number[][])=>A.map(row=>row.map(x=>Q.of(x)));
const mv=(A:QMat,v:Q[])=>A.map(row=>row.reduce((s,x,j)=>s.add(x.mul(v[j])),Q.ZERO));
const add=(a:Q[],b:Q[])=>a.map((x,i)=>x.add(b[i]));
const sub=(a:Q[],b:Q[])=>a.map((x,i)=>x.sub(b[i]));
const key=(o:SOp)=>o.W.map(row=>row.join(',')).join(';')+'|'+o.w.map(x=>x.mod1().toString()).join(',');
const identity=(o:SOp)=>o.W.every((row,i)=>row.every((x,j)=>x===+(i===j)));
export function mappedQ(q:string[],P:QMat,E:number[][],M:QMat):string[] {
 const names=paramsOf(q), rows=q.map(row=>row.split(',').map(x=>parseLinear(x,names))), Pi=matInv(P);
 return E.map((row,i)=>formatQ([0,1,2].map(j=>({
  coef:names.map((_,k)=>rows.reduce((s,r,l)=>s.add(r.reduce((v,x,h)=>v.add(x.coef[k].mul(Pi[h][j]).mulInt(row[l])),Q.ZERO)),Q.ZERO)),
  c:rows.reduce((s,r,l)=>s.add(r.reduce((v,x,h)=>v.add(x.c.mul(Pi[h][j]).mulInt(row[l])),Q.ZERO)),Q.ZERO).add(M[i].reduce((s,x,h)=>s.add(x.mul(Pi[h][j])),Q.ZERO))
 })),names));
}
export function randomSetting(r:GroupRecord,sg:number,seed:number):{ops:SOp[];q:string[]} {
 const rnd=rng(seed), pick=<T>(a:T[])=>a[Math.floor(rnd()*a.length)];
 const perms=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
 const signed=perms.flatMap(p=>Array.from({length:8},(_,s)=>p.map((j,i)=>[0,1,2].map(k=>k===j?((s>>i)&1?-1:1):0)))).filter(A=>imatDet(A.map(row=>row.map(BigInt)))===1n);
 const mul=(A:number[][],B:number[][])=>A.map(row=>B[0].map((_,j)=>row.reduce((s,x,k)=>s+x*B[k][j],0)));
 const sys=crystalSystemOf(sg);let P=I(3);
 if(sys==='triclinic'||sys==='monoclinic') {
  for(let k=0;k<2;k++){const A=I(3),[i,j]=pick(sys==='monoclinic'?[[0,1],[1,0]]:[[0,1],[1,0],[0,2],[2,0],[1,2],[2,1]]);A[i][j]=pick([-1,1]);P=mul(A,P);}P=mul(pick(signed),P);
 } else if(sys==='orthorhombic'||sys==='cubic')P=pick(signed);
 else if(sys==='tetragonal')P=pick([I(3),[[0,-1,0],[1,0,0],[0,0,1]],[[0,1,0],[1,0,0],[0,0,-1]],[[0,1,0],[-1,0,0],[0,0,1]]]);
 else {for(let k=Math.floor(rnd()*6);k>0;k--)P=mul([[1,-1,0],[1,0,0],[0,0,1]],P);if(rnd()<0.5)P=mul([[0,1,0],[1,0,0],[0,0,-1]],P);}
 let E:number[][];do {E=I(r.d).map(row=>row.map(()=>Math.floor(rnd()*3)-1));}while(![1n,-1n].includes(imatDet(E.map(row=>row.map(BigInt)))));
 const frac=()=>Q.of(Math.floor(rnd()*12),12),p=[frac(),frac(),frac()],p4=E.map(frac),M=E.map(()=>[Q.ZERO,Q.ZERO,Q.ZERO]);
 const ops=transformOps(buildOps(r),qm(P),p,E,M,p4);if(!ops)throw Error('invalid test setting');
 return {ops,q:mappedQ(r.q,qm(P),E,M)};
}
/** TS analogue of check_map: unimodular superspace lattice, symbolic affine
 * family, every conjugated operation, and equal operation/lattice multiplicity.
 * The target is the reported setting, reconstructed and checked against its symbol.
 */
export function checkMap(ops:SOp[],q:string[],r:Standardised):boolean {
 const d=q.length,n=3+d,P=r.P.map(row=>row.map(Q.parse)),M=r.M.map(row=>row.map(Q.parse));
 const E=qm(r.E),S=[...P.map(row=>[...row,...Array(d).fill(Q.ZERO)]),...M.map((row,i)=>[...row,...E[i]])];
 const Si=matInv(S),shift=[...r.p,...r.p4].map(Q.parse);
 if(![1n,-1n].includes(imatDet(r.E.map(row=>row.map(BigInt)))))return false;
 const aq=mappedQ(q,P,r.E,M),names=paramsOf([...aq,...r.qActual]).sort();
 const flat=(rows:string[])=>rows.flatMap(row=>row.split(',').map(x=>parseLinear(x,names)));
 const aa=flat(aq),bb=flat(r.qActual);
 if(aa.some((x,i)=>!x.c.eq(bb[i].c)||x.coef.some((v,k)=>!v.eq(bb[i].coef[k]))))return false;
 const target=transformOps(ops,P,r.p.map(Q.parse),r.E,M,r.p4.map(Q.parse));if(!target)return false;
 const sourceC=ops.filter(identity),targetC=target.filter(identity);
 const inLattice=(v:Q[],cs:SOp[])=>cs.some(c=>sub(v,c.w).every(x=>x.isInt()));
 for(let j=0;j<n;j++)if(!inLattice(S.map(row=>row[j]),targetC)||!inLattice(Si.map(row=>row[j]),sourceC))return false;
 for(const c of sourceC)if(!inLattice(mv(S,c.w),targetC))return false;
 for(const c of targetC)if(!inLattice(mv(Si,c.w),sourceC))return false;
 const targetKeys=new Set(target.map(key));
 for(const o of ops){const W=matMul(matMul(S,qm(o.W)),Si);if(W.some(row=>row.some(x=>!x.isInt())))return false;
  const w=sub(add(mv(S,o.w),shift),mv(W,shift));if(!targetKeys.has(key({W:W.map(row=>row.map(x=>Number(x.n))),w})))return false;}
 if(ops.length*targetC.length!==target.length*sourceC.length)return false;
 const form=[...r.symbol.q.matchAll(/\(([^()]*)\)/g)].map(m=>m[1]);
 if(ssgSymbolInAxes(target,form,crystalSystemOf(r.sgNumber))?.full!==r.symbol.full)return false;
 // Equality of symbolic affine spans, including constants (parameter names may change).
 const sourceNames=paramsOf(r.qActual),targetNames=paramsOf(form);
 const parse=(rows:string[],names:string[])=>rows.flatMap(row=>row.split(',').map(x=>parseLinear(x,names)));
 const a=parse(r.qActual,sourceNames),b=parse(form,targetNames);
 const rank=(rows:Q[][])=>{if(!rows.length)return 0;let c=0;for(let j=0;j<rows[0].length&&c<rows.length;j++){let k=rows.findIndex((row,i)=>i>=c&&!row[j].isZero());if(k<0)continue;[rows[c],rows[k]]=[rows[k],rows[c]];const v=rows[c][j];rows[c]=rows[c].map(x=>x.div(v));for(let i=0;i<rows.length;i++)if(i!==c){const f=rows[i][j];rows[i]=rows[i].map((x,h)=>x.sub(f.mul(rows[c][h])));}c++;}return c;};
 const A=sourceNames.map((_,k)=>a.map(x=>x.coef[k])),B=targetNames.map((_,k)=>b.map(x=>x.coef[k]));
 const rA=rank(A.map(row=>[...row])),rB=rank(B.map(row=>[...row]));
 return rA===rB&&rank([...A,...B,a.map((x,i)=>x.c.sub(b[i].c))].map(row=>[...row]))===rA;
}

const DATA=join(new URL(".", import.meta.url).pathname,"..","local","data");
export function records(d:number):(GroupRecord & {label:string})[]{
 const all:(GroupRecord & {label:string})[]=[];
 for(const file of readdirSync(DATA).filter(f=>/^sg\d{3}\.json$/.test(f)).sort())for(const [label,r] of Object.entries<GroupRecord>(JSON.parse(readFileSync(join(DATA,file),'utf8')).groups))if(r.d===d)all.push({...r,label});
 if(d===2)return all;
 const cls=new Map(),sg=new Map();for(const r of all){if(!cls.has(r.cls))cls.set(r.cls,r);const k=r.label.split('.')[0];if(!sg.has(k))sg.set(k,r);}
 return [...new Map([...cls.values(),...sg.values()].map(r=>[r.label,r])).values()].sort((a,b)=>a.label<b.label?-1:a.label>b.label?1:0);
}