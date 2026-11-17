import {beforeAll,expect,it} from 'vitest';
import {loadMoyo,getMoyo,fromColMajor} from '../src/moyo';
import {numToQ} from '../src/identify';
import {Q,matInv,matMul} from '../src/rational';
import {buildOps,formatOp,paramsOf,type GroupRecord} from '../src/ssg';
import {standardise} from '../src/standardise';
import {tSubgroups} from '../src/subgroups';
import {checkMap} from './standardise-multi-helpers';
beforeAll(()=>loadMoyo());
it('names a specialised P1 family using its full plain type family',()=>{
 const q=['g1,g1,2g1','g2,g2,2g2'];
 const r:GroupRecord={sym:'hand P1',cls:'',d:2,q,params:paramsOf(q),cen:[],gen:[]};
 const ops=buildOps(r),st=standardise(ops,q);expect(st.ok).toBe(true);if(!st.ok)return;
 expect(st.symbol.full).toBe('P1(a1,b1,g1)0(a2,b2,g2)0');
 expect(checkMap(ops,q,st)).toBe(true);
 // Membership must still reject an actual component excluded by the type.
 expect(checkMap(ops,q,{...st,symbol:{...st.symbol,q:'(a1,b1,0)(a2,b2,0)',full:'P1(a1,b1,0)0(a2,b2,0)0'}})).toBe(false);
});
it('accepts a hand affine P1 specialisation with rational offsets',()=>{
 // Independently chosen rational coefficients; saved catalogue regressions stay local.
 const q=['2a+1/5,a-1/7,0','0,b+2/5,3b'];
 const r:GroupRecord={sym:'hand P1',cls:'',d:2,q,params:paramsOf(q),cen:[],gen:[]};
 const ops=buildOps(r),st=standardise(ops,q);expect(st.ok).toBe(true);if(st.ok){expect(st.symbol.full).toBe('P1(a1,b1,g1)0(a2,b2,g2)0');expect(checkMap(ops,q,st)).toBe(true);}
});
it('names every subgroup of an algebraic Ia-3d cubic-star lift',()=>{
 // Public 3D Ia-3d operations from moyo, lifted algebraically. Q=a*A,
 // E=A*R*A^-1, M=0 and zero internal translations. No catalogue input.
 const A=[[1,1,1],[1,-1,-1],[-1,1,-1]].map(r=>r.map(x=>Q.of(x))),Ai=matInv(A);
 const spatial=getMoyo().operations_from_number(230,{type:'Standard'} as never,false) as {rotation:number[];translation:number[]}[];
 const gen=spatial.map(o=>{
  const R=fromColMajor(o.rotation),E=matMul(matMul(A,R.map(r=>r.map(x=>Q.of(x)))),Ai).map(r=>r.map(x=>{expect(x.isInt()).toBe(true);return Number(x.n);}));
  return formatOp({W:[...R.map(r=>[...r,0,0,0]),...E.map(r=>[0,0,0,...r])],w:[...o.translation.map(numToQ),Q.ZERO,Q.ZERO,Q.ZERO]});
 });
 const cen=spatial.filter(o=>fromColMajor(o.rotation).every((r,i)=>r.every((x,j)=>x===+(i===j))))
  .map(o=>[...o.translation.map(numToQ),Q.ZERO,Q.ZERO,Q.ZERO]).filter(w=>w.some(x=>!x.isZero())).map(w=>w.map(String).join(','));
 const r:GroupRecord={sym:'hand cubic star',cls:'',d:3,q:['a,a,a','a,-a,-a','-a,a,-a'],params:['a'],cen,gen};
 const classes=tSubgroups(r);expect(classes).toHaveLength(33);
 for(const [i,c] of classes.entries()){const st=standardise(c.ops,r.q);expect(st.ok,`class ${i+1}`).toBe(true);if(st.ok)expect(checkMap(c.ops,r.q,st),`class ${i+1} exact map`).toBe(true);}
},600_000);
