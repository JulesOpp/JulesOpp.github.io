import {beforeAll,describe,expect,it} from 'vitest';
import {Q} from '../src/rational';
import {buildOps,type GroupRecord} from '../src/ssg';
import {loadMoyo} from '../src/moyo';
import {reciprocalOffsetRepresentative,standardise} from '../src/standardise';
import {checkMap} from './standardise-multi-helpers';
const q=(xs:string[])=>xs.map(Q.parse);
const centres={C:[q(['1/2','1/2','0'])],A:[q(['0','1/2','1/2'])],I:[q(['1/2','1/2','1/2'])],F:[q(['0','1/2','1/2']),q(['1/2','0','1/2']),q(['1/2','1/2','0'])]};
const sameCoset=(a:Q[],b:Q[],cs:Q[][])=>{
 const delta=a.map((v,j)=>v.sub(b[j]));
 return delta.every(v=>v.isInt())&&cs.every(c=>c.reduce((v,x,j)=>v.add(x.mul(delta[j])),Q.ZERO).isInt());
};
describe('offset cosets of centred basic reciprocal lattices',()=>{
 for(const [name,allowed,forbidden] of [
  ['C',['1','1','0'],['1','0','0']],
  ['A',['0','1','1'],['0','1','0']],
  ['I',['1','1','0'],['1','0','0']],
  ['F',['1','1','1'],['1','0','0']],
 ] as const)it(`${name}: removes reciprocal offsets without removing a nonreciprocal integer vector`,()=>{
  const cs=centres[name];
  expect(reciprocalOffsetRepresentative(q([...allowed]),cs).map(String)).toEqual(['0','0','0']);
  const x=q([...forbidden]),r=reciprocalOffsetRepresentative(x,cs);
  expect(r.some(v=>!v.isZero())).toBe(true);
  expect(sameCoset(x,r,cs)).toBe(true);
 });
 for(const [name,cs] of Object.entries(centres))it(`${name}: exact, deterministic and idempotent across reciprocal translates`,()=>{
  const x=q(['7/4','-5/6','13/8']),a=reciprocalOffsetRepresentative(x,cs);
  expect(sameCoset(x,a,cs)).toBe(true);
  expect(reciprocalOffsetRepresentative(a,cs).map(String)).toEqual(a.map(String));
  for(let i=-2;i<=2;i++)for(let j=-2;j<=2;j++)for(let k=-2;k<=2;k++){
   const v=[i,j,k].map(x=>Q.of(x));if(!sameCoset(v,[Q.ZERO,Q.ZERO,Q.ZERO],cs))continue;
   expect(reciprocalOffsetRepresentative(x.map((a,h)=>a.add(v[h])),cs).map(String)).toEqual(a.map(String));
  }
 });
 it('minimises absolute offsets and chooses the nonnegative half at an absolute-value tie',()=>{
  expect(reciprocalOffsetRepresentative(q(['3/4','3/4','0']),centres.C).map(String)).toEqual(['-1/4','-1/4','0']);
  expect(reciprocalOffsetRepresentative(q(['-1/2','-1/2','0']),centres.C).map(String)).toEqual(['1/2','1/2','0']);
  expect(reciprocalOffsetRepresentative(q(['2','0','0']),centres.F).map(String)).toEqual(['0','0','0']);
 });
});
beforeAll(()=>loadMoyo());
it('removes a reciprocal (1,1,0) offset from a hand C-centred star and retains its exact map',()=>{
 const r:GroupRecord={sym:'',cls:'',d:2,q:['a,b,0','-a+1,b+1,0'],params:['a','b'],cen:['1/2,1/2,0,0,0'],gen:['-x,-y,z,-t,-u','-x,y,-z,u-x-y,t-x+y']};
 const ops=buildOps(r),got=standardise(ops,r.q);expect(got.ok).toBe(true);if(!got.ok)return;
 expect(got.symbol.q).not.toMatch(/[abg][0-9]*[+-]1/);
 expect(checkMap(ops,r.q,got)).toBe(true);
});
