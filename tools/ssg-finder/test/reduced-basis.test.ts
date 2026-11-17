import {beforeAll,expect,it} from 'vitest';
import {loadMoyo} from '../src/moyo';
import {buildOps,paramsOf,type GroupRecord} from '../src/ssg';
import {standardise} from '../src/standardise';
import {checkMap} from './standardise-multi-helpers';
beforeAll(()=>loadMoyo());
// Hand tetragonal families: the extra vector shares removable star components.
for(const [q,gen,maximum] of [
 [['a,b,g','-b,a,g','a,b,0'],['-y,x,z,2t-u-2v,t,t-u-v'],5],
 [['a,b,g+1/2','-b,a,g+1/2','0,0,g'],['-y,x,z,-u+2v+z+1/4,t+1/4,v+1/4'],5],
 [['a,a,g','a,a,0','a,0,g'],['-y,x,z+1/4,-t+2v,-2t+u+2v,-u+v','x,-y,-z,-3t+2u+2v,-2t+u+2v,-2t+2u+v'],3],
] as const)it(`reduces shared components in ${q.join(';')}`,()=>{
 const r:GroupRecord={sym:'',cls:'',d:3,q:[...q],params:paramsOf([...q]),cen:[],gen:[...gen]};
 const ops=buildOps(r),got=standardise(ops,r.q);expect(got.ok).toBe(true);if(!got.ok)return;
 expect((got.symbol.q.match(/[abg]\d*/g)??[]).length).toBeLessThanOrEqual(maximum);
 expect(checkMap(ops,r.q,got)).toBe(true);
});
