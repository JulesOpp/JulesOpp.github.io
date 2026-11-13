import { beforeAll, describe, expect, it } from 'vitest';
import { existsSync } from 'node:fs';
import { join } from 'node:path';
import { loadMoyo } from '../src/moyo';
import { buildOps, type GroupRecord } from '../src/ssg';
import { standardise } from '../src/standardise';
import { checkMap, randomSetting, records } from './standardise-multi-helpers';
beforeAll(async()=>{await loadMoyo();});
const hex:GroupRecord={sym:'P6/mmm(0,0,g1)0000(0,0,g2)0000',cls:'',d:2,q:['0,0,g1','0,0,g2'],params:['g1','g2'],cen:[],gen:['x-y,x,z,t,u','x,x-y,z,t,u','x,y,-z,-t,-u']};
// Hand I-centred axial pair. The integer fixed part is integral on the
// conventional axes; removing it can require a compensating free-coordinate shift.
const centred:GroupRecord={sym:'',cls:'',d:2,q:['1,0,g1','0,0,g2'],params:['g1','g2'],cen:['1/2,1/2,1/2,0,0'],gen:['-x,-y,z,t-2x,u','x,-y,-z,-t+2x,-u']};
describe('standardise higher dimensions (public)',()=>{
 for(const [name,r,sg] of [['P6/mmm',hex,191],['I-centred integer fixed component',centred,71]] as const)it(name,()=>{
  const own=standardise(buildOps(r),r.q);expect(own.ok).toBe(true);if(!own.ok)return;
  if(r===hex)expect(own.symbol.full).toBe(hex.sym);
  else expect(own.symbol.basic.startsWith('I')).toBe(true);
  expect(checkMap(buildOps(r),r.q,own)).toBe(true);
  for(let k=0;k<3;k++){const changed=randomSetting(r,sg,20261113+k),got=standardise(changed.ops,changed.q);expect(got.ok).toBe(true);if(got.ok){expect(got.symbol.full).toBe(own.symbol.full);expect(checkMap(changed.ops,changed.q,got)).toBe(true);}}
 });
});
// These stars exercise the unbounded monoclinic shear and the soft q-form rank.
const star={sym:'',cls:'',d:2,q:['a,b,g','-a,-b,g'],params:['a','b','g'],cen:[]};
describe('canonical naming regressions',()=>{
 for(const [sg,gen] of [[3,'-x,-y,z,u,t'],[4,'-x,-y,z+1/2,u,t'],[6,'x,y,-z,-u,-t']] as const)it(`monoclinic star ${sg}`,()=>{
  const r:GroupRecord={...star,gen:[gen]},own=standardise(buildOps(r),r.q);
  expect(own.ok).toBe(true);if(!own.ok)return;
  expect(checkMap(buildOps(r),r.q,own)).toBe(true);
  for(const seed of [20261158,20261159,20261181]){
   const t=randomSetting(r,sg,seed),got=standardise(t.ops,t.q);expect(got.ok).toBe(true);
   if(got.ok){expect(got.symbol.full).toBe(own.symbol.full);expect(checkMap(t.ops,t.q,got)).toBe(true);}
  }
 });
 it('compares parsed q vectors rather than term order, without accepting changed coefficients',()=>{
  const own=standardise(buildOps(hex),hex.q);expect(own.ok).toBe(true);if(!own.ok)return;
  const r={...own,P:[['1','0','0'],['0','1','0'],['0','0','1']],p:['0','0','0'],E:[[1,1],[0,1]],M:[['0','0','0'],['0','0','0']],p4:['0','0'],qActual:['0,0,g2+g1','0,0,g2']};
  expect(checkMap(buildOps(hex),hex.q,r)).toBe(true);
  expect(checkMap(buildOps(hex),hex.q,{...r,qActual:['0,0,g2-g1','0,0,g2']})).toBe(false);
 });
});
const DATA=join(__dirname,'..','local','data');

describe.skipIf(!process.env.SSG_SWEEP || !existsSync(DATA))('canonical d=2,3 invariance (local-only)',()=>{
 for(const d of [2,3])it(`d=${d}: every record in three seeded settings and every reported map`,async()=>{
  const recs=records(d);expect(recs).toHaveLength(d===2?3338:430);const failures:string[]=[];let matches=0,checks=0;
  for(const [i,r] of recs.entries()){
   const own=standardise(buildOps(r),r.q);if(own.ok&&own.symbol.full===r.sym)matches++;
   if(!own.ok)failures.push(`${r.label}\town\t${own.reason}`);else if(!checkMap(buildOps(r),r.q,own))failures.push(`${r.label}\town map`);
   for(let k=0;k<3;k++){const t=randomSetting(r,Number(r.label.split('.')[0]),20261113+i*3+k),got=standardise(t.ops,t.q);checks++;
    if(!own.ok||!got.ok||own.symbol.full!==got.symbol.full)failures.push(`${r.label}\tseed ${20261113+i*3+k}\t${own.ok?own.symbol.full:own.reason}\t${got.ok?got.symbol.full:got.reason}`);
    if(got.ok&&!checkMap(t.ops,t.q,got))failures.push(`${r.label}\tseed ${20261113+i*3+k}\tmap`);
   }
   await new Promise<void>(resolve=>setImmediate(resolve));
  }
  console.log(`d=${d} BYU exact ${matches}/${recs.length}; invariance cases ${checks}; failures ${failures.length}`);
  expect(failures,failures.slice(0,30).join('\n')).toEqual([]);
 },7_200_000);
});
