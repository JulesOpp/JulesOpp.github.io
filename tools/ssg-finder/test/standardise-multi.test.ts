import { beforeAll, describe, expect, it } from 'vitest';
import { existsSync } from 'node:fs';
import { join } from 'node:path';
import { Q } from '../src/rational';
import { loadMoyo } from '../src/moyo';
import { buildOps, paramsOf, parseLinear, type GroupRecord } from '../src/ssg';
import { compareQBlocks, simpleParameterForm, starParameterForm, standardise, transformOps } from '../src/standardise';
import { checkMap, mappedQ, randomSetting, records } from './standardise-multi-helpers';
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

// Hand families: these contain no catalogue labels or oracle-derived records.
describe('simple canonical q spelling',()=>{
 it('reparametrises independent compound components with component/block names',()=>{
  const q=['a+b,a-b,0','0,0,c+d'],names=paramsOf(q).sort();
  expect(simpleParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names)).toEqual(['a1,b1,0','0,0,g2']);
 });
 it('preserves shared parameter relations and rational offsets',()=>{
  const q=['a+b+1/3,0,0','-a-b,0,1/3'],names=paramsOf(q).sort();
  expect(simpleParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names)).toEqual(['a+1/3,0,0','-a,0,1/3']);
 });
 it('does not pretend dependent sums are independent parameters',()=>{
  const q=['a,b,0','a+b,0,0'],names=paramsOf(q).sort();
  expect(simpleParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names)).toBeNull();
 });
 it('orders free, rational descending, then zero by component',()=>{
  const q=['0,0,g1','1/3,0,g2','1/2,0,g3','a4,0,0'];
  expect(q.sort(compareQBlocks)).toEqual(['a4,0,0','1/2,0,g3','1/3,0,g2','0,0,g1']);
 });
 it('keeps a simple fallback through setting-cost pruning',()=>{
  const r:GroupRecord={...star,gen:['-x,-y,z,u,t']};
  const got=standardise(buildOps(r),r.q);expect(got.ok).toBe(true);if(!got.ok)return;
  for(const block of got.symbol.q.matchAll(/\(([^)]*)\)/g))for(const x of block[1].split(','))
   expect((x.match(/[abg]\d*/g)??[]).length).toBeLessThanOrEqual(1);
  expect(checkMap(buildOps(r),r.q,got)).toBe(true);
 });
});

describe('exact star compound exception',()=>{
 const R=[[-1,1,0],[-1,0,0],[0,0,1]];
 const form=(q:string[],rotations=[R],cents:Q[][]=[])=>{const names=paramsOf(q).sort();return starParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names,rotations,cents);};
 it('allows a compound only in a later exact star image',()=>{
  expect(form(['a,b,0','-a-b,a,0'])).toEqual(['a,b,0','-a-b,a,0']);
  expect(form(['a,b,0','-a-b,a,0'],[])).toBeNull();
  expect(form(['a,b,0','-a-2b,a,0'])).toBeNull();
  expect(form(['a+b,a,0','b,0,0'])).toBeNull();
 });
 it('checks reciprocal offsets exactly, including centred lattices',()=>{
  expect(form(['a,b,0','-a-b+1,a,0'])).toEqual(['a,b,0','-a-b+1,a,0']);
  expect(form(['a,b,0','-a-b+1/2,a,0'])).toBeNull();
  expect(form(['a,b,0','-a-b+1,a,0'],[R],[[Q.of(1,2),Q.ZERO,Q.ZERO]])).toBeNull();
 });
 it('names a hand threefold star and verifies its exact map',()=>{
  const r:GroupRecord={sym:'',cls:'',d:2,q:['a,b,0','-a-b,a,0'],params:['a','b'],cen:[],gen:['-x+y,-x,z,u,-t-u']};
  const own=standardise(buildOps(r),r.q);expect(own.ok).toBe(true);if(!own.ok)return;
  expect(checkMap(buildOps(r),r.q,own)).toBe(true);
  expect(own.symbol.q).not.toMatch(/[abg][0-9]/);
 });
});

it('retains the existing rational-setting cost and carries phases with reordered blocks',()=>{
 const r:GroupRecord={sym:'',cls:'',d:2,q:['1/2,0,g1','0,0,g2'],params:['g1','g2'],cen:[],gen:['-x,-y,z,t-x,u+1/2']};
 const ops=buildOps(r),own=standardise(ops,r.q);expect(own.ok).toBe(true);if(!own.ok)return;
 const form=[...own.symbol.q.matchAll(/\(([^)]*)\)/g)].map(m=>m[1]);
 expect(form).toEqual(['1/2,0,g1','0,0,g2']);
 expect(checkMap(ops,r.q,own)).toBe(true);
 const P=[[1,0,0],[0,1,0],[0,0,1]].map(r=>r.map(x=>Q.of(x))),E=[[0,1],[1,0]],M=[[Q.ZERO,Q.ZERO,Q.ZERO],[Q.ZERO,Q.ZERO,Q.ZERO]];
 const changed=transformOps(ops,P,[Q.ZERO,Q.ZERO,Q.ZERO],E,M,[Q.ZERO,Q.ZERO])!;
 const q=mappedQ(r.q,P,E,M),swapped=standardise(changed,q);expect(swapped.ok).toBe(true);if(!swapped.ok)return;
 expect(swapped.symbol.full).toBe(own.symbol.full);
 expect(checkMap(changed,q,swapped)).toBe(true);
});
