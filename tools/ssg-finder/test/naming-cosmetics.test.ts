import {beforeAll,describe,expect,it} from 'vitest';
import {simpleParameterForm,starParameterForm,phaseOrderCost,standardise} from '../src/standardise';
import {paramsOf,parseLinear,buildOps,type GroupRecord} from '../src/ssg';
import type {SsgSymbol} from '../src/ssgsymbol';
import {loadMoyo} from '../src/moyo';
import {checkMap} from './standardise-multi-helpers';
const simple=(q:string[])=>{const names=paramsOf(q);return simpleParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names);};
describe('cosmetic naming conventions',()=>{
 it('uses unindexed component names for one shared set, including opposite q pairs',()=>{
  expect(simple(['x,0,0','-x,x,0'])).toEqual(['a,0,0','-a,a,0']);
  expect(simple(['x,y,0','-x,y+1/2,0'])).toEqual(['a,b,0','-a,b+1/2,0']);
 });
 it('keeps independent sets distinct and numbers a star plus an independent row by set',()=>{
  expect(simple(['x,y,0','u,v,0'])).toEqual(['a1,b1,0','a2,b2,0']);
  expect(simple(['x,y,0','-x,y,0','0,0,z'])).toEqual(['a1,b1,0','-a1,b1,0','0,0,g2']);
 });
 it('uses unindexed names for a verified compound star',()=>{
  const q=['x,y,0','-x-y,x,0'],names=paramsOf(q);
  expect(starParameterForm(q.map(r=>r.split(',').map(x=>parseLinear(x,names))),names,[[[-1,1,0],[-1,0,0],[0,0,1]]],[])).toEqual(['a,b,0','-a-b,a,0']);
 });
 const symbol=(q:string[],shifts:string[]):SsgSymbol=>({basic:'P1',q:q.map(r=>`(${r})`).join(''),shifts,full:'P1'+q.map((r,i)=>`(${r})${shifts[i]}`).join('')});
 it('puts a nontrivial phase before a trivial phase only within identical q-order keys',()=>{
  expect(phaseOrderCost(symbol(['a1,b1,0','a2,b2,0'],['0','s']))).toBeGreaterThan(phaseOrderCost(symbol(['a1,b1,0','a2,b2,0'],['s','0'])));
  expect(phaseOrderCost(symbol(['a1,b1,0','0,0,g2'],['0','s']))).toBe(0);
  expect(phaseOrderCost(symbol(['a1,b1,0','a2,b2,0'],['0','0']))).toBe(0);
  expect(phaseOrderCost(symbol(['0,0,g1','0,0,g2','0,0,g3'],['s','-t','0']))).toBe(0);
 });
});
beforeAll(()=>loadMoyo());
it('carries the nontrivial mirror phase with the first independent planar block',()=>{
 const r:GroupRecord={sym:'',cls:'',d:2,q:['a1,b1,0','a2,b2,0'],params:['a1','b1','a2','b2'],cen:['1/2,0,1/2,0,0'],gen:['x,y,-z,t+1/2,u']};
 const ops=buildOps(r),got=standardise(ops,r.q);expect(got.ok).toBe(true);if(!got.ok)return;
 expect(got.symbol.shifts).toEqual(['s','0']);
 expect(checkMap(ops,r.q,got)).toBe(true);
});
