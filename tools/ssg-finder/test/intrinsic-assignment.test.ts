import { expect, it } from 'vitest';
import { Q } from '../src/rational';
import { bestIntrinsicAssignment, compareIntrinsicAssignments } from '../src/ssgsymbol';

it('minimizes the phase product exactly like exhaustive whole-assignment comparison',()=>{
 let state=20261113;
 const pick=(n:number)=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state%n;};
 const values=['0','1/2','1/3','-1/3','1/4','-1/4','1/6','-1/6'].map(Q.parse);
 for(let test=0;test<400;test++){
  const d=2+test%2,g=1+test%4;
  const options=Array.from({length:g},()=>Array.from({length:1+pick(4)},()=>Array.from({length:d},()=>values[pick(values.length)])));
  const flat=(c:Q[][])=>Array.from({length:d},(_,i)=>c.map(v=>v[i])).flat();let best:Q[][]|undefined;
  const visit=(chosen:Q[][])=>{if(chosen.length===g){if(!best||compareIntrinsicAssignments(flat(chosen),flat(best))<0)best=[...chosen];return;}for(const v of options[chosen.length])visit([...chosen,v]);};
  visit([]);const got=bestIntrinsicAssignment(options);expect(got).not.toBeNull();
  expect(compareIntrinsicAssignments(flat(got!),flat(best!))).toBe(0);
 }
 expect(bestIntrinsicAssignment([[]])).toBeNull();
 expect(bestIntrinsicAssignment([])).toEqual([]);
});
