import { existsSync, readdirSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { beforeAll, describe, expect, it, vi } from "vitest";

vi.setConfig({ testTimeout: 300_000 });
import { Q } from "../src/rational";
import { fromColMajor, getMoyo, loadMoyo } from "../src/moyo";
import { identify, numToQ, type Op3 } from "../src/identify";
import { crystalSystemOf, itShortForm, settingSymbol } from "../src/symbol";

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");
const FIX = join(ROOT, "local", "fixtures", "identify");
const haveFixtures = existsSync(FIX) && readdirSync(FIX).some((f) => f.endsWith(".json"));

beforeAll(async () => {
  await loadMoyo();
});

// ------------------------------------------------------------------ helpers

type MoyoOp = { rotation: number[]; translation: number[] };
const toOp3 = (o: MoyoOp): Op3 => ({ R: fromColMajor(o.rotation), t: o.translation.map(numToQ) });

function standardOps(setting: { type: "Standard" } | { type: "HallNumber"; value: number }, n: number): Op3[] {
  return (getMoyo().operations_from_number(n, setting as never, false) as MoyoOp[]).map(toOp3);
}
function hallOps(hall: number): Op3[] {
  const e = getMoyo().hall_symbol_entry(hall);
  return standardOps({ type: "HallNumber", value: hall }, e.number);
}

const ORTHO = [[5.1, 0, 0], [0, 6.3, 0], [0, 0, 7.7]];
const HEX = [[4, 0, 0], [-2, 2 * Math.sqrt(3), 0], [0, 0, 6.5]];
const CUBIC = [[5, 0, 0], [0, 5, 0], [0, 0, 5]];

const inv3 = (A: number[][]): number[][] => {
  const [a, b, c] = A;
  const det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
  const cof = [
    [b[1] * c[2] - b[2] * c[1], a[2] * c[1] - a[1] * c[2], a[1] * b[2] - a[2] * b[1]],
    [b[2] * c[0] - b[0] * c[2], a[0] * c[2] - a[2] * c[0], a[2] * b[0] - a[0] * b[2]],
    [b[0] * c[1] - b[1] * c[0], a[1] * c[0] - a[0] * c[1], a[0] * b[1] - a[1] * b[0]],
  ];
  return cof.map((r) => r.map((x) => x / det));
};
const mm = (A: number[][], B: number[][]) => A.map((r) => B[0].map((_, j) => r.reduce((s, a, k) => s + a * B[k][j], 0)));
const mv = (A: number[][], v: number[]) => A.map((r) => r.reduce((s, a, k) => s + a * v[k], 0));
const near = (x: number) => Math.abs(x - Math.round(x)) < 1e-6;

/** Check x_conv = P x_std + p maps our ops onto ops of the reported Hall setting (mod lattice). */
function checkStdTransform(ops: Op3[], res: Extract<ReturnType<typeof identify>, { ok: true }>) {
  const P = res.stdLinear, p = res.stdOriginShift, Pi = inv3(P);
  const std = standardOps({ type: "HallNumber", value: res.hallNumber }, res.number);
  for (const o of ops) {
    const W = mm(mm(Pi, o.R), P);
    const w = mv(Pi, mv(o.R, p).map((x, i) => x + o.t[i].toNumber() - p[i]));
    expect(W.every((r) => r.every(near))).toBe(true);
    const hit = std.some(
      (s) =>
        s.R.every((r, i) => r.every((x, j) => Math.abs(x - W[i][j]) < 1e-6)) &&
        s.t.every((x, i) => near(x.toNumber() - w[i])),
    );
    expect(hit).toBe(true);
  }
}

function basisFor(n: number): number[][] {
  if (n >= 143 && n <= 194) return HEX;
  if (n >= 195) return CUBIC;
  return ORTHO;
}

// ------------------------------------------------------------------ synthetic (no BYU data)

describe("identify: synthetic groups", () => {
  it("all 230 standard settings (operations_from_number) are identified as themselves", () => {
    const bad: string[] = [];
    for (let n = 1; n <= 230; n++) {
      const ops = standardOps({ type: "Standard" }, n);
      const r = identify(ops, basisFor(n));
      if (!r.ok || r.number !== n) bad.push(`${n}: ${r.ok ? r.number : r.reason}`);
      else checkStdTransform(ops, r);
    }
    expect(bad).toEqual([]);
  });

  it("all 230 standard settings with a non-invariant (unit) basis: metric averaging", () => {
    const bad: string[] = [];
    for (let n = 1; n <= 230; n++) {
      const r = identify(standardOps({ type: "Standard" }, n), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
      if (!r.ok || r.number !== n) bad.push(`${n}: ${r.ok ? r.number : r.reason}`);
    }
    expect(bad).toEqual([]);
  });

  it("all 530 Hall settings are identified with the right number", () => {
    const bad: string[] = [];
    for (let hall = 1; hall <= 530; hall++) {
      const n = getMoyo().hall_symbol_entry(hall).number;
      const r = identify(hallOps(hall), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
      if (!r.ok || r.number !== n) bad.push(`hall ${hall} (#${n}): ${r.ok ? r.number : r.reason}`);
    }
    expect(bad).toEqual([]);
  });

  it("hexagonal P6_3/mmc in hexagonal axes", () => {
    const ops = standardOps({ type: "Standard" }, 194);
    const r = identify(ops, HEX);
    expect(r.ok && r.number).toBe(194);
    expect(r.ok && r.nOps).toBe(24);
    expect(settingSymbol(ops, "hexagonal")).toBe("P 6_3/m 2/m 2/c");
    expect(itShortForm("P 6_3/m 2/m 2/c")).toBe("P63/mmc");
  });

  it("trigonal R-3m in hexagonal (obverse) and rhombohedral axes", () => {
    const hexOps = standardOps({ type: "Standard" }, 166);
    const r = identify(hexOps, HEX);
    expect(r.ok && r.number).toBe(166);
    expect(r.ok && r.nCentring).toBe(3);
    expect(settingSymbol(hexOps, "trigonal")).toBe("R -3 2/m 1");
    expect(itShortForm("R -3 2/m 1")).toBe("R-3m");
    expect(itShortForm("F 4/m -3 2/m")).toBe("Fm-3m");
    expect(itShortForm("P 4_2/m 2_1/n 2/m")).toBe("P42/mnm");
    const rhOps = (getMoyo().operations_from_number(166, { type: "Spglib" }, true) as MoyoOp[]).map(toOp3);
    const rh = identify(rhOps, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    expect(rh.ok && rh.number).toBe(166);
  });

  it("trigonal P3_1 21 / P3_2 21 enantiomorphs, and a left-handed basis swaps them", () => {
    const o152 = standardOps({ type: "Standard" }, 152);
    expect(identify(o152, HEX)).toMatchObject({ ok: true, number: 152 });
    expect(settingSymbol(o152, "trigonal")).toBe("P 3_1 2 1");
    expect(settingSymbol(standardOps({ type: "Standard" }, 154), "trigonal")).toBe("P 3_2 2 1");
    const leftHex = [HEX[0], HEX[1], HEX[2].map((x) => -x)];
    expect(identify(o152, leftHex)).toMatchObject({ ok: true, number: 154 });
  });

  it("tetragonal and cubic setting symbols", () => {
    expect(settingSymbol(standardOps({ type: "Standard" }, 136), "tetragonal")).toBe("P 4_2/m 2_1/n 2/m");
    expect(settingSymbol(standardOps({ type: "Standard" }, 122), "tetragonal")).toBe("I -4 2 d");
    expect(settingSymbol(standardOps({ type: "Standard" }, 225), "cubic")).toBe("F 4/m -3 2/m");
    expect(settingSymbol(standardOps({ type: "Standard" }, 230), "cubic")).toBe("I 4_1/a -3 2/d");
    expect(settingSymbol(standardOps({ type: "Standard" }, 2), "triclinic")).toBe("P -1");
    expect(settingSymbol(standardOps({ type: "Standard" }, 62), "orthorhombic")).toBe("P 2_1/n 2_1/m 2_1/a");
    expect(settingSymbol(standardOps({ type: "Standard" }, 64), "orthorhombic")).toBe("C 2/m 2/c 2_1/e");
    // mm2 subgroup of P-42m keeping [110]: the mirrors sit in the tertiary (diagonal) position
    const sub = standardOps({ type: "Standard" }, 111).filter((o) => {
      const v = [o.R[0][0] + o.R[0][1], o.R[1][0] + o.R[1][1], o.R[2][0] + o.R[2][1]];
      return v[2] === 0 && Math.abs(v[0]) === 1 && v[0] === v[1];
    });
    expect(sub.length).toBe(4);
    expect(identify(sub, [[4, 0, 0], [0, 4, 0], [0, 0, 6]])).toMatchObject({ ok: true, number: 35 }); // Cmm2: diagonal mirrors in a P square lattice
    expect(settingSymbol(sub, "tetragonal")).toBe("P 2 1 m");
    // orthorhombic subgroup of P4/mmm with x and y no longer equivalent
    const ortho = standardOps({ type: "Standard" }, 127).filter((o) => o.R[0][1] === 0);
    expect(settingSymbol(ortho, "tetragonal")).toBe("P 2/m [2_1/b;2_1/a] 1");
  });

  it("returns null when an axis is not a parental direction", () => {
    // tetragonal 4-fold along z given to an orthorhombic parent is fine (z is a direction) ...
    expect(settingSymbol(standardOps({ type: "Standard" }, 75), "orthorhombic")).not.toBeNull();
    // ... but a 3-fold along [111] is not
    expect(settingSymbol(standardOps({ type: "Standard" }, 195), "orthorhombic")).toBeNull();
  });

  function supercell(n: number, N: [number, number, number]): Op3[] {
    // ops of the standard conventional cell expressed in an N1 x N2 x N3 supercell:
    // every old lattice/centring translation becomes a centring of the supercell
    const base = standardOps({ type: "Standard" }, n);
    const out: Op3[] = [];
    for (const o of base)
      for (let i = 0; i < N[0]; i++)
        for (let j = 0; j < N[1]; j++)
          for (let k = 0; k < N[2]; k++) {
            const c = [i, j, k];
            out.push({ R: o.R, t: o.t.map((x, a) => x.add(Q.of(BigInt(c[a]))).div(Q.of(BigInt(N[a]))).mod1()) });
          }
    return out;
  }

  it("supercell of index 4096 (Pnma, 16x16x16, 32768 operations)", () => {
    const ops = supercell(62, [16, 16, 16]);
    expect(ops.length).toBe(8 * 4096);
    const r = identify(ops, [[16 * 5.1, 0, 0], [0, 16 * 6.3, 0], [0, 0, 16 * 7.7]]);
    expect(r).toMatchObject({ ok: true, number: 62, nOps: 8, nCentring: 4096 });
    const sym = settingSymbol(ops, "orthorhombic")!;
    expect(sym.startsWith("{")).toBe(true); // 4095 centrings: explicit list
    expect(sym.split(" ").slice(1)).toEqual(["2_1/n", "2_1/m", "2_1/a"]);
  });

  it("supercell of index 4096 with a centred lattice (Cmce in 8x8x32)", () => {
    const ops = supercell(64, [8, 8, 32]);
    const r = identify(ops, [[8 * 5.1, 0, 0], [0, 8 * 6.3, 0], [0, 0, 32 * 7.7]]);
    expect(r).toMatchObject({ ok: true, number: 64, nOps: 8, nCentring: 4096 });
    if (r.ok) checkStdTransform(ops.slice(0, 64), r);
  });

  it("fails loudly on an inconsistent operation list", () => {
    const ops = standardOps({ type: "Standard" }, 62);
    const broken = ops.map((o, i) => (i === 3 ? { R: o.R, t: [Q.of(1n, 3n), o.t[1], o.t[2]] } : o));
    expect(identify(broken, ORTHO).ok).toBe(false);
    expect(identify(ops.slice(0, 3), ORTHO).ok).toBe(false); // not closed
  });
});

// ------------------------------------------------------------------ BYU-derived fixtures

interface FixClass {
  sg_number: number;
  sg_symbol: string;
  setting: string;
  ops3: { R: string[][]; t: string[] }[];
}
interface Fixture {
  label: string;
  params: Record<string, string>;
  N: number[];
  cell_len: number[];
  classes: FixClass[];
  it?: { name: string; expected: string[] };
}

const fixtures: Fixture[] = haveFixtures
  ? readdirSync(FIX)
      .filter((f) => f.endsWith(".json"))
      .sort()
      .map((f) => JSON.parse(readFileSync(join(FIX, f), "utf8")))
  : [];
const fixOps = (c: FixClass): Op3[] =>
  c.ops3.map((o) => ({ R: o.R.map((r) => r.map((x) => Q.parse(x).toNumber())), t: o.t.map(Q.parse) }));

/** Normalise a setting symbol for comparison with the prototype's hm_in_setting. */
function norm(s: string): string {
  const [cen, ...ax] = s.split(/\s+/);
  // the prototype prints "X 1 1 1" for both P1 and P-1 (it does not detect -1)
  if (ax.every((a) => a === "1" || a === "-1")) return `${cen} 1`;
  return [cen, ...ax].join(" ");
}

describe.skipIf(!haveFixtures)("identify: BYU-derived fixtures (local/fixtures/identify; skipped if absent, see reference/make_identify_fixtures.py)", () => {
  for (const fx of fixtures) {
    const parent = crystalSystemOf(Number(fx.label.split(".")[0]));
    const tag = `${fx.label} ${JSON.stringify(fx.params)}`;
    it(`${tag}: numbers match spglib, setting symbols match the prototype`, () => {
      const basis = [0, 1, 2].map((i) => [0, 1, 2].map((j) => (i === j ? fx.N[i] * fx.cell_len[i] : 0)));
      for (const c of fx.classes) {
        const ops = fixOps(c);
        const r = identify(ops, basis);
        expect(r.ok ? r.number : r.reason, `${tag} class ${c.sg_symbol}`).toBe(c.sg_number);
        const sym = settingSymbol(ops, parent);
        expect(sym && norm(sym), `${tag} class ${c.sg_symbol}`).toBe(norm(c.setting));
      }
    });
    if (fx.it) {
      it(`${tag}: IT finder short forms (${fx.it.name})`, () => {
        const got = new Set(fx.classes.map((c) => itShortForm(settingSymbol(fixOps(c), parent) ?? "?")));
        expect([...got].sort()).toEqual([...fx.it!.expected].sort());
      });
    }
  }
  it("covers the 19 IT-finder rows", () => {
    expect(fixtures.filter((f) => f.it).length).toBe(19);
  });
});
