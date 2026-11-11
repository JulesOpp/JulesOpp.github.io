# Superspace Subgroups (k-index 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A CELLSUB-style page: paste a (3+1)D superspace group, give k-index = 1, get every subgroup (up to conjugacy in G) with its index, t-index, k-index, point group, basic space group and superspace symbol, e.g. `P-3c1(0,0,g)000` → `P321(0,0,g)000`, index 2.

**Architecture:** A new pure-TS module set in `tools/ssg-finder/src/` (subgroup enumeration, phase-shift symbol, standard setting), a new UI entry bundled next to `ssgf.js`, and a new static page `blog/superspace-subgroups.html`. Same paste-mode input and licence rule as the finder. The gitignored BYU table in `local/data/` (775 groups with d = 1) is the test oracle, through a Python reference that decides superspace-group equivalence by brute force.

**Tech Stack:** TypeScript + vitest + esbuild + moyo-wasm (existing); Python 3 (`/opt/anaconda3/bin/python3 -I`, spglib, numpy, sympy) for the reference.

**Spec:** this plan (no separate spec). Background: `CONTRACTS.md`, `~/dials/In8IrSe9_Xtal4_FullRun/ssg_finder_plan/PLAN.md`.

## Global Constraints

- Licence rule D1 (CONTRACTS.md): nothing derived from BYU ISO(3+d)D is committed or published. Oracle fixtures go in `local/fixtures/subgroups/`; tests that need them skip when `local/` is absent (pattern: `test/golden.test.ts`).
- Scope: d = 1 and k-index = 1 only. The code must refuse d > 1 and k > 1 with a clear message, not give wrong answers.
- Exact arithmetic only (`Q`, bigint `intmat`); no floating point in group theory.
- Reuse `ssg.ts` (`GroupRecord`, `SOp`, `buildOps`, `splitOp`, `opMul`), `paste.ts`, `identify.ts`, `symbol.ts`. Do not change their behaviour; the finder's tests must stay green.
- Phase-shift letters: `0`=0, `s`=1/2, `t`=1/3, `q`=1/4, `h`=1/6 (BYU / IT C convention).
- Commands: `npm test`, `npm run typecheck`, `npm run build` in `tools/ssg-finder/`; Python tests with `/opt/anaconda3/bin/python3 -I -m pytest`.

## The group theory (k-index 1)

Let G be a (3+1)D superspace group with superspace lattice Σ (3D lattice L plus q) and point group K_s = G/Σ of pairs (R, ε), with Rᵀq ≡ εq mod L*.

1. **K_s ≅ K.** For incommensurate q, (1, −1) ∉ K_s (it needs 2q ∈ L*), so (R, ε) ↦ R is injective and K_s is isomorphic to the point group K of the basic space group. ε is a function of R.
2. **k-index 1 = translationengleiche.** A subgroup H ≤ G with the same Σ is the full preimage of a subgroup P ≤ K: H = {g ∈ G : R(g) ∈ P}. So the subgroups with k-index 1 correspond one-to-one to subgroups of K. Index = t-index = |K| / |P|, k-index = 1.
3. **Conjugacy.** Σ is normal and lies in every such H, so H₁, H₂ are conjugate in G iff P₁, P₂ are conjugate in K. List one row per K-conjugacy class, with the class size.
4. **q is unchanged.** H has the same q. Its *type* is the BYU group whose q family contains this q: e.g. the P-1 subgroup of P-3c1(0,0,g)000 has type P-1(a,b,g)0 with q = (0,0,g) as a special member. Show both the type symbol and the actual q in the type's setting.
5. **Phase shifts.** The symbol letter for a position is the internal translation of the representative operation with ε = +1, after the best origin choice. In H fewer operations constrain the origin, so letters can become 0 that were not in G (and the same type can arise twice from different conjugacy classes). The exact rule is fixed in Task 3 against the oracle, not derived by hand.
6. **Sanity check, P-3c1(0,0,g)000 (165.1.24.1).** K = -3m1 has 10 subgroup classes. Expected basic groups: -3m1→165, 3m1→158, 321→150, -3→147, 3→143, 2/m→15, 2→5, m→9, -1→2, 1→1. Index 1,2,2,2,4,3,6,6,6,12.

**Equivalence (the oracle).** Two (3+1)D groups A, B with the same basic space-group number are equivalent iff there is a superspace affine map x' = Px + p, x₄' = m·x + ηx₄ + p₄ (P rational with P·L_A = L_B, η = ±1, m integer row) that sends the operations of A onto those of B modulo Σ_B and sends q_A into the q family of B (mod L_B*). Search: standardise both basic groups with spglib, enumerate P from the linear parts of the affine normaliser of the standard setting (integer matrices, entries in {−2..2}, det ±1, normalising K and preserving the centred lattice), and for each P solve η, m from q and then (p, p₄) from the translation congruences with Smith normal form. This is slow but only runs offline.

## File structure

- Create `reference/ssg_equiv.py`: `equivalent(recA, qA, recB) -> Transform | None`.
- Create `reference/subgroups.py`: t-subgroup enumeration (Python twin of `src/subgroups.ts`).
- Create `reference/make_subgroup_fixtures.py`: writes `local/fixtures/subgroups/<label>.json`.
- Create `reference/test_subgroups.py`.
- Create `src/subgroups.ts`: point group, subgroup classes, preimages.
- Create `src/ssgsymbol.ts`: superspace symbol in given axes (basic HM + q + phase letters).
- Create `src/standardise.ts`: move H to the BYU-convention setting and name it.
- Create `src/subfinder.ts`: orchestration and filters (pure, no DOM).
- Create `src/ui/subgroups.ts`: page UI.
- Create `blog/superspace-subgroups.html`.
- Create tests `test/subgroups.test.ts`, `test/ssgsymbol.test.ts`, `test/standardise.test.ts`, `test/subfinder.test.ts`.
- Modify `build.mjs` (second entry → `assets/blog/ssg-finder/ssgsub.js`), `CONTRACTS.md` (§5 fixture format), `blog/index.html` (link).

## Review Focus

1. Groups with ε = −1 operations and q with a rational part, e.g. `P2_1/b(1/2,0,g)00` (14.1.6.3): the internal origin and M ≠ 0 rows change τ₄ under 3D origin shifts. Covered by the all-775 symbol test in Task 3.
2. Centred basic groups and BYU's non-standard centrings (`B2/b`, `X` lattices): subgroups whose centring is not P must keep the centring vectors in their operator lists. Test in Task 2 on 15.1.7.3 `B2/b(0,0,g)s0`.
3. Subgroups that fall to a lower crystal family (hexagonal → monoclinic): the parent-axes symbol is not the type symbol. Covered by Task 4's "hide the setting" test.
4. Two conjugacy classes of the same type (same symbol, different classes): must stay two rows, as in CELLSUB. Test in Task 5.
5. Pasted page with d = 2 or 3: refuse with "only d = 1 for now". Test in Task 5.

---

### Task 1: Python oracle — equivalence and subgroup fixtures (local only)

**Files:**
- Create: `reference/ssg_equiv.py`, `reference/subgroups.py`, `reference/make_subgroup_fixtures.py`
- Test: `reference/test_subgroups.py`
- Modify: `CONTRACTS.md` (add §5)

**Interfaces:**
- Produces: `local/fixtures/subgroups/<label>.json`:
  ```jsonc
  { "label": "165.1.24.1", "sym": "P-3c1(0,0,g)000", "order": 12,
    "classes": [  // sorted by (index, byu label, first op indices)
      { "pg_ops": [0, 3, 7], // indices into the parent's canonical op list (CONTRACTS §1) of one op per coset of Σ
        "index": 2, "class_size": 1,
        "byu": "150.1.24.1", "byu_sym": "P321(0,0,g)000",
        "P": [["1","0","0"],...], "p": ["0","0","0"], "eta": 1, "m": [0,0,0], "p4": "0" } ] }
  ```
- Produces (Python): `equivalent(rec_a: dict, q_a: list[list[Fraction]], rec_b: dict) -> dict | None`, `t_subgroups(rec: dict) -> list[dict]`.

- [ ] **Step 1: Write the failing tests** in `reference/test_subgroups.py`, all marked skip if `local/data` is missing:
  - `test_self_equivalent`: for every d = 1 record r, `equivalent(r, q_generic(r), r)` is not None.
  - `test_distinct_not_equivalent`: for every pair of different d = 1 records with the same basic SG number, `equivalent(a, q_generic(a), b)` is None.
  - `test_p3c1_classes`: `t_subgroups(rec("165.1.24.1"))` has 10 classes; basic SG numbers sorted equal `[1,2,5,9,15,143,147,150,158,165]`; the class with basic SG 150 matches `"150.1.24.1"`.
  - `test_every_subgroup_identified`: for every d = 1 record, every class matches exactly one BYU label.
  (`q_generic` substitutes 3001/10007, 4723/10007, 7177/10007 for the free parameters, as in CONTRACTS §3.1.)
- [ ] **Step 2: Run** `/opt/anaconda3/bin/python3 -I -m pytest reference/test_subgroups.py -v` — expect import errors.
- [ ] **Step 3: Implement `equivalent`** as in "Equivalence (the oracle)" above, reusing `reference/intmat.py` (SNF) and `reference/exact.py` (op closure).
- [ ] **Step 4: Implement `t_subgroups`**: subgroups of K by closing every set of ≤ 3 point operations, dedupe, group by K-conjugacy, take preimages; identify each against the candidates with the same basic SG number.
- [ ] **Step 5: Run the tests.** Expected: all pass. If `test_distinct_not_equivalent` fails, the P search is too loose; if `test_self_equivalent` fails, it is too narrow. Record runtime of the all-775 test in the commit message.
- [ ] **Step 6: Write fixtures** with `make_subgroup_fixtures.py` for all 775 groups; add CONTRACTS.md §5 describing the format above.
- [ ] **Step 7: Commit** (code and CONTRACTS only; `local/` is ignored).

### Task 2: Subgroup enumeration in TS

**Files:**
- Create: `src/subgroups.ts`
- Test: `test/subgroups.test.ts`

**Interfaces:**
- Consumes: `buildOps`, `splitOp`, `opMul`, `GroupRecord`, `SOp` from `ssg.ts`.
- Produces:
  ```ts
  export interface PointOp { R: number[][]; eps: number }           // eps = ±1
  export interface SubgroupClass {
    pointOps: PointOp[];      // P, identity first
    ops: SOp[];               // canonical list of H (CONTRACTS §1 order), centrings included
    index: number;            // |K| / |P|
    classSize: number;        // number of conjugates in K
    maximal: boolean;         // no subgroup class strictly between P and K
  }
  export function pointGroup(ops: SOp[], d: number): PointOp[];
  export function tSubgroups(rec: GroupRecord): SubgroupClass[];  // sorted by (index, |P| desc, ops of P)
  ```

- [ ] **Step 1: Write failing tests:** (a) a hand-written P-3c1(0,0,g)000 record (operators from IT A No. 165 plus ε, written in the test, not taken from BYU) gives 10 classes with indices `[1,2,2,2,3,4,6,6,6,12]` sorted, class sizes summing to 16 (2, m and 2/m have 3 conjugates each, every other class 1); (b) every `ops` list is closed under `opMul`; (c) the order of `ops` is `index` times smaller than the parent's; (d) local-only: for all 775 records, the number of classes and their `index` multiset equal the Task 1 fixture; (e) local-only: 15.1.7.3 `B2/b(0,0,g)s0` — every class keeps the centring translation `1/2,0,1/2,0`.
- [ ] **Step 2: Run** `npm test -- subgroups` — FAIL (module missing).
- [ ] **Step 3: Implement.** Subgroups of K by closure of all ≤ 3-element generating sets (|K| ≤ 48, so this is cheap); conjugacy by acting with every element of K; preimage = parent ops whose point part lies in P.
- [ ] **Step 4: Run** `npm test` — all pass, finder tests included.
- [ ] **Step 5: Commit.**

### Task 3: Superspace symbol in given axes (phase letters)

**Files:**
- Create: `src/ssgsymbol.ts`
- Test: `test/ssgsymbol.test.ts`

**Interfaces:**
- Consumes: `settingSymbol`, `CrystalSystem` from `symbol.ts`; `Op3` from `identify.ts`.
- Produces:
  ```ts
  export interface SsgSymbol { basic: string; q: string; shifts: string; full: string } // full e.g. "P321(0,0,g)t00"
  export function ssgSymbolInAxes(ops: SOp[], qFamily: string[], system: CrystalSystem): SsgSymbol | null;
  ```
  `basic` is the BYU-style compact HM (no spaces, `P2_1/b`, `P-3c1`); `shifts` has one letter per symbol position (BYU omits positions that hold only 1, as in `P321` → 3 letters, `P-3` → 1 letter).

- [ ] **Step 1: Write the failing test (local-only, the oracle):** for every one of the 775 d = 1 records r, `ssgSymbolInAxes(buildOps(r), r.q, systemOf(r)).full === r.sym`. Print the first 20 mismatches grouped by basic SG.
- [ ] **Step 2: Run** — FAIL (module missing).
- [ ] **Step 3: Implement.** For each symbol position take the representative operation `settingSymbol` chose; letter = its τ₄ reduced mod 1 after choosing the internal origin p₄ (and 3D origin shifts that keep the basic HM symbol) to make the letter string smallest in the order 0 < s < t < q < h, position by position. Operations with ε = −1 always give 0. If the oracle disagrees on a family, fix the rule, not the test; any truly irreducible exception goes in a commented table keyed by basic SG number and q form (no BYU labels in source).
- [ ] **Step 4: Run** until 775/775 pass. Also a public test: the hand-written P-3c1 record → `"P-3c1(0,0,g)000"`.
- [ ] **Step 5: Commit.**

### Task 4: Standard (BYU-convention) setting and type name

**Files:**
- Create: `src/standardise.ts`
- Test: `test/standardise.test.ts`

**Interfaces:**
- Consumes: `identify` (moyo standardisation), `ssgSymbolInAxes`.
- Produces:
  ```ts
  export interface Standardised {
    sgNumber: number;
    symbol: SsgSymbol;            // type symbol, q written as the family (a, b, g letters)
    qActual: string;              // the subgroup's own q in that setting, e.g. "0,0,g" or "-g,g,0"
    P: string[][]; p: string[];   // 3D: x_std = P x_parent + p
    eta: 1 | -1; m: number[]; p4: string;
  }
  export function standardise(ops: SOp[], Q: string[]): Standardised | { ok: false; reason: string };
  ```
- [ ] **Step 1: Write failing tests (local-only):** (a) for each of the 775 records, apply 3 seeded random setting changes (axis permutation allowed by the crystal system, integer origin-preserving shear for triclinic/monoclinic, rational origin shift, internal shift p₄, η = −1) to the operators and q, then `standardise` returns the record's own `sym`; (b) for every Task 1 fixture class, `standardise(cls.ops).symbol.full === byu_sym`. Public test: the hand-written P-3c1 → its C2/c class gives a monoclinic `B2/b(...)`-style symbol with unique axis c.
- [ ] **Step 2: Run** — FAIL.
- [ ] **Step 3: Implement.** Go to moyo's standard setting of the basic group; enumerate the finite list of setting changes of that crystal system that BYU uses (monoclinic unique axis c with cell choices and B/A/I centrings; the 6 orthorhombic axis permutations; trigonal obverse/hexagonal; identity otherwise) times the normaliser linear parts from Task 1; for each, transform q, reduce its rational part mod the reciprocal lattice (centring included) and match it to a BYU q form for that Bravais class (forms listed in source as data such as `"1/2,0,g"`, not labels); keep candidates with a valid q form; pick by the rule that test (a) pins.
- [ ] **Step 4: Run** until both tests pass 100%.
- [ ] **Step 5: Commit.**

### Task 5: Orchestration and filters

**Files:**
- Create: `src/subfinder.ts`
- Test: `test/subfinder.test.ts`

**Interfaces:**
- Consumes: `tSubgroups`, `standardise`, `ssgSymbolInAxes`, `parseByuPage`.
- Produces:
  ```ts
  export interface SubFilter { kIndex: number; onlyK?: boolean; pointGroup?: string; system?: CrystalSystem; centring?: string }
  export interface SubRow {
    index: number; tIndex: number; kIndex: 1; classSize: number; maximal: boolean;
    pointGroup: string;            // HM, e.g. "321"
    sgNumber: number; sgSymbol: string;
    inParentAxes: string | null;   // ssgSymbolInAxes in parent axes, null if the system drops (Review Focus 3)
    type: Standardised;
    generators: string[];          // xyzt operators of H in parent axes, centrings excluded
  }
  export function runSubgroups(raw: string, f: SubFilter): { parent: string; rows: SubRow[]; warnings: string[] };
  ```
  Errors: `SubgroupError("only d = 1 is supported for now")`, `SubgroupError("only k-index 1 is supported for now")`.
- [ ] **Step 1: Write failing tests:** P-3c1 hand record → 10 rows, row with `sgNumber 150` has `type.symbol.full === "P321(0,0,g)000"`, `index 2`, `kIndex 1`; filter `system: "trigonal"` → 6 rows; filter `pointGroup: "321"` → 1 row; `kIndex: 2` throws; a d = 2 page throws; local-only: one parent where two rows share `type.symbol.full` stays two rows (pick it from the Task 1 fixtures, record the label in a comment).
- [ ] **Step 2: Run** — FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `npm test && npm run typecheck`.
- [ ] **Step 5: Commit.**

### Task 6: Page

**Files:**
- Create: `src/ui/subgroups.ts`, `blog/superspace-subgroups.html`
- Modify: `build.mjs`, `blog/index.html`

**Interfaces:**
- Consumes: `runSubgroups`, `loadMoyo`.

- [ ] **Step 1: Build the page** in the style of `blog/superspace-group-finder.html`: paste box; k-index field (1 only, others disabled with a note); the CELLSUB optional filters (only this k-index, point group, crystal system, centring); results table with columns Index | t-index | k-index | Point group | Basic SG | Superspace group (type) | q | Conjugates | Maximal | Transformation; each row expands to show the generators and has a "Send to finder" link that opens `superspace-group-finder.html` with the subgroup's operators in parent axes (via `sessionStorage`, wrapped in try/catch; the finder needs a small reader for it).
- [ ] **Step 2: Add the second esbuild entry** → `assets/blog/ssg-finder/ssgsub.js`; `npm run build`.
- [ ] **Step 3: Verify in the browser preview:** paste the cached 165.1.24.1 page; the table shows 10 rows and the P321(0,0,g)000 row; filters change the row count as in Task 5; no console errors; phone width has no horizontal page scroll.
- [ ] **Step 4: Link** the page from `blog/index.html` and from the finder page.
- [ ] **Step 5: Commit.**

## Later (not in this plan)

- **k-index > 1.** Subgroups with a K-invariant sublattice Σ' ⊂ Σ of index k: 3D supercells (q's rational part changes relative to the new L*) and internal sublattices (index in the q direction, e.g. q → 2q only matters commensurately). Needs Σ'-enumeration (HNFs of index k invariant under K) and then t-subgroups of each; CELLSUB's "general" subgroups.
- **d = 2, 3.** Theory step 1 still holds for generic q, but symbols carry several q's and E permutes them.
- **Jana2020 cross-check** once you have a machine for it.

## Decisions for you

1. **Public naming.** With D1, the public page cannot show BYU numbers. It can show the type *symbol* only if Task 4 produces it algorithmically (the plan's route). If Task 4 stalls below 100%, the fallback is to show the parent-axes symbol and the basic space group, and keep BYU labels in a local-only mode.
2. **Separate page or a tab** on the finder page. The plan assumes a separate page sharing the code.
