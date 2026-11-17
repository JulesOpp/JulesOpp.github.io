# Superspace Subgroups — Finishing Plan (d = 2, 3; k-index > 1; minors)

## Accepted page scope — 2026-11-14

The page scope is frozen at **d = 1, k ≤ 3; d = 2, k ≤ 2; d = 3, k = 1**. All supported caps are fully validated. The completed d = 2, k = 2 reference sweep covers all 3,338 parents and 768,066 classes, each identified with exactly one reference type, with zero unidentified or multiply identified classes.

Future work, **not acceptance blockers**:

- Full d = 1, k = 4 validation: approximately 10 hours projected wall time (10.22 h from the saved 30-parent sample); the sample is not full coverage.
- d = 3, k > 1: the pre-fix, censored k = 2 projection was at least 350 hours (1,273,943 seconds on four workers). It is infeasible for the reference oracle in the present scope; post-fix target timings do not establish a replacement full projection.
- Higher-dimension k = 3 validation: unmeasured. Earlier extrapolations are not measured validation evidence.

The 2026-11-14 Phase D decision supersedes the former naming gate: BYU string agreement is not required. Exact identification remains authoritative. The current work implements that decision with TDD and no commits before review. Nothing under local/ may be staged or published. Do not merge or push.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The finished Superspace Subgroups page. It covers (3+1)D, (3+2)D and (3+3)D parents and k-index 1 up to a small maximum k, with one row per G-conjugacy class and a BYU-style type for each row. The page carries no "in progress" note and shows its proper index entry.

**Architecture:** The same module set as the k-index-1 work (`src/subgroups.ts`, `ssgsymbol.ts`, `standardise.ts`, `subfinder.ts`, `ui/subgroups.ts`, `ui/subflow.ts`), generalised from d = 1 to any d ≤ 3. It gains two new modules, `src/sublattices.ts` and `src/ksubgroups.ts`, for k > 1. A subgroup with k > 1 is named by rebasing it onto its own lattice Σ′ and then calling the k = 1 `standardise`. The Python oracle (`reference/ssg_equiv.py`, `subgroups.py`, a new `ksubgroups.py`) is generalised the same way. It is validated against the gitignored BYU table in `local/data`: 775 (d = 1), **3338 (d = 2)** and **12584 (d = 3)** records, counted 2026-11-11.

**Tech Stack:** TypeScript + vitest + esbuild + moyo-wasm; Python 3 (`/opt/anaconda3/bin/python3 -I`, spglib, numpy, sympy).

**Spec:** this plan, together with the first plan `plans/2026-11-10-ssg-subgroups-k1.md` (its group theory §1–6 still holds) and `CONTRACTS.md` §1, §2 and §5.

## Global Constraints

- **Licence rule D1** (CONTRACTS.md): nothing derived from BYU ISO(3+d)D is committed or published. Fixtures go under `local/fixtures/`. Every test that reads `local/` uses `describe.skipIf(!have)` (TS) or `pytest.mark.skipif` (Python). Source may contain *rules* fitted to BYU, as `rank` in `standardise.ts` already does. It may not contain BYU labels or tables keyed by them (see Decision 1 and Decision 2).
- **CPU:** never use more than 4 cores. Python pools use `min(requested, 4)`, and vitest keeps `maxWorkers: 4`. Every worker sets `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = MKL_NUM_THREADS = VECLIB_MAXIMUM_THREADS = 1` before numpy is imported.
- **Python:** always run it as `/opt/anaconda3/bin/python3 -I`.
- **Arithmetic:** exact only (`Q`, bigint `intmat`, `fractions.Fraction`). Floats may appear only as a documented prefilter whose result is re-checked exactly.
- **Don't touch** `assets/stl/` (untracked) or the finder's behaviour. The finder tests must stay green after every task.
- **Commands** (run in `tools/ssg-finder/`): `npm test`, `npm run typecheck`, `npm run build`, and `/opt/anaconda3/bin/python3 -I -m pytest reference/...`. Run any sweep that may exceed 2 minutes in the background, with `--resume`.
- **Phase letters:** `0`=0, `s`=1/2, `t`=1/3, `q`=1/4, `h`=1/6. At d ≥ 2 BYU also writes signed letters (`-h`, `-t`, `-q`; 28 records, e.g. `P622(0,0,g1)-h00(0,0,g2)000`), so the sign is part of the symbol there.
- **Fixture paths:** k = 1 fixtures go in `local/fixtures/subgroups/d<d>/<label>.json`; the existing d = 1 files move into `d1/`. k > 1 fixtures go in `local/fixtures/ksubgroups/d<d>/k<k>/<label>.json`.

## Group theory for d = 2, 3 (k-index 1)

Notation: n = 3 + d. The superspace lattice is Σ = Z^n plus the record's centrings. An operation is g = (W, w) with W = [[R, 0], [M, E]], where E is d×d. Records satisfy M = Q R − E Q exactly, with Q the d×3 matrix of q's at generic parameters (CONTRACTS §2).

1. **K_s ≅ K still holds for every d.** Suppose (1, E) ∈ K_s with E ≠ 1. Then (E − 1)Q ≡ 0 mod L*, so Σ_j c_j q_j ∈ L* with c ≠ 0 an integer row. That contradicts rank M* = 3 + d, which is the definition of a (3+d)D group. So E is a function of R, and §1–3 of the first plan carry over word for word: the k = 1 subgroups are the preimages of subgroups P ≤ K, and they are classified up to K-conjugacy.
   *Caveat:* the oracle evaluates at rational "generic" parameters, so a rank drop is possible in principle. Task 4 checks that rank M* = 3 + d at the generic point for every record.
2. **How E acts.** R ↦ E(R) is the integral representation of K on M*/L* ≅ Z^d, written in the record's q-basis. The data shows three patterns:
   - **stars:** E permutes the q's, possibly with signs or with integer sums. Examples are the hexagonal pair `(a,0,1/3)(-a,a,1/3)` and the cubic triple `(a,a,a)(a,-a,-a)(-a,a,-a)`.
   - **parallel q's:** the q's are parallel and independent, e.g. `(0,0,g1)(0,0,g2)`, and E = ±I.
   - **mixtures** of the two.
   The subgroup H keeps M*, so it keeps every q. Restricted to P, though, the representation can split: a star can break into orbits, and a basis that is natural for H need not be the parent's. **Consequence:** the equivalence map must allow E ∈ GL(d, Z), not just η = ±1. That means S = [[P, 0], [M, E]] with q′ = (E Q + M) P⁻¹ (rows). The type's family has at least as many free parameters as the parent's, and the parent's q is a special member of it.
3. **E search in the oracle.** E is not enumerated blindly. For a fixed P it must satisfy two linear conditions over Q:
   - (i) **intertwining:** E · E_A(R) = E_B(P R P⁻¹) · E for every R in K_A;
   - (ii) **family:** each E·Q_k·P⁻¹ (Q_k = coefficient matrix of parameter k) must lie in the span of B's coefficient matrices.

   The integer solutions form a lattice: take its basis with `intmat.hnf` and enumerate combinations with coefficients in −2..2 and det E = ±1. M, p and p4 then follow exactly as for d = 1, by Smith normal form, now with d rows.
4. **Invariant prefilter.** Two cheap invariants are kept equal by every equivalence:
   - the multiset over K of (order R, det R, tr R, tr E(R));
   - nops.

   Candidates must match both before any search (the same idea as the basic SG number at d = 1). This matters because d = 3 has up to 809 records per space group, against 18 at d = 1.
5. **Symbols.** These have the form basic + Σ_i "(q_i)" + letters_i: one letter string per q, with the same positions as at d = 1. A letter for q_i is read from row i of the internal translation after the best internal origin (p4 ∈ Q^d). An operation whose E is not diagonal in the q-basis (a star permutation) can still carry a letter, because the internal origin removes (1 − E)p4, a sublattice of the translations. The rule, including signs, is fitted against all 3338 + 12584 symbols, not derived by hand (Task 7).
6. **Settings.** BYU uses more basic settings at d ≥ 2 than at d = 1. Monoclinic has `Pa`/`Pb`, `P2/a`/`P2/b` and `P2_1/a`/`P2_1/b`; orthorhombic has up to three of the six settings (`P222_1`, `P22_12`, `P2_122`; `A2_122`, `B22_12`, `C222_1`). The q-forms of the (3+2)D and (3+3)D Bravais classes take the place of Q_FORMS (83 and 215 Bravais classes appear in the data; see Decision 2 for where the forms come from). Some q's have integer components that are not reduced, e.g. `1,0,g` in centred lattices.

## Group theory for k-index > 1

A subgroup H ≤ G has point group P ≤ K and lattice Σ′ = H ∩ Σ. Its indices are k = [Σ : Σ′] and t = |K|/|P|, with index = t·k.

1. **Σ′** is a W(P)-invariant sublattice of index k. Since kΣ ⊆ Σ′, the candidates are the index-k subgroups of Σ/kΣ ≅ (Z/k)^n. Enumerate them as lower-HNF matrices of determinant k in the Σ basis (the centrings are folded into the basis first). Counts:
   - n = 4, k = 2: 15
   - n = 4, k = 3: 40
   - n = 6, k = 2: 63
   - n = 6, k = 3: 364

   For each Σ′, keep its stabiliser in K.
2. **Kinds of Σ′.** Each Σ′ is labelled 3D, internal or mixed, and the page shows the label:
   - **3D supercell:** Σ′ ∩ V_I = Σ ∩ V_I. The q's are unchanged, but they are re-read modulo the larger L′*, so a rational part can vanish. Example: q = (1/2,0,g) with a doubled gives the (0,0,g) form.
   - **internal:** the projection of Σ′ onto V_E is all of L. Here q′ = q/k or similar, a different (3+d)D type with the same 3D lattice.
   - **mixed:** neither of the above.
3. **Subgroups for a given (P, Σ′).** These are the complements of Σ/Σ′ in G_P/Σ′, where G_P is the k = 1 preimage of P. For generators g_1..g_r of P (r ≤ 3), try every choice of translations t_i ∈ Σ/Σ′ (k^r choices). Close the set modulo Σ′. Keep the result when the closure has exactly |P| cosets of Σ′, which is the same as H ∩ Σ = Σ′.
4. **Conjugacy in G.** G acts on the pairs (Σ′, H) through its generators: the operations plus a basis of Σ. kΣ acts trivially, so every orbit is finite. The canonical key is (lower HNF of Σ′, sorted list over W ∈ P of (W, w reduced modulo Σ′)). Each row is one orbit, and its class size is the orbit length.
5. **Naming by rebasing.** Choose a basis of Σ′ in block form: d vectors spanning Σ′ ∩ V_I, completed by 3 more. The basis matrix B = [[A, 0], [C, D]] gives a superspace transformation S = B⁻¹. In the new coordinates H is an ordinary record with Σ′ = Z^n, no centrings, and new q′. The k = 1 `standardise` then names it. Every class must be identified with exactly one BYU record of the same d, since BYU's tables are complete.
6. **Checks without BYU** (public tests on hand-written records):
   - P1(a,b,g)0 at k = 2 gives 15 classes, and at k = 3 gives 40.
   - P-1(a,b,g)0 at k = 2 gives 30 classes: 15 lattices × 2 inversion centres (t vs t + Σ/Σ′ representative), with no merging, since 2t ∈ Σ′.
   - The 3D analogue is IT A1, where P-1 has 14 classes of index-2 k-subgroups.
   - The k = 1 rows of the k > 1 code equal `tSubgroups` exactly.

## Cost and runtime (honest estimates; Task 4 and Task 5 measure them first)

| | records | k = 1 classes (est.) | records per SG (max) | measured / estimated |
|---|---|---|---|---|
| d = 1 | 775 | 8 204 | 18 | Python fixtures **96 s** (4 cores). TS sweep (b) **≈ 17 min** under load. |
| d = 2 | 3 338 | ≈ 37 900 | 152 | Python ≈ 0.5–3 h. TS sweeps ≈ 1–3 h. |
| d = 3 | 12 584 | ≈ 145 000 + ≈ 10 000 cubic | 809 | Python ≈ 10–40 h without the invariant prefilter, ≈ 3–10 h with it. TS sweeps ≈ 6–20 h. |

Class counts are |subgroup classes of K| summed per record. The d = 1 counts per space group are reused, and the cubic count is estimated, since cubic parents exist only at d = 3: 244 records with up to 98 classes for m-3m.

The cost per class grows with the candidate list and the E search, which is why the d = 3 range is wide. For k > 1, each parent produces about (number of invariant Σ′) × (complements) / (orbit size) rows. Measure it on a sample before running anything in full.

**Gates.**
- Task 5 times a 1 % stratified sample (every Bravais class and every space group present) and then extrapolates. Above 12 h, ask the user (Decision 4).
- Full sweeps are opt-in: `SSG_FULL=d2` or `SSG_FULL=d3` for Python, `SSG_SWEEP=d2` or `SSG_SWEEP=d3` for vitest. The default local test runs the stratified sample in under 10 minutes.

### Task5 measured cost checkpoint (2026-11-12)

Three workers plus coordinator, BLAS one thread, aggregate RSS watchdog9GB.

| Validation | Parents/classes | Measured wall | Peak aggregate RSS |
| --- | --- | --- | --- |
| d1 full |775 /8204|49.79s|759MB|
| d2 stratified |265 /2293|41.14s|745MB|
| d2 full resumed |3338 /37891|522.23s additional;563.37s cumulative|1170MB|
| d3 stratified only |430 /3499|198.19s|702MB|

All classes uniquely identified and emitted maps independently checked. The d2 sample projected313s at3workers; actual full generation took longer. The d3 sample projects2958s(~49min) ideal full walltime; this is descriptive, not a guarantee. **Full d3 remains deferred by explicit user decision.** Older source/overnight timings were invalidated after exact affine-sieve and monotonic-clock corrections.

## Review Focus

1. **Parallel q's** `(0,0,g1)(0,0,g2)`: any unimodular E mixes them, so a search limited to permutations or signs wrongly reports "not equivalent". Test in Task 4: the 177.2.81.x family (local), and a hand record mapped by E = [[1,1],[0,1]].
2. **Signed letters** at d ≥ 2: `P622(0,0,g1)-h00(0,0,g2)000` must not come out as `h00`. Covered by Task 7's all-record test, plus an explicit hand test of the sign under E = −I on one q.
3. **Unreduced integer q components** (`1,0,g`, `0,1,g` in I, F and C lattices): the type's q text must keep BYU's form. Test in Task 8 on one I-centred and one F-centred d = 2 record.
4. **Cubic d = 3 parents** (m-3m: 98 subgroup classes), the slowest page case. Task 9 adds a local timing test, and the page runs in a Web Worker with progress and Cancel.
5. **k > 1 supercells that remove q's rational part** (q = (1/2,0,g), a doubled): the row's type must move to the (0,0,g) Bravais class, not fail. Test in Task 12 on a hand record.

---

## Session checkpoint — 2026-11-11

Phase A (Tasks 1–3) is complete on `codex/ssg-subgroups-phase-a`; stop for user feedback.
Resume with Phase B / Task 4 only in a later approved session. The original checkout
and its untracked plan were preserved; this worktree contains the implementation.

User constraints override the original runtime gates: at most four cores, aggregate
computation memory below 10 GB, and a checkpoint before a full d = 3 calculation
estimated above four hours. Prefer the published Stokes, Campbell & van Smaalen
(2011) Bravais-class tables for d = 2, 3 q-forms, with citation.

Phase A retains the minimum-letter rule: two nearest-copy disagreements in the d = 1
symbol sweep (basic SG 126 and 133, q-form `(1/2,1/2,g)`). Higher-d checks remain deferred.
Middle-click handover explicitly initializes the new tab's session storage because
default new-tab storage cloning failed in browser verification. Finder receive logic
is unchanged. Long validation tests yield between parents to keep Vitest RPC responsive.

## Phase A — the deferred minors (fix now; they touch files Phase B rewrites)

All of them are fixed now except two: the "only this k-index" checkbox is fixed in Task 13, where it gets a meaning, and the "nearest copies disagree" behaviour is measured first and then decided in Task 2.

### Task 1: Python oracle hygiene

**Files:**
- Modify: `reference/make_subgroup_fixtures.py`, `reference/ssg_equiv.py`, `reference/subgroups.py`, `reference/test_subgroups.py`, `CONTRACTS.md` §5

**Interfaces:**
- Produces: `ssg_equiv.load_records(d: int, data_dir=DATA) -> dict[str, dict]` replaces `load_d1()` and `test_subgroups.all_d1()`. `load_d1()` stays as a one-line alias until Task 5 removes it.
- Produces: `make_subgroup_fixtures.cap_jobs(n: int) -> int` = `max(1, min(n, 4))`.

- [x] **Step 1: Write failing tests** in `reference/test_subgroups.py` (these do not need `local/`):
  - `test_cap_jobs`: `cap_jobs(16) == 4`, `cap_jobs(0) == 1`, `cap_jobs(3) == 3`.
  - `test_worker_env`: after `import make_subgroup_fixtures`, `os.environ["OPENBLAS_NUM_THREADS"] == "1"` (likewise for OMP, MKL and VECLIB).
  - `test_d_guard_raises`: `t_subgroups({"d": 4, ...})` raises `ValueError`, not `AssertionError`.
  - `test_summary_seconds_after_resume`: run `main(["1.1.1.1"], out=tmp)` twice, the second time with `--resume`. `_summary.json["seconds"]` equals `prior_seconds + this_run_seconds`, both are ≥ 0, and a new `"wall_seconds_this_run"` field is present.
- [x] **Step 2: Run** `/opt/anaconda3/bin/python3 -I -m pytest reference/test_subgroups.py -k "cap_jobs or worker_env or guard or resume" -v`. Expect FAIL.
- [x] **Step 3: Implement the fixes:**
  - Set the thread variables with `os.environ.setdefault` at the top of `make_subgroup_fixtures.py`, before `import numpy`, and in the pool `initializer`.
  - Use `cap_jobs` for both `--jobs` and `SSG_JOBS`.
  - Replace the d-guards in `ssg_equiv.py:172` and `subgroups.py` with `ValueError`. Keep the asserts that check internal invariants (`ssg_equiv.py:155,165,182`), but add a message to each.
  - Make `--resume` read the old summary's `seconds` and add to it.
  - Write a comment block on each float prefilter in `ssg_equiv.py` saying what it filters and where the exact check is done.
  - Delete the duplicate loader in the tests.
- [x] **Step 4: Run** the full `reference/test_subgroups.py` (local subset). Expect PASS.
- [x] **Step 5: Commit** with the message "Oracle: cap jobs at 4, single-thread BLAS in workers, ValueError guards, one loader, resume timing".

### Task 2: TS core minors

**Files:**
- Modify: `src/subgroups.ts:120-128`, `src/ssgsymbol.ts` (header and letter choice), `src/standardise.ts` (lines with `Number(` at 194, 311, 428, 479 and 574; header), `src/subfinder.ts:117-124`
- Test: `test/subgroups.test.ts`, `test/ssgsymbol.test.ts`, `test/standardise.test.ts`, `test/subfinder.test.ts`

**Interfaces:**
- Produces: `standardise.ts`: `export function isFail(x: Standardised | Fail): x is Fail` and `export interface Fail { ok: false; reason: string }`.
- Produces: `rational.ts` or `standardise.ts`: `export function toSafeInt(x: Q | bigint, what: string): number`. It throws `RangeError(what)` unless x is an integer with |x| ≤ 2^53 − 1.
- Produces: `ssgsymbol.ts`: `export const _stats = { disagree: 0 }`. This is a test-only counter, incremented when the nearest copies disagree.

- [x] **Step 1: Write failing tests:**
  - `subgroups.test.ts` "maximal, all classes": for the hand P-3c1 record and three more hand records (P2/m, Pmmm, P4/mmm; operators from IT A, written in the test), `maximal` equals a brute-force recomputation: no class Q with P ⊊ Q ⊊ K up to conjugacy. Local-only: the same check for all 775 parents.
  - `standardise.test.ts`: `toSafeInt(Q.of(2n**60n), "x")` throws. `isFail(standardise(badOps, ["0,0,g"]))` narrows the type, checked by a `// @ts-expect-error` line on the non-narrowed access.
  - `subfinder.test.ts`: `runSubgroups` called twice on the same paste calls `standardise` for the parent only once. Spy with `vi.spyOn` on a module-level `parentCache` hit counter that is exported as `_cacheHits()`.
  - `ssgsymbol.test.ts` (local-only): after the 775 sweep, log `_stats.disagree`.
- [x] **Step 2: Run** `npm test -- subgroups standardise subfinder ssgsymbol`. Expect FAIL.
- [x] **Step 3: Implement:**
  - **Sort:** sort an array of `{cls, rep}` pairs, then map to `cls`. No `as any`.
  - **Integer conversions:** replace every unchecked conversion of an exact integer with `toSafeInt`.
  - **`isFail`** replaces every `"ok" in st && st.ok === false` pattern.
  - **Parent cache:** a `Map` keyed by the sorted canonical op strings plus the q text, capped at 8 entries.
  - **Documentation:** add three header paragraphs to `standardise.ts`:
    - why the normaliser entries are {−1, 0, 1};
    - why p4 = 0 is enough for k = 1 (the internal origin is absorbed by m's residues);
    - that the letter-position order is empirical and fixed by the 775 sweep.
  - **Comment fix:** in `ssgsymbol.ts`, correct the "distance" comment. The squared fractional offset is measured in lattice coordinates, not Cartesian ones. For hexagonal 3/6-fold axes and for diagonal planes it is a *tie-break order*, not a metric distance. Say so, and name the cases where it differs.
- [x] **Step 4: Decide the disagreeing-copies rule from data.** If `_stats.disagree === 0` over the 775 sweep (and, after Task 7, over the d = 2, 3 sweeps), return `null` on disagreement and add a test with a hand-made op list that triggers it. If it is > 0, keep the minimum, document the parents it affects as basic SG numbers plus q forms (no labels), and add a test that pins the current behaviour.
- [x] **Step 5: Run** `npm test && npm run typecheck`. Expect PASS, with the finder tests unchanged.
- [x] **Step 6: Commit.**

### Task 3: Page minors

**Files:**
- Modify: `src/ui/subgroups.ts`, `src/ui/subflow.ts`
- Test: `test/subflow.test.ts`

**Interfaces:**
- Produces: `subflow.ts`: `export function displayMap(t: Pick<Standardised, "P"|"p"|"eta"|"m"|"p4">): typeof t`, with p and p4 reduced to [0, 1). This is valid because Σ_B = Z^n contains every integer shift. `affineText` is called on `displayMap(t)`. Task 8 renames the fields to E, M and p4[], and this function follows.
- Produces: `subflow.ts`: `export function prepareSend(parent, row, sg): { ok: true; json: string } | { ok: false; reason: string }`. This separates building the record from writing to storage.

- [x] **Step 1: Write failing tests** in `subflow.test.ts`:
  - `displayMap` on p = `["0","0","1"]` gives `["0","0","0"]`, so the P-3c1 C2/c row's text no longer contains `z+1`.
  - `prepareSend` with a row whose generators do not parse gives `{ok:false}` with a reason that does not mention storage.
- [x] **Step 2: Run** `npm test -- subflow`. Expect FAIL.
- [x] **Step 3: Implement** in `ui/subgroups.ts`:
  - **Parse once:** `run()` reuses the record parsed in `onPaste` (store it as `pasted: GroupRecord | null`). `runSubgroups` takes the record as well as the raw text: `runSubgroups(raw: string | GroupRecord, f)`.
  - **Expanded rows:** keep a `Set<string>` of expanded row keys (key = sorted generators) and re-open those rows after `render()`.
  - **Run button:** `done()` sets `elRun.disabled = !(pasted && supported(pasted))`, and runs only when `id === runId`. A superseded run returns without touching the button.
  - **Send to finder:** handle both `click` and `auxclick` (`button === 1`) on `a[data-send]`. Call `prepareSend` first. On failure show "Could not build the subgroup record: …". Only a storage exception shows the storage message.
- [x] **Step 4: Verify in the browser preview:**
  - Paste the cached 165.1.24.1 page and expand two rows.
  - Change the crystal-system filter and then reset it. Both rows are still expanded.
  - Middle-click "Send to finder". A new tab opens with the subgroup loaded.
  - The C2/c row shows no `+1`.
  - There are no console errors.
- [x] **Step 5: Run** `npm test && npm run build`. Commit.

---

## Phase B — d = 2, 3, k-index 1

### Task 4: Oracle equivalence for any d

**Files:**
- Modify: `reference/ssg_equiv.py`
- Create: `reference/test_equiv_d23.py`, `reference/count_d23.py`

**Interfaces:**
- Produces: `equivalent(rec_a, Q_a: list[list[Fraction]], rec_b) -> dict | None`. The returned map is `{"P": 3×3, "p": [3], "E": d×d int, "M": d×3, "p4": [d]}` (all rationals as Fractions). It works for d = 1, 2, 3. For d = 1, `E = [[eta]]`, `M = [m]` and `p4 = [p4]`.
- Produces: `invariant(rec, Q) -> tuple` (theory step 4) and `rank_Mstar(rec, Q) -> int`.
- Produces: `q_generic(rec) -> list[list[Fraction]]`. Free parameters take the §3.1 values in order and cycle with `+1009·a` when there are more than 3 (up to 9 parameters at d = 3).

- [ ] **Step 1: Count and record.** `count_d23.py` prints:
  - the records per d and per (d, Bravais class);
  - the records per invariant bucket;
  - the number of same-bucket pairs per d;
  - the median and maximum bucket size.

  Write the output to `local/fixtures/subgroups/_counts.json`. These numbers set the size of the distinct-pairs test below.
- [ ] **Step 2: Write failing tests** (local-only unless marked):
  - `test_rank_Mstar`: `rank_Mstar(r, q_generic(r)) == 3 + r["d"]` for every d = 2, 3 record.
  - `test_self_equivalent_d2`: every d = 2 record is equivalent to itself, and `check_map` passes.
  - `test_self_equivalent_d3`: the stratified sample (every Bravais class and every space group present; `SSG_FULL=d3` for all) is equivalent to itself.
  - `test_distinct_not_equivalent`: for every same-bucket pair (all at d = 2; the d = 3 sample's buckets), the result is `None`.
  - `test_parallel_q_mixing` (public, hand records): `P1(a1,b1,g1)0(a2,b2,g2)0`, and the same group with q-basis (q1+q2, q2). The result is not None, and the map has `E == [[1,1],[0,1]]` or another unimodular E.
  - `test_star_permutation` (local): 221.3.213.10 transformed by a 3-fold axis along [111] is equivalent to itself, with E a permutation matrix.
  - `test_d1_unchanged`: the existing d = 1 tests still pass with the new map format.
- [ ] **Step 3: Run.** Expect FAIL.
- [ ] **Step 4: Implement** theory steps 3–4:
  - prefilter by `invariant`;
  - for each P, solve the intertwining and family conditions for E with `intmat`, then enumerate small unimodular E;
  - solve M and (p, p4) with the existing SNF congruence solver, generalised to d rows.
- [ ] **Step 5: Run.** Expect all PASS. Record the per-test runtime in the commit message.
- [ ] **Step 6:** Update CONTRACTS.md §5 with the new map fields. Commit.

### Task 5: Oracle t-subgroups for any d, fixtures, cost gate

**Files:**
- Modify: `reference/subgroups.py`, `reference/make_subgroup_fixtures.py`, `reference/test_subgroups.py`, `CONTRACTS.md` §5

**Interfaces:**
- Produces: `t_subgroups(rec, db=None)` for any d. `db` defaults to `load_records(rec["d"])`. Class dicts keep the d = 1 keys, and `map` uses the Task 4 format.
- Produces: `make_subgroup_fixtures.py --d {1,2,3} [--sample] [--jobs N] [--resume]`. Output goes to `local/fixtures/subgroups/d<d>/`, and `_summary.json` has per-label seconds.

- [ ] **Step 1: Write failing tests:**
  - `test_every_subgroup_identified_d2_sample` (local): every class matches exactly one d = 2 record.
  - `test_star_split` (local): the subgroups of the cubic parent 221.3.213.10 with point group 4/mmm or lower are identified, and their type family has more than one free parameter.
  - `test_p1_d2` (public): the hand `P1(...)0(...)0` gives 1 class.
- [ ] **Step 2: Implement.** Move the existing d = 1 fixtures to `d1/` and regenerate them (≈ 96 s) in the new map format.
- [ ] **Step 3: Cost gate.** Run `--d 3 --sample` and `--d 2 --sample` in the background. Extrapolate the full runtimes from the per-label seconds and write them to `_summary.json`, then into this plan's cost table. **If the d = 3 estimate exceeds 12 h, stop and ask the user (Decision 4).**
- [ ] **Step 4:** Run the full d = 2 sweep (background, `--resume`). If the user approved it, also run the full d = 3 sweep (overnight, background, `--resume`). Required result: `unidentified == []` and every class matched exactly once.
- [ ] **Step 5: Commit** the code and CONTRACTS only.

### Task 6: TS subgroup enumeration for any d

**Files:**
- Modify: `src/subgroups.ts`
- Test: `test/subgroups.test.ts`

**Interfaces:**
- Produces: `PointOp { R: number[][]; E: number[][] }`. This replaces `eps`; the callers in `ssgsymbol.ts` and `standardise.ts` use `E[0][0]` at d = 1.
- Produces: `pointGroup(ops, d)` and `tSubgroups(rec)` for d ≤ 3. They throw `Error("d > 3 is not supported")`.

- [ ] **Step 1: Write failing tests:**
  - hand `P6/mmm(0,0,g1)0000(0,0,g2)0000` (operators from IT A No. 191 with E = ±I; written in the test): its class count equals the d = 1 count for 191 (local value read from `d1/` when present, or the literal number recorded in the test from the IT A subgroup lattice of 6/mmm);
  - local: for every d = 2 record and the d = 3 sample, the class count and the multiset of `index` values equal the Task 5 fixtures;
  - local: for cubic parents, `classSize` sums to the number of subgroups.
- [ ] **Step 2: Implement.** Only `splitOp` and the E bookkeeping change; the closure and conjugacy code is already independent of d.
- [ ] **Step 3: Run** `npm test`. Commit.

### Task 7: Symbols with several q's

**Files:**
- Modify: `src/ssgsymbol.ts`
- Test: `test/ssgsymbol.test.ts`

**Interfaces:**
- Produces: `ssgSymbolInAxes(ops, qFamily: string[], system)` for any d, with `SsgSymbol.shifts: string[]` (one string per q). `full` concatenates `(q_i)shifts_i` in order. d = 1 callers use `shifts[0]`.

- [ ] **Step 1: Write failing tests** (local-only): `ssgSymbolInAxes(buildOps(r), r.q, systemOf(r)).full === r.sym` for every d = 2 record, and for the d = 3 sample (`SSG_SWEEP=d3` for all). On failure, print the first 20 mismatches grouped by (basic SG, q-form pattern).
- [ ] **Step 2: Write the public tests:**
  - the hand P6/mmm d = 2 record gives `P6/mmm(0,0,g1)0000(0,0,g2)0000`;
  - a hand P622 record with τ = (−1/6, 0) for the 6-fold gives `-h00` for q1. **This one pins the sign rule found in Step 3; write its expected value after Step 3.**
- [ ] **Step 3: Implement.** Apply the per-q letter rule of theory step 5 with p4 ∈ Q^d. Then fit the sign rule from the 28 signed records: test "BYU keeps the sign that survives every internal origin and the orientation choice fixed by E" first, and fall back to an order such as −h < h if the data say so. If one family stays irreducible, give it a commented rule keyed by basic SG and q form (no labels). Re-run the Task 2 Step 4 disagreement count on d = 2, 3.
- [ ] **Step 4: Run** until 3338/3338 pass and the d = 3 sample passes 100 %, with d = 1 still 775/775. Commit.

### Task 8: Standard setting for any d

**Files:**
- Modify: `src/standardise.ts`, `src/ui/subflow.ts` (`affineText`, `displayMap`)
- Test: `test/standardise.test.ts`, `test/subflow.test.ts`

**Interfaces:**
- Produces: `Standardised { ok: true; sgNumber; symbol; qActual: string[]; P; p; E: number[][]; M: string[][]; p4: string[]; alternates: string[] }`. This replaces `eta`, `m` and the scalar `p4`. The d = 1 fixtures map `eta → E[0][0]`, `m → M[0]`.
- Produces: `affineText` writes `(x', y', z'; t', u')`, using the variables t, u, v (CONTRACTS §1).

- [ ] **Step 1: Write failing tests** (local-only, as for d = 1):
  - (a) every d = 2 record and the d = 3 sample, in 3 seeded random settings, give back their own `sym`. A random setting change now also applies a random unimodular E with entries in −1..1 and a random p4 ∈ Q^d.
  - (b) every Task 5 fixture class gives back `byu_sym`, or a member of `alternates` (see Task 14).

  Public tests:
  - the hand P6/mmm d = 2 record in a rotated setting comes back as its own symbol;
  - an I-centred d = 2 record keeps `1,0,g` (Review Focus 3).
- [ ] **Step 2: Implement:**
  - **Settings:** extend the BYU setting list with the d ≥ 2 settings in theory step 6.
  - **E basis:** choose E's basis among the unimodular E from the Task 4 method. Prefer E = I. Otherwise use an order that the sweep pins: fewest non-zero off-diagonal entries, then positive diagonal, then lexicographic.
  - **q-forms:** check q′ against the d = 2, 3 Bravais-class q-forms from the source chosen in Decision 2. Until the user decides, the tests run with forms derived from rule (a) of that decision.
  - **Ties:** measure how often d ≥ 2 ties occur and report the count.
- [ ] **Step 3: Run.** Fix the rules until (a) and (b) pass 100 % on d = 1, d = 2 and the d = 3 sample. Then run the full d = 3 sweep in the background if Decision 4 allows it.
- [ ] **Step 4: Commit.**

### Task 9: Orchestration and page for d = 2, 3 (Web Worker)

**Files:**
- Modify: `src/subfinder.ts`, `src/ui/subgroups.ts`, `blog/superspace-subgroups.html`, `build.mjs`
- Create: `src/ui/subworker.ts`, bundled to `assets/blog/ssg-finder/ssgsubw.js`
- Test: `test/subfinder.test.ts`

**Interfaces:**
- Produces: `runSubgroups(raw | rec, f, onProgress?: (done: number, total: number) => void)`. The worker protocol is `{type:"run", raw, filter} → {type:"progress", done, total} | {type:"rows", ...} | {type:"error", message}`, and cancelling calls `worker.terminate()`.

- [ ] **Step 1: Write failing tests:**
  - the hand P6/mmm d = 2 record gives rows; the row with basic SG 177 has the type `P622(0,0,g1)000(0,0,g2)000`;
  - local timing for the largest cubic d = 3 parent: record the total, and fail above 120 s in Node (sets the expectation for the browser);
  - the d = 1 rows do not change.
- [ ] **Step 2: Implement:**
  - remove the d = 1 refusals;
  - name the types incrementally, so rows appear as they are named;
  - add the worker, with a fallback to the main thread if worker construction throws;
  - write the paste-status text for d = 2, 3.
- [ ] **Step 3: Verify in the browser preview** with one d = 2 and one cubic d = 3 page from `local/cache` (local only, never committed):
  - progress is shown;
  - Cancel works;
  - there are no console errors;
  - there is no horizontal scroll at phone width.
- [ ] **Step 4: Commit.**

---

## Phase C — k-index > 1

Work in this order:
1. d = 1 up to k ≤ 4.
2. d = 2, 3 up to k ≤ 3, or as in Decision 3.

### Task 10: Oracle k-subgroups (Python)

**Files:**
- Create: `reference/ksubgroups.py`, `reference/make_ksubgroup_fixtures.py`, `reference/test_ksubgroups.py`
- Modify: `CONTRACTS.md` (new §6)

**Interfaces:**
- Produces:
  - `sublattices(n: int, k: int) -> list[list[list[int]]]`, returning lower HNFs;
  - `k_subgroups(rec, k, db=None) -> list[dict]`, with keys `pg_ops`, `sigma` (HNF in the Σ basis), `kind` ∈ {"3d","internal","mixed"}, `tw` (translation part of each point op mod Σ′, as rationals), `index`, `t_index`, `k_index`, `class_size`, `byu`, `byu_sym` and `map`;
  - `brute_subgroups(rec, k) -> set[key]`, an independent enumeration of all subgroups of G/kΣ by joining cyclic subgroups, then filtered by H ∩ Σ of index k;
  - `rebase(rec, ops, sigma) -> (rec2, S)`, as in theory step 5.
- Fixture §6: `local/fixtures/ksubgroups/d<d>/k<k>/<label>.json`, with the §5 keys plus `sigma`, `kind`, `tw`, `t_index` and `k_index`.

- [ ] **Step 1: Write failing tests:**
  - public: `len(sublattices(4,2)) == 15`, `len(sublattices(4,3)) == 40`, `len(sublattices(6,2)) == 63`;
  - public: hand P1(a,b,g)0 at k = 2 gives 15 classes, and at k = 3 gives 40;
  - public: hand P-1(a,b,g)0 at k = 2 gives 30;
  - public: the key sets of `brute_subgroups` and `k_subgroups` are equal for the hand P2/m(0,0,g)00 at k = 2;
  - public: no two classes of k_subgroups are conjugate, checked over all of G/kΣ;
  - local: for every d = 1 parent with k = 2, every class is identified with exactly one d = 1 record;
  - local: brute force equals k_subgroups on 20 random d = 1 parents with |G/2Σ| ≤ 2048.
- [ ] **Step 2: Implement** theory steps 1–5.
- [ ] **Step 3: Measure** the d = 1 runtime per k on a sample, then write the fixtures in the background: d = 1 for k = 2, 3, 4, and the d = 2, 3 samples for k = 2, 3. Update the cost table.
- [ ] **Step 4: Commit** the code and CONTRACTS.

### Task 11: TS sublattices and k-subgroups

**Files:**
- Create: `src/sublattices.ts`, `src/ksubgroups.ts`
- Test: `test/sublattices.test.ts`, `test/ksubgroups.test.ts`

**Interfaces:**
- Produces:
  ```ts
  export function indexKSublattices(n: number, k: number): IMat[];                 // lower HNF, Σ basis
  export function isInvariant(L: IMat, Ws: IMat[]): boolean;
  export type LatticeKind = "3d" | "internal" | "mixed";
  export interface KSubgroupClass extends SubgroupClass { k: number; sigma: IMat; kind: LatticeKind }
  export function kSubgroups(rec: GroupRecord, k: number): KSubgroupClass[]; // sorted (index, k, |P| desc, sigma, ops)
  ```
- [ ] **Step 1: Write failing tests:** the same public numbers as in Task 10, Step 1; `kSubgroups(rec, 1)` equals `tSubgroups(rec)` up to the extra fields; local: the counts and the multiset of (index, kind) equal the Task 10 fixtures.
- [ ] **Step 2: Implement.** Precompute the K-stabiliser of each sublattice, so every P only scans its own invariant list.
- [ ] **Step 3: Run** `npm test`. Commit.

### Task 12: Naming k > 1 subgroups by rebasing

**Files:**
- Modify: `src/ksubgroups.ts` (adds `rebase`), `src/subfinder.ts`
- Test: `test/ksubgroups.test.ts`

**Interfaces:**
- Produces: `rebase(ops: SOp[], q: string[], sigma: IMat, d: number): { ops: SOp[]; q: string[]; S: QMat }`. The returned q′ is the q text in the new basis, with the rational part written exactly. `standardise(rebase(...).ops, rebase(...).q)` gives the type, and the row's transformation is the composite of `S` and the standardise map.

- [ ] **Step 1: Write failing tests:**
  - public: the hand `P2/m(1/2,0,g)00` with Σ′ = a doubled gives a row whose type q-form is `(0,0,g)` (Review Focus 5);
  - public: rebasing preserves the group order modulo Σ′;
  - local: every class of the Task 10 fixtures gives `byu_sym`, or a member of `alternates`.
- [ ] **Step 2: Implement.** Make `SubRow.kIndex` a `number`. Add `lattice: { kind: LatticeKind; basis: string }`, where basis is written like `"(2a, b, c; t)"`.
- [ ] **Step 3: Run.** Commit.

### Task 13: Page for k > 1

**Files:**
- Modify: `src/ui/subgroups.ts`, `src/subfinder.ts` (`SubFilter`, `keep`), `blog/superspace-subgroups.html`
- Test: `test/subfinder.test.ts`

**Interfaces:**
- Produces: `SubFilter { kIndex: number; onlyK: boolean; kinds?: LatticeKind[]; ... }`. With `onlyK` false, the page lists k′ ≤ kIndex; with `onlyK` true, it lists only k′ = kIndex. Check this against CELLSUB's help text before implementing, and follow CELLSUB if it says otherwise.

- [ ] **Step 1: Write failing tests:**
  - hand P-1 with kIndex 2 and onlyK false gives 2 + 30 rows (k = 1 adds P-1 and P1);
  - onlyK true gives 30 rows;
  - `kinds: ["3d"]` drops the internal and mixed rows.
- [ ] **Step 2: Implement:**
  - a k selector from 1 to the maximum in Decision 3;
  - a Lattice column (kind and basis);
  - a kinds filter;
  - "Send to finder" passes Σ′ through the rebased record.
- [ ] **Step 3: Verify in the browser preview** as in Task 9. Commit.

---

## Phase D — canonical naming decision (2026-11-14)

The user supersedes the former BYU-agreement acceptance gate: all groups are already oracle-verified. BYU spelling is a measurement, not acceptance. No commits until review; nothing under local/ staged or published.

- [x] TDD: add public component-parametrisation and semantic q-block-order regressions before implementation.
- [x] TDD: add hand-family reduced-basis regressions. Among valid plain/star candidates rank total parameter occurrences, nonzero rational offsets and nontrivial phase letters ahead of conventional/setting cost. Generate partial-support lattice seeds and preserve admissible pruning.
- [x] Record witness exit code, signal, timeout and stderr; retry a failed witness subprocess once before stopping its row. Report non-gating more/fewer parameter-occurrence counts by crystal system and list remaining more rows.
- [x] Prefer plain component parameters with rational offsets via unimodular q-basis selection and invertible reparametrisation. Complexity outranks setting cost. Remove raw-prefix pruning that could discard the new minimum.
- [x] Name a, b, g by component; unindexed for one shared set, indexed only for distinct sets by first appearance; preserve shared relations. Normalise rational components to [0,1) subject to exact lattice validity and the existing sign rule.
- [x] Keep basic-axis orientation unchanged. Within that preference, order blocks by free < nonzero rational < 0, larger rationals first, with phases attached and deterministic final ties.
- [x] Add the own-setting / BYU axis- or q-permutation note next to d ≥ 2 symbols.
- [x] Measure **d2 k=1 only**, all **3,338 parents**, **≤4 workers**, ~1 hour. Report five buckets: exact / pure q-permutation / basic-axis setting / phase-after-permutation / other.
- [x] Each parameter first appears plainly in its first block. Allow a later compound block only after exact verification as a basic-point-operation image of an earlier block modulo the basic reciprocal lattice. Count allowed star compounds separately by crystal system; expect only trigonal, hexagonal, cubic. Stop and report any other compound or any other case without a setting-permutation witness.
- [x] Count **basic-cell shear** separately by crystal system only with exact integral det ±1 basic-lattice automorphism, unchanged lattice symbol, point group mapped onto itself, monoclinic unique axis fixed, and an exact group map (origin shifts allowed). Stop otherwise.
- [x] Complete the 430-parent d3 k=1 sample and conditional full-run gate: no stops and a 2.63-hour stratified four-worker projection.
- [x] Verify the d2 phase cohort and retained witnesses under the general unimodular rule; preserve the 120/9 writer split and 19 registered-type rows, apart from the new mixing bucket.
- [x] Finish the full d3 k=1 parent run and independently audit every saved map, label and star compound. Final buckets, timings (including restarts) and stops remain in local/naming-reduced-basis-2026-11-15/full-d3/final-report.md. The later subclass naming gate is separate.
- [ ] After d3 completes, run the authorised 2,000-class stratified d2 k=2 pilot. **2026-11-16 attempt stopped at 48/2,000 on a production naming failure; no population projection or full run.** Local report: local/naming-reduced-basis-2026-11-15/k2-pilot/pilot-report.md. Time actual rebasing, production naming and exact composed-map verification; report buckets, stops and the weighted projection for all 768,066 classes. Ask before the full k=2 run. No full-catalogue type lookup substituted for class-level naming evidence.

**2026-11-14 reciprocal-offset correction:**
- [x] TDD: C, A, I and F reciprocal vectors vanish; nonreciprocal integer vectors do not. Test idempotence, invariance under reciprocal translates, absolute-offset ties, and an exact map for a hand centred star.
- [x] Canonicalise whole offset vectors modulo the centred basic reciprocal lattice by exact L1 size, componentwise absolute size, [0,1) membership, nonnegative signs and numeric order. Transform M and the operations and recompute phases.
- [x] Permit a lattice-symbol change only with an exact signed basic-axis-permutation group map (det ±1); count it in basic-axis setting. Non-permutation shears keep the stricter unchanged-symbol rule.
- [x] Fresh d2 k=1 run from parent 1, ≤4 workers; all earlier spelling rows invalid for this run. Compute each old spelling afresh with a local frozen baseline and report the symbols changed by this fix.

**2026-11-14 cosmetic naming follow-up:**
- [x] TDD: one shared set uses unindexed a,b,g; independent sets keep indices (including a mixed star plus an independent row).
- [x] Prefer the nontrivial phase first within identical q-order keys before Unicode fallback; verify a hand mirror's exact map.
- [x] Fresh d2 k=1 rerun, same buckets/stops. Compare every row with the last complete report; count movements into and out of literal exact agreement, and identify residual disagreements beyond these cosmetic causes. The superseded seven-row expectation is withdrawn; list the 19 s00/0s0 rows and compare Task 7 mechanisms, and audit all 129 phase-after-permutation rows. No current exact row may regress unnoticed.
- [x] Then run exactly the 430 previously validated d3 k=1 parents, ≤4 workers, same buckets/stops. Report crystal-system-stratified wall-time projection for all 12,584 d3 parents and apply the conditional full-run gate above; no k=2 class-naming run.

**Acceptance — 2026-11-15:** d2 k=1 naming is accepted over all 3,338 parents. The pure-origin-only stop is withdrawn. Every one of the 129 phase-after-permutation rows has an exact allowed signed basic-axis/q-permutation and origin map (det ±1); 120 reproduce the reference phase string and nine retain registered generator-assignment differences. Fold the 19 residual phase-order rows into registered-type phase assignment, with no new writer rule. Cosmetics moved 680 rows into literal exact agreement and zero out. Finish the seven missing cubic parents in the validated 430-parent d3 sample, then apply the no-stop / at-most-four-hour stratified-projection gate to the detached, resumable full d3 k=1 run. The full run has the same buckets and stop rules. Propose, but do not start, a 2,000-class stratified k=2 timing pilot. All derived row data remain local; nothing staged or committed.

**d3 measurement — 2026-11-15:** The validated sample is complete at 430/430, with no stops. Its crystal-system-stratified four-worker projection for all 12,584 parents is 2.63 hours. The conditional gate passed and the full d3 k=1 run was launched detached and resumable with four single-thread workers. The full run stopped after 4328 recorded parents under its stop rules. Current counts, crystal-system timing strata and exact witnesses are retained in the local naming-accepted-2026-11-15 report. A 2,000-class stratified k=2 timing pilot is proposed there; it has not started. Nothing staged or committed.

**Unimodular q-basis continuation — 2026-11-15:** General E ∈ GL(d,Z), det ±1, is a legitimate setting when exact checks establish onto q-module equality modulo the centred basic reciprocal lattice, agreement of the printed affine q-family after invertible reparametrisation, and an exact group map with basic P det ±1, rational origin shifts and permitted centred-reciprocal M offsets. Preserve plain/star rules and existing basic-axis/shear restrictions. Prefer a signed q-permutation when one exists; otherwise count **q-basis mixing (unimodular)** by crystal system. General GL(d,Z) maps are allowed for equivalence checking only. The canonical writer must choose the reduced q-basis: minimise total parameter occurrences, then nonzero rational offsets (including free-component offsets), then nontrivial phase letters, ahead of conventional preference and existing setting cost. Generate independent vectors with shared star parameters removed and keep admissible pruning; retain plain/star and q-order rules. The former signed-q-only stop is superseded; all other failures still stop. Verify the previous stop, then resume d3. The 2,000-class pilot is now authorised after d3 completes, superseding the earlier proposal-only wording. The full k=2 naming run still requires approval. All derived data remain local; nothing staged or committed.

## Phase E — finishing

### Task 15: Page text, note removal, index

**Files:**
- Modify: `blog/superspace-subgroups.html`, `blog/index.html`, `blog/superspace-group-finder.html` (link text only)

- [ ] **Step 1:** Delete the `<p class="wip-note" …>Superspace Subgroups (in progress)</p>` at `blog/superspace-subgroups.html:120`.
- [ ] **Step 2: Rewrite the intro and method text for the new scope:**
  - (3+d)D, d = 1, 2, 3;
  - k-index up to the maximum, with the three lattice kinds;
  - how the types are named (algorithmic, validated against the BYU tables locally, nothing BYU-derived shipped);
  - the tie note from Task 14;
  - the expected runtime for cubic parents.
- [ ] **Step 3:** Update the `blog/index.html:39` blurb to "The subgroups of a (3+d)D superspace group (d = 1, 2, 3), with k-index up to N, one row per conjugacy class, named in BYU style". Add a dedicated thumbnail only if the user supplies one.
- [ ] **Step 4:** Run `npm test && npm run typecheck && npm run build`, plus the Python tests (local subset). Then verify in the browser preview:
  - one d = 1, one d = 2 and one d = 3 page, each at k = 1 and k = 2;
  - the index link works;
  - no console errors;
  - phone width.
- [ ] **Step 5: Commit.** Do not push without the user.

---

## Decisions for the user

1. **The orthorhombic naming tie — historical options, superseded by Phase D (2026-11-14).** The former options were:
   - **(a)** Keep the IT A-order main name and show "≡ alternate". This is licence-safe and is the current behaviour. *Recommended unless a rule turns up.*
   - **(b)** Ship a preference that matches BYU. This is only clean if it is a general rule; a fitted table of setting pairs is BYU-derived data and may fall under D1.
   - **(c)** Keep the BYU-matching table in `local/` for tests and personal use only.
2. **Where the d = 2, 3 Bravais-class q-forms come from.** There are 83 classes at d = 2 and 215 at d = 3. The choices are:
   - **(a)** Derive them by rule (the invariant q-forms of each holohedry and centring under E, written in a canonical form fitted to the sweep). This is licence-safe, but more work.
   - **(b)** Transcribe them from the published Stokes, Campbell & van Smaalen (2011) Bravais-class tables, with a citation. The authors are the BYU group, so decide whether D1 covers the paper.
3. **Maximum k on the page, and the internal/mixed lattices.** The suggestion is k ≤ 4 for d = 1 and k ≤ 3 for d = 2, 3. Internal and mixed Σ′ (q → q/k-type subgroups, which CELLSUB has no analogue of) would be listed but off by default, or hidden entirely.
4. **The d = 3 validation budget.** Full sweeps may take 3–40 h of 4-core time (Task 5 measures this first). The choices are full overnight runs with `--resume`, or a stratified sample of every Bravais class and every space group (about 1–5 %) as the acceptance criterion for d = 3.
5. **Order of work.** As written, the plan does Phase B (d = 2, 3) before Phase C (k > 1). Doing k > 1 for d = 1 first (Tasks 10–13 at d = 1) would give a useful page sooner, at the cost of touching `standardise.ts` twice.
6. **Nearest copies disagreeing in `ssgsymbol.ts`.** The plan decides from data (Task 2, Step 4): return null if it never happens, otherwise keep the minimum and document it. Say if you would rather have null in both cases, with a warning on the row.
