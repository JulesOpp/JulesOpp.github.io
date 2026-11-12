# Shared formats (contract between build_data/, reference/ and src/)

The full plan is in `~/dials/In8IrSe9_Xtal4_FullRun/ssg_finder_plan/PLAN.md`. This file fixes the formats that more than one part reads or writes. If you need to change one, update this file.

**Licence rule (decision D1):** anything derived from BYU ISO(3+d)D goes under `local/` (gitignored), never anywhere else in this public repo. That covers operators, cached pages, JSON and fixtures that contain BYU groups.

## 1. Rationals and operators

* A rational is a string: `"0"`, `"1/4"`, `"-3/8"`. It is always reduced, with a positive denominator.
* An operator is in xyzt notation over the variables `x,y,z,t,u,v` (the first 3+d are used). Write `"-x+1/2,y,z,-t+u,..."`: no spaces, terms ordered by variable, constant last, translation reduced to [0,1).
* The **canonical operator list** of a group is the closure of its generators and the centring translations, modulo Z^{3+d}. Translations are reduced to [0,1). The list is sorted lexicographically by `(W flattened row-major, w)`, comparing rationals numerically. Operator indices in fixtures refer to this order. The identity always comes first, as an exception to that sort; every other operator, centring translations included, follows in sorted order.

## 2. Group record (`local/data/sgNNN.json`, also the parsed form of paste mode)

```jsonc
{ "sg": 55, "groups": {
  "55.3.122.96": {
    "sym": "Pbam(a1,b1,0)000(-a1,b1,0)00s(0,0,g2)000",
    "cls": "3.122",
    "d": 3,
    "q":   ["a1,b1,0", "-a1,b1,0", "0,0,g2"],   // basic-setting q's, components as written in the label
    "params": ["a1", "b1", "g2"],                // free parameters, first-appearance order
    "cen": ["0,0,0,1/2,1/2,0"],                  // non-zero basic-setting centring vectors (list of rationals, comma-joined)
    "gen": ["<op>", "..."],                           // basic-setting generators (not including centrings or identity)
    "nops": 8                                    // size of the canonical operator list (incl. centring combos), for checking
  } } }
```

Paste mode parses the BYU `ssginfoform.php` page text into the same record, with `gen` = all non-lattice operators and `cls` taken from "Bravais class:" if present.

## 3. Golden fixture (`local/fixtures/golden/<label>__<params>.json`, written by `reference/make_golden.py`)

```jsonc
{
  "label": "55.3.122.96", "d": 3,
  "params": {"a1": "1/4", "b1": "1/4", "g2": "1/4"},
  "Q": [["1/4","1/4","0"], ...],                 // d rows × 3, basic setting
  "ops": ["x,y,z,t,u,v", ...],                   // canonical operator list (§1)
  "supercell": {
    "hnf": [[...],[...],[...]],                  // integer, rows = basis of L_super in basic coords, lower-triangular HNF
    "conv": [4,4,4],                             // N_j = lcm of denominators of column j of Q
    "centring": [["0","1/2","1/2"], ...],        // in conventional-supercell fractions, excluding 0
    "index": 32
  },
  "Lambda": [["1/4","0","0"], ...],              // d×d basis of Λ = Z^d + Q Z^3 (columns = basis vectors), HNF-derived
  "classes": [                                   // sorted by (-dim, sg_number, setting, generic_t)
    {
      "dim": 2,
      "flat": {"C": [[0,0,1]], "f": ["1/2"]},    // canonical connected flat in s-coords (t = Λ·s): C s ≡ f (mod 1), C saturated HNF
      "locus": "t3 ≡ 1/8 (mod 1/4)",             // human text in t (basic setting)
      "generic_t": ["...","...","..."],          // exact point used for identification
      "example_t": ["0","0","1/8"],              // nice point in the same class
      "surviving": [[0,[0,0,0]], [5,[0,1,0]]],   // [op index, matched n reduced to 0<=n_j<N_j, lexicographically smallest]
      "sg_number": 7,
      "sg_symbol": "Pc",                         // standard short HM from spglib
      "setting": "A e 1 1",                      // symbol in SSG axes, conventional supercell
      "orbit": 0                                 // id of the SSG-orbit merge group (§4.4 of the plan)
    } ]
}
```

The TS port in `src/` must reproduce every field except `sg_symbol` and `setting`, which may differ only in spacing or formatting. `test/golden.test.ts` normalises those.

### 3.1 Field rules (fixed by `reference/exact.py`, M2)

These make every field reproducible byte for byte. `exact.py` is the executable definition; the rules below are what the TS port must follow.

* **File name.** `<label>__<k>=<v>,<k>=<v>.json` with params in record order and `/` written as `_`, e.g. `55.3.122.96__a1=1_4,b1=1_4,g2=1_4.json`. `_summary.json` in the same folder is a run log, not a fixture.
* **HNF conventions.** *Row HNF* (`intmat.hnf`): rows in echelon form, zero rows dropped, pivot positive, entries **above** each pivot in `[0, pivot)`. *Lower HNF* (`intmat.hnf_lower`): the same lattice with pivots on the diagonal and entries **below** each pivot in `[0, pivot)` (reverse the columns, take the row HNF, reverse columns and rows back).
* **`supercell.hnf`** = lower HNF of the generators `{n in the box ∏[0,N_j) : Qn ∈ Z^d, n ≠ 0} ∪ {N_j e_j}`. `index` = det. `centring` = those box points n ≠ 0, as `n_j/N_j`, in lexicographic order of n.
* **`Lambda`** = B, columns = basis vectors. With D = lcm of all denominators in Q, Bᵀ = (lower HNF of the rows `D·e_i` (i<d) and `D·(column j of Q)` (j<3)) / D. So B is upper-triangular. s-coordinates are `s = B⁻¹ t`.
* **`flat`.** C is the row HNF of a saturated basis (pivots positive, full row rank r = d − dim); `f` is transformed with the same unimodular U (`f ← U f`) and reduced to [0,1). The full torus is `{"C": [], "f": []}`.
* **`generic_t`.** In s: the non-pivot coordinates of C are free; free coordinate m (0-based) at attempt a = 0,1,… is `k/10007` with `k = (S[m] + 1009·a) mod 10007` (k = 0 → 1), `S = (3001, 4723, 7177)`. Pivot coordinates are solved from `C s = f` exactly, last row first. Take the first attempt whose point lies in no proper sub-flat of the closure. Then `t = B s`, each component reduced to [0,1).
* **`example_t`.** For m = 1..24: run over the free coordinates `k/m`, k ∈ [0,m)^dim in lexicographic order, pivots solved as above, keeping points in F's open part. For the first m with any such point, form all `(B s + λ) mod 1` for λ in the Λ/Z^d representatives `{Qn mod 1}`, and pick the minimum of (lcm of denominators, number of non-zero components, t as a tuple of rationals). If no m works, `generic_t` is used.
* **`surviving`** is evaluated **at `generic_t`**, sorted by op index. The matched n depends on the Λ-translate of t: moving t by Qm changes n to `n + (R − 1)m` (mod L_super, then reduced to the lexicographically smallest representative in the box). So n identifies a connected piece of the class in t-space, not the class. Two different classes can even share all keys: in 55.1.9.4 with g = 2/5, t = 0 and t = 1/2 are different classes (IT columns n/N and 1/2N+n/N, both Pbam) with the same surviving set and n's.
* **Class order.** `(-dim, sg_number, setting, generic_t)`, with `generic_t` compared as a tuple of rationals. `orbit` ids are numbered in order of first appearance in that list. Orbits are the SSG's own action `s ↦ Ẽ s + c̃` (Ẽ = B⁻¹EB, c̃ = B⁻¹(τ_I − Qτ_E)), with no normaliser merge.
* **`sg_symbol`** = spglib `international_short` (e.g. `"P2_1/c"`). The spglib call gets the lattice `L` with `L Lᵀ = diag(N) (Σ_R RᵀR) diag(N)`, an invariant metric of the basic point group in the conventional supercell. **`setting`** = `ssg_descent.hm_in_setting` on the 3D operators in the conventional N₁×N₂×N₃ cell (orthogonal-axis symbol only, so it is unreliable for hexagonal and trigonal groups).
* **`locus`** is the same for every class of an orbit: the union over the orbit's flats (Λ-translates are already included, since s lives on R^d/Λ). Grammar:
  * `general` for dim = d. Variables are `t` for d = 1 and `t1..td` otherwise.
  * Flats are grouped by the row lattice W of `K = C B⁻¹`; its basis is the row HNF of `den·K`, divided by den. For each group the values `φ = W t mod 1` form a finite set R. Let H be the stabiliser of R (`{h : R + h = R}`); the coarser rows are `z·W` for z running over the row HNF of the dual of `Z^r + H`.
  * Each row w is written `c·t ≡ value (mod μ)`, where c is primitive integer with its first non-zero entry positive, `w = ±c/μ`, and the value is reduced to `[0, μ)`. A linear form prints as `t1 + t2`, `t1 - 2t3`, `-t2`. Rows whose value is the same over R are printed alone, joined by `, `. Rows that take several values are printed together as `lhs ≡ a or b (mod μ)` (one row), or `(lhs1, lhs2) ≡ (a1, a2) or (b1, b2) (mod μ1, μ2)` (several rows).
  * Several W-groups are joined by ` or `, sorted, each one bracketed `[…]` when it has more than one condition.
* **d = 1 only:** top-level `"it_N"` = |Λ/Z| (Λ = (1/it_N)Z), and per class `"it_label"`: `general` for dim 1, otherwise from a = frac(t·it_N): `n/N`, `1/4N+n/N`, `1/2N+n/N`, `3/4N+n/N` (other values print as `a/N+n/N`). Whether this N equals the IT finder's N has **not** been checked against its column labels yet (needs plan T1).

## 4. Commands

* Python: `/opt/anaconda3/bin/python3 -I` (spglib 2.7, numpy, sympy). Tests: `... -m pytest`.
* Node 24: `npm test` (vitest) and `npm run build` (esbuild → `../../assets/blog/ssg-finder/`).

## 5. Subgroup fixture (`local/fixtures/subgroups/<label>.json`, written by `reference/make_subgroup_fixtures.py`)

The k-index-1 (translationengleiche) subgroups of a (3+1)D group (d = 1 only), one row per conjugacy class, each identified with its BYU type. These are the full preimages of the subgroups of the basic point group K, so classes correspond to K-conjugacy classes of subgroups of K.

```jsonc
{ "label": "165.1.24.1", "sym": "P-3c1(0,0,g)000", "order": 12,   // order = |K|
  "classes": [                    // sorted by (index, byu label as a tuple of ints, pg_ops)
    { "pg_ops": [0, 1, 3, 6, 8, 9], // canonical-list indices (§1) of the parent: for each point operation of
                                    // the class representative, the first op with that linear part; sorted
      "index": 2,                   // |K| / |subgroup|  (t-index; k-index is 1)
      "class_size": 1,              // number of K-conjugates
      "byu": "150.1.24.1", "byu_sym": "P321(0,0,g)000",
      // the map from the subgroup (parent setting, parent's generic q) onto the BYU record:
      // x' = P x + p,  x4' = m·x + eta x4 + p4
      "P": [["1","0","0"],["0","1","0"],["0","0","1"]], "p": ["0","0","-1/4"],
      "eta": 1, "m": [0,0,0], "p4": "0" } ] }
```

* **Representative.** Among the conjugates in a class, the one whose `pg_ops` list is lexicographically smallest.
* **Identification** (`reference/ssg_equiv.py`, `equivalent`). The subgroup, at the parent's generic Q (the §3.1 values `3001/10007, 4723/10007, 7177/10007` for the parameters in order), is tested against every d = 1 record with the same basic space-group number. Exactly one must match; the type's q is then `q' = (eta q + m) P⁻¹`, a member of the BYU q family (not reduced modulo L*). The map is one valid choice, not a canonical one: `P`, `p`, `m` are not unique. `m` is an integer row; `p` is not reduced into [0,1).
* **Equivalence.** `S = [[P, 0], [m, eta]]` maps Σ_A onto Σ_B, `S g S⁻¹` lies in B modulo Σ_B for every g in A, and q' lies in B's q family exactly. P is searched as `T_B⁻¹ N T_A`, where `T` is spglib's standardisation and N runs over integer matrices with entries in −2..2, det +1, that normalise the standard point group and keep the standard centred lattice (one N per coset N·K). Then m and (p, p4) are solved from mixed real/integer congruences with the Smith normal form.
* `_summary.json` in the same folder is a run log (per-group class counts, timings, any unidentified classes), not a fixture. `seconds` is cumulative wall time across resumed runs; `wall_seconds_this_run` is wall time for this invocation. `jobs` records the effective capped worker count.
* Tests: `reference/test_subgroups.py` runs a representative subset (every Bravais class, plus 14, 15, 62, 165, 194); `SSG_FULL=1` runs the sweep over all 775 (3+1)D records.
* Generation: `make_subgroup_fixtures.py [--resume]` uses at most 4 worker processes (`--jobs` or `SSG_JOBS` selects 1–4, clamped to that range), each with BLAS/OMP limited to one thread; `--resume` keeps complete fixtures and adds this run’s wall time to the prior summary’s `seconds`. All 775 take roughly 4–5 minutes on 4 cores.

* Oracle record loading: `ssg_equiv.load_records(d, data_dir=DATA)` caches records by dimension and data directory. `load_d1()` remains a compatibility alias for `load_records(1)`. Public hygiene tests use synthetic data and run without `local/`.

### 5.1 Exact oracle matrix API (d = 1, 2, 3)

`q_generic(rec)` returns a d×3 matrix of Fractions. Parameter i takes
`(S[i % 3] + 1009*(i // 3))/10007`, S = (3001,4723,7177).
`equivalent(rec_a, Q_a, rec_b)` returns `None` or a map with:

* `P`: 3×3 Fractions; `p`: three Fractions.
* `E`: d×d unimodular integers; `M`: d×3 Fractions; `p4`: d Fractions.

The map is `x' = P x + p`, `t' = M x + E t + p4`, with
`Q' = (E Q_a + M) P^-1` exactly in B's affine parameter family. It preserves
external orientation: `det(P) > 0` and `det(P)*covol(Σ_A) = covol(Σ_B)`;
`det(E)` may be either sign. This follows van Smaalen, Campbell and Stokes
(2013), §2.1 equation (7) and §6
(<https://journals.iucr.org/a/issues/2013/01/00/pc5018/index.html>).
Thus enantiomorphic types remain distinct. Positive normaliser candidates are
selected before point-group coset representatives are taken.

For d = 1, `E=[[eta]]`, `M=[m]`, `p4=[p4]`. The existing d1 fixture writer
currently serialises the legacy scalar `eta`, row `m`, scalar `p4` fields in
§5; it adapts matrix oracle results explicitly. `check_map` accepts those legacy
fixture maps while fixtures and subgroup generation are migrated separately.
New oracle callers use only matrix fields.

`invariant(rec,Q)` returns dimension, point-group order, symbolic module rank,
and the sorted multiset `(order R, det R, tr R, tr E(R))`. Raw canonical `nops`
is not a conjugacy invariant between primitive and centred settings; use
`|K| = nops*covol(Σ)` instead. Corpus pair buckets additionally use the basic
space-group number (known directly from the local catalogue label).

`rank_Mstar(rec,Q)` uses independent *formal family coefficients*, stacking
the constant and every parameter coefficient of each q row alongside the
three ordinary reciprocal generators. A rational numeric Q always has
Q-linear rank at most three; it cannot represent incommensurate rank 3+d.
The helper verifies exact Q shape and affine-family membership; exact operator
faithfulness and map validity are checked by the oracle/checker. It does not
claim that the supplied rational evaluation is itself irrational or Q-independent.

For each external candidate P, intertwining and symbolic-family constraints
are solved exactly; their integer solution lattice is saturated and put in
row HNF. E is searched with lattice coefficients in −2..2 and determinant ±1
(identity is tried first when allowed). The bounded external normaliser and
bounded E search can miss equivalences in arbitrary distant settings; a
returned map is independently checked exactly. Exact determinant-polynomial
content and affine-family torsion obstructions reduce negative searches.
The torsion obstruction applies only to product lattices L3⊕Z^d and equal
family dimensions, so a narrower subgroup family can still match a canonical
family with more free parameters. M and all origin coordinates follow from
Smith-normal-form mixed congruence solves.

A further product-lattice obstruction keeps the finite abelian group of
constant q-row offsets transverse to all free physical directions. In primitive
external coordinates, let V span every parameter-coefficient row, N span its
annihilator, L0=N Z^3, and L1=L0 plus the integer span of N Q0[i]^T. Rational
HNF of exactly these generators (without adding an ambient unit lattice) and
SNF of B1^-1 B0 give the nonunit invariant factors of L1/L0. Compare only when
the physical row ranks agree and both lattices have no internal centring:
family inclusion then identifies V under unimodular external coordinates;
unimodular E preserves the projected-row subgroup and integral M contributes
only L0. Unequal physical row ranks deliberately skip this obstruction for
narrower subgroup families. It distinguishes a single order-two offset from
two independent order-two offsets, which scalar family torsion cannot.

For a fixed external P, a stronger embedded-lattice check runs before E search.
On product lattices require U=L3B^-1 P L3A integral unimodular. In B's cached
annihilator N, all source coefficient rows transformed by P^-1 L3B must project
to zero. Compare the exact rational-HNF lattices generated by N Z^3 and the
projected target constants against those generated by N Z^3 and source constants
transformed by P^-1 L3B. Primitive M contributes only N Z^3 and unimodular E
preserves the generated subgroup, so embedded lattice equality is necessary.
Unlike the abstract quotient comparison, this per-P condition needs no equality
of physical row ranks or parameter-family dimensions. It is skipped for
internal/mixed centrings, and the projection condition is vacuous when N is empty.
It distinguishes offset directions that are abstractly isomorphic but cannot
be related by an allowed external normalizer, while accepting an external shear
when that shear is actually available.

A threefold acting on three chosen arms of a four-arm tetrahedral star can
have integer sums in E: the omitted arm is minus the sum of the chosen arms.
A permutation-matrix-only expectation is invalid in that basis.

`count_d23.py` writes local `_counts.json`, including total and within-basic-SG
pair counts. `check_equiv_d23.py --d 2 --jobs 4` validates all d2 pairs with
source-and-corpus-hashed per-case JSONL resume checkpoints under local/.
`--d 3` uses a stable sample covering every Bravais class and every basic SG;
`--full` is explicit and requires a prior runtime estimate and user approval
above four hours. The user has deferred all full d3 sweeps for Phase B.
Default pytest pair tests use up to 100 evenly spaced pairs from stratified records; `SSG_FULL=d2` requests all
d2 pairs, and `SSG_FULL=d3` requests all d3 records. Full sweeps should use the
resumable harness in the background. All local corpus tests skip without local/.
