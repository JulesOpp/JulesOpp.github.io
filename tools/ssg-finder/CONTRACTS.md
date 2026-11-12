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

## 5. Subgroup fixture (`local/fixtures/subgroups/d<d>/<label>.json`)

`reference/subgroups.py:t_subgroups(rec, db=None)` enumerates k-index-1
(translationengleiche) subgroups for d=1,2,3. `db` defaults to
`load_records(rec["d"])`. Classes are full preimages of K-subgroups, modulo
K-conjugacy; each representative has the lexicographically smallest canonical
parent `pg_ops` list in its conjugacy class. Classes sort by `(index, BYU label
as an integer tuple, pg_ops)` and retain `pg_ops`, `index`, `class_size`, `sg`,
`matches`, and, on exactly one match, `byu`, `byu_sym`, `map`.

A fixture has `label`, `d`, `sym`, `order=|K|`, and `classes`. Each class has
`pg_ops`, `index=|K|/|H|`, `class_size`, `byu`, `byu_sym`, and the flattened map
fields `P`, `p`, `E`, `M`, `p4` from §5.1. All rational entries are reduced
strings; `E` entries are integers. Thus d1 now uses `E=[[eta]]`, `M=[m]`, and
`p4=[origin]`; new fixtures contain no scalar `eta` or row `m`. Unidentified
or multiply identified classes have `byu=null`, `matches`, and `sg` instead.
Every emitted unique map passes the independent exact `check_map` checker.

Identification considers every candidate with the same basic SG and exact
point-group characters (including internal traces), dimension, order, and the
exact affine fixed-torus profile described below.
Raw canonical `nops` and free-parameter counts are not bucket keys: centred
settings change the former, and a narrower parent family may match a canonical
family with more parameters. The oracle enforces external positive orientation.

`make_subgroup_fixtures.py --d {1,2,3} [--sample] [--jobs N] [--resume]`
writes only under ignored `local/`. Default dimension is 1. `--sample` uses
`check_equiv_d23.sample`, covering every Bravais class and every basic SG.
Full d3 sweeps are deferred and the CLI rejects an implicit full d3 request;
explicit d3 labels remain available for regressions. Pools clamp to 1–4 workers;
BLAS/OMP/MKL/VECLIB threads are set to one before numerical imports. The operator
must coordinate concurrent jobs to keep the total at four cores and use a 9 GB
aggregate-RSS watchdog for background generation.

Fixtures are atomic per-label checkpoints. `_signature` binds writer, subgroup,
oracle, exact-arithmetic, integer-matrix sources and the entire dimension's
corpus. `_content_sha256` binds the artifact to its validated contents, and
`_seconds` retains that label's measured worker wall time. `--resume` reuses only
complete, checksum-valid fixtures with identical signatures; failed classes and
stale/tampered fixtures are recomputed. No legacy unbound fixture is reused.

`_summary.json` is atomically updated after every completed label. It records
per-label seconds, class counts and unidentified classes, cumulative invocation
wall `seconds`, `wall_seconds_this_run`, effective `jobs`, signature and completion.
The cost estimate includes every label with its measured/projected seconds,
supporting sample labels and projection method: measured, same SG/class mean,
class mean, SG mean, or overall sample mean in that order. Summing projected
worker times and dividing by jobs gives an ideal balanced wall estimate; startup,
stragglers and sample bias are explicitly excluded. It is not a confidence bound.
Run stratified d2/d3 timing samples before full d2 generation. The user deferred
all full d3 sweeps; any future full d3 estimate above four hours requires a prior
user checkpoint. Never infer approval from completion of a sample.

Tests in `reference/test_subgroups.py` include public synthetic P1 d2 and fixture/
checkpoint tests, plus individually marked local corpus tests (skip without
local data): d1 regression subset, stratified d2 exact-one matches and checked
maps, and cubic `221.3.213.10` star splitting. TS d1 tests read `d1/` and preserve
Phase A event-loop yields and application cache/storage behavior. Existing legacy
root d1 fixtures may be archived locally after regeneration; never stage them.

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

For d = 1, `E=[[eta]]`, `M=[m]`, `p4=[p4]`. All new fixtures and oracle
callers use matrix fields. `check_map` continues accepting legacy d1 maps for
compatibility with locally archived fixtures.

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

### 5.2 Exact affine fixed-torus profile for subgroup candidates

For each linear part W choose one affine lift (W,w), and fold the full
superspace centring into a lattice basis L for Σ. In lattice coordinates the
integer matrix is A=L⁻¹(I−W)L and translation is b=L⁻¹w. An affine action fixes
a torus point exactly when A y=b modulo Zⁿ. Let U A V=D be the exact Smith
normal form. The transformed equations are D z=U b modulo Zⁿ. Each nonzero
Smith diagonal permits a real solution; each zero row requires its U b entry
to be an integer. If nonempty, the fixed set has dimension n−rank(A) and
product of the nonzero absolute Smith factors connected components.

The profile is the multiset `(nonzero absolute Smith factors, nonempty)` over
linear parts. Replacing a lift changes w by Σ and cannot change this result.
Any valid affine lattice map induces an integral unimodular conjugacy between
lattice-coordinate torus actions. It preserves Smith factors and bijects fixed
sets, so profile equality is necessary even for coupled centring and narrower
parent q families. This is a pruning invariant, not a sufficient equivalence
test. It uses no q-family comparison. Origin shifts and M terms are accounted
for by affine conjugacy; a basic screw with empty fixed set is not incorrectly
separated by a compensable internal half-phase. The existing exact oracle and
map checker remain authoritative after this candidate filter.

The pure exact caches are bounded: up to 8,192 `(W,L)` Smith decompositions and
20,000 `(generators,centring,d)` profiles per worker. No BYU table is embedded
or written outside ignored local data. Public tests cover P2 d3 phase pruning
before any E search, a known valid screw-phase M map, internal/coupled centring
lifts, and rational external settings.

### 5.3 Scalar-action phase-character quotient

This additional necessary filter applies only when both records have product
translation lattice L3⊕Z^d, zero mixed point-operation blocks, and scalar internal
point actions ε_g I_d with ε_g=±1. Otherwise comparison is skipped, never rejected.
Product eligibility is checked exactly from integer internal centring coordinates.
A conservative order≤24 guard also skips larger point actions. This filter makes
no comparison of q-family widths or parameter values.

In primitive external coordinates R_g=L3⁻¹ R_g L3, an admissible map has integer
rows j=M L3 satisfying j R_g=ε_g j. Compute the complete integral row lattice U
by rational nullspace followed by integer saturation. For all point elements g,
write h_u[g]=u·L3⁻¹w3_g, c[g]=1−ε_g and F_i[g]=wI_g[i]. Let N be any rational
full-row-rank basis of the annihilator of c. Form the full-rank rational lattices
L0=N Z^|K| + Σ_(u in basis U) Z N h_u and L1=L0 + Σ_i Z N F_i. Since L0⊂L1,
B1⁻¹B0 is an integer matrix for column bases B0,B1. Its nonunit Smith factors
are exactly the finite quotient L1/L0; these factors constitute the profile.

Under affine equivalence, point elements are permuted. External lattice changes
identify the integer covariant row lattices. Internal rows transform by an
integer unimodular E. Allowed M adds only combinations of h_u; changing affine
lifts adds integer vectors. Internal origins add multiples of c. Physical-origin
changes alter h_u by u(I−R_g)p3=(1−ε_g)u p3. N kills both origin contributions.
Thus projected L0 and L1 undergo the same invertible rational change of basis,
and their finite quotient is preserved. Permuting the point elements or changing
annihilator basis does not affect its Smith factors. Allowing every covariant M
rather than only family-compatible M weakens pruning and cannot cause a false
rejection for a narrower parent family. The existing exact equivalence/map check
remains authoritative after this filter.

The profile cache is bounded to 20,000 structural records; no data table is
embedded. Candidate indexing is lazy by basic SG: the catalogue is grouped cheaply
by label once, and expensive invariant/fixed-torus buckets are built only for SGs
actually requested by subgroup classes. At most three database identities are
retained. Scalar-phase profiles are evaluated only for surviving candidate pairs.
Public controls include C4v screw/glide internal phases invisible to fixed-torus
nonemptiness, compensable P21 screw phases, primitive-to-centred rational P maps,
negative scalar origin removal, and mixed/nonscalar/coupled skip cases.


### 5.4 Simultaneous affine feasibility sieve

Subgroup matching supplies an optional `e_candidates` callback to the exact
oracle. Other callers retain the original bounded search. For product lattices
ΣA=L3A⊕Z^d and ΣB=L3B⊕Z^d, fixed positive rational P determines the target point
operation WB for each source point operation WA through faithful external
projection. The sieve uses a smallest generating set of the finite point group.
It skips presentations with nonzero internal centring components.

The integer unknowns are E, J=M L3A and target-lattice lift vectors. Real unknowns
are the target q parameters and a single common external/internal origin s.
The system imposes E EA=EB E, M RA+E MA=MB P+EB M,
(B0+Σ yk Bk)P=E A0+M, source-direction inclusion into the target family, and
LB^-1(S wA+(I-WB)s-wB) integral for every selected generator. It first requires
L3B^-1 P L3A to be integral unimodular. Using source A0 is valid even for a
narrower source family: the transformed source parameter terms can be absorbed
into target parameters by the homogeneous direction conditions, without changing
M or s. All arithmetic is exact rational/integer arithmetic.

Eliminating real variables by rational row reduction leaves integer linear
Diophantine equations. Smith form gives their full affine integer solution
lattice; projection onto E followed by HNF retains all congruences. The projected
lattice is **not saturated**. Empty solutions, or an affine residue space whose
matrices are all singular modulo 2 or modulo 3, reject the fixed P. Otherwise,
exact coset membership filters the original bounded E stream before solving M
and origins. The callback preserves the original candidate order and bound; it
does not enumerate the projected affine basis using the old coefficient bound.
Every accepted map still passes the independent exact checker.

These are necessary conditions: any originally accepted map supplies integer
E,J and affine lifts satisfying the system, and unimodular E is nonsingular over
both residue fields. Eliminating and projecting therefore cannot reject a map
in the previous bounded search. Source hashes include `affine_search.py`.
Per-label and incremental/final elapsed times all use the same monotonic clock,
so wall-clock changes do not corrupt timing or full-corpus estimates.
