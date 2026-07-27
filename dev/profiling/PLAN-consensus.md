# Consensus optimisation & Jansson comparison — plan

Date opened: 2026-06-02. Author: Claude (Opus 4.8), at MRS's request (user unavailable; no plan approval possible — proceeding after advisor consult).

## Goal (user's words)

> Aggressively /profile the new Consensus implementation. Fast throughput at
> high-k/low-n, high-n/low-k, and high-n/high-k. We have an O(kn) impl that
> doesn't mirror Jansson's; implement Jansson's algo in parallel; optimise both;
> use data to decide which wins where. Crossover ⇒ switch near it; one [almost]
> always wins (within a couple %) ⇒ keep just that, drop redundant code.
> Correctness paramount: **no systematic error**. Micro-probability of split
> collision is acceptable.

## What the code does today (`src/consensus.cpp`, `R/Consensus.R`)

- **Strict (p = 1)**: single-reference O(kn) pass over tree 0 (`CLUSTONL/R` +
  contiguity test). Already optimal. *Out of scope.*
- **Majority / threshold (0.5 ≤ p < 1)**: count every split's frequency in one
  pass, keep those with `count >= thresh` where
  `thresh = min(floor(k·p)+1, k)` (i.e. strictly > k·p). Correct *by
  construction*: every kept split occurs in > k/2 trees ⇒ pairwise- (hence
  globally-) compatible ⇒ `splits_to_edge` forms a tree with no merge step.
  Two counting cores share one postorder primitive (`for_each_internal_node`):
  - `count_splits_hashed` — **default**. O(kn). 128-bit splitmix64 sum-hash of
    leaf ids, accumulated incrementally up the tree. Probabilistic (≈1e-30
    collision). Tiny constant: one postorder pass, O(1)/node.
  - `count_splits_exact` — `hash = FALSE`. O(k·n·h). Builds an exact packed
    leaf-set bitmask key by iterating L..R per cluster (cost = Σ cluster sizes =
    O(n·h)). Deterministic. Slow for tall trees (h ~ n).

## Jansson's algorithm (arXiv 1307.7821 = SODA/RECOMB 2013; J.ACM 2016)

- Deterministic O(kn) **standard majority** lives in ref [18] (RECOMB 2013).
  The arXiv paper's `Maj_Rule_Plus` (Fig. 2) is the same Phase-1 candidate
  generation; only the keep-test differs.
- **Phase 1** (candidate generation, ⊇ all majority clusters): maintain a
  candidate *tree* T (start = T_1, all counts 1). For each later tree T_j:
  Boyer–Moore-style update — `+1` if a candidate cluster occurs in T_j, `-1` if
  incompatible with T_j, unchanged if compatible-but-absent; delete count-0
  nodes (top-down); then insert T_j's clusters that are compatible with T but
  absent. Uses **Day's algorithm** (= our `ClusterTable`, O(1) "occurs in?"
  after O(n) prep), **`One_Way_Compatible`** and **`Merge_Trees`** (ref [17],
  SODA/J.ACM — still need to fetch their exact procedures).
- **Phase 2 (standard majority)**: recount K(v) = #trees each surviving
  candidate occurs in; keep iff K(v) > k·p (same threshold as today). (The
  arXiv paper's Maj+ Phase 2 instead keeps K(v) > Q(v); we do **not** want
  Maj+ — it is a *different, more resolved* tree and would be a systematic
  error vs `ape::consensus`.)

Derivation soundness: majority ⊆ majority(+), so the Phase-1 candidate set
(⊇ all Maj+ clusters, Lemma 4) ⊇ all majority clusters for any p ≥ 0.5. Phase 2
keeps exactly the > k·p subset; survivors live in tree T so are mutually
compatible ⇒ valid tree. ✓

## Why this is subtle / the key prior

- **FACT (Jansson's own C++ reference impl) majority code is BROKEN** (verified
  by build 2026-05; see memory `fact-majority-broken.md` + this dir). Crashes /
  drops majority splits. ⇒ Implement from the **paper with proofs**, never port
  FACT. The bug is in FACT's *code*, not Jansson's *algorithm*.
- **The current method is already O(kn) (hashed) and correct-by-construction**,
  needing no merge — a genuine simplification that Jansson does not have.
- **Early evidence hashed dominates**: the paper's own C++ `Maj_Rule_Plus`
  timings (their 2.2 GHz box) are ≈10× slower than TreeTools' historical hashed
  numbers (e.g. (k1000,n100): FACT 1.29 s vs TreeTools ≈0.13 s). Jansson does
  several O(n) passes/tree (Day prep, 2× One_Way_Compatible, Merge_Trees, BST);
  hashed does one O(1)/node pass. So Jansson's realistic edge is **determinism**
  + **beating `exact` on tall trees** (O(kn) vs O(knh)), NOT beating hashed.

## Framing the comparison

Jansson competes with **`exact`** (the deterministic path), not with `hashed`.
- vs `hashed`: user has waived collision risk; Jansson must also be *faster* to
  win, which the analysis says it won't be. Expect hashed = default everywhere.
- vs `exact`: Jansson O(kn) should beat exact O(knh) for **tall, high-n** trees;
  exact's smaller constant may win for **short/balanced or small** trees. A real
  crossover may exist *within the deterministic option*.

## Regimes & covariates to benchmark

(k, n): high-k/low-n (k≈1–5k, n≈20–50); low-k/high-n (k≈5–20, n≈2–10k);
high-k/high-n (k≈500–2k, n≈1–5k). **Covariates that matter as much as (k,n):**
- **Tree shape / height h**: balanced (h~log n) vs pectinate (h~n). Drives
  exact & Jansson. *Must vary.* (Existing benches use random n=100 only.)
- **Concordance**: identical / highly-concordant / random-independent. Drives
  #distinct splits (hashmap size; Jansson insert-delete churn).
- **p**: 0.5 and 0.75.

## Phases

- **P0 — Scaffold + harness.** dev/profiling/ skeleton; deterministic
  cross-regime benchmark grid (hashed vs exact first); upgrade
  verify-consensus.R into the correctness gate (hashed == exact == ape across a
  big grid + adversarial fixtures). Every new method must pass this gate.
- **P1 — Profile & optimise existing.** profvis the R wrapper (RenumberTips,
  Preorder, metadata strip, splits_to_edge, RootTree — likely dominate at
  low-n/high-k); VTune the C++ core. Verify each fix with std::chrono micro-bench.
  Prime candidates: incremental bitmask-OR for `exact` (O(knh/W)); ClusterTable
  ctor cost (k constructions); hashmap reserve/flat-map for hashed; trim R wrapper.
- **P2 — Implement Jansson.** Fetch [17] for One_Way_Compatible/Merge_Trees (or
  reconstruct + prove). Add internal third path. Gate through the harness — must
  equal hashed/exact/ape on the full grid + adversarial fixtures.
- **P3 — Benchmark all three across grid; analyse crossovers.** Data table +
  plots; identify winners per cell.
- **P4 — Integrate decision + document.** Wire dispatch (expected: hashed
  default; deterministic fallback = winner of {exact, Jansson}, maybe
  shape-switched). Update NEWS, roxygen, memory note.

## Open questions for advisor

1. Scope: full Jansson is a large, error-prone C++ effort whose likely payoff
   (per the paper's own timings) is "confirm hashed wins." Stage it behind the
   harness + exact-optimisation + benchmark, and make *shipping* it conditional
   on data showing exact is a real deterministic-path bottleneck? (Still
   implement+verify it to GET the data, per the explicit request.)
2. Agree Jansson competes with `exact`, not `hashed`?
3. Any correctness trap in the standard-majority Phase-2 derivation above?
4. Anything in the benchmark design that would let a wrong conclusion slip
   through (e.g. measuring only random trees, or core-only timing hiding R
   overhead)?
