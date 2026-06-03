# Consensus profiling findings

Task-row format. P = priority (1 highest). Status: OPEN / DONE / WONTFIX.

| ID | P | Status | Verified Δ | Finding | Evidence |
|----|---|--------|-----------|---------|----------|
| C-001 | 1 | DONE | hashed median **x1.23**; tall(10000,5) **x13.0**, tall(5000,10) x8.6, rand(10000,5) x2.5, tall(2000,500) x2.9; nsplits identical in every cell. | [Optimise] hashed: DEFER split-pattern materialisation until thresholding — keep a `(tree,L,R)` witness per distinct split, build the packed pattern only for splits reaching `thresh`. | At high n the dominant hashed cost is per-distinct-split work, most wasted (few/no splits survive on low-concordance input). Strict path (no materialisation) is 5–30× faster than hashed-majority at high n. Implemented in src/consensus.cpp (SplitWitness + materialise_witness). Gate: 590 checks PASS (hashed==exact==ape@0.5; incl. consensus content n≤3000 and SplitFrequency split+count n=2000/5000); test-consensus 8/8, test-Support 6/6, verify-consensus green. Speedup verified: compare-grids.R. |
| C-002 | 2 | DONE | RenumberTips **1.8–9.0×** (30,5000: 0.335→0.037s 9.0×; 50,2000 6.7×; 100,1000 6.3×; 500,300 2.6×; 2000,100 1.8×). Results identical. | [Optimise] `RenumberTips()` on an unlabelled multiPhylo/list now relabels the whole forest in ONE C++ pass (`renumber_tips_to`), replacing the per-tree R loop; no-op pass-through for trees already in target order. Shared by Consensus() and all RenumberTips callers. Returns NULL → existing R fallback for numeric tipOrder / differing label sets / non-phylo elements, so edge-case behaviour is preserved exactly. | **"Tell without looking?"** Only a *labelled* multiPhylo (TipLabel attr) guarantees a shared order without inspecting each tree; otherwise inspection is unavoidable, but now in C++ (and a matching tree is passed through untouched). Gate: full suite GREEN (296); correctness-gate 590; new equivalence test (C++ == per-tree across mixed orderings). LESSON: the existing RenumberTips test caught that the C++ list initially dropped `names` (the lapply path preserves them) — fixed by copying input names onto the result. Bench: drivers/renumber-bench.R. |
| C-004 | 2 | DONE | strip→`c(trees)`: forest201.80 3.38→2.54ms (**-25%**), forestMaj 2.60→2.10 (-19%), forest21.260 1.23→1.12 (-9%), forest1k.888 72.1→68.0 (-6%, strict core is 92%). Output byte-identical. | [Optimise] Consensus() metadata-strip replaced by bare `trees <- c(trees)`. The old `lapply(c(trees), \\(tr){tr$edge.length<-NULL; tr$node.label<-NULL; tr})` copied EVERY tree; the consensus core reads only edge+tip.label and the output is rebuilt from scratch, so neither field can change the answer. `c()` alone still materialises a labelled multiPhylo's shared labels onto each tree (the load-bearing job). **Floor finding (key):** the naive presence-guard (`if(!is.null(...))`) was a near-wash — the `is.null` reads cost ~as much as the assignments they skip (0.82→0.75ms), and REGRESSED the metadata-present case; the win comes from skipping the lapply entirely (strip 0.82→c() 0.11ms, ~7x). | gha-bench-diagnose.R: strip was 19–29% of the SMALL forests, 6% of forest1k.888. Equivalence: c004-equiv.R 108 cells (plain/unaligned/edge_len/node_lab/both/labelled/lab_el/list/diff_size × p × check × hash) byte-identical pre/post; core RawMatrix metadata-independent 36 cells; new test-Consensus.R block (edge.length/node.label/labelled, check T/F, p∈{.5,⅔,1}); full suite 4192 PASS; gate 590. |
| C-003 | 3 | OPEN | — | [Optimise] hashed: large `unordered_map<Hash128>` churn dominates the LOW-height high-distinct (`rand`) case (deferred materialisation does NOT help there — clusters are tiny). Try a flat open-addressing map / better reserve. Asymptotic floor is Ω(#distinct) ops; only Jansson's O(n) candidate tree removes it. | jansson-bound.R: rand high-n strict ≪ hashed; the gap is hashmap ops, not materialisation. Realistic (concordant) inputs have few distinct splits ⇒ small map ⇒ low priority. |
| C-005 | 3 | OPEN | — | [Optimise] `NTip(trees)` in the Consensus() repeat-loop is now the next-largest wrapper cost after the C-004 strip removal: forest201.80 0.73ms (22% of total), forest1k.888 3.4ms (5%). It builds a full k-vector of tip counts every loop iteration just to test `length(unique(.)) > 1`. Could early-exit / cache / vectorise via the materialised list. Consensus()-local, low risk. | gha-bench-diagnose.R breakdown. |

## Why the GHA `bench-Consensus.R` cells show NSD after C-001/C-002 (2026-06-03)

The shipped wins don't appear in the GHA suite because **none of its 5 cells
exercise the regimes the optimisations target** (gha-bench-diagnose.R):

- All four forests (`as.phylo(int, n)`) are **metadata-free, already in matching
  tip order, and high-concordance** (consecutive tree-integers ⇒ few distinct
  splits).
- **C-002 (RenumberTips):** 4/5 cells pass `check = FALSE` (partial-matches
  `check.labels = FALSE`) ⇒ RenumberTips never runs. The one cell that runs it
  (`forest21.260`, the −6.96% blip) is k=21 and already-ordered, so the C++
  path's fixed setup doesn't amortise vs the old R loop — real, but NSD.
- **C-001 (deferred materialisation):** 3/5 cells use default `p = 1` (strict)
  ⇒ the single-reference path, which never calls `count_splits_hashed`. The 2
  majority cells (`p = 0.5`) are 1–3 ms and high-concordance ⇒ nothing to defer.
- The strict C++ core dominates the one big cell (`forest1k.888`: 66 of 72 ms).

To make GHA *demonstrate* the wins, add cells that hit the regimes (OFFERED to
user, not yet added): a high-k unlabelled forest with **shuffled** tip order
(real relabel → C-002), and a tall **low-concordance `p < 1`** forest (→ C-001).
C-004 does move the small cells (strip was 19–29% of them).

## Threshold semantics at `p > 0.5` — open question for maintainer (2026-06-03)

Prompted by "how are p = 0.5 ties resolved?". Two distinct findings (repro:
`drivers/tie-check.R`):

- **p = 0.5 ties are handled correctly.** `thresh = floor(n·p) + 1`, so a split
  in exactly n/2 trees (`count = n/2 < thresh`) is **dropped**. Two conflicting
  50% splits are both dropped ⇒ valid tree. Provably correct (any two retained
  splits each occur in > n/2 trees ⇒ counts sum > n ⇒ co-occur ⇒ compatible) and
  **matches `ape::consensus`** (which bumps `p <- 0.5000001` to the same effect).
  Now pinned by `ApeTest(tie, 0.5)` in test-Consensus.R (the even-n tie fixture
  was previously only in the hashed==exact test — internal counters only, no
  oracle).

- **p > 0.5 is off-by-one *stricter* than the docstring and `ape`.** The
  docstring promises "a proportion `p` or more"; `ape` keeps `count >= p·n`
  (`bs >= p * ntree`). TreeTools keeps `count >= floor(n·p) + 1`, i.e. *strictly
  more than* `floor(n·p)`. These **diverge only when `n·p` computes to an exact
  integer k**: ape keeps `count = k`, TreeTools requires `count >= k+1`. E.g.
  n=3, p=2/3: a clade in exactly 2 of 3 trees is kept by ape, dropped by
  TreeTools. Both choices are safe (valid tree) and both are FP-fragile at the
  boundary. **Convention + doc-consistency decision is the maintainer's:** either
  align the code to "p or more" (a tolerance-aware integer threshold, *not* a
  copy of ape's fragile `>=`), or correct the docstring to "more than `p`".
  Not pre-built. No code/threshold change made.

## Jansson decision (tracking)

Open question: does deterministic Jansson (Maj_Rule_Plus, O(kn)) beat OPTIMISED
hashed in any *realistic* regime? Realistic consensus inputs are concordant (few
distinct splits) where hashed is O(kn) tiny-constant and wins. Jansson's only
edge is high-n LOW-concordance (≈ independent random trees — not a real consensus
input). Lower bound: Jansson ≥ 2×strict. Decision pending grid-v2 + a realistic
high-concordance scenario. See log.md / PLAN-consensus.md.
