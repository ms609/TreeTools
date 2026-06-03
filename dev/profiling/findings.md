# Consensus profiling findings

Task-row format. P = priority (1 highest). Status: OPEN / DONE / WONTFIX.

| ID | P | Status | Verified Δ | Finding | Evidence |
|----|---|--------|-----------|---------|----------|
| C-001 | 1 | DONE | hashed median **x1.23**; tall(10000,5) **x13.0**, tall(5000,10) x8.6, rand(10000,5) x2.5, tall(2000,500) x2.9; nsplits identical in every cell. | [Optimise] hashed: DEFER split-pattern materialisation until thresholding — keep a `(tree,L,R)` witness per distinct split, build the packed pattern only for splits reaching `thresh`. | At high n the dominant hashed cost is per-distinct-split work, most wasted (few/no splits survive on low-concordance input). Strict path (no materialisation) is 5–30× faster than hashed-majority at high n. Implemented in src/consensus.cpp (SplitWitness + materialise_witness). Gate: 590 checks PASS (hashed==exact==ape@0.5; incl. consensus content n≤3000 and SplitFrequency split+count n=2000/5000); test-consensus 8/8, test-Support 6/6, verify-consensus green. Speedup verified: compare-grids.R. |
| C-002 | 2 | OPEN | — | [Port/Optimise] Consensus() R wrapper: `RenumberTips` is 54% of the wrapper at high-k/low-n; the metadata-strip `lapply` copies every tree (pure overhead when no edge.length/node.label) and downgrades the forest to a plain list (forces per-tree RenumberTips). | Rprof at (n=30,k=5000): RenumberTips.phylo self 22% / total 54%; wrapper = 57% of end-to-end (baselines.md). **Two fixes, different scope:** (a) SAFE/localised — a short-circuit consistency check in Consensus() that SKIPS RenumberTips when all trees already share tree1's tip.label order (a genuine no-op then); ~0s overhead on mismatch, saves ~0.135s/call (~2× end-to-end) when it triggers, BUT only helps already-consistent input (ape::rtree shuffles label ORDER, so the bench's `rand`/`tall` don't trigger it; real bootstrap/posterior sets often do). (b) GENERAL/risky — route genuine relabelling through the C++ batch path `_TreeTools_renumber_tips_batch` (needs a labelled multiPhylo; helps inconsistent input too). Note: passing a precomputed char order is *slower* (TipLabels.character does more). Shared-code risk for (b) — gate with full test suite. |
| C-003 | 3 | OPEN | — | [Optimise] hashed: large `unordered_map<Hash128>` churn dominates the LOW-height high-distinct (`rand`) case (deferred materialisation does NOT help there — clusters are tiny). Try a flat open-addressing map / better reserve. Asymptotic floor is Ω(#distinct) ops; only Jansson's O(n) candidate tree removes it. | jansson-bound.R: rand high-n strict ≪ hashed; the gap is hashmap ops, not materialisation. Realistic (concordant) inputs have few distinct splits ⇒ small map ⇒ low priority. |

## Jansson decision (tracking)

Open question: does deterministic Jansson (Maj_Rule_Plus, O(kn)) beat OPTIMISED
hashed in any *realistic* regime? Realistic consensus inputs are concordant (few
distinct splits) where hashed is O(kn) tiny-constant and wins. Jansson's only
edge is high-n LOW-concordance (≈ independent random trees — not a real consensus
input). Lower bound: Jansson ≥ 2×strict. Decision pending grid-v2 + a realistic
high-concordance scenario. See log.md / PLAN-consensus.md.
