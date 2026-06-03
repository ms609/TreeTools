# Profiling log — Consensus

last_focus: 1

---

## Round 1 — deferred materialisation + Jansson decision (2026-06-02)

**Optimisation C-001 (DONE, gated, measured).** Hashed counter now defers
split-pattern materialisation: each distinct split keeps a `(tree,L,R)` witness;
only splits reaching `thresh` are materialised (consensus), or all (split_freq).
Gate: 576 checks PASS. Measured (compare-grids.R): hashed median x1.23, with the
wins where it was weakest — **tall(10000,5) x13.0**, tall(5000,10) x8.6,
rand(10000,5) x2.5, tall/rand high-k/high-n x1.5–2.9. nsplits identical in every
cell (no behaviour change). exact untouched (x1.02).

**Jansson decision (bound re-run vs OPTIMISED hashed, + concord70).**
- Concordant (realistic) high-k/high-n: 2*strict ≥ hashed (x1.24–1.32) ⇒
  **Jansson can't win**. Same for all high-k/low-n.
- **Moderate conflict** (realistic gene-tree / bootstrap-with-rogues: base + ~5%
  rogue taxa relocated/replicate) at high-k/high-n: 2*strict ≥ hashed
  (2000,500: 1.54; 1000,1000: 1.42) ⇒ **Jansson can't win**. low-k/high-n
  (5000,10): 1.09 ⇒ loses.
- Only NOT-proven cells: (10000,5) (any concordance) — but <60ms absolute, a 2×
  gap is irrelevant; and high-k/high-n EXTREME conflict (rand/tall) — degenerate
  near-star consensus, not a real input.
- **DECISION (advisor-confirmed): do NOT implement full Jansson.** Proven to lose
  (2*strict ≥ hashed) in every realistic, time-meaningful regime. Rationale:
  user's "avoid redundant code"; high correctness risk (FACT ships broken
  majority); collision risk waived. Honest scope: this is a rigorous *lower-bound*
  result where it concludes — not a direct Jansson measurement; the un-proven
  cells are negligible-time or degenerate. Re-openable (algorithm + subroutine
  specs in PLAN-consensus.md / DECISION.md). See **DECISION.md** (user-facing).

**Advisor passes (2):** (1) flagged compare-grids' nsplits check is COUNT not
CONTENT ⇒ added high-n (n=1000/3000) hashed==exact split-SET content cells;
confirmed verify-consensus.R green; added moderate-conflict generator. (2) flagged
the symmetric gap — the witness change also rewired SplitFrequency's hashed
assembly, only checked at n≤30 ⇒ added high-n (n=2000/5000) SplitFrequency
hashed==exact split-SET+COUNT cells. Gate now 590 checks PASS. Every modified
code path gated at shipping scale (+ test-consensus 8/8, test-Support 6/6).
Advisor also corrected its earlier "commit" → do NOT commit (on main, not asked).

**Next (if continuing):** C-002 wrapper fast-path (high-k/low-n, 57% overhead) —
the remaining real throughput win, in a target regime. Shared-code (RenumberTips)
risk ⇒ gate with full test suite.

---

## Round 0 — scaffold + baseline + Jansson lower bound (2026-06-02)

**Setup.** Installed -O2 -g build into scratch lib (`.libpath.txt`); all drivers
load from it (not load_all). Built the correctness gate
(`correctness-gate.R`, 576 checks: hashed==exact everywhere; ape agrees at p=0.5)
and the benchmark grid (`drivers/consensus-grid.R` → `results-grid.csv`,
summarised in `baselines.md`).

**Correctness fixes made en route (real bugs):**
- gate `which.min(<character>)` → `integer(0)` silently skipped EVERY ape check;
  fixed to tip-index 1. Hard hashed==exact gate was still valid (both used the
  same broken-but-consistent canonicaliser).
- `dev/red-team/verify-consensus.R` called `consensus_tree(forest, p, hash=FALSE)`
  but the export arg is `exact` → "unused argument" → that adversarial gate had
  been DEAD since the rename. Fixed to `exact = TRUE` (lines 49, 77).
- Threshold-convention finding (surfaced to MRS, not changed): TreeTools keeps
  `count > k*p`; ape keeps `>= k*p` for p>0.5; roxygen says "p or more" (>=) but
  code is strict (>). Out of scope here.

**Baseline result.** hashed beats exact in all 80 cells (1.1x–3x). R-wrapper is
57% of end-to-end at high-k/low-n (30,5000). See `baselines.md`.

**Jansson lower bound** (`drivers/jansson-bound.R`; advisor-validated: Jansson ≥
2×strict-path). Verdict per cell:
- low-n/high-k: 2*strict ≥ hashed (ratio 1.1–1.4) ⇒ **Jansson can't win** there.
- HIGH-n (all shapes): 2*strict ≪ hashed (ratio 0.07–0.51) ⇒ **INCONCLUSIVE**.
  Strict is 5–30× faster than hashed-majority at high-n because hashed handles
  ~k·n DISTINCT splits (big hash map + per-split materialisation) while strict
  matches against ONE reference's O(n) clusters. Jansson's candidate tree is also
  O(n) clusters ⇒ it plausibly shares strict's advantage and could BEAT hashed at
  high-n. **Must implement Jansson to settle this** (user asked; advisor concurs).

**Caveat noticed:** cost structure is subtler than "materialisation dominates"
(identical is *slower* than rand at (10000,5) despite fewer distinct splits) ⇒
must VTune the core, not guess, before optimising it.

**Next:** (R1) profvis the R wrapper at high-k/low-n → fast-path. (R2) VTune the
hashed/exact core at high-n to find the true hotspot. (R3) implement Jansson
(R-prototype-first, gated) and benchmark vs optimised hashed.
