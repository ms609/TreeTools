# Consensus baselines (installed -O2 -g build; R 4.7.0; this machine)

Core timing = `TreeTools:::consensus_tree(forest, p=0.5)` on RenumberTips+Preorder'd
forests (min of 4 reps after warmup). Full grid in `results-grid.csv`. Driver:
`drivers/consensus-grid.R`. Run 1: 2026-06-02.

## Headline: hashed dominates exact in EVERY cell (1.1x–3x). No crossover.

| regime (n,k)        | scenario  | hashed min(s) | exact min(s) | exact/hashed |
|---------------------|-----------|---------------|--------------|--------------|
| highK_lowN (30,5000)| rand      | 0.0647        | 0.0723       | 1.12 |
| highK_lowN (50,2000)| tall      | 0.0450        | 0.0495       | 1.10 |
| lowK_highN (10000,5)| rand      | 0.0422        | 0.1269       | 3.01 |
| lowK_highN (10000,5)| tall      | 0.1885        | 0.2885       | 1.53 |
| lowK_highN (5000,10)| rand      | 0.0291        | 0.0743       | 2.55 |
| highK_highN(2000,500)| rand     | 0.6648        | 1.2550       | 1.89 |
| highK_highN(2000,500)| tall     | 1.2358        | 1.9502       | 1.58 |
| highK_highN(1000,1000)| concord70| 0.2308       | 0.3917       | 1.70 |

## Cost decomposition (what actually dominates)

- **R-wrapper overhead** (`Consensus()` end-to-end minus core), by regime:
  | (n,k)        | Consensus() | core   | wrapper | wrapper % |
  |--------------|-------------|--------|---------|-----------|
  | (30, 5000)   | 0.1480      | 0.0636 | 0.0844  | **57%**   |
  | (100, 100)   | 0.0061      | 0.0038 | 0.0023  | 37%       |
  | (2000, 10)   | 0.0112      | 0.0086 | 0.0026  | 23%       |
  | (1000, 1000) | 0.6081      | 0.5814 | 0.0267  | 4%        |
  ⇒ At **high-k/low-n** the R wrapper (RenumberTips / metadata strip / Preorder)
    is the bottleneck, NOT the counter. This is a named target regime.

- **ClusterTable construction** dominates at high-n: tall-vs-rand gap at fixed
  (n,k) for hashed (e.g. n=10000,k=5: rand 0.042 vs tall 0.189) reflects the
  O(n·h) split-materialisation + construction work shared by both counters; the
  hashed/exact gap is secondary to this shared floor.

## Implication for Jansson

Jansson (deterministic O(kn)) competes with `exact` (deterministic), not hashed.
But (1) hashed already beats exact everywhere by ≤3x; (2) exact never blows up
(shared construction floor dominates); (3) Jansson structurally does ~2k
ClusterTable constructions + per-tree One_Way_Compatible + Merge_Trees + Day
match + delete/insert — several× hashed's single pass — and the authors' own C++
(Maj_Rule_Plus) benchmarks ~10x slower than TreeTools hashed. So Jansson is very
unlikely to win any regime. The high-leverage throughput wins are the **R
wrapper** (high-k/low-n) and **ClusterTable construction** (high-n).

## Threshold-convention note (separate finding, not a bug I introduced)

For 0.5 < p < 1, TreeTools keeps `count > k*p` (`thresh=floor(k*p)+1`); ape keeps
`count >= k*p`. They agree at p=0.5 and p=1. The roxygen says "p or more" (>=) but
the code is strict (>). Surfacing for MRS; out of scope to change here.
