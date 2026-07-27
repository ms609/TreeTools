# Profiling focus areas — Consensus

Ranked by (per-call cost) × (call frequency on the consensus hot path).
See PLAN-consensus.md for the broader project.

| # | Area | Files | Why hot | Baseline cost | Last profiled | Status |
|---|------|-------|---------|---------------|---------------|--------|
| 1 | `count_splits_hashed` (majority default) | src/consensus.cpp:267 | One postorder pass/tree; the default majority counter. O(kn). | tbd | — | NEW |
| 2 | `count_splits_exact` (deterministic) | src/consensus.cpp:86 | Builds packed bitkey by L..R loop/cluster ⇒ O(k·n·h); slow for tall trees. | tbd | — | NEW |
| 3 | `ClusterTable` constructor | inst/include/TreeTools/ClusterTable.h:352 | Built once per tree per call ⇒ k constructions; `root_on_node` + allocs. | tbd | — | NEW |
| 4 | `Consensus()` R wrapper | R/Consensus.R:46 | Per-call R overhead: metadata strip (lapply), NTip, RenumberTips, Preorder, splits_to_edge, RootTree. May dominate at low-n/high-k. | tbd | — | NEW |
| 5 | strict single-ref path | src/consensus.cpp:165 | p=1 path; already O(kn) optimal. Profile only to confirm AT-LIMIT. | tbd | — | NEW |
| 6 | `for_each_internal_node` primitive | src/consensus.cpp:44 | Shared postorder/stack arithmetic backing #1,#2,#5. | tbd | — | NEW |
| 7 | Jansson core (to be written) | src/consensus.cpp (new) | Deterministic O(kn) candidate-tree path; compare vs #1,#2. | tbd | — | NEW |

Parked (called once/session; not in rotation): IO, label checking.
