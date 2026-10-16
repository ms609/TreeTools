# Red-team focus areas — TreeTools

Rotation table for the `/red-team` skill. One area per invocation. Each area
pairs the R surface with the C++ that backs it, so the finder sees the
`.Call`/Rcpp interface contract in one pass.

Definition only — no tier and no status column. Round records live in GitHub
Discussions (one category per area); findings live as issues labelled
`red-team` + `sev:*` + `area:N` in `agent-issues/TreeTools`.

`area:N` labels are per row. **Adding a row means creating its label**, or
filing against it fails.

Area 13 audits the rotation itself rather than the package.

TreeTools is the foundation of the stack: PlotTools, TreeDist, TreeSearch,
Quartet and Rogue consume its exported R functions, and several link against
`inst/include/TreeTools/*.h` directly. A defect here has the widest blast
radius of anywhere in the stack — weight severity accordingly.

| # | Area | Files | Key questions |
|---|------|-------|---------------|
| 1 | **Downstream C++ headers (ABI surface)** | `inst/include/TreeTools/*.h` (57 kB, 11 files) | These are compiled into TreeDist, TreeSearch and Quartet, so a silent semantic change ships as a defect in packages whose tests never run here. Bounds and ownership in `SplitList.h`, `ClusterTable.h`, `keep_tip.h`. Does `assert.h` degrade safely under `NDEBUG`? Which headers have *no* direct test in this package? |
| 2 | **Tree numbering / Preorder** | `R/tree_numbering.R` (23 kB), `inst/include/TreeTools/renumber_tree.h`, `src/ape_reorder.{cpp,h}` | The Preorder guarantee is load-bearing for the whole stack. Is the postorder/preorder assertion actually checked (#266)? Behaviour on already-ordered, malformed and non-binary input. Classification of "Preorder" (#93, #229). Edge-matrix-only input (#27). |
| 3 | **Splits** | `R/Splits.R`, `R/SplitFunctions.R`, `src/splits.cpp` (21 kB), `inst/include/TreeTools/{SplitList,edge_to_splits}.h`, `src/{tips_in_splits,first_matching_split}.cpp` | `raw` bitset packing: the multi-word path at >64 tips (and the 8-bit boundary). Trailing-bit hygiene in the final word — do unused bits leak into comparisons? `in.Splits()` deprecation (#34). Duplicate and trivial splits. |
| 4 | **Consensus** | `R/Consensus.R`, `src/consensus.cpp` | Partly audited (C-001…C-005 in the legacy profiling ledger; C-003 and C-005 believed open — verify before refiling). Threshold semantics at `p > 0.5` were fixed once; are they pinned? Hash collision handling in `count_splits_hashed`. |
| 5 | **Tip manipulation** | `R/DropTip.R`, `R/AddTip.R`, `src/kept_vertices.cpp`, `inst/include/TreeTools/keep_tip.h`, `R/KeptVerts.R`, `R/KeptPaths.R` | Attribute preservation is the open cross-cutting question: which of `edge.length`, `node.label`, `root.edge` survive, and is dropping them ever correct here? Dropping all/all-but-one/zero tips. Tips absent from the tree. Polytomy collapse when an internal node loses degree. |
| 6 | **Tree generation** | `R/tree_generation.R`, `R/ImposeConstraint.R`, `R/ConsistentSplits.R` | Uniformity of `RandomTree` across the topology space; RNG reproducibility under `set.seed`. Constraint imposition when constraints are contradictory or already satisfied. `n = 0, 1, 2, 3` boundaries. |
| 7 | **Tree numbers / combinatorics** | `R/TreeNumber.R` (17 kB), `R/Combinatorics.R`, `R/BigInteger.R`, `src/int_to_tree.cpp`, `inst/include/TreeTools/tree_number.h` | Integer overflow and the documented leaf ceiling — "This many leaves cannot be supported" is a live external report (#141). Is the limit the right one, correctly stated, and reached with a clear error rather than a wrong answer? Round-tripping `as.phylo`/`as.integer` at the boundary. |
| 8 | **Parsers and IO** | `R/parse_files.R` (43 kB — the largest file in the package), `R/ReadTntTree.R`, `R/ReadMrBayes.R`, `R/PhyToString.R`, `src/as_newick.cpp`, `src/fast_paste.cpp` | The only untrusted-input surface in the package. Malformed Newick/NEXUS: unbalanced parens, missing semicolon, comments, quoted labels with punctuation, CRLF, BOM, empty file, truncated mid-token. Multiline TNT (#99). Unbounded recursion or allocation driven by file content. |
| 9 | **Tree properties and information** | `R/tree_properties.R`, `R/Information.R`, `R/tree_information.R`, `src/node_depth.cpp`, `R/Treeness.R`, `R/EdgeRatio.R` | Functions that *should* consume `edge.length` — do they, and do they error clearly when it is absent rather than returning a topological answer silently? Numerical stability of information measures at large n. |
| 10 | **Tree shape and balance** | `R/tree_shape.R`, `R/RUtreebalance.R`, `R/TotalCopheneticIndex.R`, `R/Cherries.R`, `R/Stemwardness.R`, `src/{tree_shape,n_cherries}.cpp` | Third-party-derived code (RUtreebalance) — does it share this package's input conventions? Behaviour on polytomies and on unrooted input. Balance index definitions against their published sources. |
| 11 | **Rearrangement, paths and distances** | `R/tree_rearrangement.R` (19 kB), `R/PathLengths.R`, `R/mst.R`, `R/Decompose.R`, `src/{path_lengths,descendant_edges,minimum_spanning_tree}.cpp` | `CollapseNode()` performance (#113). `MatchNodes()` optimisation (#161). Merging `AllDescendantEdges()` into `DescendantEdges()` (#31). Do path-length routines handle missing or negative edge lengths? |
| 12 | **Display, support and annotation** | `R/RoguePlot.R`, `R/Support.R`, `R/PaintTree.R`, `R/tree_display.R`, `R/sort.R`, `R/MatchNodes.R` | Node/edge index alignment between a tree and its plot — an off-by-one here mislabels support values silently. Snapshot tests exist (`_snaps/`); do they pin behaviour or just pixels? |

| 13 | **Red-team process meta-review** | `dev/red-team/focus-areas.md`, `dev/next-issue/config.md`, `.github/workflows/sync-upstream.yml`, `AGENTS.md`, the issue and Discussion record itself | Is the rotation finding defects or performing diligence? Per-area yield: which areas have produced nothing across several rounds, and should they be merged, retired or re-scoped? Are `sev:*` labels calibrated, or has everything drifted to `sev:med`? Are issues actually closing via `Fixes #N` on merge into `agent`, or accumulating as stale `in-progress`? Is `sync-upstream.yml` still running, or silently disabled by the 60-day inactivity rule? Are findings being re-hunted because an earlier round's reopening condition was never recorded? This area audits the apparatus, not the package. |

## Not an area

**Attribute preservation across the API** (`edge.length`, `node.label`,
`root.edge`, tip-label order) is cross-cutting, not an area. It is tracked as a
`task` issue with a generated support matrix over every exported function that
returns a `phylo`/`multiPhylo`, and findings from it are filed against whichever
area owns the offending function.
