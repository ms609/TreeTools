# Is tree rooted?

`TreeIsRooted()` is a fast alternative to
[`ape::is.rooted()`](https://rdrr.io/pkg/ape/man/root.html).

## Usage

``` r
TreeIsRooted(tree)
```

## Arguments

- tree:

  A phylogenetic tree of class `phylo`.

## Value

`TreeIsRooted()` returns a logical specifying whether a root node is
resolved.

## See also

Other tree properties:
[`Cherries()`](https://ms609.github.io/TreeTools/reference/Cherries.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/reference/ConsensusWithout.md),
[`LongBranch()`](https://ms609.github.io/TreeTools/reference/LongBranch.md),
[`MatchEdges()`](https://ms609.github.io/TreeTools/reference/MatchEdges.md),
[`NSplits()`](https://ms609.github.io/TreeTools/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/reference/NTip.md),
[`NodeNumbers()`](https://ms609.github.io/TreeTools/reference/NodeNumbers.md),
[`PathLengths()`](https://ms609.github.io/TreeTools/reference/PathLengths.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/reference/TipLabels.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
TreeIsRooted(BalancedTree(6))
#> [1] TRUE
TreeIsRooted(UnrootTree(BalancedTree(6)))
#> [1] FALSE
```
