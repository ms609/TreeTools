# Ratio of external:internal edge length

Reports the ratio of tree length associated with external edges (i.e.
edges whose child is a leaf) and internal edges. Where tree length is
dominated by internal edges, variation between tips is dominantly
controlled by phylogenetic history.

## Usage

``` r
EdgeRatio(x)

# S3 method for class 'phylo'
EdgeRatio(x)
```

## Arguments

- x:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`EdgeRatio()` returns a numeric specifying the ratio of external to
internal edge length (\> 1 means the length of a tree is predominantly
in external edges), with attributes `external`, `internal`, and `total`
specifying the total length associated with edges of that nature.

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
[`TreeIsRooted()`](https://ms609.github.io/TreeTools/reference/TreeIsRooted.md),
[`Treeness()`](https://ms609.github.io/TreeTools/reference/Treeness.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
