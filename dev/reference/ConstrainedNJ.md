# Constrained neighbour-joining tree

Constructs an approximation to a neighbour-joining tree, modified in
order to be consistent with a constraint. Zero-length branches are
collapsed at random.

## Usage

``` r
ConstrainedNJ(dataset, constraint, weight = 1L, ratio = TRUE, ambig = "mean")
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree.

- constraint:

  Either an object of class `phyDat`, in which case returned trees will
  be perfectly compatible with each character in `constraint`; or a tree
  of class `phylo`, in which each node in `constraint` will occur in the
  returned tree. See
  [vignette](https://ms609.github.io/TreeSearch/articles/tree-search.html)
  for further examples.

- weight:

  Numeric specifying degree to up-weight characters in `constraint`.

- ambig, ratio:

  Settings of `ambig` and `ratio` to be used when computing
  [`Hamming()`](https://ms609.github.io/TreeTools/dev/reference/Hamming.md)
  distances between sequences.

## Value

`ConstrainedNJ()` returns a tree of class `phylo`.

## See also

Other tree generation functions:
[`GenerateTree`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md),
[`NJTree()`](https://ms609.github.io/TreeTools/dev/reference/NJTree.md),
[`TreeNumber`](https://ms609.github.io/TreeTools/dev/reference/TreeNumber.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
dataset <- MatrixToPhyDat(matrix(
  c(0, 1, 1, 1, 0, 1,
    0, 1, 1, 0, 0, 1), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
constraint <- MatrixToPhyDat(
  c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))
plot(ConstrainedNJ(dataset, constraint))
```
