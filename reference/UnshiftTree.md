# Add tree to start of list

`UnshiftTree()` adds a phylogenetic tree to the start of a list of
trees. This is useful where the class of a list of trees is unknown, or
where names of trees should be retained.

## Usage

``` r
UnshiftTree(add, treeList)
```

## Arguments

- add:

  Tree to add to the list, of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- treeList:

  A list of trees, of class `list`,
  [`multiPhylo`](https://rdrr.io/pkg/ape/man/multiphylo.html), or, if a
  single tree, [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`UnshiftTree()` returns a list of class `list` or `multiPhylo`
(following the original class of `treeList`), whose first element is the
tree specified as \`add.

## Details

Caution: adding a tree to a `multiPhylo` object whose own attributes
apply to all trees, for example trees read from a Nexus file, causes
data to be lost.

## See also

[`c()`](https://rdrr.io/r/base/c.html) joins a tree or series of trees
to a `multiPhylo` object, but loses names and does not handle lists of
trees.

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/reference/as.multiPhylo.md),
[`match,phylo,phylo-method`](https://ms609.github.io/TreeTools/reference/match.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
forest <- as.phylo(0:5, 6)
tree <- BalancedTree(6)

UnshiftTree(tree, forest)
#> 7 phylogenetic trees
UnshiftTree(tree, tree)
#> 2 phylogenetic trees
```
