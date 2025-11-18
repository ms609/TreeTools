# Write a phylogenetic tree in Newick format

`as.Newick()` creates a character string representation of a
phylogenetic tree, in the Newick format, using R's internal tip
numbering. Use
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md)
to ensure that the internal numbering follows the order you expect.

## Usage

``` r
as.Newick(x)

# S3 method for class 'phylo'
as.Newick(x)

# S3 method for class 'list'
as.Newick(x)

# S3 method for class 'multiPhylo'
as.Newick(x)
```

## Arguments

- x:

  Object to convert to Newick format. See Usage section for supported
  classes.

## Value

`as.Newick()` returns a character string representing `tree` in Newick
format.

## See also

- Retain leaf labels:
  [`NewickTree()`](https://ms609.github.io/TreeTools/dev/reference/NewickTree.md)

- Change R's internal numbering of leaves:
  [`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md)

- Write tree to text or file:
  [`ape::write.tree()`](https://rdrr.io/pkg/ape/man/write.tree.html)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
trees <- list(BalancedTree(1:8), PectinateTree(8:1))
trees <- lapply(trees, RenumberTips, 1:8)
as.Newick(trees)
#> [1] "(((0,1),(2,3)),((4,5),(6,7)));" "(((((((0,1),2),3),4),5),6),7);"
```
