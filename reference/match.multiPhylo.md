# Tree matching

[`match()`](https://ms609.github.io/TreeTools/reference/match.Splits.md)
returns a vector of the positions of (first) matches of trees in its
first argument in its second. `%in%` is a more intuitive interface as a
binary operator, which returns a logical vector indicating whether there
is a match or not for each tree in its left operand.

## Usage

``` r
# S4 method for class 'phylo,phylo'
match(x, table, nomatch = NA_integer_, incomparables = NULL)

# S4 method for class 'multiPhylo,phylo'
match(x, table, nomatch = NA_integer_, incomparables = NULL)

# S4 method for class 'phylo,multiPhylo'
match(x, table, nomatch = NA_integer_, incomparables = NULL)

# S4 method for class 'multiPhylo,multiPhylo'
match(x, table, nomatch = NA_integer_, incomparables = NULL)

# S4 method for class 'multiPhylo,multiPhylo'
x %in% table

# S4 method for class 'multiPhylo,phylo'
x %in% table

# S4 method for class 'phylo,multiPhylo'
x %in% table

# S4 method for class 'phylo,phylo'
x %in% table
```

## Arguments

- x, table:

  Object of class `phylo` or `multiPhylo`.

- nomatch:

  Integer value that will be used in place of `NA` in the case where no
  match is found.

- incomparables:

  Ignored. (Included for consistency with generic.)

## Value

[`match()`](https://ms609.github.io/TreeTools/reference/match.Splits.md)
returns an integer vector specifying the position in `table` that
matches each element in `x`, or `nomatch` if no match is found.

## See also

Corresponding base functions are documented in
[`match()`](https://rdrr.io/r/base/match.html).

Other utility functions:
[`ClusterTable`](https://ms609.github.io/TreeTools/reference/ClusterTable.md),
[`ClusterTable-methods`](https://ms609.github.io/TreeTools/reference/ClusterTable-methods.md),
[`Hamming()`](https://ms609.github.io/TreeTools/reference/Hamming.md),
[`MSTEdges()`](https://ms609.github.io/TreeTools/reference/MSTEdges.md),
[`SampleOne()`](https://ms609.github.io/TreeTools/reference/SampleOne.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/reference/TipTimedTree.md),
[`UnshiftTree()`](https://ms609.github.io/TreeTools/reference/UnshiftTree.md),
[`as.multiPhylo()`](https://ms609.github.io/TreeTools/reference/as.multiPhylo.md),
[`sapply64()`](https://ms609.github.io/TreeTools/reference/sapply64.md),
[`sort.multiPhylo()`](https://ms609.github.io/TreeTools/reference/sort.multiPhylo.md)

## Examples

``` r
tree1 <- BalancedTree(7)
trees <- c(PectinateTree(7), BalancedTree(7))

match(tree1, trees)
#> [1] 2
```
