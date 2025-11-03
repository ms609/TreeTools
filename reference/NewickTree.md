# Write Newick Tree

`NewickTree()` encodes a tree as a Newick-format string. This differs
from [`write.tree()`](https://rdrr.io/pkg/ape/man/write.tree.html) in
the encoding of spaces as spaces, rather than underscores.

## Usage

``` r
NewickTree(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

`NewickTree()` returns a character string denoting `tree` in Newick
format.

## See also

Use tip numbers, rather than leaf labels:
[`as.Newick`](https://ms609.github.io/TreeTools/reference/as.Newick.md)

## Examples

``` r
NewickTree(BalancedTree(LETTERS[4:9]))
#> [1] "(((D,E),F),((G,H),I));"
```
