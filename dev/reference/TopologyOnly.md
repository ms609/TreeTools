# Remove metadata from trees

`TopologyOnly()` removes all information from trees except for their
topologies and leaf labels.

## Usage

``` r
TopologyOnly(tree)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

## Value

Returns `tree`, with each tree in
[`Preorder`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md),
with edge lengths, node labels and other attributes removed.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
