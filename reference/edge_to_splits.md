# Efficiently convert edge matrix to splits

Wrapper for internal C++ function for maximum efficiency. Improper input
may crash R. Behaviour not guaranteed. It is advisable to contact the
package maintainers before relying on this function.

## Usage

``` r
edge_to_splits(
  edge,
  edgeOrder,
  tipLabels = NULL,
  asSplits = TRUE,
  nTip = NTip(edge),
  ...
)
```

## Arguments

- edge:

  A matrix with two columns, with each row listing the parent and child
  node of an edge in a phylogenetic tree. Property `edge` of objects of
  class `phylo`.

- edgeOrder:

  Integer vector such that `edge[edgeOrder, ]` returns a postorder
  ordering of edges.

- tipLabels:

  Character vector specifying sequence in which to order tip labels.
  Label order must (currently) match to combine or compare separate
  `Splits` objects.

- asSplits:

  Logical specifying whether to return a `Splits` object, or an
  unannotated two-dimensional array (useful where performance is
  paramount).

- nTip:

  Integer specifying number of leaves in tree.

- ...:

  Presently unused.

## Value

`edge_to_splits()` uses the same return format as
[`as.Splits()`](https://ms609.github.io/TreeTools/reference/Splits.md).

## See also

[`as.Splits()`](https://ms609.github.io/TreeTools/reference/Splits.md)
offers a safe access point to this function that should be suitable for
most users.
