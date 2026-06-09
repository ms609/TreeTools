# Wrapper for internal C function `root_on_node()`

Direct entry point to `root_on_node()`; recommended for expert use only.
[`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md)
checks that input is properly formatted and is recommended for general
use.

## Usage

``` r
root_on_node(phy, outgroup)
```

## Arguments

- phy:

  Minimally, a named list with entries `edge` and `Nnode`, in the format
  of equivalent entries in a tree of class `phylo`. `edge.length` will
  also be considered if supplied. The root node must be numbered
  `n_tip + 1`, per the `phylo` convention.

- outgroup:

  Integer specifying index of leaf or node to set as the outgroup.

## Value

`root_on_node()` returns `phy` rooted on the specified node, in
preorder.

## Details

`phy` may be supplied in any valid edge order; it is preordered
internally, so there is no need to call
[`Preorder()`](https://ms609.github.io/TreeTools/dev/reference/Reorder.md)
first. The returned tree is always in preorder, with node numbers
reassigned accordingly.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
