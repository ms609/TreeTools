# Generate a neighbour joining tree

`NJTree()` generates a rooted neighbour joining tree from a phylogenetic
dataset.

## Usage

``` r
NJTree(dataset, edgeLengths = FALSE, ratio = TRUE, ambig = "mean")
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree.

- edgeLengths:

  Logical specifying whether to include edge lengths.

- ambig, ratio:

  Settings of `ambig` and `ratio` to be used when computing
  [`Hamming()`](https://ms609.github.io/TreeTools/dev/reference/Hamming.md)
  distances between sequences.

## Value

`NJTree` returns an object of class `phylo`.

## See also

Other tree generation functions:
[`ConstrainedNJ()`](https://ms609.github.io/TreeTools/dev/reference/ConstrainedNJ.md),
[`GenerateTree`](https://ms609.github.io/TreeTools/dev/reference/GenerateTree.md),
[`TreeNumber`](https://ms609.github.io/TreeTools/dev/reference/TreeNumber.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo")
NJTree(Lobo.phy)
#> 
#> Phylogenetic tree with 48 tips and 47 internal nodes.
#> 
#> Tip labels:
#>   Tubiluchus_Priapulida, Cricocosmia, Aysheaia, Siberion, Onychodictyon_ferox, Onychodictyon_gracilis, ...
#> 
#> Rooted; no branch length.
```
