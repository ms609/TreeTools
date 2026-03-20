# Frequency of splits

`SplitFrequency()` provides a simple way to count the number of times
that bipartition splits, as defined by a reference tree, occur in a
forest of trees. May be used to calculate edge ("node") support for
majority consensus or bootstrap trees.

## Usage

``` r
SplitFrequency(reference, forest = NULL)
```

## Arguments

- reference:

  A tree of class `phylo`, a `Splits` object. If `NULL`, the frequencies
  of all splits in `forest` will be returned.

- forest:

  A list of trees of class `phylo`, or a `multiPhylo` object; or a
  `Splits` object. See
  [vignette](https://ms609.github.io/TreeTools/articles/load-trees.html)
  for possible methods of loading trees into R.

## Value

`SplitFrequency()` returns the number of trees in `forest` that contain
each split in `reference`. If `reference` is a tree of class `phylo`,
then the sequence will correspond to the order of nodes (use
[`ape::nodelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html) to
view). Note that the three nodes at the root of the tree correspond to a
single split; see the example for how these might be plotted on a tree.

If `reference` is `NULL`, then `SplitFrequency()` returns a list of
splits (in the order encountered in forest) with attribute `"count"`
stating the number of times each split occurs in `forest`.

## Details

If multiple calculations are required, some time can be saved by using
the constituent functions (see examples).

## See also

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# An example forest of 100 trees, some identical
forest <- as.phylo(c(1, rep(10, 79), rep(100, 15), rep(1000, 5)), nTip = 9)

# Generate an 80% consensus tree
cons <- ape::consensus(forest, p = 0.8)
plot(cons)


# Calculate split frequencies
splitFreqs <- SplitFrequency(cons, forest)

# Optionally, colour edges by corresponding frequency.
# Note that not all edges are associated with a unique split
# (and two root edges may be associated with one split - not handled here)
edgeSupport <- rep(1, nrow(cons$edge)) # Initialize trivial splits to 1
childNode <- cons$edge[, 2]
edgeSupport[match(names(splitFreqs), childNode)] <- splitFreqs / 100

plot(cons, edge.col = SupportColour(edgeSupport), edge.width = 3)

# Annotate nodes by frequency 
LabelSplits(cons, splitFreqs, unit = "%",
            col = SupportColor(splitFreqs / 100),
            frame = "none", pos = 3L)

```
