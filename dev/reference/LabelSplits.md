# Label splits

Labels the edges associated with each split on a plotted tree.

## Usage

``` r
LabelSplits(tree, labels = NULL, unit = "", ...)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- labels:

  Named vector listing annotations for each split. Names should
  correspond to the node associated with each split; see
  [`as.Splits()`](https://ms609.github.io/TreeTools/dev/reference/Splits.md)
  for details. If `NULL`, each splits will be labelled with its
  associated node.

- unit:

  Character specifying units of `labels`, if desired. Include a leading
  space if necessary.

- ...:

  Additional parameters to
  [`ape::edgelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html).

## Value

`LabelSplits()` returns
[`invisible()`](https://rdrr.io/r/base/invisible.html), after plotting
`labels` on each relevant edge of a plot (which should already have been
produced using `plot(tree)`).

## Details

As the two root edges of a rooted tree denote the same split, only the
rightmost (plotted at the bottom, by default) edge will be labelled. If
the position of the root is significant, add a tip at the root using
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md).

## See also

Calculate split support:
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md)

Colour labels according to value:
[`SupportColour()`](https://ms609.github.io/TreeTools/dev/reference/SupportColour.md)

Other Splits operations:
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
[`Splits`](https://ms609.github.io/TreeTools/dev/reference/Splits.md),
[`SplitsInBinaryTree()`](https://ms609.github.io/TreeTools/dev/reference/SplitsInBinaryTree.md),
[`TipLabels()`](https://ms609.github.io/TreeTools/dev/reference/TipLabels.md),
[`TipsInSplits()`](https://ms609.github.io/TreeTools/dev/reference/TipsInSplits.md),
[`match,Splits,Splits-method`](https://ms609.github.io/TreeTools/dev/reference/match.Splits.md),
[`xor()`](https://ms609.github.io/TreeTools/dev/reference/xor.md)

## Examples

``` r
tree <- BalancedTree(LETTERS[1:5])
splits <- as.Splits(tree)
plot(tree)
LabelSplits(tree, as.character(splits), frame = "none", pos = 3L)
LabelSplits(tree, TipsInSplits(splits), unit = " tips", frame = "none",
            pos = 1L)


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
