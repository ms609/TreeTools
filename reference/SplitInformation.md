# Phylogenetic information content of splitting leaves into two partitions

Calculate the phylogenetic information content (*sensu* Steel and Penny
2006 ) of a split, which reflects the probability that a uniformly
selected random tree will contain the split: a split that is consistent
with a smaller number of trees will have a higher information content.

## Usage

``` r
SplitInformation(A, B = A[1])

MultiSplitInformation(partitionSizes)
```

## Arguments

- A, B:

  Integer specifying the number of taxa in each partition.

- partitionSizes:

  Integer vector specifying the number of taxa in each partition of a
  multi-partition split.

## Value

`SplitInformation()` and `MultiSplitInformation()` return the
phylogenetic information content, in bits, of a split that subdivides
leaves into partitions of the specified sizes.

## Details

`SplitInformation()` addresses bipartition splits, which correspond to
edges in an unrooted phylogeny; `MultiSplitInformation()` supports
splits that subdivide taxa into multiple partitions, which may
correspond to multi-state characters in a phylogenetic matrix.

A simple way to characterise trees is to count the number of edges.
(Edges are almost, but not quite, equivalent to nodes.) Counting edges
(or nodes) provides a quick measure of a tree's resolution, and
underpins the Robinson-Foulds tree distance measure. Not all edges,
however, are created equal.

An edge splits the leaves of a tree into two subdivisions. The more
equal these subdivisions are in size, the more instructive this edge is.
Intuitively, the division of mammals from reptiles is a profound
revelation that underpins much of zoology; recognizing that two species
of bat are more closely related to each other than to any other mammal
or reptile is still instructive, but somewhat less fundamental.

Formally, the phylogenetic (Shannon) information content of a split *S*,
*h(S)*, corresponds to the probability that a uniformly selected random
tree will contain the split, *P(S)*: *h(S)* = -log *P(S)*. Base 2
logarithms are typically employed to yield an information content in
bits.

As an example, the split `AB|CDEF` occurs in 15 of the 105 six-leaf
trees; *h*(`AB|CDEF`) = -log *P*(`AB|CDEF`) = -log(15/105) ~ 2.81 bits.
The split `ABC|DEF` subdivides the leaves more evenly, and is thus more
instructive: it occurs in just nine of the 105 six-leaf trees, and
*h*(`ABC|DEF`) = -log(9/105) ~ 3.54 bits.

As the number of leaves increases, a single even split may contain more
information than multiple uneven splits – see the examples section
below.

Summing the information content of all splits within a tree, perhaps
using the '[TreeDist](https://ms609.github.io/TreeDist/)' function
[`SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html),
arguably gives a more instructive picture of its resolution than simply
counting the number of splits that are present – though with the caveat
that splits within a tree are not independent of one another, so some
information may be double counted. (This same charge applies to simply
counting nodes, too.)

Alternatives would be to count the number of quartets that are resolved,
perhaps using the '[Quartet](https://ms609.github.io/Quartet/)' function
[`QuartetStates()`](https://ms609.github.io/Quartet/reference/QuartetState.html),
or to use a different take on the information contained within a split,
the clustering information: see the 'TreeDist' function
[`ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)
for details.

## References

Steel MA, Penny D (2006). “Maximum parsimony and the phylogenetic
information in multistate characters.” In Albert VA (ed.), *Parsimony,
Phylogeny, and Genomics*, 163–178. Oxford University Press, Oxford.

## See also

Sum the phylogenetic information content of splits within a tree:
[[`TreeDist::SplitwiseInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)](https://ms609.github.io/TreeDist/reference/TreeInfo.html)

Sum the clustering information content of splits within a tree:
[[`TreeDist::ClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeInfo.html)](https://ms609.github.io/TreeDist/reference/TreeInfo.html)

Other split information functions:
[`CharacterInformation()`](https://ms609.github.io/TreeTools/reference/CharacterInformation.md),
[`SplitMatchProbability()`](https://ms609.github.io/TreeTools/reference/SplitMatchProbability.md),
[`TreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/TreesMatchingSplit.md),
[`UnrootedTreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/UnrootedTreesMatchingSplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# Eight leaves can be split evenly:
SplitInformation(4, 4)
#> [1] 5.529821

# or unevenly, which is less informative:
SplitInformation(2, 6)
#> [1] 3.459432

# A single split that evenly subdivides 50 leaves contains more information
# that seven maximally uneven splits on the same leaves:
SplitInformation(25, 25)
#> [1] 47.50376
7 * SplitInformation(2, 48)
#> [1] 45.98899
# Three ways to split eight leaves into multiple partitions:
MultiSplitInformation(c(2, 2, 4))
#> [1] 5.97728
MultiSplitInformation(c(2, 3, 3))
#> [1] 6.714246
MultiSplitInformation(rep(2, 4))
#> [1] 6.714246

```
