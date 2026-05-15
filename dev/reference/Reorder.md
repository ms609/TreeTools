# Reorder tree edges and nodes

Functions for systematically ordering the internal edges of trees.

## Usage

``` r
RenumberTree(parent, child, weight)

RenumberEdges(parent, child)

Cladewise(tree, nTip, edge)

# S3 method for class 'phylo'
Cladewise(tree, nTip = NTip(tree), edge = tree[["edge"]])

# S3 method for class 'list'
Cladewise(tree, nTip, edge)

# S3 method for class 'multiPhylo'
Cladewise(tree, nTip, edge)

# S3 method for class 'matrix'
Cladewise(tree, nTip = min(tree[, 1]) - 1L, edge)

# S3 method for class '`NULL`'
Cladewise(tree, nTip = min(tree[, 1]) - 1L, edge)

ApePostorder(tree, nTip, edge)

# S3 method for class 'phylo'
ApePostorder(tree, nTip = NTip(tree), edge = tree[["edge"]])

# S3 method for class 'list'
ApePostorder(tree, nTip, edge)

# S3 method for class '`NULL`'
ApePostorder(tree, nTip, edge)

# S3 method for class 'multiPhylo'
ApePostorder(tree, nTip, edge)

Postorder(tree, force = FALSE)

# S3 method for class 'phylo'
Postorder(tree, force = FALSE)

# S3 method for class '`NULL`'
Postorder(tree, force = FALSE)

# S3 method for class 'list'
Postorder(tree, force = FALSE)

# S3 method for class 'multiPhylo'
Postorder(tree, force = FALSE)

# S3 method for class 'numeric'
Postorder(tree, force = FALSE)

PostorderOrder(tree)

# S3 method for class 'phylo'
PostorderOrder(tree)

# S3 method for class 'numeric'
PostorderOrder(tree)

Pruningwise(tree, nTip, edge)

# S3 method for class 'phylo'
Pruningwise(tree, nTip = NTip(tree), edge = tree[["edge"]])

# S3 method for class 'list'
Pruningwise(tree, nTip, edge)

# S3 method for class 'multiPhylo'
Pruningwise(tree, nTip, edge)

# S3 method for class '`NULL`'
Pruningwise(tree, nTip, edge)

Preorder(tree, topologyOnly = FALSE)

# S3 method for class 'phylo'
Preorder(tree, topologyOnly = FALSE)

# S3 method for class 'numeric'
Preorder(tree, topologyOnly = FALSE)

# S3 method for class 'multiPhylo'
Preorder(tree, topologyOnly = FALSE)

# S3 method for class 'list'
Preorder(tree, topologyOnly = FALSE)

# S3 method for class '`NULL`'
Preorder(tree, topologyOnly = FALSE)

TntOrder(tree)

TNTOrder(tree)

# S3 method for class 'phylo'
TntOrder(tree)

# S3 method for class 'numeric'
TntOrder(tree)

# S3 method for class 'multiPhylo'
TntOrder(tree)

# S3 method for class 'list'
TntOrder(tree)

# S3 method for class '`NULL`'
TntOrder(tree)
```

## Arguments

- parent:

  Integer vector corresponding to the first column of the edge matrix of
  a tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html),
  i.e. `tree[["edge"]][, 1]`

- child:

  Integer vector corresponding to the second column of the edge matrix
  of a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), i.e.
  `tree[["edge"]][, 2]`.

- weight:

  Optional vector specifying the weight of each edge, corresponding to
  the `edge.length` property of a `phylo` object.

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- nTip:

  Integer specifying number of tips (leaves).

- edge:

  Two-column matrix listing the parent and child of each edge in a tree,
  corresponding to `tree[["edge"]]`. Optional in `Cladewise()`.

- force:

  Logical specifying whether to rearrange trees already in postorder, in
  order to ensure edges are ordered in the "TreeTools" fashion.

- topologyOnly:

  Logical; if `TRUE`, edge weights may not be retained.

## Value

`RenumberTree()` returns an edge matrix for a tree of class `phylo`
following the preorder convention for edge and node numbering.

`RenumberEdges()` formats the output of `RenumberTree()` into a list
whose two entries correspond to the new parent and child vectors, in
preorder.

`ApePostorder()`, `Cladewise()`, `Postorder()`, `Preorder()` and
`Pruningwise()` each return a tree of class `phylo` with nodes following
the specified numbering scheme.

`Postorder.numeric` accepts a numeric matrix corresponding to the `edge`
entry of a tree of class `phylo`, and returns a two-column array
corresponding to `tree`, with edges listed in postorder

`PostorderOrder()` returns an integer vector. Visiting edges in this
order will traverse the tree in postorder.

## Details

`Reorder()` is a wrapper for `ape:::.reorder_ape`. Calling this C
function directly is approximately twice as fast as using
`ape::`[`cladewise`](https://rdrr.io/pkg/ape/man/reorder.phylo.html) or
`ape::`[`postorder`](https://rdrr.io/pkg/ape/man/reorder.phylo.html)

`Cladewise()`, `ApePostorder()` and `Pruningwise()` are convenience
functions to the corresponding functions in "ape". Single nodes may need
to be collapsed using
[ape::collapse.singles](https://rdrr.io/pkg/ape/man/collapse.singles.html)
first. "ape" functions can cause crashes if nodes are numbered
unconventionally – sometimes arising after using tree rearrangement
functions, e.g. `phangorn::SPR()`.

`Preorder()` is more robust: it supports polytomies, nodes may be
numbered in any sequence, and edges may be listed in any order in the
input tree. Its output is guaranteed to be identical for any tree of an
equivalent leaf labelling (see
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md))
and topology, allowing unique trees to be detected by comparing sorted
edge matrices alone.

Nodes and edges in a preorder tree are numbered starting from the
deepest node. Each node is numbered in the sequence in which it is
encountered, and each edge is listed in the sequence in which it is
visited.

At each node, child edges are sorted from left to right in order of the
lowest-numbered leaf in the subtree subtended by each edge; i.e. an edge
that leads eventually to tip 1 will be to the left of an edge leading to
a subtree containing tip 2.

Numbering begins by following the leftmost edge of the root node, and
sorting its descendant subtree into preorder. Then, the next edge at the
root node is followed, and its descendants sorted into preorder, until
each edge has been visited.

`RenumberTree()` and `RenumberEdges()` are wrappers for the C function
`preorder_edges_and_nodes()`; they do not perform the same checks on
input as `Preorder()` and are intended for use where performance is at a
premium.

`Postorder()` numbers nodes as in `Preorder()`, and lists edges in
descending order of parent node number, breaking ties by listing child
nodes in increasing order. If a tree is already in postorder, it will
not be rearranged unless `force = TRUE`.

Methods applied to numeric inputs do not check input for sanity, so
should be used with caution: malformed input may cause undefined
results, including crashing R.

Trees with \>8191 leaves require additional memory and are not handled
by `Postorder()` at present. If you need to process such large trees,
please contact the maintainer for advice.

## Functions

- `Cladewise()`: Reorder tree cladewise.

- `ApePostorder()`: Reorder tree in Postorder using ape's `postorder`
  function, which is robust to unconventional node numbering.

- `Pruningwise()`: Reorder tree Pruningwise.

- `Preorder()`: Reorder tree in Preorder (special case of cladewise).

- `TntOrder()`: Reorder tree in postorder, numbering internal nodes
  according to [TNT's
  rules](https://stackoverflow.com/a/54296100/3438001), which number the
  root node as `nTip + 1`, then the remaining nodes in the sequence
  encountered when traversing the tree in postorder, starting from each
  tip in sequence.

## See also

Rotate each node into a consistent orientation with
[`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md).

Other tree manipulation:
[`AddTip()`](https://ms609.github.io/TreeTools/dev/reference/AddTip.md),
[`CollapseNode()`](https://ms609.github.io/TreeTools/dev/reference/CollapseNode.md),
[`ConsensusWithout()`](https://ms609.github.io/TreeTools/dev/reference/ConsensusWithout.md),
[`DropTip()`](https://ms609.github.io/TreeTools/dev/reference/DropTip.md),
[`ImposeConstraint()`](https://ms609.github.io/TreeTools/dev/reference/ImposeConstraint.md),
[`KeptPaths()`](https://ms609.github.io/TreeTools/dev/reference/KeptPaths.md),
[`KeptVerts()`](https://ms609.github.io/TreeTools/dev/reference/KeptVerts.md),
[`LeafLabelInterchange()`](https://ms609.github.io/TreeTools/dev/reference/LeafLabelInterchange.md),
[`MakeTreeBinary()`](https://ms609.github.io/TreeTools/dev/reference/MakeTreeBinary.md),
[`Renumber()`](https://ms609.github.io/TreeTools/dev/reference/Renumber.md),
[`RenumberTips()`](https://ms609.github.io/TreeTools/dev/reference/RenumberTips.md),
[`RootTree()`](https://ms609.github.io/TreeTools/dev/reference/RootTree.md),
[`SortTree()`](https://ms609.github.io/TreeTools/dev/reference/SortTree.md),
[`Subtree()`](https://ms609.github.io/TreeTools/dev/reference/Subtree.md),
[`TipTimedTree()`](https://ms609.github.io/TreeTools/dev/reference/TipTimedTree.md),
[`TrivialTree`](https://ms609.github.io/TreeTools/dev/reference/TrivialTree.md)

Other C wrappers:
[`Neworder`](https://ms609.github.io/TreeTools/dev/reference/Neworder.md)

Other C wrappers:
[`Neworder`](https://ms609.github.io/TreeTools/dev/reference/Neworder.md)

## Author

`Preorder()` and `Postorder()`: Martin R. Smith.

`Cladewise()`, `ApePostorder()` and `Pruningwise()`: modified by Martin
R. Smith from `.reorder_ape()` in ape (Emmanuel Paradis).
