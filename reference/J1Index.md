# Robust universal tree balance index

Calculate tree balance index J¹ (when `nonRootDominance = FALSE`) or
J^(1c) (when `nonRootDominance = TRUE`) from Lemant J, Le Sueur C,
Manojlović V, Noble R (2022). “Robust, Universal Tree Balance Indices.”
*Systematic Biology*, **71**(5), 1210–1224.
[doi:10.1093/sysbio/syac027](https://doi.org/10.1093/sysbio/syac027) . .

## Usage

``` r
J1Index(tree, q = 1, nonRootDominance = FALSE)

JQIndex(tree, q = 1, nonRootDominance = FALSE)
```

## Arguments

- tree:

  Either an object of class 'phylo', or a dataframe with column names
  Parent, Identity and (optionally) Population. The latter is similar to
  `tree$edge`, where tree is an object of class 'phylo'; the differences
  are in class (data.frame versus matrix) and column names. The
  dataframe may (but isn't required to) include a row for the root node.
  If population sizes are omitted then internal nodes will be assigned
  population size zero and leaves will be assigned population size one.

- q:

  Numeric between zero and one specifying sensitivity to type
  frequencies. If `q < 1`, the J^(q) index - based on generalized
  entropy - will be returned; see Lemant et al. (2022) , page 1223.

- nonRootDominance:

  Logical specifying whether to use non-root dominance factor.

## Details

If population sizes are not provided, then the function assigns size 0
to internal nodes, and size 1 to leaves.

## References

Lemant J, Le Sueur C, Manojlović V, Noble R (2022). “Robust, Universal
Tree Balance Indices.” *Systematic Biology*, **71**(5), 1210–1224.
[doi:10.1093/sysbio/syac027](https://doi.org/10.1093/sysbio/syac027) .

## See also

Other tree characterization functions:
[`CladisticInfo()`](https://ms609.github.io/TreeTools/reference/CladisticInfo.md),
[`Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.md),
[`Stemwardness`](https://ms609.github.io/TreeTools/reference/Stemwardness.md),
[`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/reference/TotalCopheneticIndex.md)

## Author

Rob Noble, adapted by Martin R. Smith

## Examples

``` r
# Using phylo object as input:
phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
J1Index(phylo_tree)
#> [1] 0.7924813
phylo_tree2 <- read.tree(text='((A, B), ((C, D), (E, F)));')
J1Index(phylo_tree2)
#> [1] 0.9693609

# Using edges lists as input:
tree1 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
                    Identity = 1:7,
                    Population = c(1, rep(5, 6)))
J1Index(tree1)
#> [1] 0.6666667
tree2 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
                    Identity = 1:7,
                    Population = c(rep(0, 4), rep(1, 3)))
J1Index(tree2)
#> [1] 0.5
tree3 <- data.frame(Parent = c(1, 1, 1, 1, 2, 3, 4),
                    Identity = 1:7,
                    Population = c(0, rep(1, 3), rep(0, 3)))
J1Index(tree3)
#> [1] 1
cat_tree <- data.frame(Parent = c(1, 1:14, 1:15, 15),
                       Identity = 1:31,
                       Population = c(rep(0, 15), rep(1, 16)))
J1Index(cat_tree)
#> [1] 0.4740741

# If population sizes are omitted then internal nodes are assigned population
# size zero and leaves are assigned population size one:
sym_tree1 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
                       Identity = 1:31,
                       Population = c(rep(0, 15), rep(1, 16)))
# Equivalently:                        
sym_tree2 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
                       Identity = 1:31)
J1Index(sym_tree1)
#> [1] 1
J1Index(sym_tree2)
#> [1] 1
```
