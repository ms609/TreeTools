# Convert object to `Splits`

`as.Splits()` converts a phylogenetic tree to a `Splits` object
representing its constituent bipartition splits.

## Usage

``` r
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'phylo'
as.Splits(x, tipLabels = NULL, asSplits = TRUE, ...)

# S3 method for class 'multiPhylo'
as.Splits(x, tipLabels = unique(unlist(TipLabels(x))), asSplits = TRUE, ...)

# S3 method for class 'Splits'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'list'
as.Splits(x, tipLabels = NULL, asSplits = TRUE, ...)

# S3 method for class 'matrix'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'integer'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'numeric'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'logical'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'character'
as.Splits(x, tipLabels = NULL, ...)

# S3 method for class 'Splits'
as.logical(x, tipLabels = attr(x, "tip.label"), ...)
```

## Arguments

- x:

  Object to convert into splits: perhaps a tree of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html). If a logical
  matrix is provided, each row will be considered as a separate split.
  If an integer matrix is provided, entries assigned the same integer
  will be assigned to the same split.

- tipLabels:

  Character vector specifying sequence in which to order tip labels.
  Label order must (currently) match to combine or compare separate
  `Splits` objects.

- ...:

  Presently unused.

- asSplits:

  Logical specifying whether to return a `Splits` object, or an
  unannotated two-dimensional array (useful where performance is
  paramount).

## Value

`as.Splits()` returns an object of class `Splits`, or (if
`asSplits = FALSE`) a two-dimensional array of `raw` objects, with each
bit specifying whether or not the leaf corresponding to the respective
bit position is a member of the split. Splits are named according to the
node at the non-root end of the edge that defines them. In rooted trees,
the child of the rightmost root edge names the split.

## See also

Other Splits operations:
[`LabelSplits()`](https://ms609.github.io/TreeTools/dev/reference/LabelSplits.md),
[`NSplits()`](https://ms609.github.io/TreeTools/dev/reference/NSplits.md),
[`NTip()`](https://ms609.github.io/TreeTools/dev/reference/NTip.md),
[`PolarizeSplits()`](https://ms609.github.io/TreeTools/dev/reference/PolarizeSplits.md),
[`SplitFrequency()`](https://ms609.github.io/TreeTools/dev/reference/SplitFrequency.md),
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
splits <- as.Splits(BalancedTree(letters[1:6]))
summary(splits)
#> 3 bipartition splits dividing 6 tips, a .. f
#>      123456
#>  8   ***...
#>  9   **....
#>  11  ...**.
#> 
#>  Tip 1: a     Tip 2: b    Tip 3: c    Tip 4: d    Tip 5: e   
#>  Tip 6: f    
TipsInSplits(splits)
#>  8  9 11 
#>  3  2  2 
summary(!splits)
#> 3 bipartition splits dividing 6 tips, a .. f
#>      123456
#>  8   ...***
#>  9   ..****
#>  11  ***..*
#> 
#>  Tip 1: a     Tip 2: b    Tip 3: c    Tip 4: d    Tip 5: e   
#>  Tip 6: f    
TipsInSplits(!splits)
#>  8  9 11 
#>  3  4  4 

length(splits + !splits)
#> [1] 6
length(unique(splits + !splits))
#> [1] 3

summary(c(splits[[2:3]], !splits[[1:2]]))
#> 4 bipartition splits dividing 6 tips, a .. f
#>      123456
#>  9   **....
#>  11  ...**.
#>  8   ...***
#>  9   ..****
#> 
#>  Tip 1: a     Tip 2: b    Tip 3: c    Tip 4: d    Tip 5: e   
#>  Tip 6: f    

moreSplits <- as.Splits(PectinateTree(letters[6:1]), tipLabel = splits)
print(moreSplits, details = TRUE)
#> 3 bipartition splits dividing 6 tips, a .. f
#>      123456
#>  9   ****..
#>  10  ***...
#>  11  **....
match(splits, moreSplits)
#> [1]  2  3 NA
moreSplits %in% splits
#>     9    10    11 
#> FALSE  TRUE  TRUE 

as.Splits("....**", letters[1:6])
#> 1 bipartition split dividing 6 tips, a .. f
```
