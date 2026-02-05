# Reweight phylogenetic characters

`Reweight()` allows the weights of specific characters in phylogenetic
datasets to be arbitrarily adjusted.

## Usage

``` r
Reweight(dataset, weights)
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, or as a matrix
  in the format produced by
  [`PhyDatToMatrix()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md).

- weights:

  Unnamed integer vector specifying desired weight of each character in
  turn; or named integer vector specifying weights of each character;
  unnamed entries will be assigned weight 1.

## Value

`Reweight()` returns `dataset` after adjusting the weights of the
specified characters. For a matrix, this is attained by repeating each
column the `weights` times. For a `phyDat` object, the "weight"
attribute will be modified.

## Details

This functionality should be employed with care. The underlying
principle of parsimony is that all evolutionary steps are equivalent.
Setting different weights to different characters is at odds with that
principle, so analysis of a re-weighted matrix using a parsimony-based
framework is arguably no longer parsimony analysis; on the most
permissive view, the criteria used to determine a weighting scheme will
always be arbitrary.

It can be useful to relax the criterion that all evolutionary steps are
equivalent – for example, implied weighting (Goloboff 1997) typically
recovers better trees than equal-weights parsimony (Smith 2019) . This
said, assigning different weights to different characters tacitly
imposes a model of evolution that differs from that implicit in
equal-weights parsimony. Whereas probabilistic models can be evaluated
by various methods (e.g. fit, marginal likelihood, posterior predictive
power), there are no principled methods of comparing different models
under a parsimony framework.

As such, `Reweight()` is likely to be useful for a narrow set of uses.
Examples may include:

- informal robustness testing, to explore whether certain characters are
  more or less influential on the resulting tree;

- Imposing constraints on a dataset, by adding each constraint as a
  column in a dataset whose weight exceeds the total amount of data.

## References

Goloboff PA (1997). “Self-Weighted Optimization: Tree Searches and
Character State Reconstructions under Implied Transformation Costs.”
*Cladistics*, **13**(3), 225–245.
[doi:10.1111/j.1096-0031.1997.tb00317.x](https://doi.org/10.1111/j.1096-0031.1997.tb00317.x)
.  
  
Smith MR (2019). “Bayesian and Parsimony Approaches Reconstruct
Informative Trees from Simulated Morphological Datasets.” *Biology
Letters*, **15**(2), 20180632.
[doi:10.1098/rsbl.2018.0632](https://doi.org/10.1098/rsbl.2018.0632) .

## See also

Other phylogenetic matrix conversion functions:
[`Decompose()`](https://ms609.github.io/TreeTools/reference/Decompose.md),
[`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md),
[`StringToPhyDat()`](https://ms609.github.io/TreeTools/reference/PhyToString.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
mat <- rbind(a = c(0, 2, 0), b = c(0, 2, 0), c = c(1, 3, 0), d = c(1, 3, 0))
dat <- MatrixToPhyDat(mat)

 # Set character 1 to weight 1, character 2 to weight 2; omit character 3
Reweight(mat, c(1, 2, 0))
#>   _1 _1 _2
#> a  0  2  2
#> b  0  2  2
#> c  1  3  3
#> d  1  3  3
# Equivalently:
Reweight(dat, c("3" = 0, "2" = 2))
#> $a
#> [1] 1 3 1
#> 
#> $b
#> [1] 1 3 1
#> 
#> $c
#> [1] 2 4 1
#> 
#> $d
#> [1] 2 4 1
#> 
#> attr(,"weight")
#> [1] 1 2 0
#> attr(,"nr")
#> [1] 3
#> attr(,"nc")
#> [1] 4
#> attr(,"index")
#> [1] 1 2 3
#> attr(,"levels")
#> [1] "0" "1" "2" "3"
#> attr(,"allLevels")
#> [1] "0" "1" "2" "3"
#> attr(,"type")
#> [1] "USER"
#> attr(,"contrast")
#>      0 1 2 3
#> [1,] 1 0 0 0
#> [2,] 0 1 0 0
#> [3,] 0 0 1 0
#> [4,] 0 0 0 1
#> attr(,"class")
#> [1] "phyDat"
```
