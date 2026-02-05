# Decompose additive (ordered) phylogenetic characters

`Decompose()` decomposes additive characters into a series of binary
characters, which is mathematically equivalent when analysed under equal
weights parsimony. (This equivalence is not exact under implied weights
or under probabilistic tree inference methods.)

## Usage

``` r
Decompose(dataset, indices)
```

## Arguments

- dataset:

  A phylogenetic data matrix of phangorn class `phyDat`, whose names
  correspond to the labels of any accompanying tree.

- indices:

  Integer or logical vector specifying indices of characters that should
  be decomposed

## Value

`Decompose()` returns a `phyDat` object in which the specified ordered
characters have been decomposed into binary characters. The attribute
`originalIndex` lists the index of the character in `dataset` to which
each element corresponds.

## Details

An ordered (additive) character can be rewritten as a mathematically
equivalent hierarchy of binary neomorphic characters (Farris et al.
1970) . Two reasons to prefer the latter approach are:

- It makes explicit the evolutionary assumptions underlying an ordered
  character, whether the underlying ordering is linear, reticulate or
  branched (Mabee 1989) .

- It avoids having to identify characters requiring special treatment to
  phylogenetic software, which requires the maintenance of an up-to-date
  log of which characters are treated as additive and which sequence
  their states occur in, a step that may be overlooked by re-users of
  the data.

Careful consideration is warranted when evaluating whether a group of
related characteristics ought to be treated as ordered (Wilkinson 1992)
. On the one hand, the 'principle of indifference' states that we should
treat all transformations as equally probable (/ surprising /
informative); ordered characters fail this test, as larger changes are
treated as less probable than smaller ones. On the other hand, ordered
characters allow more opportunities for homology of different character
states, and might thus be defended under the auspices of Hennig’s
Auxiliary Principle (Wilkinson 1992) .

For a case study of how ordering phylogenetic characters can affect
phylogenetic outcomes in practice, see Brady et al. (2024) .

## References

Brady PL, Castrellon Arteaga A, López-Torres S, Springer MS (2024). “The
Effects of Ordered Multistate Morphological Characters on Phylogenetic
Analyses of Eutherian Mammals.” *Journal of Mammalian Evolution*,
**31**(3), 28.
[doi:10.1007/s10914-024-09727-2](https://doi.org/10.1007/s10914-024-09727-2)
.  
  
Farris JS, Kluge AG, Eckardt MJ (1970). “A Numerical Approach to
Phylogenetic Systematics.” *Systematic Biology*, **19**(2), 172–189.
[doi:10.2307/2412452](https://doi.org/10.2307/2412452) .  
  
Mabee PM (1989). “Assumptions Underlying the Use of Ontogenetic
Sequences for Determining Character State Order.” *Transactions of the
American Fisheries Society*, **118**(2), 151–158.
[doi:10.1577/1548-8659(1989)118\<0151:AUTUOO\>2.3.CO;2](https://doi.org/10.1577/1548-8659%281989%29118%3C0151%3AAUTUOO%3E2.3.CO%3B2)
.  
  
Wilkinson M (1992). “Ordered versus Unordered Characters.” *Cladistics*,
**8**(4), 375–385.
[doi:10.1111/j.1096-0031.1992.tb00079.x](https://doi.org/10.1111/j.1096-0031.1992.tb00079.x)
.

## See also

Other phylogenetic matrix conversion functions:
[`MatrixToPhyDat()`](https://ms609.github.io/TreeTools/reference/MatrixToPhyDat.md),
[`Reweight()`](https://ms609.github.io/TreeTools/reference/Reweight.md),
[`StringToPhyDat()`](https://ms609.github.io/TreeTools/reference/PhyToString.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo")

# Identify character 11 as additive
# Character 11 will be replaced with two characters
# The present codings 0, 1 and 2 will be replaced with 00, 10, and 11.
decomposed <- Decompose(Lobo.phy, 11)

NumberOfChars <- function(x) sum(attr(x, "weight"))
NumberOfChars(Lobo.phy)   # 115 characters in original
#> [1] 115
NumberOfChars(decomposed) # 116 characters in decomposed
#> [1] 116
```
