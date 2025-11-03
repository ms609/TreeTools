# Cladistic information content of a tree

`CladisticInfo()` calculates the cladistic (phylogenetic) information
content of a phylogenetic object, *sensu* Thorley *et al.* (1998).

## Usage

``` r
CladisticInfo(x)

# S3 method for class 'phylo'
CladisticInfo(x)

# S3 method for class 'Splits'
CladisticInfo(x)

# S3 method for class 'list'
CladisticInfo(x)

# S3 method for class 'multiPhylo'
CladisticInfo(x)

CladisticInformation(x)
```

## Arguments

- x:

  Tree of class `phylo`, or a list thereof.

## Value

`CladisticInfo()` returns a numeric giving the cladistic information
content of the input tree(s), in bits. If passed a `Splits` object, it
returns the information content of each split in turn.

## Details

The CIC is the logarithm of the number of binary trees that include the
specified topology. A base two logarithm gives an information content in
bits.

The CIC was originally proposed by Rohlf (1982) , and formalised, with
an information-theoretic justification, by Thorley et al. (1998) . Steel
and Penny (2006) term the equivalent quantity "phylogenetic information
content" in the context of individual characters.

The number of binary trees consistent with a cladogram provides a more
satisfactory measure of the resolution of a tree than simply counting
the number of edges resolved (Page 1992) .

## References

Page RD (1992). “Comments on the information content of
classifications.” *Cladistics*, **8**(1), 87–95.
[doi:10.1111/j.1096-0031.1992.tb00054.x](https://doi.org/10.1111/j.1096-0031.1992.tb00054.x)
.  
  
Rohlf FJ (1982). “Consensus indices for comparing classifications.”
*Mathematical Biosciences*, **59**(1), 131–144.
[doi:10.1016/0025-5564(82)90112-2](https://doi.org/10.1016/0025-5564%2882%2990112-2)
.  
  
Steel MA, Penny D (2006). “Maximum parsimony and the phylogenetic
information in multistate characters.” In Albert VA (ed.), *Parsimony,
Phylogeny, and Genomics*, 163–178. Oxford University Press, Oxford.  
  
Thorley JL, Wilkinson M, Charleston M (1998). “The information content
of consensus trees.” In Rizzi A, Vichi M, Bock H (eds.), *Advances in
Data Science and Classification*, 91–98. Springer, Berlin. ISBN
978-3-540-64641-9,
[doi:10.1007/978-3-642-72253-0](https://doi.org/10.1007/978-3-642-72253-0)
.

## See also

Other tree information functions:
[`NRooted()`](https://ms609.github.io/TreeTools/reference/NRooted.md),
[`TreesMatchingTree()`](https://ms609.github.io/TreeTools/reference/TreesMatchingTree.md)

Other tree characterization functions:
[`Consensus()`](https://ms609.github.io/TreeTools/reference/Consensus.md),
[`J1Index()`](https://ms609.github.io/TreeTools/reference/J1Index.md),
[`Stemwardness`](https://ms609.github.io/TreeTools/reference/Stemwardness.md),
[`TotalCopheneticIndex()`](https://ms609.github.io/TreeTools/reference/TotalCopheneticIndex.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
