# Character information content

`CharacterInformation()` calculates the cladistic information content
(Steel and Penny 2006) of a given character, in bits. The total
information in all characters gives a measure of the potential utility
of a dataset (Cotton and Wilkinson 2008) , which can be compared with a
profile parsimony score (Faith and Trueman 2001) to evaluate the degree
of homoplasy within a dataset.

## Usage

``` r
CharacterInformation(tokens)
```

## Arguments

- tokens:

  Character vector specifying the tokens assigned to each taxon for a
  character. Example: `c(0, 0, 0, 1, 1, 1, "?", "-")`.

  Note that ambiguous tokens such as `(01)` are not supported, and
  should be replaced with `?`.

## Value

`CharacterInformation()` returns a numeric specifying the phylogenetic
information content of the character (*sensu* Steel and Penny 2006 ), in
bits.

## References

Cotton JA, Wilkinson M (2008). “Quantifying the potential utility of
phylogenetic characters.” *Taxon*, **57**(1), 131–136.  
  
Faith DP, Trueman JWH (2001). “Towards an inclusive philosophy for
phylogenetic inference.” *Systematic Biology*, **50**(3), 331–350.
[doi:10.1080/10635150118627](https://doi.org/10.1080/10635150118627) .  
  
Steel MA, Penny D (2006). “Maximum parsimony and the phylogenetic
information in multistate characters.” In Albert VA (ed.), *Parsimony,
Phylogeny, and Genomics*, 163–178. Oxford University Press, Oxford.

## See also

Other split information functions:
[`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.md),
[`SplitMatchProbability()`](https://ms609.github.io/TreeTools/reference/SplitMatchProbability.md),
[`TreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/TreesMatchingSplit.md),
[`UnrootedTreesMatchingSplit()`](https://ms609.github.io/TreeTools/reference/UnrootedTreesMatchingSplit.md)

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
