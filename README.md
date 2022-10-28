# TreeTools

[![codecov](https://codecov.io/gh/ms609/TreeTools/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeTools)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeTools)](https://cran.r-project.org/package=TreeTools)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeTools)](https://cran.r-project.org/package=TreeTools)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3522726.svg)](http://doi.org/10.5281/zenodo.3522725)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

'TreeTools' is an R package that provides efficient implementations of 
functions for the creation, modification and analysis of phylogenetic trees.

Applications include:
generation of trees with specified shapes;
analysis of tree shape;
rooting of trees and extraction of subtrees;
calculation and depiction of node support;
calculation of ancestor-descendant relationships;
import and export of trees from Newick, Nexus and 
[TNT](https://www.lillo.org.ar/phylogeny/tnt/) formats;
and analysis of partitions and cladistic information.

It complements packages such as 
['ape'](https://cran.r-project.org/package=ape),
['phangorn'](https://cran.r-project.org/package=phangorn) and
['phytools'](https://cran.r-project.org/package=phytools),
aiming for efficient and robust implementations of functions, typically
applied to unweighted trees (i.e. those without edge lengths).

# Installation

Install and load the library from CRAN as follows:
```r
install.packages("TreeTools")
library("TreeTools")
```

Install the very latest version, which may be under development, with:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeTools')
```

Please note that the 'TreeTools' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
