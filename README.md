# TreeTrunk

[![Build Status](https://travis-ci.org/ms609/TreeTrunk.svg?branch=master)](https://travis-ci.org/ms609/TreeTrunk)
[![codecov](https://codecov.io/gh/ms609/TreeTrunk/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeTrunk)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeTrunk)](https://cran.r-project.org/package=TreeTrunk)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeTrunk)](https://cran.r-project.org/package=TreeTrunk)
[![DOI](https://zenodo.org/badge/98171642.svg)](https://zenodo.org/badge/latestdoi/98171642)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

This package provides functions for creating, modifying and analysing 
phylogenetic trees.  It complements packages such as 
ape,
phangorn and
phytools,
aiming for efficient and robust implementations of functions, typically
applied to unweighted trees (i.e. those without edge lengths).

# Installation

Install and load the library from CRAN as follows:
```
install.packages('TreeTrunk')
library('TreeTrunk')
```

If you're feeling brave, you can install the development version thus:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeTrunk')
```

Please note that the 'TreeTrunk' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
