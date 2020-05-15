## Test environments

* Local Windows 10 installation, R 4.0.0
* Ubuntu 16.04.6 LTS, R 3.4.0, release and devel, via 
  [Travis CI](https://travis-ci.org/ms609/TreeTools)
* Mac OS X 10.13.6, R release, via Travis
* R-hub, with `check_for_cran()`

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

> Maintainer: 'Martin R. Smith <martin.smith@durham.ac.uk>'
> 
> Possibly mis-spelled words in DESCRIPTION:
>   cladistic (26:30)

Spelling confirmed as correct

> 
> Found the following (possibly) invalid URLs:
>   URL: http://epubs.siam.org/doi/abs/10.1137/0403005
>     From: man/NRooted.Rd
>           man/UnrootedTreesMatchingSplit.Rd
>     Message: libcurl error code 56:
>       	Recv failure: Connection was reset
>     Status: Error
>   URL: https://doi.org/10.1093/sysbio/46.4.590
>     From: man/ReadCharacters.Rd
>     Status: Error
>     Message: libcurl error code 56:
>       	Recv failure: Connection was reset
>   URL: https://doi.org/10.1137/0403005
>     From: man/NRooted.Rd
>           man/UnrootedTreesMatchingSplit.Rd
>     Status: Error
>     Message: libcurl error code 56:
>       	Recv failure: Connection was reset
> 
> Found the following (possibly) invalid DOIs:
>   DOI: 10.1093/sysbio/46.4.590
>     From: DESCRIPTION
>     Status: libcurl error code 56:
>     	Recv failure: Connection was reset
>     Message: Error

The DOI and URLs are correct.

## Downstream dependencies

revdepcheck::revdep_check() found no changes to worse in downstream dependencies

√ Quartet 1.1.0                          -- E: 0     | W: 0     | N: 0
√ TreeSearch 0.4.0                       -- E: 0     | W: 0     | N: 0

