## Test environments

* Local Windows 10 installation, R 4.0.0
* Ubuntu 16.04.6 LTS, R 3.4.0, release and devel, via 
  [Travis CI](https://travis-ci.org/ms609/TreeTools)
* Mac OS X 10.13.6, R release, via Travis
* R-hub, with `check_for_cran()` and `check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

> Maintainer: 'Martin R. Smith <martin.smith@durham.ac.uk>'
> 
> Possibly mis-spelled words in DESCRIPTION:
>   cladistic (26:30)

Spelling confirmed as correct

> Found the following (possibly) invalid URLs:
>  URL: https://doi.org/10.1137/0403005
>    From: man/NRooted.Rd
>          man/UnrootedTreesMatchingSplit.Rd
>    Status: Error
>    Message: libcurl error code 56:
>      	Recv failure: Connection was reset

The DOI and URL are correct, and work in my browser.

## Downstream dependencies

Downstream dependencies R CMD CHECKed locally: No changes to worse in

√ Quartet 1.1.0                          -- E: 0     | W: 0     | N: 0
√ TreeSearch 0.4.0                       -- E: 0     | W: 0     | N: 0

