## Test environments

* Local Windows 10 installation, R 3.6.3
* Windows 10 via check_win_devel(quiet = TRUE), R devel
* Ubuntu 16.04.6 LTS, R 3.4.0, release and devel, via 
  [Travis CI](https://travis-ci.org/ms609/TreeTools)
* Mac OS X 10.13.6, R release, via Travis
* R-hub, with `check_for_cran()`

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

Apologies for overlooking this in the initial submission.
 
> Found the following (possibly) invalid URLs:
>   URL: https://doi.org/10.1093/sysbio/46.4.590
>     From: man/ReadCharacters.Rd
>     Status: Error
>     Message: libcurl error code 56:
>       	Recv failure: Operation timed out
> 
> Found the following (possibly) invalid DOIs:
>   DOI: 10.1093/sysbio/46.4.590
>     From: DESCRIPTION
>     Status: libcurl error code 56:
>     	Recv failure: Operation timed out
>     Message: Error

The DOI and URL are correct (though the page seems a little slow to load and
is presumably causing a timeout).

## Downstream dependencies

revdepcheck::revdep_check() found no changes to worse in downstream dependencies

√ Quartet 1.1.0                          -- E: 0     | W: 0     | N: 0
√ TreeSearch 0.4.0                       -- E: 0     | W: 0     | N: 0

