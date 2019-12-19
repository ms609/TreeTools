This re-submission catches the overloading of sqrt on Solaris identified by 
the CRAN package check results, and highlighted by Brian Ripley.
I have re-read the guidance for portable C/C++ code in 'Writing R Extensions'
and made concominant changes; I have also added Solaris to my test environments.

## Test environments
* Windows 10 on local machine, R 3.6.1
* Windows 10 via check_win_devel(quiet = TRUE), R devel
* ubuntu 16.04.6 LTS (on travis-ci), R 3.4.0, release and devel
* Mac OS X 10.13.3 (on travis-ci), R devel
* Using check_rhub(platforms = rhub::platforms()[[1]])

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
