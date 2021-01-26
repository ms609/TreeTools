## Test environments

* Local PC:
  - Windows 10, R 4.0.3

* [GitHub Actions](https://github.com/ms609/TreeTools/actions)
  - Ubuntu 20.04
    - R 3.4.4
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.
There was one NOTE:

>  Possibly mis-spelled words in DESCRIPTION:
>    al (22:62)
>    cladistic (25:30)
>    et (22:59)
>    Maddison (22:50)
>    Newick (22:35)

These spellings have been verified.

v Quartet 1.1.0                          -- E: 0     | W: 0     | N: 1          
v TBRDist 1.0.2                          -- E: 0     | W: 0     | N: 1          
v TreeDist 1.2.1                         -- E: 0     | W: 0     | N: 2          
v TreeSearch 0.4.3                       -- E: 0     | W: 0     | N: 0      

