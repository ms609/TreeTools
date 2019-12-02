## Test environments
* Windows 10 on local machine, R 3.6.1
* Windows 10 via check_win_devel(), R devel
* ubuntu 16.04.6 LTS (on travis-ci), R 3.4.0, release and devel
* Mac OS X 10.13.3 (on travis-ci), R release
* Using check_rhub()

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

This is a new submission, maintained by Martin R. Smith <martin.smith@durham.ac.uk>.

## Response to CRAN comments
### Original submission (Martina Schmirl, comments dated 2019-11-25)
If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file [...]

> This is a workhorse package and does not pertain to a published study;
> no references are relevant.

Please always make sure to reset to user's options, wd or par after you
changed it in examples and vignettes.

> Fixed.

\dontrun{} should be only used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Examples of ReadTntTree.Rd are wrapped in \dontrun{}.
Does not seem necessary.
Please add small files needed for the examples in the inst/extdata
subfolder of your package and use system.file() to get the correct
package path.
Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.

> I have included files in inst/extdata to reduce the amount of code within
> \dontrun{}; now only examples that rely on a local system path that will
> not be present on a user's machine are \dontrun{}ed.

### Resubmission (comments dated 2019-11-28)
>> Note that the included small files are in the format expected for TNT output
>> files.  As such, they include hard-coded example paths that will not exist 
>> on the user's machine.  But the examples are intended to show a user how they
>> can easily work with files they have generated on their own machine using
>> their own copy of TNT.  Thus it is necessary to retain the `\dontrun{}` 
>> section, unfortunately.

Please improve documentation of C_node_depth.Rd and similar ones, or do
not export those functions.

>> I have now either documented or removed all internal functions.


## Downstream dependencies
There are currently no downstream dependencies for this package.
