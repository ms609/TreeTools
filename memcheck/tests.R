# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < tests/thisfile.R
# First build and install the package.
library('TreeTools')
devtools::test()
