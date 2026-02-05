# File system navigation in R

Before you can open a file, you need to tell R where to find it. You can
do this by providing the full path to the file on your system. Be
careful to use forward slashes (`/`, not `\`, which you’ll get if you
copy file paths in Windows).

``` r
filename <- "C:/nexus/matrix.nex"
```

You can save typing by giving R a working directory. You can think of R
as having a file explorer window open invisibly in the background.

You can see the folder that’s open at the moment by typing
[`getwd()`](https://rdrr.io/r/base/getwd.html) at the console.

[`setwd()`](https://rdrr.io/r/base/getwd.html) tells R to open a
different folder instead.

`setwd('../')` tells R to go up to a parent directory. (You can [do this
using the Graphical User
Interface](https://www.princeton.edu/~otorres/RStudio101.pdf) in
RStudio).

By setting the directory that your files are in as the working
directory, you only need to specify the filename:

``` r
setwd("C:/nexus/") # You only need to do this once
filename <- "matrix.nex"
# Do something with this file
#

filename <- "tree.nex"
# Do something with this file
#
```

## What next?

Now you know how to locate files, you might want to load a
[dataset](https://ms609.github.io/TreeTools/articles/load-data.md) or
[phylogenetic
tree](https://ms609.github.io/TreeTools/articles/load-trees.md) into R.
