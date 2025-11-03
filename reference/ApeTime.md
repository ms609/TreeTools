# Read modification time from "ape" Nexus file

`ApeTime()` reads the time that a tree written with "ape" was modified,
based on the comment in the Nexus file.

## Usage

``` r
ApeTime(filepath, format = "double")
```

## Arguments

- filepath:

  Character string specifying path to the file.

- format:

  Format in which to return the time: "double" as a sortable numeric;
  any other value to return a string in the format
  `YYYY-MM-DD hh:mm:ss`.

## Value

`ApeTime()` returns the time that the specified file was created by ape,
in the format specified by `format`.

## Author

[Martin R. Smith](https://orcid.org/0000-0001-5660-1727)
(<martin.smith@durham.ac.uk>)
