# Brewer palettes

A list of eleven Brewer palettes containing one to eleven colours that
are readily distinguished by colourblind viewers, followed by a twelfth
12-colour palette adapted for colour blindness.

## Usage

``` r
brewer
```

## Format

An object of class `list` of length 12.

## Source

- [ColourBrewer2.org](https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3)

- [Martin Krzywinski](https://mk.bcgsc.ca/colorblind/)

## Examples

``` r
data("brewer", package = "TreeTools")
plot(0, type = "n", xlim = c(1, 12), ylim = c(12, 1),
     xlab = "Colour", ylab="Palette")
for (i in seq_along(brewer)) text(seq_len(i), i, col = brewer[[i]])

```
