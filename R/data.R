#' Brewer palettes
#'
#' A list of eleven Brewer palettes containing one to eleven colours that
#' are readily distinguished by colourblind viewers, followed by a twelfth
#' 12-colour palette adapted for colour blindness.
#'
#' @source{
#'
#' * [ColourBrewer2.org](https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3)
#'
#' * [Martin Krzywinski](https://mk.bcgsc.ca/colorblind/)
#' }
#'
#' @examples
#' data("brewer", package = "TreeTools")
#' plot(0, type = "n", xlim = c(1, 12), ylim = c(12, 1),
#'      xlab = "Colour", ylab="Palette")
#' for (i in seq_along(brewer)) text(seq_len(i), i, col = brewer[[i]])
#'
#' @keywords datasets
#' @name brewer
"brewer"

#' Double factorials
#'
#' A vector with pre-calculated values of double factorials up to 300!!,
#' and the logarithms of double factorials up to 50 000!!.
#'
#' 301!! is too large to store as an integer; use `logDoubleFactorials` instead.
#'
#' @family double factorials
#' @keywords datasets
#' @name doubleFactorials
#' @export
"doubleFactorials"

#' Natural logarithms of double factorials
#'
#' `logDoubleFactorials` is a numeric vector with pre-calculated values of
#' double factorials up to 50 000!!.
#'
#' @family double factorials
#' @keywords datasets
#' @name logDoubleFactorials
#' @export
"logDoubleFactorials"

#' @rdname TreeShape
#' @format `unrootedKeys` is a list of length `r length(unrootedKeys)`; each
#' entry is a vector of integers corresponding to they keys (not shape numbers)
#' of the different unrooted tree shapes with `nTip` leaves.
"unrootedKeys"

#' Number of rooted / unrooted tree shapes
#'
#' `nRootedShapes` and `nUnrootedShapes` give the number of (un)rooted binary
#' trees on _n_ unlabelled leaves.
#'
#' @source
#' `nRootedShapes` corresponds to the Wedderburn-Etherington numbers,
#' [\acronym{OEIS} A001190](https://oeis.org/A001190)
#'
#' `nUnrootedShapes` is [\acronym{OEIS} A000672](https://oeis.org/A000672)
#'
#' @keywords datasets
"nRootedShapes"

#' @rdname nRootedShapes
"nUnrootedShapes"

#' Data from Zhang et al. 2016
#'
#' Phylogenetic data from \insertCite{Zhang2016;textual}{TreeTools} in raw
#' (`Lobo.data`) and `phyDat` (`Lobo.phy`) formats.
#'
#' @template LoboMods
#'
#' @examples
#' data("Lobo", package = "TreeTools")
#' Lobo.data
#' Lobo.phy
#' @source \insertCite{Zhang2016;textual}{TreeTools}
#' 
#' @references \insertAllCited{}
#'
#' @keywords datasets
"Lobo.data"

#' @rdname Lobo.data
"Lobo.phy"
