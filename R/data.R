#' Brewer palettes
#'
#' A list of eleven Brewer palettes containing one to eleven colours that
#' are readily distinguished by colourblind viewers, followed by a twelfth
#' 12-colour palette adapted for colour blindness.
#'
#' @source {
#'
#' * [ColourBrewer2.org](http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3)
#' 
#' * [Martin Krzywinski](http://mkweb.bcgsc.ca/colorblind/)
#' }
#'
#' @examples 
#' data("brewer", package="TreeTools")
#' plot(0, type='n', xlim=c(1, 12), ylim=c(12, 1), 
#'      xlab = 'Colour', ylab='Palette')
#' for (i in seq_along(brewer)) text(seq_len(i), i, col=brewer[[i]])
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
#' @family Double factorial
#' @keywords datasets
#' @name doubleFactorials
#' @export
"doubleFactorials"

#' Natural logarithms of double factorials
#' 
#' A vector with pre-calculated values of double factorials up to 50 000!!.
#' 
#' @family Double factorial
#' @keywords datasets
#' @name logDoubleFactorials
#' @export
"logDoubleFactorials"

#' Raw data from Zhang et al. 2016
#'
#' @template LoboMods
#'
#' @source 
#'  \insertRef{Zhang2016}{TreeTools}
#' 
#' @keywords datasets
"Lobo.data"

#' Data from Zhang et al. 2016 in phyDat format
#'
#' @template LoboMods
#'
#' @source
#'  \insertRef{Zhang2016}{TreeTools}
#' 
#' @keywords datasets
"Lobo.phy"
