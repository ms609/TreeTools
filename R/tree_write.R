#' Write a phylogenetic tree in Newick format
#'
#' `as.Newick()` creates a character string representation of a phylogenetic
#' tree, in the Newick format, using R's internal tip numbering.
#' Use [`RenumberTips()`] to ensure that the internal numbering follows the
#' order you expect.
#'
#'
#' @param x Object to convert to Newick format.
#' See Usage section for supported classes.
#'
#' @return `as.Newick()` returns a character string representing `tree` in Newick
#' format.
#'
#' @examples
#' trees <- list(BalancedTree(1:8), PectinateTree(8:1))
#' trees <- lapply(trees, RenumberTips, 1:8)
#' as.Newick(trees)
#'
#' @seealso
#' - Retain leaf labels: [`NewickTree()`]
#'
#' - Change R's internal numbering of leaves: [`RenumberTips()`]
#'
#' - Write tree to text or file: [`ape::write.tree()`]
#'
#' @template MRS
#' @export
as.Newick <- function (x) UseMethod('as.Newick')

#' @rdname as.Newick
#' @export
as.Newick.phylo <- function (x) {
  as_newick(Preorder(x)$edge - 1L)
}

#' @rdname as.Newick
#' @export
as.Newick.list <- function (x) {
  vapply(x, as.Newick, character(1L))
}

#' @rdname as.Newick
#' @export
as.Newick.multiPhylo <- as.Newick.list

#' Write morphological character matrix to TNT file
#' @param dataset Morphological dataset of class `phyDat` or `matrix`.
#' @param filepath Path to file; if `NULL`, returns a character vector.
#' @param comment Optional comment with which to entitle matrix.
#'
#' @template MRS
#' @examples
#' data('Lobo', package = 'TreeTools')
#' WriteTntCharacters(Lobo.phy)
#' # Read with extended implied weighting
#' WriteTntCharacters(Lobo.phy, 'TEST.tnt', pre = 'piwe=10;', post = 'xpiwe=;')
#' @export
WriteTntCharacters <- function (dataset, filepath = NULL,
                                comment = 'Dataset written by `TreeTools::WriteTntCharacters()`',
                                pre = '', post = '') {
  UseMethod('WriteTntCharacters')
}

#' @rdname WriteTntCharacters
#' @export
WriteTntCharacters.phyDat <- function (dataset, filepath = NULL,
                                       comment = 'Dataset written by `TreeTools::WriteTntCharacters()`',
                                       pre = '', post = '') {
  WriteTntCharacters(PhyDatToMatrix(dataset), filepath, comment, pre, post)
}

#' @rdname WriteTntCharacters
#' @export
WriteTntCharacters.matrix <- function (dataset, filepath = NULL,
                                       comment = 'Dataset written by `TreeTools::WriteTntCharacters()`',
                                       pre = '', post = '') {
  EOL <- '\n'
  ret <- paste(
    pre,
    paste0("xread '", comment, "'"),
    paste(rev(dim(dataset)), collapse = ' '),
    '&[num]',
    paste(rownames(dataset), apply(dataset, 1, paste0, collapse = ''), collapse = EOL),
    ';',
    post,
    sep = EOL)
  if (is.null(filepath)) {
    ret
  } else {
    writeLines(ret, filepath)
  }
}
