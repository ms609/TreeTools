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
as.Newick <- function(x) UseMethod("as.Newick")

#' @rdname as.Newick
#' @export
as.Newick.phylo <- function(x) {
  as_newick(x[["edge"]] - 1L)
}

#' @rdname as.Newick
#' @export
as.Newick.list <- function(x) {
  vapply(x, as.Newick, character(1L))
}

#' @rdname as.Newick
#' @export
as.Newick.multiPhylo <- as.Newick.list

#' Write morphological character matrix to TNT file
#' @param dataset Morphological dataset of class `phyDat` or `matrix`.
#' @param filepath Path to file; if `NULL`, returns a character vector.
#' @param comment Optional comment with which to entitle matrix.
#' @param types Optional list specifying where different data types begin.
#' `c(num = 1, dna = 10)` sets characters 1..9 as numeric, 10..end as DNA.
#' @param pre,post Character vector listing text to print before and after the
#' character matrix.  Specify `pre = 'piwe=;` if the matrix is to be analysed
#' using extended implied weighting (`xpiwe=`).
#'
#' @seealso [`ReadTntCharacters()`]
#' @examples
#' data("Lobo", package = "TreeTools")
#'
#' WriteTntCharacters(Lobo.phy)
#'
#' # Read with extended implied weighting
#' WriteTntCharacters(Lobo.phy, pre = "piwe=10;", post = "xpiwe=;")
#'
#' # Write to a file with:
#' # WriteTntCharacters(Lobo.phy, "example_file.tnt")
#' @template MRS
#' @export
WriteTntCharacters <- function(dataset, filepath = NULL,
                                comment = "Dataset written by `TreeTools::WriteTntCharacters()`",
                                types = NULL,
                                pre = "", post = "") {
  UseMethod("WriteTntCharacters")
}

#' @rdname WriteTntCharacters
#' @export
WriteTntCharacters.phyDat <- function(dataset, filepath = NULL,
                                       comment = "Dataset written by `TreeTools::WriteTntCharacters()`",
                                       types = NULL,
                                       pre = "", post = "") {
  WriteTntCharacters(PhyDatToMatrix(dataset), filepath, comment, types,
                     pre, post)
}

#' @rdname WriteTntCharacters
#' @export
WriteTntCharacters.matrix <- function(dataset, filepath = NULL,
                                       comment = "Dataset written by `TreeTools::WriteTntCharacters()`",
                                       types = NULL,
                                       pre = "", post = "") {
  EOL <- "\n"
  dataset <- gsub("(", "[", fixed = TRUE, dataset)
  dataset <- gsub(")", "]", fixed = TRUE, dataset)

  ret <- paste(
    paste(pre, collapse = "\n"),
    paste0("xread '", paste(comment, collapse = " "), "'"),
    paste(rev(dim(dataset)), collapse = " "),
    if (is.null(types)) {
      paste(rownames(dataset),
            apply(dataset, 1, paste0, collapse = ""),
            collapse = EOL)
    } else {
      typeEnds <- c(unname(types[-1]) - 1L, ncol(dataset))
      paste(paste0("&[", names(types), "]\n"),
            vapply(seq_along(types), function(i)
              paste(rownames(dataset),
                    apply(dataset[, types[i]:typeEnds[i]], 1, paste0, collapse = ""),
                    collapse = EOL),
              character(1)), collapse = EOL)
    },
    ";",
    paste(post, collapse = "\n"),
    sep = EOL)
  if (is.null(filepath)) {
    ret
  } else {
    writeLines(ret, filepath)
  }
}
