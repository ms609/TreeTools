#' Colour a tree by topology
#'
#' `PaintTree()` assigns a colour to every edge, leaf, and internal node of a
#' tree such that sister clades occupy adjacent hue bands proportional to their
#' tip counts, and saturation grows from zero at the root to one at every tip.
#' The result is a list of vectors that correspond to 
#' [`plot.phylo()`][ape::plot.phylo]'s
#' `edge.color`, `tip.color`, and `node.color` arguments.
#'
#' Hue is allocated recursively from a 0–360° budget at the root.  At each
#' internal node, the parent's budget is partitioned across its descendant
#' edges in proportion to the number of leaves descended from each child; the
#' colour reported for an edge or node is the midpoint of its assigned range.
#' Saturation `s = (nTip - nDesc) / (nTip - 1)`, where `nDesc` is the number of
#' leaves descended from the node (including itself for tips), so the root is
#' achromatic (`s = 0`) and every leaf is fully saturated (`s = 1`).
#'
#' @template treeParam
#' @param palette Either a character string naming one of the built-in
#' palettes (`"default"`, `"protanopia"`, `"tritanopia"`,
#' or a function `function(h, s)` taking vectors of hues in
#' `[0, 360]` and saturations in `[0, 1]` and returning a character vector of
#' colours.  The `"protanopia"` and `"tritanopia"` palettes simulate
#' dichromatic perception of the default HCL spectrum using the
#' Viénot–Brettel–Mollon projection.
#'
#' @return A list with three character vectors of hex colours:
#' \describe{
#'   \item{`edgeCol`}{One colour per edge in `tree$edge`, taken from the
#'   child node's `(hue, saturation)`.}
#'   \item{`tipCol`}{One colour per leaf, indexed `1:NTip(tree)`.}
#'   \item{`nodeCol`}{One colour per internal node, indexed so that entry `i`
#'   corresponds to node `NTip(tree) + i`; the root (entry 1) is grey
#'   because its saturation is zero.}
#' }
#'
#' @examples
#' tree <- BalancedTree(1:8) + PectinateTree(9:14)
#' tree <- ape::bind.tree(BalancedTree(1:8), PectinateTree(9:12), 8, 4)
#' paint <- PaintTree(tree)
#' plot(tree,
#'      edge.color = paint$edgeCol,
#'      tip.color = paint$tipCol,
#'      edge.width = 3)
#'
#' # Colour-blind-safe variants
#' paintP <- PaintTree(tree, "protanopia")
#' plot(tree,
#'      edge.color = paintP$edgeCol,
#'      tip.color = paintP$tipCol,
#'      edge.width = 3)
#'
#' paintT <- PaintTree(tree, "tritanopia")
#' plot(tree,
#'      edge.color = paintT$edgeCol,
#'      tip.color = paintT$tipCol,
#'      edge.width = 3)
#'
#' # User-supplied palette function (greyscale)
#' grey_pal <- function(h, s) grey(1 - s * 0.8)
#' paintG <- PaintTree(tree, grey_pal)
#' plot(tree, edge.color = paintG$edgeCol, edge.width = 3)
#'
#' @seealso [`CladeSizes()`], [`DescendantEdges()`]
#' @template MRS
#' @family tree navigation
#' @importFrom grDevices hcl col2rgb rgb
#' @export
PaintTree <- function(tree, palette = "default") {
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  nTip <- NTip(tree)
  nNode <- tree[["Nnode"]]
  maxNode <- nTip + nNode

  Pal <- .ResolvePalette(palette)

  if (nTip < 2L) {
    return(list(
      edgeCol = if (length(parent)) Pal(rep_len(180, length(parent)),
                                        rep_len(0, length(parent)))
                else character(0),
      tipCol = Pal(rep_len(180, nTip), rep_len(1, nTip)),
      nodeCol = Pal(rep_len(180, nNode), rep_len(0, nNode))
    ))
  }

  nDesc <- CladeSizes(tree, internal = FALSE)

  loH <- numeric(maxNode)
  hiH <- numeric(maxNode)
  root <- setdiff(parent, child)
  loH[root] <- 0
  hiH[root] <- 360

  # Order in which to process parents: preorder (root first, descendants later).
  # Within each parent, children are taken in their original `tree$edge` order
  # so the result does not depend on whether the tree was passed in preorder,
  # postorder, or any other valid edge ordering.
  parentOrder <- unique(parent[rev(PostorderOrder(tree))])
  groups <- split(seq_along(parent),
                  factor(parent, levels = parentOrder))

  for (g in groups) {
    p <- parent[g[1]]
    lo <- loH[p]
    hi <- hiH[p]
    kids <- child[g]
    w <- nDesc[kids]
    cuts <- lo + (hi - lo) * cumsum(c(0, w)) / sum(w)
    loH[kids] <- cuts[-length(cuts)]
    hiH[kids] <- cuts[-1]
  }

  hue <- (loH + hiH) / 2
  sat <- (nTip - nDesc) / (nTip - 1)

  list(
    edgeCol = Pal(hue[child], sat[child]),
    tipCol  = Pal(hue[seq_len(nTip)], rep_len(1, nTip)),
    nodeCol = Pal(hue[nTip + seq_len(nNode)], sat[nTip + seq_len(nNode)])
  )
}

.ResolvePalette <- function(palette) {
  if (is.function(palette)) {
    return(palette)
  }
  if (!is.character(palette) || length(palette) != 1L) {
    stop("`palette` must be a single string or a `function(h, s)`.")
  }
  choices <- c("default", "protanopia", "tritanopia")
  idx <- pmatch(tolower(palette), choices)
  if (is.na(idx)) {
    stop("`palette` must match one of: ",
         paste(choices, collapse = ", "),
         " (or be a function).")
  }
  switch(choices[idx],
         default    = .PaletteDefault,
         protanopia = .PaletteProtanopia,
         tritanopia = .PaletteTritanopia)
}

.PaletteDefault <- function(h, s) {
  hcl(h = h %% 360, c = s * 70, l = 70, fixup = TRUE)
}

.srgb2linear <- function(x) {
  ifelse(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ^ 2.4)
}

.linear2srgb <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  ifelse(x <= 0.0031308, 12.92 * x, 1.055 * x ^ (1 / 2.4) - 0.055)
}

# Viénot, Brettel & Mollon (1999) sRGB projection for protanopia.
.protanopiaMatrix <- matrix(c(
  0.567, 0.433, 0.000,
  0.558, 0.442, 0.000,
  0.000, 0.242, 0.758
), nrow = 3, byrow = TRUE)

.tritanopiaMatrix <- matrix(c(
  0.950, 0.050, 0.000,
  0.000, 0.433, 0.567,
  0.000, 0.475, 0.525
), nrow = 3, byrow = TRUE)

.ApplyDichromat <- function(h, s, M) {
  base <- .PaletteDefault(h, s)
  rgbMat <- col2rgb(base) / 255
  lin <- matrix(.srgb2linear(as.numeric(rgbMat)), nrow = 3)
  proj <- M %*% lin
  out <- matrix(.linear2srgb(as.numeric(proj)), nrow = 3)
  rgb(out[1, ], out[2, ], out[3, ])
}

.PaletteProtanopia <- function(h, s) {
  .ApplyDichromat(h, s, .protanopiaMatrix)
}

.PaletteTritanopia <- function(h, s) {
  .ApplyDichromat(h, s, .tritanopiaMatrix)
}
