#' Visualise position of rogue taxa
#'
#' Plots a consensus of trees with a rogue taxon omitted, with edges coloured
#' according to the proportion of trees in which the taxon attaches to that
#' edge, after Klopfstein &amp; Spasojevic (2019).
#'
#' @param trees List or `multiPhylo` object containing phylogenetic trees
#' of class `phylo` to be summarized.
#' @param tip Numeric or character identifying rogue leaf, in format accepted
#' by `drop.tip()`.
#' @param p A numeric value between 0.5 and 1 giving the proportion for a clade
#' to be represented in the consensus tree (see `consensus()`).
#' @param \dots Additional parameters to `plot.phylo()`.
#' @param plot Logical specifying whether to plot the tree.
#' @param Palette Function that takes a parameter `n` and generates a colour
#' palette with `n` entries.
#' @param thin,fat Numeric specifying width to plot edges if the rogue tip
#' never / sometimes does attach to them.
#' @return `RoguePlot()` returns a vector of integers specifying the number of
#' trees in `trees` in which the rogue leaf is attached to each edge in turn
#' of the consensus tree.
#' @references
#' \insertRef{Klopfstein2019}{TreeTools}
#' @examples
#' trees <- list(read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
#'               read.tree(text = '(rogue, (a, (b, (c, (d, (e, f))))));'),
#'               read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),
#'               read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'),
#'               read.tree(text = '(a, (b, ((c, (rogue, d)), (e, f))));'),
#'               read.tree(text = '(a, (b, (c, (d, (rogue, (e, f))))));'))
#' RoguePlot(trees, 'rogue')
#' @template MRS
#' @importFrom ape consensus drop.tip
#' @importFrom phangorn allDescendants
#' @importFrom grDevices colorRampPalette
#' @export
RoguePlot <- function (trees, tip, p = 1, plot = TRUE,
                       Palette = colorRampPalette(c(par('fg'), 'red'), space = 'Lab'),
                       thin = par('lwd'), fat = thin + 1L,
                       ...) {
  trees <- RenumberTips(trees, trees[[1]])
  if (is.character(tip)) {
    tip <- match(tip, trees[[1]]$tip.label)
  }

  noRogue <- trees
  noRogue[] <- lapply(trees, drop.tip, tip)
  cons <- Preorder(consensus(noRogue, p = p))
  consTip <- NTip(cons)

  nTip <- NTip(trees[[1]])
  #allTips <- logical(ceiling((nTip) / 8) * 8L + 2L) # Multiple of 8 for packBits
  allTips <- logical(nTip)
  index <- vapply(trees, function (tr) {
    edge <- tr$edge
    parent <- edge[, 1]
    child <- edge[, 2]
    rogueEdge <- match(tip, child)
    rogueParent <- parent[rogueEdge]
    aboveRogue <- match(rogueParent, child)
    splitTips <- allTips
    if (is.na(aboveRogue)) {
      splitKids <- child
    } else {
      edgeInSplit <- DescendantEdges(edge = aboveRogue, parent = parent, child = child)
      splitKids <- child[edgeInSplit]
    }
    splitTips[splitKids[splitKids <= nTip]] <- TRUE
    splitTips <- splitTips[-tip]

    #ret <- packBits(splitTips[-1]) # TODO can do something clever like this?
    if (splitTips[1]) {
      !splitTips[-1]
    } else {
      splitTips[-1]
    }
  #}, raw((length(allTips) - 2L)  / 8L))
  }, logical((length(allTips) - 2L)))
  if (is.null(dim(index))) {
    index <- matrix(index, nrow = 1)
  }

  indexSums <- colSums(index)
  nAtTip <- c(sum(indexSums == 0L), double(nTip - 1L))
  atTip <- indexSums == 1L
  tipMatches <- 1L + apply(index[, atTip, drop = FALSE], 2, which)
  tab <- table(as.integer(tipMatches))
  nAtTip[as.integer(names(tab))] <- tab
  nAtTip <- nAtTip[-tip]

  unmatched <- indexSums > 1L
  consSplits <- as.Splits(cons)
  splits <- as.logical(consSplits)
  nSplits <- nrow(splits)
  edgeMatches <- match(data.frame(index[, unmatched, drop = FALSE]),
                       data.frame(t(splits[, -1])))


  nAtSplit <- double(nSplits)
  tab <- table(edgeMatches)
  nAtSplit[as.integer(names(tab))] <- tab
  decipher <- c(nAtTip, rep.int(NA, cons$Nnode))
  decipher[as.integer(names(consSplits))] <- nAtSplit
  nOnEdge <- decipher[cons$edge[, 2]]
  nOnEdge[is.na(nOnEdge)] <- nOnEdge[1] # second root edge


  nAtNode <- double(cons$Nnode)
  unmatched[unmatched] <- is.na(edgeMatches)
  if (any(unmatched)) {
    nodeDescs <- vapply(allDescendants(cons)[-seq_len(nTip)], function (tips) {
      ret <- logical(nTip)
      ret[tips[tips <= nTip]] <- TRUE
      if (ret[1]) !ret[-1] else ret[-1]
    }, logical(nTip - 1L))
    index <- rbind(FALSE, index)

    for (i in seq_len(ncol(nodeDescs))) {
      overlap <- index[nodeDescs[, i], unmatched]
      active <- apply(overlap, 2, any)
      within <- apply(overlap, 2, all)
      unmatched[active & within] <- FALSE
      nAtNode[i] <- sum(active & within)
      if (!any(unmatched)) {
        break
      }
    }
  }

  pal <- Palette(length(trees) + 1L)
  if (plot) {
    .plot.phylo(cons, edge.color = pal[nOnEdge + 1L],
                node.color = pal[c(double(consTip), nAtNode) + 1L],
                edge.width = ifelse(nOnEdge > 0, fat, thin), ...)
  }

  # Return:
  list(onEdge = nOnEdge, atNode = nAtNode)
}

# Modified from ape::plot.phylo
.plot.phylo <- function (x, type = "phylogram", use.edge.length = TRUE,
                         node.pos = NULL, show.tip.label = TRUE,
                         show.node.label = FALSE, edge.color = "black",
                         node.color = par('fg'),
                         edge.width = 1, edge.lty = 1, font = 3,
                         cex = par("cex"), adj = NULL, srt = 0,
                         no.margin = FALSE, root.edge = FALSE, label.offset = 0,
                         underscore = FALSE, x.lim = NULL, y.lim = NULL,
                         direction = "rightwards", lab4ut = NULL,
                         tip.color = "black", plot = TRUE, rotate.tree = 0,
                         open.angle = 0, node.depth = 1,
                         align.tip.label = FALSE, ...) {
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height,
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge),
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth,
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[,
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge,
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[,
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length),
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan",
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards",
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label)
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length ||
        is.ultrametric(x))
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length ||
      is.null(x$root.edge) || !x$root.edge)
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos))
      node.pos <- if (type == "cladogram" && !use.edge.length)
        2
    else 1
    if (node.pos == 1)
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[,
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge),
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2)
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge,
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge,
                                 Nedge, z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360,
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge,
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge,
                        node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge,
                          node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode,
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip,
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards")
        xx <- xx + x$root.edge
      if (direction == "upwards")
        yy <- yy + x$root.edge
    }
  }
  if (no.margin)
    par(mai = rep(0, 4))
  if (show.tip.label)
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin))
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1,
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal)
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label)
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
                         cex)
    if (type == "radial")
      x.lim[1] <- if (show.tip.label)
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards")
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal)
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2,
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy *
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal)
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label)
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
                         cex)
    if (type == "radial")
      y.lim[1] <- if (show.tip.label)
        -1 - max(nchar.tip.label * 0.018 * max.yy *
                   cex)
    else -1
  }
  if (phyloORclado && direction == "downwards")
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards")
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards")
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted"))
    1
  else NA
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim,
               xlab = "", ylab = "", axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj))
      adj <- if (phyloORclado && direction == "leftwards")
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 -
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] -
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      .phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal,
                      edge.color, edge.width, edge.lty,
                      color.v = if (is.null(node.color)) rep(par('fg'), Nnode)
                                    else node.color[-seq_len(Ntip)])
    } else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta,
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color,
                          edge.width, edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) edge.color else "black"
      rootw <- if (length(edge.width) == 1) edge.width else 1
      rootlty <- if (length(edge.lty) == 1) edge.lty else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol,
                 lwd = rootw, lty = rootlty)
      } else {
        switch(direction, rightwards = segments(0, yy[ROOT],
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw,
                                                lty = rootlty), leftwards = segments(xx[ROOT],
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT],
                                                                                     col = rootcol, lwd = rootw, lty = rootlty),
               upwards = segments(xx[ROOT], 0, xx[ROOT],
                                  x$root.edge, col = rootcol, lwd = rootw,
                                  lty = rootlty), downwards = segments(xx[ROOT],
                                                                       yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge,
                                                                       col = rootcol, lwd = rootw, lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label))
        underscore <- TRUE
      if (!underscore)
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]),
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip],
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip],
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]),
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp,
                   lty = align.tip.label.lty)
        } else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label,
             adj = adj, font = font, srt = srt, cex = cex,
             col = tip.color)
      } else {
        angle <- if (type == "unrooted") {
          XY$axe
        } else {
          atan2(yy[1:Ntip], xx[1:Ntip])
        }
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") "horizontal" else "axial"
        } else {
          match.arg(lab4ut, c("horizontal", "axial"))
        }
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] *
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 *
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] *
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] *
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj *
                 cex, x$tip.label, adj = c(adj, 0), font = font,
               srt = srt, cex = cex, col = tip.color)
        } else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips,
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i],
                                 x$tip.label[i], font = font[i], cex = cex[i],
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) {
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
           x$node.label, adj = adj, font = font, srt = srt,
           cex = cex)
    }
  }
  L <- list(type = type, use.edge.length = use.edge.length,
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label,
            show.node.label = show.node.label, font = font, cex = cex,
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
            x.lim = x.lim, y.lim = y.lim, direction = direction,
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time,
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx,
                                      yy = yy)), envir = .PlotPhyloEnv)
  invisible(L)
}

# Modified from ape::phylogram.plot
.phylogram.plot <- function (edge, Ntip, Nnode, xx, yy, horizontal, edge.color,
          edge.width, edge.lty, color.v = rep(par('fg'), Nnode)) {
  nodes <- (Ntip + 1):(Ntip + Nnode)
  if (!horizontal) {
    tmp <- yy
    yy <- xx
    xx <- tmp
  }
  x0v <- xx[nodes]
  y0v <- y1v <- numeric(Nnode)
  NodeInEdge1 <- vector("list", Nnode)
  e1 <- edge[, 1]
  for (i in seq_along(e1)) {
    j <- e1[i] - Ntip
    NodeInEdge1[[j]] <- c(NodeInEdge1[[j]], i)
  }
  for (i in 1:Nnode) {
    j <- NodeInEdge1[[i]]
    tmp <- range(yy[edge[j, 2]])
    y0v[i] <- tmp[1]
    y1v[i] <- tmp[2]
  }
  x0h <- xx[edge[, 1]]
  x1h <- xx[edge[, 2]]
  y0h <- yy[edge[, 2]]
  nc <- length(edge.color)
  nw <- length(edge.width)
  nl <- length(edge.lty)
  if (nc + nw + nl == 3) {
    color.v <- edge.color
    width.v <- edge.width
    lty.v <- edge.lty
  } else {
    Nedge <- dim(edge)[1]
    edge.color <- rep(edge.color, length.out = Nedge)
    edge.width <- rep(edge.width, length.out = Nedge)
    edge.lty <- rep(edge.lty, length.out = Nedge)
    DF <- data.frame(edge.color, edge.width, edge.lty, stringsAsFactors = FALSE)
    width.v <- rep(1, Nnode)
    lty.v <- rep(1, Nnode)
    for (i in 1:Nnode) {
      br <- NodeInEdge1[[i]]
      if (length(br) == 1) {
        A <- br[1]
        color.v[i] <- edge.color[A]
        width.v[i] <- edge.width[A]
        lty.v[i] <- edge.lty[A]
      } else if (length(br) > 2) {
        x <- unique(DF[br, 1])
        if (length(x) == 1) {
          color.v[i] <- x
        }
        x <- unique(DF[br, 2])
        if (length(x) == 1) {
          width.v[i] <- x
        }
        x <- unique(DF[br, 3])
        if (length(x) == 1) {
          lty.v[i] <- x
        }
      } else {
        A <- br[1]
        B <- br[2]
        if (any(DF[A, ] != DF[B, ])) {
          color.v[i] <- edge.color[B]
          width.v[i] <- edge.width[B]
          lty.v[i] <- edge.lty[B]
          y0v <- c(y0v, y0v[i])
          y1v <- c(y1v, yy[i + Ntip])
          x0v <- c(x0v, x0v[i])
          color.v <- c(color.v, edge.color[A])
          width.v <- c(width.v, edge.width[A])
          lty.v <- c(lty.v, edge.lty[A])
          y0v[i] <- yy[i + Ntip]
        } else {
          color.v[i] <- edge.color[A]
          width.v[i] <- edge.width[A]
          lty.v[i] <- edge.lty[A]
        }
      }
    }
  }
  if (horizontal) {
    segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width,
             lty = edge.lty)
    segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v,
             lty = lty.v)
  }
  else {
    segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width,
             lty = edge.lty)
    segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v,
             lty = lty.v)
  }
}
