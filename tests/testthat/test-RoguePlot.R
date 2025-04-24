AllTreesCounted <- function(trees, rogue) {
  x <- RoguePlot(trees, rogue, plot = FALSE)
  expect_equal(length(trees), sum(x$onEdge, x$atNode))
}

ExpectedCons <- function(text) {
  Preorder(read.tree(text = text))
}

LegendLabels <- function(n) {
  paste0(n:1, " tree", c(rep("s", n - 1), ""))
}

test_that(".apply() helper function", {
  Splat <- function (...) paste(..., collapse = " ")
  mat <- matrix(1:12, 3, 4)
  expect_equal(.apply(mat, 2, paste), apply(mat, 2, paste))
  expect_equal(.apply(mat, 2, Splat), matrix(apply(mat, 2, Splat), 1, 4))
  mat1 <- mat[1, , drop = FALSE]
  expect_equal(.apply(mat1, 2, paste), matrix(apply(mat1, 2, paste), 1, 4))
})

test_that("Simple rogue plot", {
  trees <- list(read.tree(text = "(a, (b, (c, (rogue, (d, e)))));"),
                read.tree(text = "(a, (b, (c, ((d, e), rogue))));"),
                read.tree(text = "(a, (b, (c, (rogue, (d, e)))));"),
                read.tree(text = "(a, (b, (c, (d, (rogue, e)))));"))
  expect_equal(
    RoguePlot(trees, tip = "rogue", plot = FALSE),
    list(cons = ExpectedCons("(a, (b, (c, (d, e))));"),
         onEdge = c(0, 0, 0, 0, 0, 3, 0, 1),
         atNode = double(4),
         legendLabels = LegendLabels(4)
         )
  )
  expect_equal(
    RoguePlot(trees, tip = "rogue", plot = FALSE, sort = TRUE)[["cons"]],
    Preorder(SortTree(RoguePlot(trees, tip = "rogue", plot = FALSE)[["cons"]]))
  )

  skip_if_not_installed("vdiffr", "1.0")
  RoguePlotTest <- function() {
    par(mar = rep(0, 4))
    RoguePlot(trees, "rogue",
              Palette = function(...) hcl.colors(..., palette = "inferno"),
              thin = 2, fat = 4, edge.lty = 2,
              legend = "topleft", legend.inset = 0.02
              )
  }
  vdiffr::expect_doppelganger("RoguePlot(simple)", RoguePlotTest)
})

test_that("RoguePlot(sort = TRUE)", {
  trees <- c(PectinateTree(7), PectinateTree(7))
  rp <- RoguePlot(trees, "t5", sort = TRUE, plot = FALSE)
  expect_equal(rp[["atNode"]], rep(0, 5))
  expect_equal(rp[["onEdge"]], `[<-`(double(10), 4, 2))
})

test_that("polytomy id", {
  trees <- list(read.tree(text = "(a,(((b,(c,d)),e),rogue));"),
                read.tree(text = "(a,((((c,d),e),rogue),b));"),
                read.tree(text = "(a,(rogue,(e,((b,d),c))));"))
  AllTreesCounted(trees, "rogue")

  expect_equal(
    RoguePlot(trees, "rogue", plot = FALSE),
    list(cons = ExpectedCons("(a, (b, c, d, e));"),
         onEdge = c(0, 2, 0, 0, 0, 0),
         atNode = c(0, 1),
         legendLabels = LegendLabels(3))
  )

  skip_if_not_installed("vdiffr", "1.0")
  RoguePlotTest <- function() {
    par(mar = rep(0, 4))
    RoguePlot(trees, "rogue", legend = "bottomleft")
  }
  vdiffr::expect_doppelganger("RoguePlot(poly)", RoguePlotTest)
})

# TODO test tree with rogue to left and right of balanced root
test_that("Complex rogue plot", {

  trees1 <- list(read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),  # node 9
                 read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),  # node 9
                 read.tree(text = "(a, (b, (c, (rogue, (d, (f, e))))));"),  # node 9
                 read.tree(text = "(a, (b, (c, (rogue, ((e, f), d)))));"),  # node 9
                 
                 read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),  # node 9 x 5
                 read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),  # node 7 x 1
                 read.tree(text = "(a, (rogue, (b, (c, (d, (e, f))))));"),  # edge 2 x 1
                 read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));"),  # edge 1
                 
                 read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));"),  # edge 1 x 2
                 read.tree(text = "(a, (b, ((c, d), (rogue, (f, e)))));"),  # edge 7
                 read.tree(text = "(a, (b, (((rogue, d), c), (e, f))));"),  # edge 6 x 1
                 read.tree(text = "(a, (b, (c, (d, (rogue, (e, f))))));"))  # edge 7 x 2
  expect_equal(
    RoguePlot(trees = trees1, tip = "rogue", plot = FALSE),
    list(cons = ExpectedCons("(a, (b, (c, d, (e, f))));"),
         onEdge = c(2, 1, 0, 0, 0, 1, 2, 0, 0),
         atNode = c(1, 0, 5, 0),
         legendLabels = LegendLabels(6))
    )

  expectedCons <- Preorder(
    RenumberTips(read.tree(text = "(f, (e, (d, c, (b, a))));"), letters[1:6])
  )
  
  expect_equal(
    RoguePlot(trees1, "rogue", outgroupTips = "f", plot = FALSE),
    list(cons = expectedCons,
         onEdge = c(0, 2, 0, 4, 0, 0, 1, 0, 0),
         atNode = c(0, 0, 5, 0),
         legendLabels = LegendLabels(6)
         )
  )

  trees2 <- list(read.tree(text = "(a, (b, (rogue, ((d, c), (e, f)))));"),
                 read.tree(text = "(a, (b, ((rogue, (d, c)), (e, f))));"),
                 read.tree(text = "(a, (b, ((rogue, (d, c)), (e, f))));"),
                 read.tree(text = "(a, ((b, (c, d)), (rogue, (e, f))));"),
                 read.tree(text = "(a, (b, ((c, d), (rogue, (e, f)))));"),
                 read.tree(text = "(a, (b, ((c, d), (rogue, (e, f)))));"))
  expected <- list(cons = Preorder(read.tree(text = "(a, (b, (c, d), (e, f)));")),
                   onEdge = c(0, 0, 0, 2, 0, 0, 3, 0, 0),
                   atNode = c(0, 1, 0, 0),
                   legendLabels = LegendLabels(5))
  actual <- RoguePlot(trees = trees2, tip = "rogue", plot = FALSE)
  expect_equal(names(actual), names(expected))
  expect_true(all.equal(actual$cons, expected$cons))
  expect_equal(actual$onEdge, expected$onEdge)
  expect_equal(actual$atNode, expected$atNode)

  trees3 <- lapply(c(9, 10, 13), AddTip, tree = BalancedTree(8), label = "rogue")
  AllTreesCounted(trees3, "rogue")

  skip_if_not_installed("vdiffr", "1.0")
  RoguePlotTest <- function() {
    par(mar = rep(0, 4))
    RoguePlot(trees1, "rogue",
              Palette = function(...) hcl.colors(..., palette = "inferno"),
              edgeLen = 1, thin = 2, fat = 4,
              legend = "bottomright"
              )
  }
  vdiffr::expect_doppelganger("RoguePlot(trees1)", RoguePlotTest)

  RoguePlotTest <- function() {
    par(mar = rep(0, 4))
    RoguePlot(trees2, "rogue", edgeLength = 1:7, thin = 2, fat = 4)
  }
  vdiffr::expect_doppelganger("RoguePlot(trees2)", RoguePlotTest)
})
