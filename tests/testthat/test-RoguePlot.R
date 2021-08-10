devtools::load_all('c:/research/r/ape')

AllTreesCounted <- function (trees, rogue) {
  x <- RoguePlot(trees, rogue, plot = FALSE)
  expect_equal(length(trees), sum(x$onEdge, x$atNode))
}

test_that('Simple rogue plot', {
  skip_if_not_installed('ape', '5.5.2')
  trees <- list(read.tree(text = '(a, (b, (c, (rogue, (d, e)))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, e)))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, e)))));'),
                read.tree(text = '(a, (b, (c, (d, (rogue, e)))));'))
  expect_equal(list(cons = Preorder(read.tree(text = '(a, (b, (c, (d, e))));')),
                    onEdge = c(0, 0, 0, 0, 0, 3, 0, 1),
                    atNode = double(4)),
               RoguePlot(trees, 'rogue', plot = FALSE))

  skip_if_not_installed('vdiffr')
  RoguePlotTest <- function () {
    par(mar = rep(0, 4))
    RoguePlot(trees, 'rogue',
              Palette = function(...) hcl.colors(..., palette = 'inferno'),
              thin = 2, fat = 4, edge.lty = 2)
  }
  vdiffr::expect_doppelganger('RoguePlot(simple)', RoguePlotTest)
})

test_that("polytomy id", {
  trees <- list(read.tree(text = "(a,(((b,(c,d)),e),f));"),
                read.tree(text = "(a,((((c,d),e),f),b));"),
                read.tree(text = "(a,(f,(e,((b,d),c))));"))
  AllTreesCounted(trees, 'f')
})

# TODO test tree with rogue to left and right of balanced root
test_that('Complex rogue plot', {

  skip_if_not_installed('ape', '5.5.2')
  trees1 <- list(read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),  # node 9
                 read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),  # node 9
                 read.tree(text = '(a, (b, (c, (rogue, (d, (f, e))))));'),  # node 9
                 read.tree(text = '(a, (b, (c, (rogue, ((e, f), d)))));'),  # node 9
                 read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),  # node 9 x 5
                 read.tree(text = '(rogue, (a, (b, (c, (d, (e, f))))));'),  # node 7 x 1
                 read.tree(text = '(a, (rogue, (b, (c, (d, (e, f))))));'),  # edge 2 x 1
                 read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),  # edge 1
                 read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),  # edge 1 x 2
                 read.tree(text = '(a, (b, ((c, d), (rogue, (f, e)))));'),  # edge 7
                 read.tree(text = '(a, (b, (((rogue, d), c), (e, f))));'),  # edge 6 x 1
                 read.tree(text = '(a, (b, (c, (d, (rogue, (e, f))))));'))  # edge 7 x 2
  expect_equal(list(cons = Preorder(read.tree(text = '(a, (b, (c, d, (e, f))));')),
                    onEdge = c(2, 1, 0, 0, 0, 1, 2, 0, 0),
                    atNode = c(1, 0, 5, 0)),
               RoguePlot(trees1, 'rogue', plot = FALSE))

  trees2 <- list(read.tree(text = '(a, (b, (rogue, ((d, c), (e, f)))));'),
                 read.tree(text = '(a, (b, ((rogue, (d, c)), (e, f))));'),
                 read.tree(text = '(a, (b, ((rogue, (d, c)), (e, f))));'),
                 read.tree(text = '(a, ((b, (c, d)), (rogue, (e, f))));'),
                 read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'),
                 read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'))
  expect_equal(list(cons = Preorder(read.tree(text = '(a, (b, (c, d, (e, f))));')),
                    onEdge = c(0, 0, 0, 2, 0, 0, 3, 0, 0),
                    atNode = c(0, 1, 0, 0)),
               RoguePlot(trees2, 'rogue', plot = FALSE))

  trees3 <- lapply(c(9, 10, 13), AddTip, tree = BalancedTree(8), label = 'rogue')
  AllTreesCounted(trees3, 'rogue')

  skip_if_not_installed('vdiffr')
  RoguePlotTest <- function () {
    par(mar = rep(0, 4))
    RoguePlot(trees1, 'rogue',
              Palette = function(...) hcl.colors(..., palette = 'inferno'),
              edgeLen = 1, thin = 2, fat = 4)
  }
  vdiffr::expect_doppelganger('RoguePlot(trees1)', RoguePlotTest)

  RoguePlotTest <- function () {
    par(mar = rep(0, 4))
    RoguePlot(trees2, 'rogue', edgeLength = 1:7, thin = 2, fat = 4)
  }
  vdiffr::expect_doppelganger('RoguePlot(trees2)', RoguePlotTest)
})

