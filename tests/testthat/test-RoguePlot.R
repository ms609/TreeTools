test_that('Simple rogue plot', {
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
    RoguePlot(trees, 'rogue', Palette = viridisLite::inferno,
              thin = 2, edge.lty = 2)
  }
  vdiffr::expect_doppelganger('RoguePlot()', RoguePlotTest)
})

# TODO test tree with rogue to left and right of balanced root
test_that('Complex rogue plot', {

  trees <- list(read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
                read.tree(text = '(a, (b, (c, (rogue, (d, (e, f))))));'),
                read.tree(text = '(rogue, (a, (b, (c, (d, (e, f))))));'),
                read.tree(text = '(a, (rogue, (b, (c, (d, (e, f))))));'),
                read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),
                read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),
                read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'),
                read.tree(text = '(a, (b, ((c, (rogue, d)), (e, f))));'),
                read.tree(text = '(a, (b, (c, (d, (rogue, (e, f))))));'))
  RoguePlot(trees, 'rogue', p = 1)
  expect_equal(list(cons = Preorder(read.tree(text = '(a, (b, (c, d, (e, f))));')),
                    onEdge = c(2, 1, 0, 0, 0, 1, 2, 0, 0),
                    atNode = c(1, 0, 5, 0)),
               RoguePlot(trees, 'rogue', plot = FALSE))

  skip_if_not_installed('vdiffr')
  RoguePlotTest <- function () {
    par(mar = rep(0, 4))
    RoguePlot(trees, 'rogue', Palette = viridisLite::inferno,
              thin = 2, edge.lty = 2)
  }
  vdiffr::expect_doppelganger('RoguePlot()', RoguePlotTest)
})

