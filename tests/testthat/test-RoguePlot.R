test_that('Simple rogue plot', {
  trees <- list(read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
                read.tree(text = '(a1, (b1, (c1, (d1, (rogue, e1)))));'))
  expect_equal(c(0, 0, 3, 0, 1, 0, 0, 0),
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
                read.tree(text = '(rogue, (a, (b, (c, (d, (e, f))))));'),
                read.tree(text = '((rogue, a), (b, (c, (d, (e, f)))));'),
                read.tree(text = '(a, (b, ((c, d), (rogue, (e, f)))));'),
                read.tree(text = '(a, (b, ((c, (rogue, d)), (e, f))));'),
                read.tree(text = '(a, (b, (c, (d, (rogue, (e, f))))));'))
  RoguePlot(trees, 'rogue', p = 1)
  expect_equal(c(0, 0, 3, 0, 1, 0, 0, 0),
               RoguePlot(trees, 'rogue', plot = FALSE))

  skip_if_not_installed('vdiffr')
  RoguePlotTest <- function () {
    par(mar = rep(0, 4))
    RoguePlot(trees, 'rogue', Palette = viridisLite::inferno,
              thin = 2, edge.lty = 2)
  }
  vdiffr::expect_doppelganger('RoguePlot()', RoguePlotTest)
})

