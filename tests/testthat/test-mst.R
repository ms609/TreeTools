context('mst.R')
test_that('MST edges calculated correctly', {
  set.seed(0)
  points <- matrix(c(0.1, 0, 1.9, 2, 1.1, 1,
                     0.1, 2,   0, 2, 1, 1.1,
                       0, 0,   0, 0, 1, -1), 6)
  distances <- dist(points)
  expect_equal(matrix(c(5, 6, 6, 5, 5, 1, 1:4), 5),
               MSTEdges(distances, FALSE))
  MSTPlot <- function () {
    plot(points, asp = 1, ann = FALSE)
    expect_equal(MSTEdges(distances, FALSE),
                 MSTEdges(distances, TRUE, points[, 1], points[, 2]))
  }
  vdiffr::expect_doppelganger('MST plotting', MSTPlot)
})
