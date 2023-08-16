test_that("MSTEdges() handles bad input", {
  expect_error(MSTEdges(matrix(1:12, 6, 2)), "distance")
})

test_that("minimum_spanning_tree.cpp handles bad input", {
  expect_equal(minimum_spanning_tree(numeric(0)), matrix(0, 0, 0))
  expect_error(minimum_spanning_tree(c(1:-1)),
               "`order` contains entries < 0")
  expect_error(minimum_spanning_tree(c(3, 100, 1)),
               "`order` contains entries > `length.order.`")
  expect_error(minimum_spanning_tree(c(3, 1, NA_real_)),
               "`order` contains NA")
  expect_error(minimum_spanning_tree(0:13),
               "`length.order.`.* not.* triangular")
})

test_that("MST edges calculated correctly", {
  set.seed(0)
  points <- matrix(c(0.1, 0, 1.9, 2, 1.1, 1,
                     0.1, 2,   0, 2, 1, 1.1,
                       0, 0,   0, 0, 1, -1), 6)
  distances <- dist(points)
  apeMst <- matrix(c(5, 6, 6, 5, 5, 1, 1:4), 5)
  distMat <- as.matrix(distances)
  expect_equal(MSTLength(distances, apeMst),
               MSTLength(distances))
  expect_equal(MSTLength(distances),
               MSTLength(distMat))
  MSTPlot <- function() {
    plot(points, asp = 1, ann = FALSE)
    expect_equal(MSTEdges(distances, FALSE),
                 MSTEdges(distances, TRUE, points[, 1], points[, 2]))
  }
  skip_if_not_installed("vdiffr", minimum_version = "1.0.0")
  skip_if(packageVersion("graphics") < "4.1.0")
  vdiffr::expect_doppelganger("MST plotting", MSTPlot)
})

test_that("MST handles large distance matrices", {
  x <- dist(0:300)
  expect_equal(c(300, 2), dim(MSTEdges(x)))
})
