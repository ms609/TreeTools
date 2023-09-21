test_that("Trees are sorted", {
  trees <- setNames(as.phylo(0:8, 8), paste0("as.phy", 0:8))
  expect_true(trees[[1]] == trees[[1]])
  expect_false(trees[[1]] == trees[[2]])
  
  for (i in 1:8) {
    expect_true(trees[[i]] == trees[[i]])
    expect_false(trees[[i]] > trees[[i]])
    expect_false(trees[[i]] < trees[[i]])
    for (j in (i + 1):9) {
      expect_true(trees[[i]] < trees[[j]])
      expect_true(trees[[j]] > trees[[i]])
      expect_false(trees[[j]] == trees[[i]])
      expect_false(trees[[i]] == trees[[j]])
      expect_false(trees[[i]] > trees[[j]])
      expect_false(trees[[j]] < trees[[i]])
    }
  }
  
  expect_equal(sort(trees), sort(rev(trees)))
  expect_equal(names(sort(rev(trees))), names(trees))
})
