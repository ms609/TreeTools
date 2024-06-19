test_that("Trees from Mir et al. 2013 are scored correctly", {
  Tree <- function(text) ape::read.tree(text = text)
  expect_identical(0L,  TotalCopheneticIndex(Tree("(1,2,3,4,5);")))
  expect_identical(1L,  TotalCopheneticIndex(Tree("((1,2),3,4,5);")))
  expect_identical(2L,  TotalCopheneticIndex(Tree("((1,2),(3,4),5);")))
  expect_identical(3L,  TotalCopheneticIndex(Tree("((1,2,3),4,5);")))
  expect_identical(4L,  TotalCopheneticIndex(Tree("(((1,2),3),4,5);")))
  expect_identical(4L,  TotalCopheneticIndex(Tree("((1,2,3),(4,5));")))
  expect_identical(5L,  TotalCopheneticIndex(Tree("(((1,2),3),(4,5));")))
  expect_identical(6L,  TotalCopheneticIndex(Tree("((1,2,3,4),5);")))
  expect_identical(7L,  TotalCopheneticIndex(Tree("(((1,2),3,4),5);")))
  expect_identical(8L,  TotalCopheneticIndex(Tree("(((1,2),(3,4)),5);")))
  expect_identical(9L,  TotalCopheneticIndex(Tree("(((1,2,3),4),5);")))
  expect_identical(10L, TotalCopheneticIndex(PectinateTree(5)))
  expect_identical(c(10L, 8L, 10L), TotalCopheneticIndex(as.phylo(0:2, 5)))

  expect_equal(TotalCopheneticIndex(BalancedTree(13)),
               TCIContext(13)[1, "minimum"])
  expect_equal(TotalCopheneticIndex(PectinateTree(13)),
               TCIContext(13)[1, "maximum"])

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = letters[1:8]),
    class = "phylo")
  expect_equal(28L, TotalCopheneticIndex(nasty))

  expect_equal(TCIContext(BalancedTree(5)), TCIContext(5L))
})

test_that("Expectations are correct", {
  nTip <- 5
  uniform <- DropTip(
    unlist(recursive = FALSE,
           lapply(
             lapply(as.phylo(1:NUnrooted(nTip) - 1, nTip),
                    AddTipEverywhere, "ROOT"),
             RootTree, "ROOT")
           ),
    "ROOT")
  tci <- TotalCopheneticIndex(uniform)
  context <- TCIContext(nTip)
  expect_equal(range(tci), range(context[1:2]))
  expect_equal(mean(tci), context[["uniform.expected"]])
  
  set.seed(1)
  unif <- replicate(100, TotalCopheneticIndex(RandomTree(10, root = TRUE)))
  expect_equal(
    mean(unif),
    TCIContext(10)$uniform.expected, # 76
    tolerance = 0.075
  )
})

test_that("Large trees work", {
  expect_equal(TotalCopheneticIndex(BalancedTree(164)), .MCI(164))
  H  <- function(n) sum(1 / seq_len(n))
  H2 <- function(n) sum(1 / (seq_len(n) ^ 2))
  
  x <- 1024
  
  expect_equal(unlist(TCIContext(x)),
               c(maximum = choose(x, 3), minimum = .MCI(x),
                 uniform.expected = exp(lchoose(x, 2) - log(2) +
                                          LnDoubleFactorial((x + x) - 2L) - 
                                          LnDoubleFactorial((x + x) - 3L)
                                        ) - choose(x, 2),
                 yule.expected =  (x * (x + 1)) - (2 * x * H(x)),
                 yule.variance = (
                   (1 / 12) * (x^4 - (10 * x^3) + (131 * x^2) - (2 * x))
                   ) - (4 * x^2 * H2(x)) - (6 * x * H(x))
               ))
})
