test_that("Write is successful", {
  Test <- function (tree) expect_equal(ape::write.tree(tree), as.Newick(tree))
  Test(BalancedTree(0:7))
  Test(UnrootTree(BalancedTree(0:7)))
  Test(PectinateTree(0:7))
  Test(UnrootTree(PectinateTree(0:7)))
  Test(ape::read.tree(text="(0,1,2,3,4);"))
  Test(ape::read.tree(text="((0,1,2,3,4), ((5, 6), (7, 8, 9)), 10, 11);"))
  Test(as.phylo(1:4, 5, 0:4))

  nasty <- structure(list(edge = structure(
    c(9, 12, 10, 13, 11, 10, 11, 13, 10, 13, 12, 9,
      5, 10,  1,  2,  3, 13,  9,  4, 11,  7,  8, 6),
    .Dim = c(12, 2)),
    Nnode = 5L,
    tip.label = 0:7),
    class = "phylo") # Danger: Do not plot!
  Test(nasty)
})

test_that("WriteTntCharacters()", {
  nTax <- 4L
  dataset <- matrix(seq_len(nTax), nrow = nTax, ncol = 6,
                    dimnames = list(LETTERS[1:nTax], 1:6))
  dataset[4, 4:6] <- "-"
  dataset[3, 3] <- "(23)"
  dataset[3:2, 4:6] <- "?"

  expect_equal(
    WriteTntCharacters(dataset, comment = c("COM", "MENT"),
                       pre = c("PRE", "FIX"),
                       post  = c("POST", "SCRIPT")),
    "PRE\nFIX\nxread 'COM MENT'\n6 4\nA 1 1 1 1 1 1\nB 2 2 2 ? ? ?\nC 3 3 [23] ? ? ?\nD 4 4 4 - - -\n;\nPOST\nSCRIPT"
    )
  expect_equal(WriteTntCharacters(dataset),
               WriteTntCharacters(MatrixToPhyDat(dataset)))

  expect_equal(
    WriteTntCharacters(dataset, types = c(num = 1, dna = 3)),
    "\nxread 'Dataset written by `TreeTools::WriteTntCharacters()`'\n6 4\n&[num]\n A 1 1\nB 2 2\nC 3 3\nD 4 4\n&[dna]\n A 1 1 1 1\nB 2 ? ? ?\nC [23] ? ? ?\nD 4 - - -\n;\n"
  )

  written <- tempfile()
  WriteTntCharacters(dataset, written)
  on.exit(file.remove(written))
  colnames(dataset) <- NULL
  expect_equal(ReadTntCharacters(written),
               gsub("(23)", "[23]", fixed = TRUE, dataset))
})
