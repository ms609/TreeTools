TestFile <- function(filename = "") {
  system.file("extdata", "tests", filename, package = "TreeTools")
}

test_that("File time is read correctly", {
  fileName <- TestFile("ape-tree.nex")
  expect_equal("2018-07-18 13:47:46", ApeTime(fileName, "string"))
  expect_error(ApeTime(rep(fileName, 2)))
})

test_that("Missing files fail gracefully", {
  expect_error(ReadAsPhyDat("non-existent.file"),
               "'non-existent.file' not found")
})

test_that("Nexus file can be parsed", {
  # Errors as lists:
  expect_equal("MATRIX block not found in Nexus file.",
               ReadCharacters(TestFile("ape-tree.nex"))[[1]])

  filename <- TestFile("parse-nexus.nexus")
  read <- ReadCharacters(filename)
  expect_equal(192, ncol(read))
  expect_equal(80, nrow(read))
  expect_equal("Wiwaxia", rownames(read)[4])
  expect_equal("(01)", as.character(read[1, 27]))

  filename <- TestFile("continuous.nex")
  read <- ReadCharacters(filename)
  expect_equal(1L, unique(as.integer(read[1, ])))
  expect_equal(setNames("?", "B_alienus"), read["B_alienus", 4])
  expect_equal(3L, unique(as.integer(read[3, ])))
})

test_that("NexusTokens() fails gracefully", {
  expect_error(NexusTokens("0123012301230123", integer(0)))
  expect_equal("Character number must be between 1 and 16.",
               NexusTokens("0123012301230123", 0)[[1]])
})

test_that("Matrix converts to phyDat", {
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,2,2,2,2,2,2,2,"?"),
                nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  expect_equal(mat, as.matrix(MatrixToPhyDat(mat)))
})

test_that(".PhyDatWithContrast() fails gracefully", {
  expect_error(.PhyDatWithContrast(matrix(0, 2, 2),
                                   matrix(c(1, 0, 0, 1), 2, 2, FALSE,
                                          list(0:1, 0:1))))
})

test_that("PhyDatToMatrix() with ambigs", {
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,2,2,2,2,"{12}","(01)","-","?"),
                nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  expectation <- mat
  expectation[3, 5:8] <- NA_character_
  expect_equal(PhyDatToMatrix(MatrixToPhyDat(mat), TRUE, TRUE), expectation)
  expect_equal(PhyDatToMatrix(MatrixToPhyDat(mat), TRUE, FALSE)[3, 5:8],
               c(NA, NA, "-", NA))
  expect_equal(PhyDatToMatrix(MatrixToPhyDat(mat), FALSE, TRUE)[3, 5:8],
               c("{12}", "{01}", NA_character_, "?"))
  expect_equal(
    PhyDatToMatrix(MatrixToPhyDat(mat),
                   ambigNA = FALSE,
                   inappNA = TRUE,
                   parentheses = NULL
                   )[3, 5:8],
    c("{12}", "(01)", NA_character_, "?"))
  expect_equal(
    PhyDatToMatrix(MatrixToPhyDat(mat),
                   ambigNA = FALSE,
                   inappNA = TRUE,
                   parentheses = c("<>"),
                   sep = "/"
    )[3, 5:8],
    c("<1/2>", "<0/1>", NA_character_, "?"))
})

test_that("Modified phyDat objects can be converted", {
  # Obtained by subsetting, i.e. dataset <- biggerDataset[1:4]
  dataset <- structure(list(a = c(1L, 1L, 1L, 1L), c = c(2L, 1L, 1L, 2L),
                            d = c(2L, 1L, 2L, 1L), e = c(2L, 2L, 2L, 2L)),
                       weight = c(3L, 3L, 0L, 0L), nr = 4L, nc = 2L,
                       index = c(2L, 1L, 1L, 1L, 2L, 2L),
                       .Label = c("0", "1"), allLevels = c("0", "1"),
                       type = "USER",
                       contrast = structure(c(1, 0, 0, 1), .Dim = c(2L, 2L),
                                            .Dimnames = list(NULL, c("0", "1"))),
                       class = "phyDat")
  expect_equal(c(4, 6), dim(PhyDatToMatrix(dataset)))
})

test_that("MatrixToPhyDat() warns when characters blank", {
  # May occur when loading an excel file with empty cells
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,"","","",0,1,0,1,2,2,2,2,2,2,2,"?"),
                nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  expect_warning(MatrixToPhyDat(mat))
})

test_that("MatrixToPhyDat() returns phyDat if passed", {
  expect_warning(expect_equal(Lobo.phy, MatrixToPhyDat(Lobo.phy)),
                 "phyDat")
})

test_that("StringToPhyDat()", {
  expect_equal(as.integer(StringToPhyDat("1111????", letters[1:8])),
               rep(1:2, each = 4))
  expect_equal(as.integer(StringToPhyDat("----????", letters[1:8])),
               rep(1:2, each = 4))
  expect_equal(as.integer(StringToPhyDat("----????")), rep(1:2, each = 4))
  expect_equal(names(StringToPhyDat("----????")), paste0("t", 1:8))
})

test_that("EndSentence() works correctly", {
  expect_equal(EndSentence("Hi"), "Hi.")
  expect_equal(EndSentence("Hi."), "Hi.")
  expect_equal(EndSentence("Hi?"), "Hi?")
  expect_equal(EndSentence("Hi!"), "Hi!")
  expect_equal(EndSentence(character(0)), character(0))
})

test_that("Unquote() unquotes", {
  expect_equal(Unquote("'Unquoted'"), "Unquoted")
  expect_equal(Unquote("\"Unquoted\""), "Unquoted")
  expect_equal(Unquote("'Unquoted '"), "Unquoted")
  expect_equal(Unquote("\" Unquoted \""), "Unquoted")
  expect_equal( Unquote("'Unquoted's '"), "Unquoted's")
  expect_equal(Unquote(.UnescapeQuotes("'Unquoted''s '")), "Unquoted's")
  expect_equal(Unquote("\"\""), "")
  expect_equal(Unquote("''"), "")
})

test_that("ReadNotes() reads notes", {
  notes <- ReadNotes(system.file("extdata/input/notes.nex",
                                 package = "TreeTools"))
  expect_equal(length(unlist(notes$`1`)), 0)
  expect_equal(notes[[2]][[2]], setNames("Taxon 2, char 2.", "taxon_b"))
  expect_equal(notes[[3]][[1]], "Three's a crowd.")
  expect_equal(notes[[3]][[2]], setNames("Tax1-Char3.", "taxon_a"))
})

test_that("ReadNotes() handles absence of character-taxon notes", {
  expect_equal(ReadNotes(system.file("extdata/tests/taxon-notes.nex",
                                     package = "TreeTools")),
               structure(list(), names = character(0)))
})

test_that("ReadNotes() handles misspecified encoding", {
  expect_message(
    expect_equal(ReadNotes(system.file("extdata/tests/encoding.nex",
                                     package = "TreeTools"))[[1]][[2]],
                 setNames("\u0080ncoding.", "Two")),
    "trying latin1 .*encoding")
})

test_that("ReadCharacters() reads CHARSTATELABELS", {
  labels <- ReadCharacters(system.file("extdata/input/dataset.nex",
                                       package = "TreeTools"))


  expect_equal(colnames(labels), c("Character one",
                                   "Character two",
                                   "lots-of-punctuation, and \"so on\"!",
                                   "Character n", "Character 5",  "Character 6",
                                   "final character"))

  ap <- c("absent", "present")
  expect_equal(attr(labels, "state.labels"),
               list(ap, ap,
                    c("here", "there", "everywhere"),
                    c("a long description", "present"),
                    c("simple", "more complex", "with (parentheses)",
                      "more complex, 6 still"),
                    c("this one has", "multiple lines"), ap));

  labels3 <- ReadCharacters(system.file("extdata/input/dataset.nex",
                                        package = "TreeTools"), 3)
  expect_equal(labels3, labels[, 3, drop = FALSE], ignore_attr = TRUE)
  expect_equal(attr(labels3, "state.labels"),
               attr(labels, "state.labels")[3])

})

test_that("MorphoBankDecode() decodes", {
  expect_equal("' -- x  \n 1--2", MorphoBankDecode("'' - x^n 1-2"))
})

test_that("NewickTree() works", {
  expect_equal("((Test taxon,Another test),(What's this?,Number 12.3));",
               NewickTree(BalancedTree(c("Test taxon", "Another_test",
                                             "What's this?", "Number 12.3"))))
})

test_that("as_newick() fails gracefully", {
  expect_equal(as_newick(matrix(0L, 0L, 2L)), ";")
  expect_equal(as_newick(matrix(1:0, 1L, 2L)), "(0);")
  expect_equal(as_newick(Postorder(BalancedTree(4)$edge) - 1L),
               as_newick(BalancedTree(4)$edge - 1L))
  expect_error(as_newick(matrix(0L, 8192 * 2L, 2L)),
               "Too many nodes")
  expect_error(as_newick(matrix(0L, 3, 3)),
               "`edge` must have two columns")
  expect_error(as_newick(matrix(c(4, 4, 4, 1:3), 3, 2)),
               "`min.edge.` must be zero")
  expect_error(as_newick(matrix(c(3, NA, 3, 0:2), 3, 2)),
               "`edge`.* NA")
  expect_error(as_newick(matrix(c(4, 4, 3, 0:2), 3, 2)),
               "`edge` is malformed")
})
