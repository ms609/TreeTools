TestFile <- function(filename = "") {
  system.file("extdata", "tests", filename, package = "TreeTools")
}

test_that("ReadTntCharacter()", {
  testFile <- TestFile("tnt-trees-and-matrix.tnt")
  expect_equal(
    structure(c("0", "0", "1", "1", "1", "1", "-", "-", "0", "0",
                "0", "0", "-", "-", "0", "0", "1", "1", "-", "-", "-", "-", "0",
                "1", "-", "-", "1", "1", "1", "0", "-", "-", "0", "0", "1", "-",
                "-", "-", "0", "1", "0", "-"), .Dim = 6:7,
              .Dimnames = list(c("a", "b", "c", "d", "e", "f"), NULL)),
    ReadTntCharacters(testFile)
  )
  
  dnaTest <- TestFile("tnt-dna.tnt")
  expect_equal(ReadTntCharacters(dnaTest),
               cbind(ReadTntCharacters(dnaTest, type = "num"),
                     ReadTntCharacters(dnaTest, type = "dna")))
  expect_equal(ReadTntCharacters(dnaTest),
               ReadTntCharacters(dnaTest, type = c("NUM", "Dna")))
  expect_message(expect_null(ReadTntCharacters(dnaTest, type = "NONE")))
  
  skip_if_not_installed("phangorn")
  expect_equal(
    phangorn::as.phyDat(
      ReadTntCharacters(testFile),
      type = "USER",
      contrast = structure(c(0, 0, 1, 1, 0, 0, 0, 1, 0), .Dim = c(3L, 3L),
                           .Dimnames = list(c("0", "1", "-"), c("-", "0", "1")))
    ),
    ReadTntAsPhyDat(testFile))
})

test_that("TntTextToTree()", {
  expect_equal(TNTText2Tree("(A (B (C (D E ))));"),
               ape::read.tree(text = "(A, (B, (C, (D, E))));"))
})

test_that("ReadTntTree() NULL return", {
  expect_null(ReadTntTree(TestFile("ape-tree.nex")))
})

test_that("TNT trees parsed correctly", {
  expect_warning(expect_warning(
    trees <- ReadTntTree(
      filepath = TestFile("tnt-tree.tre"),
      relativePath = TestFile()
    ),
    "Multiple tree blocks"),
    "Expected `ttags \\*N`; applying tags to first tree"
  )
  expect_equal(length(trees), 3)
  expect_equal(trees[[3]]$tip.label[2], "Novocrania")
  expect_equal(trees[[3]][["node.label"]][37:41 - NTip(trees[[3]])],
               c("Nov+1", "Nov+2", "Nov+3",
                 "Nov +4", # Concatenated from Nov & +4
                 "Nov+5"))
  expect_equal(trees[[3]][["node.label"]][69 - NTip(trees[[3]])], "Nis+1")
  trees[[2]][["node.label"]] <- NULL
  trees[[3]][["node.label"]] <- NULL
  expect_equal(trees[[2]], trees[[3]])
  expect_equal(ConsensusWithout(trees, "Paterimitra")$Nnode, 32)
  
  tipLabels <- c("Dailyatia", "Novocrania", "Craniops", "Ussunia", "Gasconsia",
                 "Heliomedusa_orienta", "Micrina", "Mickwitzia_muralensis",
                 "Micromitra", "Askepasma_toddense", "Pelagodiscus_atlanticus",
                 "Lingula", "Eoobolus", "Clupeafumosus_socialis", "Phoronis",
                 "Eccentrotheca", "Yuganotheca_elegans",
                 "Longtancunella_chengjiangensis", "Paterimitra",
                 "Lingulellotreta_malongensis", "Acanthotretella",
                 "Lingulosacculus", "Pedunculotheca_diania",
                 "Haplophrentis_carinatus", "Tomteluva_perturbata",
                 "Salanygolina", "Mummpikia_nuda", "Alisina", "Coolinia_pecten",
                 "Antigonambonites_planus", "Kutorgina_chengjiangensis",
                 "Nisusia_sulcata", "Glyptoria", "Orthis", "Terebratulina")
  fromLabels <- expect_warning(
    ReadTntTree(TestFile("tnt-tree.tre"), tipLabels = tipLabels),
    "Expected `ttags .N`;")
  expect_identical(trees, fromLabels)
  
  namedLabels <- ReadTntTree(TestFile("tnt-namedtree.tre"))[[1]]$tip.label
  expect_equal("Flustra_sp.", namedLabels[1])
  expect_equal(74L, length(namedLabels))
  
  tam <- ReadTntTree(TestFile("tnt-trees-and-matrix.tnt"))
  expect_equal(length(tam), 4L)
  expect_equal(tam[[1]],
               TNTOrder(ape::read.tree(text = "(a, (b, (c, (f, (d, e )))));")))
  expect_equal(tam[[4]],
               TntOrder(ape::read.tree(text = "(a, (b, (c, (e, d, f))));")))
  
  
  oldWD <- getwd()
  setwd(system.file(package = "TreeTools"))
  expect_equal(paste0("taxon_", letters[1:5]),
               ReadTntTree("extdata/output/numbered.tre", "./extdata", 2)[[1]]$tip.label)
  
  expect_equal(paste0("taxon_", letters[1:5]),
               ReadTntTree("extdata/output/named.tre")$tip.label)
  setwd(oldWD)
})

test_that("ReadTntTree() reads bare tree", {
  bareTree <- TestFile("tnt-bare-tree.tnt")
  expect_warning(
    expect_equal(ReadTntTree(bareTree),
                 ReadTntTree(bareTree, tipLabels = as.character(0:5))),
    "does not link to taxon names")
})

test_that("ReadTntTree() follows TNT node numbering conventions", {
  a..e <- paste("taxon", letters[1:5], sep = "_")
  namedTree <- system.file("extdata/output/named.tre", package = "TreeTools")
  expect_equal(
    ReadTntTree(namedTree),
    ReadTntTree(namedTree, tipLabels = a..e)
  )
  expect_equal(ReadTntTree(namedTree)$edge,
               Postorder(ReadTntTree(namedTree))$edge)
  
  expect_equal(
    # NB this is not a perfect test, as the order in which edges are listed
    # is not stipulated by the ReadTntTree documentation.
    ReadTntTree(namedTree, tipLabels = rev(a..e))$edge,
    cbind(rep(c(7, 8, 9, 6), each = 2), c(2, 1, 7, 3, 8, 4, 9, 5))
  )
})
