context("parse_files.R")

TestFile <- function (filename = '') {
  system.file('extdata', 'tests', filename, package = 'TreeTools')
}

test_that("File time is read correctly", {
  fileName <- TestFile('ape-tree.nex')
  expect_equal('2018-07-18 13:47:46', ApeTime(fileName, 'string'))
  expect_error(ApeTime(rep(fileName, 2)))
})

test_that("Nexus file can be parsed", {
  # Errors as lists:
  expect_equal("MATRIX block not found in Nexus file.",
               ReadCharacters(TestFile('ape-tree.nex'))[[1]])

  filename <- TestFile('parse-nexus.nexus')
  read <- ReadCharacters(filename)
  expect_equal(192, ncol(read))
  expect_equal(80, nrow(read))
  expect_equal("Wiwaxia", rownames(read)[4])
  expect_equal("(01)", as.character(read[1, 27]))

  filename <- TestFile('continuous.nex')
  read <- ReadCharacters(filename)
  expect_equal(1L, unique(as.integer(read[1, ])))
  expect_equivalent('?', read['B_alienus', 4])
  expect_equal(3L, unique(as.integer(read[3, ])))
})

test_that("ReadTntCharacter()", {
  testFile <- TestFile('tnt-trees-and-matrix.tnt')
  expect_equal(
    structure(c("0", "0", "1", "1", "1", "1", "-", "-", "0", "0",
              "0", "0", "-", "-", "0", "0", "1", "1", "-", "-", "-", "-", "0",
              "1", "-", "-", "1", "1", "1", "0", "-", "-", "0", "0", "1", "-",
              "-", "-", "0", "1", "0", "-"), .Dim = 6:7,
            .Dimnames = list(c("a", "b", "c", "d", "e", "f"), NULL)),
    ReadTntCharacters(testFile)
  )
  expect_equal(
    phangorn::as.phyDat(
      ReadTntCharacters(testFile),
      type = 'USER',
      contrast = structure(c(0, 0, 1, 1, 0, 0, 0, 1, 0), .Dim = c(3L, 3L),
                           .Dimnames = list(c("0", "1", "-"), c("-", "0", "1")))
      ),
    ReadTntAsPhyDat(testFile))

  dnaTest <- TestFile('tnt-dna.tnt')
  expect_equal(ReadTntCharacters(dnaTest),
               cbind(ReadTntCharacters(dnaTest, type = 'num'),
                     ReadTntCharacters(dnaTest, type = 'dna')))
  expect_equal(ReadTntCharacters(dnaTest),
               ReadTntCharacters(dnaTest, type = c('NUM', 'Dna')))
  expect_null(ReadTntCharacters(dnaTest, type = 'NONE'))
})

test_that("TNT trees parsed correctly", {
  trees <- ReadTntTree(TestFile('tnt-tree.tre'), relativePath = TestFile())
  expect_equal(2, length(trees))
  expect_equal(32, ConsensusWithout(trees, 'Paterimitra')$Nnode)

  tipLabels <- c('Dailyatia', 'Novocrania', 'Craniops', 'Ussunia', 'Gasconsia',
                 'Heliomedusa_orienta', 'Micrina', 'Mickwitzia_muralensis',
                 'Micromitra', 'Askepasma_toddense', 'Pelagodiscus_atlanticus',
                 'Lingula', 'Eoobolus', 'Clupeafumosus_socialis', 'Phoronis',
                 'Eccentrotheca', 'Yuganotheca_elegans',
                 'Longtancunella_chengjiangensis', 'Paterimitra',
                 'Lingulellotreta_malongensis', 'Acanthotretella',
                 'Lingulosacculus', 'Pedunculotheca_diania',
                 'Haplophrentis_carinatus', 'Tomteluva_perturbata',
                 'Salanygolina', 'Mummpikia_nuda', 'Alisina', 'Coolinia_pecten',
                 'Antigonambonites_planus', 'Kutorgina_chengjiangensis',
                 'Nisusia_sulcata', 'Glyptoria', 'Orthis', 'Terebratulina')
  fromLabels <- ReadTntTree(TestFile('tnt-tree.tre'), tipLabels = tipLabels)
  expect_identical(trees, fromLabels)

  namedLabels <- ReadTntTree(TestFile('tnt-namedtree.tre'))[[1]]$tip.label
  expect_equal('Flustra_sp.', namedLabels[1])
  expect_equal(74L, length(namedLabels))

  tam <- ReadTntTree(TestFile('tnt-trees-and-matrix.tnt'))
  expect_equal(3, length(tam))
  expect_equal(ape::read.tree(text = '(a, (b, (c, (f, (d, e )))));'), tam[[1]])


  oldWD <- getwd()
  setwd(system.file(package = 'TreeTools'))
  expect_equal(paste0('taxon_', letters[1:5]),
    ReadTntTree('extdata/output/numbered.tre', './extdata', 2)[[1]]$tip.label)

  expect_equal(paste0('taxon_', letters[1:5]),
    ReadTntTree('extdata/output/named.tre')$tip.label)
  setwd(oldWD)
})

test_that("NexusTokens() fails gracefully", {
  expect_error(NexusTokens("0123012301230123", integer(0)))
  expect_equal("Character number must be between 1 and 16.",
               NexusTokens("0123012301230123", 0)[[1]])
})

test_that("Matrix converts to phyDat", {
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,2,2,2,2,2,2,2,'?'),
                nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  expect_equal(mat, PhyDatToMatrix(MatrixToPhyDat(mat)))
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
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,'','','',0,1,0,1,2,2,2,2,2,2,2,'?'),
                nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  expect_warning(MatrixToPhyDat(mat))
})

test_that("StringToPhyDat()", {
  expect_equal(rep(1:2, each = 4),
               as.integer(StringToPhyDat('1111????', letters[1:8])))
  expect_equal(rep(1:2, each = 4),
               as.integer(StringToPhyDat('----????', letters[1:8])))
})

test_that('PhyToString() works', {
  longLevels <- phyDat(rbind(x = c('-', '?', 0:12), y = c(12:0, '-', '?')),
                       type = 'USER', levels = c(0:6, '-', 7:12))
  expect_equal("-?0123456789ABCCBA9876543210-?", PhyToString(longLevels))

  # Two -s → error
  attr(longLevels, 'allLevels')[1] <- '-'
  expect_error(PhyToString(longLevels))

  # 10 → 1
  longLevels <- phyDat(rbind(x = c('-', '?', 1:10), y = c(10:1, '-', '?')),
                       type = 'USER', levels = c(1:6, '-', 7:10))
  expect_equal("-?12345678900987654321-?", PhyToString(longLevels))

  phy <- StringToPhyDat('012[01]', letters[1:4])
  expect_equal('012{01}', PhyToString(phy))
  expect_equal('012<01>', PhyToString(phy, parentheses = '<'))
  expect_equal('012<01>', PhyToString(phy, parentheses = '>'))
  expect_equal('012(01)', PhyToString(phy, parentheses = '('))
  expect_equal('012(01)', PhyToString(phy, parentheses = ')'))
  expect_equal('012[01]', PhyToString(phy, parentheses = ']'))
  expect_equal('012[01]', PhyToString(phy, parentheses = '['))
  expect_equal('012{01}', PhyToString(phy, parentheses = '}'))
  expect_equal('012{01}', PhyToString(phy, parentheses = '{'))
  expect_equal('012{01}', PhyToString(phy, parentheses = '!'))

  str <- '012{01}0123'
  phy <- StringToPhyDat(str, letters[1:4])
  expect_equal(str, PhyToString(StringToPhyDat(str, letters[1:4])))
  expect_equal(str,
               PhyToString(StringToPhyDat(str, letters[1:4], byTaxon = TRUE),
                           byTaxon = TRUE))
})

test_that("EndSentence() works correctly", {
  expect_equal('Hi.', EndSentence('Hi'))
  expect_equal('Hi.', EndSentence('Hi.'))
  expect_equal('Hi?', EndSentence('Hi?'))
  expect_equal('Hi!', EndSentence('Hi!'))
})

test_that("Unquote() unquotes", {
  expect_equal("Unquoted", Unquote("'Unquoted'"))
  expect_equal("Unquoted", Unquote('"Unquoted"'))
  expect_equal("Unquoted", Unquote("'Unquoted '"))
  expect_equal("Unquoted", Unquote('" Unquoted "'))
  expect_equal("Unquoted's", Unquote("'Unquoted's '"))
  expect_equal("", Unquote('""'))
  expect_equal("", Unquote("''"))
})

test_that("MorphoBankDecode() decodes", {
  expect_equal("' -- x  \n 1--2", MorphoBankDecode("'' - x^n 1-2"))
})

test_that('NewickTree() works', {
  expect_equal("((Test taxon,Another test),(What's this?,Number 12.3));",
               NewickTree(BalancedTree(c('Test taxon', 'Another_test',
                                             "What's this?", "Number 12.3"))))
})

test_that('as_newick fails gracefully', {
  expect_error(as_newick(matrix(0L, 8192 * 2L, 2L)))
})
