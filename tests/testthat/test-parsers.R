context("parse_files.R")
TestFile <- function (filename = '') {
  paste0(system.file(package='TreeTools'), '/extdata/tests/', filename)
}
test_that("File time is read correctly", {
  fileName <- TestFile('ape-tree.nex')
  expect_equal('2018-07-18 13:47:46', ApeTime(fileName, 'string'))
})

test_that("Nexus file can be parsed", {
  filename <- TestFile('parse-nexus.nexus')
  read <- ReadCharacters(filename)
  expect_equal(192, ncol(read))
  expect_equal(80, nrow(read))
  expect_equal("Wiwaxia", rownames(read)[4])
  expect_equal("(01)", as.character(read[1, 27]))
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
  expect_equal('Flustra', namedLabels[1])
  expect_equal(74L, length(namedLabels))

  oldWD <- getwd()
  setwd(system.file(package = 'TreeTools'))
  expect_equal(paste0('taxon_', letters[1:5]),
    ReadTntTree('extdata/output/numbered.tre', './extdata', 2)[[1]]$tip.label)

  expect_equal(paste0('taxon_', letters[1:5]),
    ReadTntTree('extdata/output/named.tre')$tip.label)
  setwd(oldWD)
})

test_that("Matrix converts to phyDat", {
  mat <- matrix(c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,2,2,2,2,2,2,2,'?'),
                nrow=3, byrow=TRUE)
  rownames(mat) <- LETTERS[1:3]
  expect_equal(mat, PhyDatToMatrix(MatrixToPhyDat(mat)))
})

test_that('PhyToString works', {
  longLevels <- phyDat(rbind(x = c('-', '?', 0:12), y = c(12:0, '-', '?')),
                       type='USER', levels=c(0:6, '-', 7:12))
  expect_equal("-?0123456789ABCCBA9876543210-?", PhyToString(longLevels))
})

test_that('NewickTree() works', {
  expect_equal('((Test taxon,Another test),(What`sthis?,Number12.3));',
               NewickTree(BalancedTree(c('Test taxon', 'Another_test',
                                             "What`s this?", "Number 12.3"))))
})

test_that('as_newick fails gracefully', {
  expect_error(as_newick(matrix(0L, 8192 * 2L, 2L)))
})
