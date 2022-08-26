test_that("UnshiftTree() works", {
  t1 <- as.phylo(1, 8)
  t2..9 <- setNames(as.phylo(2:3, 8), letters[2:3])
  t1..9 <- setNames(as.phylo(1:3, 8), c('', letters[2:3]))
  attr(t1..9, 'tip.label') <- NULL
  expect_true(all.equal(UnshiftTree(t1, t2..9), t1..9))
  expect_equal(unclass(t1..9), UnshiftTree(t1, unclass(t2..9)))
  expectation <- as.phylo(1:2, 8)
  attr(expectation, 'tip.label') <- NULL
  expect_equal(expectation, UnshiftTree(t1, as.phylo(2, 8)))
})


test_that("SpectrumLegend()", {
  skip_if(packageVersion("graphics") < "4.1")
  skip_if(packageVersion("vdiffr") < "1.0")
  vdiffr::expect_doppelganger('SpectrumLegend', function() {
    
    # Set up blank plot
    plot(0:1, 10:11, asp = 1, type = "n", frame.plot = FALSE,
         xlab = "x", ylab = "y")
    SpectrumLegend(legend = c("Bottom", "Middle", "Top"),
                   lwd = 5,
                   palette = rev(hcl.colors(256L, "RdYlBu")),
                   text.col = c("blue", "brown", "red"),
                   title = "Blue title", title.font = 3, title.cex = 0.5)
    SpectrumLegend(0.4, 10.6, 0.2, 11, abs = TRUE,
                   legend = seq(0, 10, by = 2), palette = 0:10,
                   lty = "dotted", pos = 2,
                   lmitre = 2, lend = "round", ljoin = "round",
                   title = "Default title", xpd = NA)
  })
})
