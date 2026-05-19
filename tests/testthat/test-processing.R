test_that("gap length helpers work", {
  expect_equal(maxGapLength("AA...TT..G"), 3L)
  expect_equal(openingGapLength("...A.."), 3L)
  expect_equal(openingGapLength("A.."), 0L)
})

test_that("filterReceptorReads filters by locus and lengths", {
  airr <- data.frame(
    locus = c("IGH", "IGK", "IGH"),
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(120, 120, 50),
    j_sequence_start = c(1, 1, 1),
    j_sequence_end = c(20, 20, 20),
    junction_length = c(12, 12, 12),
    productive = c(TRUE, TRUE, TRUE)
  )
  out <- filterReceptorReads(airr, loci = "IGH", productive_only = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$locus, "IGH")
})
