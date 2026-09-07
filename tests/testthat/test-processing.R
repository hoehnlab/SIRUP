test_that("gap length helpers work", {
  expect_equal(.sirap_max_gap_length("AA...TT..G"), 3L)
  expect_equal(.sirap_opening_gap_length("...A.."), 3L)
  expect_equal(.sirap_opening_gap_length("A.."), 0L)
})

test_that("filterAirr filters by locus and lengths", {
  airr <- data.frame(
    locus = c("IGH", "IGK", "IGH"),
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(120, 120, 50),
    j_sequence_start = c(1, 1, 1),
    j_sequence_end = c(20, 20, 20),
    junction_length = c(12, 12, 12),
    productive = c(TRUE, TRUE, TRUE)
  )
  out <- filterAirr(airr, locus_val = "IGH", productive_only = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$locus, "IGH")
})
