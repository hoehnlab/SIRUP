# Tests for CIGAR + alignment helpers:
#   - find_longest_matches(cigars)
#   - analyze_cigar_patterns(cigars)
#   - has_minimum_match_length(cigar, min_length)
#   - cmv_consensusMatrix(sequences)
#   - cigar_onehot_triplets(cigars, aln_starts, aln_ends, region_start, region_end)
#   - edlib_search(queries, targets, max_distance, num_threads)
#
# CIGAR convention used by these functions: `=` denotes an exact match and
# is what gets counted by find_longest_matches / has_minimum_match_length /
# analyze_cigar_patterns. `M` (the generic align operator) is NOT counted —
# the helpers were written for the `=`/X-aware long-read CIGAR output.

# ---- find_longest_matches ------------------------------------------------

test_that("find_longest_matches returns the longest `=` run per CIGAR", {
  expect_equal(find_longest_matches(c("100=")),                 100)
  expect_equal(find_longest_matches(c("50=10X40=")),             50)   # the 50
  expect_equal(find_longest_matches(c("10=5X10=")),              10)   # tied 10s
  expect_equal(find_longest_matches(c("100=", "50=10X40=")),     c(100, 50))
})

test_that("find_longest_matches returns 0 when no `=` is present", {
  expect_equal(find_longest_matches(c("100M")),                  0)
  expect_equal(find_longest_matches(c("")),                      0)
  expect_equal(find_longest_matches(c("50M10I40M")),             0)
})

# ---- has_minimum_match_length --------------------------------------------

test_that("has_minimum_match_length compares the longest `=` run to a threshold", {
  expect_true (has_minimum_match_length("100=",     50))
  expect_true (has_minimum_match_length("100=",    100))   # boundary
  expect_false(has_minimum_match_length("100=",    101))
  expect_true (has_minimum_match_length("10=5X10=", 10))
  expect_false(has_minimum_match_length("10=5X10=", 15))
  expect_false(has_minimum_match_length("100M",     10))   # M is not counted
  expect_false(has_minimum_match_length("",          1))
})

# ---- analyze_cigar_patterns ---------------------------------------------

test_that("analyze_cigar_patterns reports longest / total / block count per CIGAR", {
  out <- analyze_cigar_patterns(c("100=10X50=", "30M", "5=5X5="))
  # Note the slightly inconsistent naming — singular `longest_match`,
  # plural `total_matches` and `match_blocks`. Locking in the names that
  # the C++ side actually returns.
  expect_named(out, c("longest_match", "total_matches", "match_blocks"))
  expect_equal(out$longest_match, c(100L,  0L,  5L))
  expect_equal(out$total_matches, c(150L,  0L, 10L))   # 100+50, 0, 5+5
  expect_equal(out$match_blocks,  c(  2L,  0L,  2L))   # two = blocks, zero, two
})

# ---- cmv_consensusMatrix -------------------------------------------------

test_that("cmv_consensusMatrix builds a 4 x ncol(seq) base-count matrix", {
  m <- cmv_consensusMatrix(c("ACGT", "ACGT", "ACAT"))
  expect_equal(rownames(m), c("A", "C", "G", "T"))
  expect_equal(dim(m), c(4L, 4L))
  expect_equal(m["A", ], c(3, 0, 1, 0))   # position 3: 2x G + 1x A
  expect_equal(m["C", ], c(0, 3, 0, 0))
  expect_equal(m["G", ], c(0, 0, 2, 0))
  expect_equal(m["T", ], c(0, 0, 0, 3))
})

test_that("cmv_consensusMatrix accepts mixed case bases", {
  # Probe: confirm case-insensitive base counting works as the C++ switch
  # suggests (covers both upper- and lowercase A/C/G/T).
  m <- cmv_consensusMatrix(c("acgt", "ACGT"))
  expect_equal(sum(m), 8)                          # 4 bases x 2 reads
  expect_equal(unname(m["A", 1]), 2)
  expect_equal(unname(m["T", 4]), 2)
})

test_that("cmv_consensusMatrix on empty input returns an empty matrix", {
  m <- cmv_consensusMatrix(character(0))
  expect_equal(nrow(m), 0)
})

# ---- cigar_onehot_triplets ----------------------------------------------

test_that("cigar_onehot_triplets emits one triplet per matched position", {
  # 5M alignment starting at ref position 10 -> covers 10..14
  out <- cigar_onehot_triplets(
    cigars       = c("5M"),
    aln_starts   = c(10L),
    aln_ends     = c(14L),
    region_start = 1L,
    region_end   = 100L
  )
  expect_named(out, c("i", "j", "nrow", "ncol"))
  expect_equal(out$i, rep(1L, 5))                # one read -> row 1
  expect_equal(out$j, 10:14)                     # column = ref position
  expect_equal(out$nrow, 1L)
  expect_equal(out$ncol, 99L)
})

test_that("cigar_onehot_triplets handles multiple reads at distinct positions", {
  out <- cigar_onehot_triplets(
    cigars       = c("3M",   "2M"),
    aln_starts   = c(10L,    20L),
    aln_ends     = c(12L,    21L),
    region_start = 1L,
    region_end   = 100L
  )
  expect_equal(out$i, c(1L, 1L, 1L, 2L, 2L))
  expect_equal(out$j, c(10L, 11L, 12L, 20L, 21L))
})

# ---- edlib_search --------------------------------------------------------

test_that("edlib_search returns distance 0 for identical query/target", {
  res <- suppressMessages(
    edlib_search(queries = "ACGTACGT", targets = "ACGTACGT", max_distance = 5)
  )
  expect_named(res, "ACGTACGT")
  expect_equal(res[["ACGTACGT"]]$target,   "ACGTACGT")
  expect_equal(res[["ACGTACGT"]]$distance, 0L)
})

test_that("edlib_search reports the edit distance and respects max_distance", {
  res <- suppressMessages(
    edlib_search(queries = "ACGT", targets = c("ACGT", "ACGA", "AAAA"),
                 max_distance = 5)
  )
  hits <- res[["ACGT"]]
  expect_setequal(hits$target, c("ACGT", "ACGA", "AAAA"))
  expect_equal(hits$distance[hits$target == "ACGT"], 0L)
  expect_equal(hits$distance[hits$target == "ACGA"], 1L)
  expect_gt   (hits$distance[hits$target == "AAAA"], 0L)

  # max_distance = 0 keeps only the exact-match target
  res0 <- suppressMessages(
    edlib_search(queries = "ACGT", targets = c("ACGT", "ACGA", "AAAA"),
                 max_distance = 0)
  )
  hits0 <- res0[["ACGT"]]
  expect_equal(hits0$target,   "ACGT")
  expect_equal(hits0$distance, 0L)
})

test_that("edlib_search handles multiple queries", {
  res <- suppressMessages(
    edlib_search(queries = c("ACGT", "TTTT"),
                 targets = c("ACGT", "TTTT", "AAAA"),
                 max_distance = 5)
  )
  expect_named(res, c("ACGT", "TTTT"))
  expect_equal(res[["ACGT"]]$distance[res[["ACGT"]]$target == "ACGT"], 0L)
  expect_equal(res[["TTTT"]]$distance[res[["TTTT"]]$target == "TTTT"], 0L)
})
