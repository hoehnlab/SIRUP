# Tests for the PAF -> splice contact matrix pipeline:
#   - read_paf()                        (R-level PAF reader; data.table out)
#   - stream_paf2mat()                  (file path -> sparse matrix)
#   - stream_paf2mat_anch()             (with anchor-overlap filter)
#   - splice_contacts_from_paf()        (synonym of stream_paf2mat)
#   - splice_contacts_from_paf_dt()     (in-memory data.frame -> matrix)
#
# Fixture rationale:
#
# PAF columns (tab-separated):
#   qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend,
#   matches, alnlen, mapq, [optional tags... cg:Z:CIGAR is the one we need]
#
# tstart in PAF is 0-based; rad converts internally to 1-based
# (ref_start1b = tstart + 1). So a read with PAF tstart = 100 and
# CIGAR `100M50N100M` has aligned blocks:
#   - block 1: ref [101, 200] (100 bases)
#   - N skip:  ref [201, 250] (50 bases — the intron)
#   - block 2: ref [251, 350] (100 bases)
# The splice junction is the pair (donor=200, acceptor=251).
#
# With bin_bp = 100, region_start = 1:
#   - bin 0: [  1, 100]
#   - bin 1: [101, 200]   <- donor lands here
#   - bin 2: [201, 300]   <- acceptor lands here
#   - bin 3: [301, 400]
# So we expect a nonzero entry at (bin 1, bin 2) (0-indexed) =
# (row 2, col 3) (1-indexed in R matrix).

make_fixture_paf <- function(rows = NULL) {
  paf <- tempfile(fileext = ".paf")
  if (is.null(rows)) {
    rows <- c(
      # junction at 200/251 -> (bin 1, bin 2)
      paste("read1", 200, 0, 200, "+", "chr1", 10000, 100, 300,
            200, 200, 60, "cg:Z:100M50N100M", sep = "\t"),
      # junction at 600/651 -> (bin 5, bin 6)
      paste("read2", 200, 0, 200, "+", "chr1", 10000, 500, 700,
            200, 200, 60, "cg:Z:100M50N100M", sep = "\t"),
      # duplicate of read1 -> doubles the (bin 1, bin 2) count
      paste("read3", 200, 0, 200, "+", "chr1", 10000, 100, 300,
            200, 200, 60, "cg:Z:100M50N100M", sep = "\t")
    )
  }
  writeLines(rows, paf)
  paf
}

# ---- read_paf ------------------------------------------------------------

test_that("read_paf parses standard PAF columns + cg tag", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  df <- read_paf(paf)
  # read_paf assigns: id, ql, qs, qe, qr, sn, sl, ss, se, mm, mt, mq + tag cols
  for (col in c("id","ql","qs","qe","qr","sn","sl","ss","se","mm","mt","mq")) {
    expect_true(col %in% names(df), info = paste("missing column:", col))
  }
  expect_equal(df$id, c("read1", "read2", "read3"))
  expect_equal(df$sn, c("chr1",  "chr1",  "chr1"))
  expect_equal(df$ss, c(100,     500,     100))
  expect_equal(df$mq, c(60,      60,      60))
  # cg tag column (named "cg" by read_paf based on prefix detection);
  # read_paf strips the leading "cg:Z:" before storing.
  expect_true("cg" %in% names(df))
  expect_equal(df$cg, c("100M50N100M", "100M50N100M", "100M50N100M"))
})

# ---- stream_paf2mat ------------------------------------------------------

test_that("stream_paf2mat returns a symmetric sparse matrix with splice contacts at the expected bins", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  mat <- stream_paf2mat(
    path         = paf,
    target       = "chr1",
    region_start = 1,
    region_end   = 1000,
    bin_bp       = 100,
    min_mapq     = 0,
    symmetric    = TRUE
  )

  expect_s4_class(mat, "dsCMatrix")
  expect_equal(nrow(mat), ncol(mat))
  expect_gt(sum(mat), 0)

  # Junction at ref 200/251 with bin_bp=100, region_start=1 lands in (bin1, bin2)
  # (0-indexed) = (row 2, col 3) in 1-indexed R.
  # Junction at 600/651 lands in (bin5, bin6) = (row 6, col 7).
  # read1 + read3 contribute to the first; read2 to the second.
  m <- as.matrix(mat)
  expect_equal(m[2, 3], 2)
  expect_equal(m[3, 2], 2)              # symmetric mirror
  expect_equal(m[6, 7], 1)
  expect_equal(m[7, 6], 1)
  # Off-junction cells should be zero
  expect_equal(m[1, 1], 0)
  expect_equal(m[2, 5], 0)
})

test_that("symmetric arg controls the matrix class, not whether content is mirrored", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  m_sym  <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  m_gen  <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = FALSE)
  expect_s4_class(m_sym, "dsCMatrix")
  expect_s4_class(m_gen, "dgCMatrix")
  # Both store the junction in both (i,j) and (j,i).
  expect_equal(as.matrix(m_sym), as.matrix(m_gen))
})

test_that("stream_paf2mat filters by min_mapq", {
  rows <- c(
    paste("read1", 200, 0, 200, "+", "chr1", 10000, 100, 300,
          200, 200,  5, "cg:Z:100M50N100M", sep = "\t"),   # mapq 5 (low)
    paste("read2", 200, 0, 200, "+", "chr1", 10000, 500, 700,
          200, 200, 60, "cg:Z:100M50N100M", sep = "\t")    # mapq 60
  )
  paf <- make_fixture_paf(rows); on.exit(unlink(paf), add = TRUE)

  mat_all <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, min_mapq = 0,  symmetric = TRUE)
  mat_q30 <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, min_mapq = 30, symmetric = TRUE)
  expect_equal(sum(mat_all), 4)         # 2 entries (symmetric mirror) × 2 junctions = 4
  expect_equal(sum(mat_q30), 2)         # only the mapq=60 junction survives, mirrored
})

test_that("stream_paf2mat filters by target chromosome", {
  rows <- c(
    paste("read1", 200, 0, 200, "+", "chr1", 10000, 100, 300,
          200, 200, 60, "cg:Z:100M50N100M", sep = "\t"),
    paste("read2", 200, 0, 200, "+", "chr2", 10000, 100, 300,
          200, 200, 60, "cg:Z:100M50N100M", sep = "\t")
  )
  paf <- make_fixture_paf(rows); on.exit(unlink(paf), add = TRUE)

  m1 <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  m2 <- stream_paf2mat(paf, "chr2", 1, 1000, bin_bp = 100, symmetric = TRUE)
  expect_equal(sum(m1), 2)              # only read1 contributes, mirrored
  expect_equal(sum(m2), 2)              # only read2 contributes, mirrored
})

test_that("stream_paf2mat with no junctions on target produces an empty matrix or errors clearly", {
  # All-on-chr1 fixture, but ask for chr_nonexistent
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  expect_error(
    stream_paf2mat(paf, "chr_nonexistent", 1, 1000, bin_bp = 100, symmetric = TRUE),
    "No edges"
  )
})

# ---- splice_contacts_from_paf == stream_paf2mat --------------------------

test_that("splice_contacts_from_paf produces the same matrix as stream_paf2mat", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  m1 <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  m2 <- splice_contacts_from_paf(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  expect_equal(as.matrix(m1), as.matrix(m2))
})

# ---- splice_contacts_from_paf_dt -----------------------------------------

test_that("splice_contacts_from_paf_dt produces the same matrix as the file-based version", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  m_file <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)

  # Equivalent in-memory table. splice_contacts_from_paf_dt's defaults expect
  # columns tname/tstart/tend/cg/mapq.
  paf_dt <- data.table::data.table(
    tname  = c("chr1", "chr1", "chr1"),
    tstart = c(100,    500,    100),
    tend   = c(300,    700,    300),
    mapq   = c(60,     60,     60),
    cg     = c("cg:Z:100M50N100M",
               "cg:Z:100M50N100M",
               "cg:Z:100M50N100M")
  )
  m_dt <- splice_contacts_from_paf_dt(
    paf_dt, target = "chr1",
    region_start = 1, region_end = 1000,
    bin_bp = 100, symmetric = TRUE
  )
  expect_equal(as.matrix(m_file), as.matrix(m_dt))
})

test_that("splice_contacts_from_paf_dt honours custom column name mappings", {
  # Use read_paf's native column names (sn/ss/cg/mq) via the col_* args.
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  paf_dt <- read_paf(paf)

  m <- splice_contacts_from_paf_dt(
    paf_dt,
    target = "chr1",
    region_start = 1, region_end = 1000,
    bin_bp = 100, symmetric = TRUE,
    col_tname = "sn", col_tstart = "ss",
    col_cg = "cg", col_mapq = "mq"
  )
  m_ref <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  expect_equal(as.matrix(m), as.matrix(m_ref))
})

# ---- stream_paf2mat_anch -------------------------------------------------

test_that("stream_paf2mat_anch with min_anchor_bp = 1 matches the unfiltered version", {
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)
  m_plain <- stream_paf2mat(paf, "chr1", 1, 1000, bin_bp = 100, symmetric = TRUE)
  m_anch  <- stream_paf2mat_anch(paf, "chr1", 1, 1000, bin_bp = 100,
                                 min_anchor_bp = 1, symmetric = TRUE)
  expect_equal(as.matrix(m_plain), as.matrix(m_anch))
})

test_that("stream_paf2mat_anch with a high min_anchor_bp drops short-anchor junctions", {
  # Anchor lengths are the M blocks flanking the N. Our fixture has 100M
  # on either side. min_anchor_bp = 50 keeps everything; min_anchor_bp =
  # 200 should drop all junctions (anchors are only 100bp).
  paf <- make_fixture_paf(); on.exit(unlink(paf), add = TRUE)

  m_50 <- stream_paf2mat_anch(paf, "chr1", 1, 1000, bin_bp = 100,
                              min_anchor_bp = 50, symmetric = TRUE)
  expect_gt(sum(m_50), 0)

  expect_error(
    stream_paf2mat_anch(paf, "chr1", 1, 1000, bin_bp = 100,
                        min_anchor_bp = 200, symmetric = TRUE),
    "No edges"
  )
})
