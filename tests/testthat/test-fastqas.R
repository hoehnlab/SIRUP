# Round-trip + contract tests for read_fastqas / write_fastqas / stream_fastqas.
# All three share a single dataframe contract documented in R/r_commands.R:
#   - `id`      always present; short id, no leading @/> prefix
#   - `full_id` present only when full_id = TRUE; full header, no prefix
#   - `seq`     always present
#   - `qual`    fastq only
#   - column order: id, full_id (if requested), seq, qual (if fastq)
#   - type accepts "fq" / "fastq" / "fa" / "fasta"

make_fixture_fq <- function() {
  tmp <- tempfile(fileext = ".fq")
  writeLines(c(
    "@read1 comment=foo",
    "ACGTACGT",
    "+",
    "IIIIIIII",
    "@read2 something else",
    "GGGGCCCC",
    "+",
    "JJJJJJJJ",
    "@read3",                # no comment / full_id == id
    "TTTTAAAA",
    "+",
    "KKKKKKKK"
  ), tmp)
  tmp
}

make_fixture_fa <- function() {
  tmp <- tempfile(fileext = ".fa")
  writeLines(c(
    ">read1 comment=foo",
    "ACGTACGT",
    ">read2 something else",
    "GGGGCCCC",
    ">read3",
    "TTTTAAAA"
  ), tmp)
  tmp
}

# ---- canonical contract --------------------------------------------------

test_that("read_fastqas strips the @ prefix from id (fastq)", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  df <- read_fastqas(fq)
  expect_equal(df$id, c("read1", "read2", "read3"))
  expect_false(any(startsWith(df$id, "@")))
})

test_that("read_fastqas strips the > prefix from id (fasta)", {
  fa <- make_fixture_fa(); on.exit(unlink(fa), add = TRUE)
  df <- read_fastqas(fa, type = "fa")
  expect_equal(df$id, c("read1", "read2", "read3"))
  expect_false(any(startsWith(df$id, ">")))
})

test_that("read_fastqas full_id captures the entire header, including the comment", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  df <- read_fastqas(fq, full_id = TRUE)
  expect_equal(df$id,      c("read1", "read2", "read3"))
  expect_equal(df$full_id, c("read1 comment=foo", "read2 something else", "read3"))
})

test_that("read_fastqas column order matches the canonical contract", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  expect_identical(names(read_fastqas(fq)),
                   c("id", "seq", "qual"))
  expect_identical(names(read_fastqas(fq, full_id = TRUE)),
                   c("id", "full_id", "seq", "qual"))
})

test_that("fasta omits qual; full_id placement is consistent", {
  fa <- make_fixture_fa(); on.exit(unlink(fa), add = TRUE)
  expect_identical(names(read_fastqas(fa, type = "fa")),
                   c("id", "seq"))
  expect_identical(names(read_fastqas(fa, type = "fa", full_id = TRUE)),
                   c("id", "full_id", "seq"))
})

# ---- read_fastqas and stream_fastqas produce the same output -------------

test_that("read_fastqas and stream_fastqas agree on column names + values (fastq)", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  r <- read_fastqas(fq, full_id = TRUE)
  s <- stream_fastqas(fq, full_id = TRUE)
  expect_identical(names(r), names(s))
  expect_identical(r$id,      s$id)
  expect_identical(r$full_id, s$full_id)
  expect_identical(r$seq,     s$seq)
  expect_identical(r$qual,    s$qual)
})

test_that("read_fastqas and stream_fastqas agree on column names + values (fasta)", {
  fa <- make_fixture_fa(); on.exit(unlink(fa), add = TRUE)
  r <- read_fastqas(fa, type = "fa", full_id = TRUE)
  s <- stream_fastqas(fa, type = "fa", full_id = TRUE)
  expect_identical(names(r), names(s))
  expect_identical(r$id,      s$id)
  expect_identical(r$full_id, s$full_id)
  expect_identical(r$seq,     s$seq)
})

# ---- type aliases --------------------------------------------------------

test_that("type aliases: fq == fastq, fa == fasta", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  fa <- make_fixture_fa(); on.exit(unlink(fa), add = TRUE)

  expect_identical(read_fastqas(fq, "fq"),   read_fastqas(fq, "fastq"))
  expect_identical(read_fastqas(fa, "fa"),   read_fastqas(fa, "fasta"))

  expect_identical(stream_fastqas(fq, "fq"), stream_fastqas(fq, "fastq"))
  expect_identical(stream_fastqas(fa, "fa"), stream_fastqas(fa, "fasta"))
})

test_that("invalid type errors clearly", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  expect_error(read_fastqas(fq, type = "sam"),    "one of")
  expect_error(write_fastqas(data.frame(id = "x", seq = "A", qual = "I"),
                             tempfile(), type = "sam"),
               "one of")
})

# ---- write_fastqas: prefix + header_col semantics ------------------------

test_that("write_fastqas auto-prepends @/> and strips existing prefix (no doubling)", {
  df <- data.frame(
    id   = c("@read1",      "read2", ">read3", "read4"),
    seq  = c("ACGT",        "GGCC",  "TTAA",   "AAAA"),
    qual = c("IIII",        "JJJJ",  "KKKK",   "LLLL"),
    stringsAsFactors = FALSE
  )
  out <- tempfile(fileext = ".fq"); on.exit(unlink(out), add = TRUE)
  write_fastqas(df, out, type = "fq")
  lines <- readLines(out)
  expect_equal(lines[1],  "@read1")     # @ stripped, single @ re-prepended
  expect_equal(lines[5],  "@read2")
  expect_equal(lines[9],  "@read3")     # > stripped, @ used for fastq
  expect_equal(lines[13], "@read4")

  out2 <- tempfile(fileext = ".fa"); on.exit(unlink(out2), add = TRUE)
  write_fastqas(df[, c("id","seq")], out2, type = "fa")
  fa_lines <- readLines(out2)
  expect_equal(fa_lines[1], ">read1")
  expect_equal(fa_lines[3], ">read2")
  expect_equal(fa_lines[5], ">read3")
  expect_equal(fa_lines[7], ">read4")
})

test_that("write_fastqas errors when required columns are missing", {
  out <- tempfile(fileext = ".fq")
  # No seq
  expect_error(
    write_fastqas(data.frame(id = "x", qual = "I"), out, type = "fq"),
    "`seq`"
  )
  # No qual, but type = fq
  expect_error(
    write_fastqas(data.frame(id = "x", seq = "A"), out, type = "fq"),
    "qual"
  )
  # Bad header_col
  expect_error(
    write_fastqas(data.frame(id = "x", seq = "A", qual = "I"), out,
                  type = "fq", header_col = "nope"),
    "header_col"
  )
})

# ---- round-trips ---------------------------------------------------------

test_that("fastq round-trip with header_col = 'full_id' is lossless", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  df1 <- read_fastqas(fq, full_id = TRUE)
  out <- tempfile(fileext = ".fq"); on.exit(unlink(out), add = TRUE)
  write_fastqas(df1, out, header_col = "full_id")
  df2 <- read_fastqas(out, full_id = TRUE)
  expect_identical(df1, df2)
})

test_that("fastq round-trip with default header_col = 'id' is lossy on the comment", {
  fq <- make_fixture_fq(); on.exit(unlink(fq), add = TRUE)
  df1 <- read_fastqas(fq, full_id = TRUE)
  out <- tempfile(fileext = ".fq"); on.exit(unlink(out), add = TRUE)
  write_fastqas(df1, out)              # default header_col = "id"
  df2 <- read_fastqas(out, full_id = TRUE)
  # id, seq, qual preserved; full_id collapses to id (comment dropped).
  expect_equal(df2$id,   df1$id)
  expect_equal(df2$seq,  df1$seq)
  expect_equal(df2$qual, df1$qual)
  expect_equal(df2$full_id, df1$id)
})

test_that("fasta round-trip with header_col = 'full_id' is lossless", {
  fa <- make_fixture_fa(); on.exit(unlink(fa), add = TRUE)
  df1 <- read_fastqas(fa, type = "fa", full_id = TRUE)
  out <- tempfile(fileext = ".fa"); on.exit(unlink(out), add = TRUE)
  write_fastqas(df1, out, type = "fa", header_col = "full_id")
  df2 <- read_fastqas(out, type = "fa", full_id = TRUE)
  expect_identical(df1, df2)
})

# ---- multi-file dispatch -------------------------------------------------

test_that("read_fastqas concatenates multiple files in order", {
  fq1 <- tempfile(fileext = ".fq"); on.exit(unlink(fq1), add = TRUE)
  fq2 <- tempfile(fileext = ".fq"); on.exit(unlink(fq2), add = TRUE)
  writeLines(c("@r1", "ACGT", "+", "IIII"), fq1)
  writeLines(c("@r2", "GGCC", "+", "JJJJ"), fq2)

  df <- read_fastqas(c(fq1, fq2))
  expect_equal(df$id, c("r1", "r2"))
  expect_equal(df$seq, c("ACGT", "GGCC"))
})
