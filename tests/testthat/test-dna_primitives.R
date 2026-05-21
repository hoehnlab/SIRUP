# Tests for the int64-packed DNA primitives + revcomp/revcomp_bits/hex.
#
# Encoding (inferred from observed behavior):
#   A = 00,  C = 01,  G = 11,  T = 10  (yes, G=11 and T=10, not the obvious
#   G=10/T=11). Two bits per base, MSB-first into an int64 — so the first
#   base of the string occupies the highest 2 bits, the second base the next
#   2, etc. One int64 (64 bits) therefore holds up to 32 bases. Longer
#   sequences require `sequence_to_bitlist()` which returns a list of int64
#   chunks (32 bases per chunk).

# ---- revcomp (string layer) ----------------------------------------------

test_that("revcomp on simple strings", {
  expect_equal(revcomp("ACGT"),  "ACGT")        # palindrome
  expect_equal(revcomp("AAAA"),  "TTTT")
  expect_equal(revcomp("TTTT"),  "AAAA")
  expect_equal(revcomp("CCCC"),  "GGGG")
  expect_equal(revcomp("GGGG"),  "CCCC")
  expect_equal(revcomp("ACGTACGT"), "ACGTACGT") # 8-base palindrome
})

test_that("revcomp preserves N and is involutive", {
  expect_equal(revcomp("NNNN"),       "NNNN")
  expect_equal(revcomp("ACNTG"),      "CANGT")
  for (s in c("ACGTACGT", "AACCGGTT", "NNNACGTNNN", "")) {
    expect_equal(revcomp(revcomp(s)), s, info = s)
  }
})

test_that("revcomp handles vector input + NA", {
  expect_equal(revcomp(c("ACGT","AAAA","NNNN","")),
               c("ACGT","TTTT","NNNN",""))
  expect_equal(as.character(revcomp(NA_character_)), NA_character_)
})

# ---- sequence_to_bits + bits_to_sequence round-trip ----------------------

test_that("sequence_to_bits + bits_to_sequence is lossless for <= 32 bases", {
  for (len in c(1, 2, 4, 8, 16, 31, 32)) {
    seqs <- c(strrep("A", len), strrep("T", len), strrep("G", len),
              strrep("C", len), substr("ACGTACGTACGTACGTACGTACGTACGTACGTACGT", 1, len))
    expect_equal(bits_to_sequence(sequence_to_bits(seqs), sequence_length = len),
                 seqs, info = paste0("len=", len))
  }
})

test_that("sequence_to_bits encodes the homopolymer constants we expect", {
  # Confirms the encoding scheme:
  #   A=00 -> 00000000 = 0
  #   C=01 -> 01010101 = 85
  #   T=10 -> 10101010 = 170
  #   G=11 -> 11111111 = 255
  expect_equal(as.numeric(sequence_to_bits(c("AAAA","CCCC","TTTT","GGGG"))),
               c(0, 85, 170, 255))
})

# ---- sequence_to_bitlist (multi-chunk encoding) --------------------------

test_that("sequence_to_bitlist returns a list parallel to inputs", {
  bl <- sequence_to_bitlist(c("AAAA", "CCCC", "ACGTACGT"))
  expect_type(bl, "list")
  expect_length(bl, 3)
  expect_equal(as.numeric(bl[[1]]), 0)
  expect_equal(as.numeric(bl[[2]]), 85)
})

test_that("sequence_to_bitlist splits into 32-base chunks for long sequences", {
  s32  <- strrep("ACGT", 8)              # 32 bases -> 1 chunk
  s64  <- strrep("ACGT", 16)             # 64 bases -> 2 chunks
  s96  <- strrep("ACGT", 24)             # 96 bases -> 3 chunks

  expect_length(sequence_to_bitlist(s32)[[1]], 1)
  expect_length(sequence_to_bitlist(s64)[[1]], 2)
  expect_length(sequence_to_bitlist(s96)[[1]], 3)

  # Both halves of s64 are ACGTACGTACGTACGTACGTACGTACGTACGT — the chunks
  # should be identical. Compare via integer64 to avoid double precision loss.
  bl64 <- sequence_to_bitlist(s64)[[1]]
  expect_true(bl64[1] == bl64[2])
})

# ---- bits_to_hex ---------------------------------------------------------

test_that("bits_to_hex agrees with hex of the int64", {
  # Use the homopolymers we already mapped:
  #   AAAA = 0   -> "00"
  #   CCCC = 85  -> "55"
  #   TTTT = 170 -> "AA"
  #   GGGG = 255 -> "FF"
  # bits_to_hex zero-pads to the byte width implied by sequence_length.
  hex_of <- bits_to_hex(sequence_to_bits(c("AAAA","CCCC","TTTT","GGGG")),
                        sequence_length = 4)
  expect_equal(toupper(hex_of), c("00", "55", "AA", "FF"))
})

# ---- revcomp_bits --------------------------------------------------------

test_that("revcomp_bits agrees with string-level revcomp through encode/decode", {
  for (len in c(4, 8, 16, 32)) {
    seqs <- c(strrep("A", len), strrep("AC", len/2), strrep("ACGT", len/4),
              "ACGTACGTACGTACGTACGTACGTACGTACGT"[1:1])
    seqs <- seqs[nchar(seqs) == len]
    if (!length(seqs)) next
    bits      <- sequence_to_bits(seqs)
    rc_bits   <- revcomp_bits(bits, sequence_length = len)
    expect_equal(bits_to_sequence(rc_bits, sequence_length = len),
                 revcomp(seqs),
                 info = paste0("len=", len))
  }
})

test_that("revcomp_bits on a palindrome returns the input unchanged", {
  acgt <- sequence_to_bits("ACGT")
  expect_equal(as.numeric(revcomp_bits(acgt, sequence_length = 4)),
               as.numeric(acgt))
})

# ---- hamming_distance_rcpp ----------------------------------------------
#
# Returns true base-level Hamming distance under the 2-bit-per-base
# encoding (A=00, C=01, G=11, T=10). Implementation collapses each 2-bit
# chunk of (a XOR b) down to a single bit before popcount, so it counts
# one per differing base regardless of whether the diff flipped one or
# two bits within the chunk. (Pre-fix it computed
# `ceil(popcount(XOR) / 2)`, which under-counted transversions —
# hamming("AAAA","CCCC") used to return 2; now returns 4.)

test_that("hamming_distance_rcpp returns 0 for identical inputs", {
  for (s in c("A","ACGT","AAAA","TTTTGGGGCCCC","NNNN")) {
    expect_equal(hamming_distance_rcpp(sequence_to_bits(s), sequence_to_bits(s)),
                 0, info = s)
  }
})

test_that("hamming_distance_rcpp is symmetric in its arguments", {
  pairs <- list(c("ACGT","TTTT"), c("AAAA","GGGG"), c("AC","CA"))
  for (p in pairs) {
    a <- sequence_to_bits(p[1]); b <- sequence_to_bits(p[2])
    expect_equal(hamming_distance_rcpp(a, b), hamming_distance_rcpp(b, a),
                 info = paste(p, collapse = " "))
  }
})

test_that("hamming_distance_rcpp counts differing bases (transversions + transitions)", {
  cases <- list(
    # Single-base
    list("A",    "A",    0),
    list("A",    "C",    1),    # transversion
    list("A",    "G",    1),    # transition
    list("A",    "T",    1),    # transversion
    # 4-base homopolymer comparisons — every base differs
    list("AAAA", "CCCC", 4),    # 4 transversions (previously returned 2)
    list("AAAA", "GGGG", 4),    # 4 transitions
    list("AAAA", "TTTT", 4),    # 4 transversions (previously returned 2)
    list("CCCC", "GGGG", 4),
    list("CCCC", "TTTT", 4),
    # Mixed
    list("ACGT", "TTTT", 3),    # A,C,G differ; T matches (previously returned 2)
    list("ACGT", "ACGA", 1),    # last base only
    list("ACGT", "ACGT", 0),    # palindrome with itself
    list("AC",   "CA",   2)     # both positions swapped (previously returned 1)
  )
  for (case in cases) {
    a_seq <- case[[1]]; b_seq <- case[[2]]; expected <- case[[3]]
    expect_equal(
      hamming_distance_rcpp(sequence_to_bits(a_seq), sequence_to_bits(b_seq)),
      expected,
      info = sprintf("hamming(%s, %s) -> %d", a_seq, b_seq, expected)
    )
  }
})

test_that("hamming_distance_rcpp matches a naive base-level reference for random sequences", {
  set.seed(42)
  base_hamming <- function(a, b) {
    sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
  }
  bases <- c("A","C","G","T")
  for (i in seq_len(40)) {
    len <- sample(c(4, 8, 16, 32), 1)
    a <- paste(sample(bases, len, replace = TRUE), collapse = "")
    b <- paste(sample(bases, len, replace = TRUE), collapse = "")
    expect_equal(
      hamming_distance_rcpp(sequence_to_bits(a), sequence_to_bits(b)),
      base_hamming(a, b),
      info = paste(a, b)
    )
  }
})
