test_that("buildSpacedCloneTable joins coordinates and aggregates clone counts", {
  clones <- data.frame(
    barcode = c("b1", "b1", "b2"),
    clone_id = c("1", "1", "2"),
    consensus_count = c(2, 3, 4),
    stringsAsFactors = FALSE
  )
  positions <- data.frame(
    barcode = c("b1", "b2"),
    x = c(10, 20),
    y = c(30, 40),
    stringsAsFactors = FALSE
  )

  out <- buildSpacedCloneTable(
    clones,
    positions = positions,
    count_col = "consensus_count",
    clone_prefix = "heavy_"
  )

  expect_equal(nrow(out), 2)
  expect_equal(out$n[out$clone_id == "heavy_1"], 5)
  expect_equal(out$x[out$barcode == "b2"], 20)
})

test_that("spacedSpatialAffinity gives nearby native spatial features partial credit", {
  heavy <- data.frame(
    clone_id = c("heavy_1", "heavy_2"),
    barcode = c("bin_a", "bin_b"),
    x = c(0, 100),
    y = c(0, 100),
    n = c(1, 1),
    stringsAsFactors = FALSE
  )
  light <- data.frame(
    clone_id = c("light_1", "light_2"),
    barcode = c("bin_a_near", "bin_b_near"),
    x = c(8, 108),
    y = c(0, 100),
    n = c(1, 1),
    stringsAsFactors = FALSE
  )

  affinity <- spacedSpatialAffinity(heavy, light, kernel_sigma = 16, kernel_radius = 48)

  expect_equal(rownames(affinity), c("heavy_1", "heavy_2"))
  expect_equal(colnames(affinity), c("light_1", "light_2"))
  expect_gt(affinity["heavy_1", "light_1"], affinity["heavy_1", "light_2"])
  expect_gt(affinity["heavy_2", "light_2"], affinity["heavy_2", "light_1"])
  expect_equal(unname(rowSums(affinity)), c(1, 1))
})

test_that("pairHeavyLightSpaced solves weighted assignments with clue", {
  spatial <- matrix(
    c(0.9, 0.1, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("heavy_1", "heavy_2"), c("light_1", "light_2"))
  )

  if (!requireNamespace("clue", quietly = TRUE)) {
    expect_error(
      pairHeavyLightSpaced(spatial, weights = 1),
      "Package 'clue' is required"
    )
    return()
  }

  pairs <- pairHeavyLightSpaced(spatial, weights = 1)

  expect_equal(pairs$heavy_id, c("heavy_1", "heavy_2"))
  expect_equal(pairs$light_id, c("light_1", "light_2"))
  expect_equal(pairs$weight, c(1, 1))
  expect_true(all(c("spatial_affinity_score", "repair_score", "combined_score", "normalized_pairing_score") %in% colnames(pairs)))
  expect_true(all(pairs$normalized_pairing_score > 0.5))
})

test_that("pairHeavyLightSpaced accepts unnamed mapping matrix paths positionally", {
  spatial <- matrix(
    c(0.1, 0.9, 0.9, 0.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("heavy_1", "heavy_2"), c("light_1", "light_2"))
  )
  mapping_path <- tempfile(fileext = ".tsv")
  utils::write.table(
    matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE),
    mapping_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  if (!requireNamespace("clue", quietly = TRUE)) {
    expect_error(
      pairHeavyLightSpaced(spatial, mapping_matrix = mapping_path, weights = 0),
      "Package 'clue' is required"
    )
    return()
  }

  pairs <- pairHeavyLightSpaced(spatial, mapping_matrix = mapping_path, weights = 0)

  expect_equal(pairs$heavy_id, c("heavy_1", "heavy_2"))
  expect_equal(pairs$light_id, c("light_1", "light_2"))
  expect_equal(pairs$weight, c(0, 0))
  expect_true(all(pairs$repair_score > 0.5))
})

test_that("pairHeavyLightSpaced returns a list for multiple weights", {
  spatial <- matrix(
    c(0.9, 0.1, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("heavy_1", "heavy_2"), c("light_1", "light_2"))
  )

  if (!requireNamespace("clue", quietly = TRUE)) {
    expect_error(
      pairHeavyLightSpaced(spatial, weights = c(0, 1)),
      "Package 'clue' is required"
    )
    return()
  }

  pairs <- pairHeavyLightSpaced(spatial, weights = c(0, 1))

  expect_type(pairs, "list")
  expect_equal(names(pairs), c("w=0", "w=1"))
  expect_equal(pairs[["w=1"]]$light_id, c("light_1", "light_2"))
})

test_that("buildRepairInputMatrix creates barcode-by-clone matrix", {
  heavy <- data.frame(
    clone_id = c("heavy_1", "heavy_2"),
    barcode = c("b1", "b2"),
    x = c(0, 10),
    y = c(0, 10),
    n = c(2, 3),
    stringsAsFactors = FALSE
  )
  light <- data.frame(
    clone_id = c("light_1", "light_1"),
    barcode = c("b1", "b2"),
    x = c(0, 10),
    y = c(0, 10),
    n = c(5, 7),
    stringsAsFactors = FALSE
  )

  repair <- buildRepairInputMatrix(heavy, light)

  expect_true(all(c("barcode", "heavy_1", "heavy_2", "light_1") %in% colnames(repair)))
  expect_equal(repair$heavy_1[repair$barcode == "b1"], 2)
  expect_equal(repair$light_1[repair$barcode == "b2"], 7)
})

test_that("repairAnalyzeCommand includes SPaCeD heavy and light patterns", {
  cmd <- repairAnalyzeCommand(
    input_path = "data/repair.tsv",
    output_dir = "repair_out",
    a_pattern = "^heavy",
    b_pattern = "^light"
  )

  expect_match(cmd, "repair")
  expect_match(cmd, "-ap '^heavy'", fixed = TRUE)
  expect_match(cmd, "-bp '^light'", fixed = TRUE)
})
