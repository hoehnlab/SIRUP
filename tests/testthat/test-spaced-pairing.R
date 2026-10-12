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
    count = "consensus_count",
    clone_prefix = "heavy_"
  )

  expect_equal(nrow(out), 2)
  expect_equal(out$n[out$clone_id == "heavy_1"], 5)
  expect_equal(out$x[out$barcode == "b2"], 20)
})

test_that("buildSpacedCloneTable uses supplied positions when data already has coordinates", {
  clones <- data.frame(
    barcode = c("b1", "b2"),
    clone_id = c("1", "2"),
    x = c(999, 999),
    y = c(999, 999),
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
    clone_prefix = "heavy_"
  )

  expect_equal(out$x, c(10, 20))
  expect_equal(out$y, c(30, 40))
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

  affinity <- spacedSpatialAffinity(
    heavy,
    light,
    kernel_sigma = 16,
    kernel_radius = 48,
    min_x_obs = 1
  )

  expect_equal(rownames(affinity), c("heavy_1", "heavy_2"))
  expect_equal(colnames(affinity), c("light_1", "light_2"))
  expect_gt(affinity["heavy_1", "light_1"], affinity["heavy_1", "light_2"])
  expect_gt(affinity["heavy_2", "light_2"], affinity["heavy_2", "light_1"])
  expect_equal(unname(rowSums(affinity)), c(1, 1))
})

test_that("spacedSpatialAffinity filters clones with REPAIR-style min_x_obs", {
  heavy <- data.frame(
    clone_id = c("heavy_1", "heavy_1", "heavy_2"),
    barcode = c("bin_a", "bin_b", "bin_c"),
    x = c(0, 1, 100),
    y = c(0, 1, 100),
    n = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  light <- data.frame(
    clone_id = c("light_1", "light_1", "light_2"),
    barcode = c("bin_a", "bin_b", "bin_c"),
    x = c(0, 1, 100),
    y = c(0, 1, 100),
    n = c(1, 1, 1),
    stringsAsFactors = FALSE
  )

  affinity <- spacedSpatialAffinity(heavy, light, kernel_sigma = 16, kernel_radius = 48)

  expect_equal(rownames(affinity), "heavy_1")
  expect_equal(colnames(affinity), "light_1")
})

test_that("pairHeavyLightSpaced returns top light candidates per heavy clone", {
  spatial <- matrix(
    c(0.9, 0.1, 0.8, 0.2, 0.1, 0.9),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("heavy_1", "heavy_2", "heavy_3"), c("light_1", "light_2"))
  )

  pairs <- pairHeavyLightSpaced(spatial, weights = 1, top_n_light_chains = 1)

  expect_equal(pairs$heavy_id, c("heavy_1", "heavy_2", "heavy_3"))
  expect_equal(pairs$light_id, c("light_1", "light_1", "light_2"))
  expect_equal(pairs$weight, c(1, 1, 1))
  expect_true(all(c(
    "spatial_affinity_score",
    "repair_score",
    "combined_score",
    "normalized_pairing_score",
    "light_specificity_score"
  ) %in% colnames(pairs)))
  expect_equal(
    pairs$light_specificity_score[pairs$heavy_id == "heavy_1" & pairs$light_id == "light_1"],
    0.9 / 1.8
  )
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

  pairs <- pairHeavyLightSpaced(spatial, mapping_matrix = mapping_path, weights = 0, top_n_light_chains = 1)

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

  pairs <- pairHeavyLightSpaced(spatial, weights = c(0, 1), top_n_light_chains = 2)

  expect_type(pairs, "list")
  expect_equal(names(pairs), c("w=0", "w=1"))
  expect_equal(nrow(pairs[["w=1"]]), 4)
  expect_equal(pairs[["w=1"]]$light_id[c(1, 3)], c("light_1", "light_2"))
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

test_that("checkSpacedRepair accepts REPAIR that exports a model matrix", {
  testthat::skip_on_os("windows")
  fake_conda <- tempfile("conda")
  writeLines(
    c(
      "#!/bin/sh",
      "echo SIRUP_REPAIR_ANALYZE_PATH=/tmp/repair/analyze.py",
      "echo SIRUP_REPAIR_EXPORTS_MODEL_MATRIX=True"
    ),
    fake_conda
  )
  Sys.chmod(fake_conda, "0755")

  repair <- checkSpacedRepair(conda_env = "repair", conda_exe = fake_conda, verbose = FALSE)

  expect_equal(repair$analyze_path, "/tmp/repair/analyze.py")
  expect_true(repair$exports_model_matrix)
})

test_that("checkSpacedRepair rejects REPAIR without a model matrix export", {
  testthat::skip_on_os("windows")
  fake_conda <- tempfile("conda")
  writeLines(
    c(
      "#!/bin/sh",
      "echo SIRUP_REPAIR_ANALYZE_PATH=/tmp/repair/analyze.py",
      "echo SIRUP_REPAIR_EXPORTS_MODEL_MATRIX=False"
    ),
    fake_conda
  )
  Sys.chmod(fake_conda, "0755")

  expect_error(
    checkSpacedRepair(conda_env = "repair", conda_exe = fake_conda, verbose = FALSE),
    "modified REPAIR fork",
    fixed = TRUE
  )
})
