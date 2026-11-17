# A SIRUP object can hold several receptor tables (slots), for example BCRs and TCRs.

make_slot_sirup <- function() {
  x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  bcr <- getReceptors(x)
  tcr <- bcr[seq_len(min(8L, nrow(bcr))), ]
  tcr$locus <- "TRB"
  createSIRUP(
    spatial = getSpatial(x),
    receptors = list(bcr = bcr, tcr = tcr),
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE,
    polygons = getPolygons(x)
  )
}

test_that("a single receptor table has the slot name only_slot", {
  x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  expect_equal(receptorSlots(x), "only_slot")
  expect_s3_class(getReceptors(x), "tbl_df")
  expect_identical(getReceptors(x, slot = "only_slot"), getReceptors(x))
  expect_error(getReceptors(x, slot = "bcr"), "Unknown receptor slot")
})

test_that("createSIRUP accepts a named list of receptor tables", {
  y <- make_slot_sirup()
  expect_s3_class(y, "sirup")
  expect_equal(receptorSlots(y), c("bcr", "tcr"))
  expect_true(all(getReceptors(y, slot = "tcr")$locus == "TRB"))
  expect_output(print(y), "2 slots")
  expect_output(print(y), "tcr:")
})

test_that("slot = NULL is an error unless there is one slot", {
  y <- make_slot_sirup()
  expect_error(getReceptors(y), "`slot` is required")
  expect_error(getReceptors(y, slot = "zzz"), "Unknown receptor slot")
  expect_error(updateReceptors(y, dplyr::filter, clone_id != "3"), "`slot` is required")
  expect_error(withPolygons(y, nrow), "`slot` is required")
  expect_error(polygonReceptorSummary(y), "`slot` is required")
  one <- createSIRUP(
    spatial = getSpatial(y), receptors = list(only = getReceptors(y, slot = "bcr")),
    spatial_id = "spatial_id", clone = "clone_id", image_scale = "hires", swap_coords = FALSE
  )
  expect_equal(nrow(getReceptors(one)), nrow(getReceptors(y, slot = "bcr")))
})

test_that("createSIRUP checks the receptor list and the shared fields", {
  y <- make_slot_sirup()
  bcr <- getReceptors(y, slot = "bcr")
  tcr <- getReceptors(y, slot = "tcr")
  spatial <- getSpatial(y)
  expect_error(createSIRUP(spatial, list(bcr, tcr), spatial_id = "spatial_id"), "unique, non-empty name")
  expect_error(createSIRUP(spatial, list(a = bcr, a = tcr), spatial_id = "spatial_id"), "unique, non-empty name")
  expect_error(createSIRUP(spatial, list(), spatial_id = "spatial_id"), "unique, non-empty name")
  expect_error(createSIRUP(spatial, list(a = 1), spatial_id = "spatial_id"), "must be a data frame")
  expect_error(
    createSIRUP(spatial, list(bcr = bcr, tcr = tcr), spatial_id = "spatial_id", clone = "no_such_column"),
    "missing from every receptor slot"
  )
  renamed <- tcr
  names(renamed)[names(renamed) == "spatial_id"] <- "spatial_id_8um"
  expect_error(createSIRUP(spatial, list(bcr = bcr, tcr = renamed)), "different spatial ID columns")
})

test_that("validateSIRUP checks every slot", {
  y <- make_slot_sirup()
  bad <- y
  bad$receptors$tcr$spatial_id <- "not-a-spot"
  expect_error(validateSIRUP(bad), "slot `tcr`")
  bad2 <- y
  bad2$receptors$tcr$spatial_id <- NA_character_
  expect_error(validateSIRUP(bad2), "x\\$receptors\\$tcr")
  bad3 <- y
  bad3$receptors$tcr <- 5
  expect_error(validateSIRUP(bad3), "tibble or a named list")
})

test_that("stored polygons annotate every slot", {
  y <- make_slot_sirup()
  for (slot in receptorSlots(y)) {
    table <- getReceptors(y, slot = slot)
    expect_true("polygon_id" %in% names(table))
  }
  # a spatial ID gets the same polygon in both slots
  bcr <- getReceptors(y, slot = "bcr")
  tcr <- getReceptors(y, slot = "tcr")
  shared <- intersect(bcr$spatial_id, tcr$spatial_id)
  expect_true(length(shared) > 0L)
  polygon_of <- function(table) stats::setNames(table$polygon_id, table$spatial_id)[shared]
  expect_equal(polygon_of(bcr), polygon_of(tcr))
})

test_that("updateReceptors changes only the chosen slot", {
  y <- make_slot_sirup()
  n_bcr <- nrow(getReceptors(y, slot = "bcr"))
  drop_clone <- getReceptors(y, slot = "tcr")$clone_id[[1L]]
  y2 <- updateReceptors(y, function(t) t[t$clone_id != drop_clone, ], slot = "tcr")
  expect_equal(nrow(getReceptors(y2, slot = "bcr")), n_bcr)
  expect_lt(nrow(getReceptors(y2, slot = "tcr")), nrow(getReceptors(y, slot = "tcr")))
  # a mapped column needs to be in one slot only
  y3 <- updateReceptors(y, function(t) dplyr::mutate(t, extra = 1), slot = "tcr", receptor_fields = list(count = "extra"))
  expect_equal(y3$metadata$receptor_fields$count, "extra")
  expect_error(
    updateReceptors(y, function(t) t, slot = "tcr", receptor_fields = list(count = "no_such_column")),
    "no_such_column"
  )
})

test_that("subsetSIRUP filters every slot, and receptor_rows only the chosen slot", {
  y <- make_slot_sirup()
  tcr_ids <- unique(getReceptors(y, slot = "tcr")$spatial_id)[1:2]
  s <- subsetSIRUP(y, spatial_id = tcr_ids)
  expect_true(all(getReceptors(s, slot = "bcr")$spatial_id %in% tcr_ids))
  expect_true(all(getReceptors(s, slot = "tcr")$spatial_id %in% tcr_ids))
  drop_clone <- getReceptors(y, slot = "tcr")$clone_id[[1L]]
  s2 <- subsetSIRUP(y, receptor_rows = function(t) t$clone_id != drop_clone, slot = "tcr")
  expect_equal(nrow(getReceptors(s2, slot = "bcr")), nrow(getReceptors(y, slot = "bcr")))
  expect_lt(nrow(getReceptors(s2, slot = "tcr")), nrow(getReceptors(y, slot = "tcr")))
  expect_error(subsetSIRUP(y, receptor_rows = function(t) rep(TRUE, nrow(t))), "`slot` is required")
})

test_that("mergeSIRUP binds each slot and needs the same slots", {
  y <- make_slot_sirup()
  m <- mergeSIRUP(y, y, samples = c("a", "b"))
  expect_equal(receptorSlots(m), c("bcr", "tcr"))
  expect_equal(nrow(getReceptors(m, slot = "bcr")), 2L * nrow(getReceptors(y, slot = "bcr")))
  expect_equal(nrow(getReceptors(m, slot = "tcr")), 2L * nrow(getReceptors(y, slot = "tcr")))
  expect_equal(nrow(getReceptors(m, sample_id = "a", slot = "tcr")), nrow(getReceptors(y, slot = "tcr")))
  single <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  expect_error(mergeSIRUP(y, single, samples = c("a", "b")), "receptor slots differ")
  s <- summary(m)
  expect_true("slot" %in% names(s$samples))
  expect_equal(sort(unique(s$samples$slot)), c("bcr", "tcr"))
})

test_that("saveSIRUP and readSIRUP keep the slots", {
  y <- make_slot_sirup()
  path <- saveSIRUP(y, tempfile(fileext = ".sirup.rds"))
  z <- readSIRUP(path)
  expect_equal(receptorSlots(z), c("bcr", "tcr"))
  expect_equal(getReceptors(z, slot = "tcr"), getReceptors(y, slot = "tcr"))
})

test_that("analysis names end with the slot for a multi-slot object", {
  y <- make_slot_sirup()
  expect_equal(.sirup_slot_analysis_name(y, "trees", "tcr"), "trees_tcr")
  one <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  expect_equal(.sirup_slot_analysis_name(one, "trees"), "trees")
  expect_error(.sirup_slot_analysis_name(y, "trees"), "`slot` is required")
})

test_that("polygon functions analyse one slot at a time", {
  y <- make_slot_sirup()
  y <- polygonReceptorSummary(y, slot = "bcr")
  y <- polygonReceptorSummary(y, slot = "tcr")
  expect_setequal(listAnalyses(y), c("polygon_receptor_summary_bcr", "polygon_receptor_summary_tcr"))
  n_bcr <- sum(getAnalysis(y, "polygon_receptor_summary_bcr")$n_receptors)
  n_tcr <- sum(getAnalysis(y, "polygon_receptor_summary_tcr")$n_receptors)
  expect_gt(n_bcr, n_tcr)

  y <- polygonDiversity(y, slot = "bcr", min_n = 2, max_q = 2, step_q = 1, nboot = 3, common_size = FALSE)
  expect_true("polygon_diversity_bcr" %in% listAnalyses(y))
  y <- summarizePolygonDiversity(y, q = 1, slot = "bcr")
  expect_true("polygon_diversity_q1_bcr" %in% listAnalyses(y))
  expect_error(summarizePolygonDiversity(y, q = 1), "`slot` is required")
  y <- polygonAbundance(y, slot = "tcr", min_n = 2, nboot = 3, common_size = FALSE)
  expect_true("polygon_abundance_tcr" %in% listAnalyses(y))

  expect_type(withPolygons(y, nrow, slot = "tcr"), "list")
  expect_error(withPolygons(y, nrow), "`slot` is required")
})

test_that("spatial statistics read one slot", {
  skip_if_not_installed("spdep")
  y <- make_slot_sirup()
  # genes and metadata need no receptors, so they need no slot
  expect_no_error(spatialMoran(y, variable = "gene1", variable_type = "feature", k = 4))
  expect_error(spatialMoran(y, variable = "1", variable_type = "clone", k = 4), "`slot` is required")
  expect_no_error(spatialMoran(y, variable = "1", variable_type = "clone", slot = "bcr", k = 4))
  spread_bcr <- spatialCloneSpread(y, slot = "bcr")
  spread_tcr <- spatialCloneSpread(y, slot = "tcr")
  expect_s3_class(spread_bcr, "data.frame")
  expect_s3_class(spread_tcr, "data.frame")
  expect_error(spatialCloneSpread(y), "`slot` is required")
  jc <- spatialCloneJoinCount(y, slot = "bcr", min_spatial_ids = 2)
  expect_true(is.list(jc))
  expect_error(spatialCloneJoinCount(y), "`slot` is required")
})

test_that("receptor summaries and plots read one slot", {
  y <- make_slot_sirup()
  s_bcr <- spatialReceptorSummary(y, slot = "bcr")
  s_tcr <- spatialReceptorSummary(y, slot = "tcr")
  expect_gt(s_bcr$n_sequences, s_tcr$n_sequences)
  expect_error(spatialReceptorSummary(y), "`slot` is required")
  expect_error(plotLengthQc(y), "`slot` is required")
  expect_s3_class(plotLengthQc(y, slot = "bcr"), "ggplot")
})

test_that("a slot with no rows is allowed", {
  y <- make_slot_sirup()
  tcr <- getReceptors(y, slot = "tcr")
  y2 <- updateReceptors(y, function(r) r[0, ], slot = "tcr")
  expect_equal(nrow(getReceptors(y2, slot = "tcr")), 0L)
  expect_no_error(validateSIRUP(y2))
  expect_gt(nrow(getReceptors(y2, slot = "bcr")), 0L)
  # an empty table is still an error for a single-table object
  single <- createSIRUP(
    spatial = getSpatial(y), receptors = getReceptors(y, slot = "bcr"),
    spatial_id = "spatial_id"
  )
  expect_error(
    subsetSIRUP(single, receptor_rows = rep(FALSE, nrow(getReceptors(single)))),
    "no non-missing spatial IDs"
  )
})

test_that("optional fields need to be in one slot only", {
  y <- make_slot_sirup()
  tcr <- getReceptors(y, slot = "tcr")
  tcr$clone_id <- NULL
  z <- createSIRUP(
    spatial = getSpatial(y),
    receptors = list(bcr = getReceptors(y, slot = "bcr"), tcr = tcr),
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(z$metadata$receptor_fields$clone, "clone_id")
  expect_no_error(validateSIRUP(z))
  expect_s3_class(spatialCloneSpread(z, slot = "bcr"), "data.frame")
  # a clone field that no slot has is still an error
  expect_error(
    createSIRUP(
      spatial = getSpatial(y), receptors = list(bcr = getReceptors(y, slot = "bcr")),
      spatial_id = "spatial_id", clone = "no_such_column"
    ),
    "no_such_column"
  )
  # updateReceptors records clone_id when one slot gets it
  w <- createSIRUP(
    spatial = getSpatial(y),
    receptors = list(bcr = dplyr::select(getReceptors(y, slot = "bcr"), -clone_id), tcr = tcr),
    spatial_id = "spatial_id", image_scale = "hires", swap_coords = FALSE
  )
  expect_null(w$metadata$receptor_fields$clone)
  w <- updateReceptors(w, function(r) { r$clone_id <- "1"; r }, slot = "bcr")
  expect_equal(w$metadata$receptor_fields$clone, "clone_id")
  # the sample column must still be in every slot
  expect_error(
    createSIRUP(
      spatial = list(a = getSpatial(y)),
      receptors = list(bcr = dplyr::mutate(getReceptors(y, slot = "bcr"), sample = "a"), tcr = tcr),
      spatial_id = "spatial_id", sample_id = "sample"
    ),
    "tcr"
  )
})
