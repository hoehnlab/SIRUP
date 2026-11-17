
#' Calculate global Moran's I
#'
#' `spatialMoran` tests one variable for spatial autocorrelation across spots using Moran's I.
#'
#' @details
#' The function builds a neighbor graph from the spot coordinates. It uses all
#' spots within `radius`. If `radius` is `NULL`, it uses the `k` nearest spots
#' (or fewer, if the image has fewer spots). Then it calls `spdep::moran.test()`
#' on the variable values. The function drops spots with non-finite values.
#'
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds one graph and runs one test for each sample or image.
#'
#' Set `variable_type = "clone"` to test where one clone occurs. The function
#' makes a 0/1 value for each spot: 1 if sthe spot contains the clone, 0 if not.
#'
#' @param object Seurat object or SIRUP object.
#' @param variable Name of the feature or the metadata column to test.
#' @param variable_type Type of `variable`: `feature` (gene expression),
#'   `metadata` (Seurat metadata column), `receptor` (receptor table column),
#'   or `clone` (clone ID).
#' @param value Metadata value to convert to a 0/1 indicator (default `NULL`).
#' @param receptors AIRR-style receptor table. You must give it if
#'   `variable_type` is `"receptor"` or `"clone"` and `object` is not a SIRUP
#'   object.
#' @param spatial_id Name of the `receptors` column that holds the spatial IDs.
#'   If `object` is a SIRUP object and this argument is `NULL`, the function uses
#'   the stored spatial ID field. Give a value if you use separate Seurat and
#'   receptor objects.
#' @param clone Name of the `receptors` column that holds the clone IDs. The
#'   function uses it only if `variable_type` is `"clone"`. If `object` is a
#'   SIRUP object and this argument is `NULL`, the function uses the stored
#'   clone field, if there is one.
#' @param assay Assay to use if `variable_type` is `"feature"`.
#' @param layer Seurat expression layer to use for features (default
#'   `"counts"`).
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param radius Distance that defines a neighborhood. If `NULL`, the function
#'   uses `k`.
#' @param k Number of nearest neighbors, used if `radius` is `NULL` (default 6).
#' @param alternative Alternative hypothesis for `spdep::moran.test()` (default
#'   `"greater"`).
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots and the function reads receptors. The function ignores `slot` if you
#'   give a plain receptor table.
#'
#' @return For one image, a list with these items: `test` (the
#'   `spdep::moran.test()` result), `data` (the values tested), `neighbors` and
#'   `weights` (from `spdep`), `radius`, `k`, and `variable`. For several
#'   images, the list has `tests`, `neighbors`, and `weights` as named lists with
#'   one entry for each image. It also has `summary`, a data frame with one row
#'   for each image (`image`, `statistic`, `p_value`, `moran_i`, `expectation`,
#'   `variance`), and the items `data`, `radius`, `k`, and `variable`.
#' @seealso [spatialLocalMoran()], [spatialLee()], [spatialCloneJoinCount()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Test gene1 for spatial autocorrelation.
#' result <- spatialMoran(x, "gene1", k = 4)
#' result$test
#'
#' @export
spatialMoran <- function(object,
                         variable,
                         variable_type = c("feature", "metadata", "receptor", "clone"),
                         value = NULL,
                         receptors = NULL,
                         spatial_id = NULL,
                         clone = NULL,
                         assay = NULL,
                         layer = "counts",
                         images = NULL,
                         image_scale = NULL,
                         swap_coords = NULL,
                         radius = NULL,
                         k = 6,
                         alternative = "greater",
                         slot = NULL) {
  dat <- .sirup_spatial_variable_values(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    spatial_id_col = spatial_id,
    clone_col = clone,
    assay = assay,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    slot = slot
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)
  groups <- .sirup_spatial_stat_groups(dat)
  if (length(groups) > 1L) {
    neighbors <- lapply(groups, function(group_dat) .sirup_spatial_neighbors(group_dat, radius = radius, k = k))
    weights <- lapply(neighbors, spdep::nb2listw, style = "W", zero.policy = TRUE)
    tests <- lapply(names(groups), function(name) {
      group_dat <- groups[[name]]
      .check_spatial_values(group_dat$value)
      spdep::moran.test(group_dat$value, weights[[name]], zero.policy = TRUE, alternative = alternative)
    })
    names(tests) <- names(groups)
    summary <- dplyr::bind_rows(lapply(names(tests), function(name) {
      estimate <- tests[[name]]$estimate
      data.frame(
        image = name,
        statistic = unname(tests[[name]]$statistic),
        p_value = tests[[name]]$p.value,
        moran_i = unname(estimate[["Moran I statistic"]]),
        expectation = unname(estimate[["Expectation"]]),
        variance = unname(estimate[["Variance"]]),
        stringsAsFactors = FALSE
      )
    }))
    return(list(
      tests = tests,
      summary = summary,
      data = dat,
      neighbors = neighbors,
      weights = weights,
      radius = radius,
      k = k,
      variable = unique(dat$variable)[[1]]
    ))
  }

  nb <- .sirup_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  test <- spdep::moran.test(dat$value, listw, zero.policy = TRUE, alternative = alternative)

  list(
    test = test,
    data = dat,
    neighbors = nb,
    weights = listw,
    radius = radius,
    k = k,
    variable = unique(dat$variable)[[1]]
  )
}

#' Calculate global Lee's L
#'
#' `spatialLee` tests the spatial association between two continuous variables.
#'
#' @details
#' The function builds one neighbor graph from all spots, with `radius` or `k`
#' as in [spatialMoran()]. For example, give two gene names to test if the expression of one gene
#' relates to the expression of the second gene in neighboring spots.
#'
#' `method = "analytical"` calls `spdep::lee.test()`. This function builds a
#' dense n by n weights matrix. The memory use grows with the square of the
#' number of spots, so be careful with larger datasets, like with Visium HD.
#' For large datasets, use `method = "mc"`. It calls `spdep::lee.mc()`, which permutes
#' the labels. It keeps the weights sparse and needs only a few GB.
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds one graph and runs one test for each sample or image.
#'
#' @param object Seurat object or SIRUP object.
#' @param x Name of the first variable. It is a feature or a numeric metadata
#'   column.
#' @param y Name of the second variable. It is a feature or a numeric metadata
#'   column.
#' @param x_type Type of `x`: `feature` or `metadata`.
#' @param y_type Type of `y`: `feature` or `metadata`.
#' @param assay Assay to use if `x_type` or `y_type` is `"feature"`.
#' @param layer Seurat expression layer to use for features (default
#'   `"counts"`).
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param radius Distance that defines a neighborhood. If `NULL`, the function
#'   uses `k`.
#' @param k Number of nearest neighbors, used if `radius` is `NULL` (default 6).
#'   Use `k = 8` for Visium HD, where the bins form a square grid. Use `k = 6`
#'   for standard Visium, where the spots form a hexagonal grid.
#' @param alternative Alternative hypothesis for the significance test (default
#'   `"greater"`).
#' @param method Significance test: `analytical` or `mc` (permutation). See
#'   Details.
#' @param nsim Number of permutations for `method = "mc"` (default 999).
#' @param progress Show a progress bar for the permutations (default `FALSE`)?
#'   The function uses it only if `method = "mc"`. Each permutation re-lags both
#'   variables over every spot. For Visium HD, a large `nsim` takes minutes for
#'   each pair.
#'
#' @return A list with these items: `test` (the test result), `L` (the observed
#'   Lee's L), `method`, `data` (the values tested), `neighbors` and `weights`
#'   (from `spdep`), `radius`, `k`, `x`, `y`, `x_type`, and `y_type`.
#' @seealso [spatialMoran()], [spatialLocalMoran()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Test if gene1 and gene2 are associated.
#' spatialLee(x, "gene1", "gene2", k = 4)$L
#'
#' # Use permutations. This method also works for large Visium HD data.
#' set.seed(1)
#' spatialLee(x, "gene1", "gene2", k = 4, method = "mc", nsim = 99)$test
#' @export
spatialLee <- function(object,
                       x,
                       y,
                       x_type = c("feature", "metadata"),
                       y_type = c("feature", "metadata"),
                       assay = NULL,
                       layer = "counts",
                       images = NULL,
                       image_scale = NULL,
                       swap_coords = NULL,
                       radius = NULL,
                       k = 6,
                       alternative = "greater",
                       method = c("analytical", "mc"),
                       nsim = 999,
                       progress = FALSE) {
  x_type <- match.arg(x_type)
  y_type <- match.arg(y_type)
  method <- match.arg(method)

  x_dat <- .sirup_spatial_variable_values(
    object = object,
    variable = x,
    variable_type = x_type,
    value = NULL,
    receptors = NULL,
    spatial_id_col = NULL,
    clone_col = NULL,
    assay = assay,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  y_dat <- .sirup_spatial_variable_values(
    object = object,
    variable = y,
    variable_type = y_type,
    value = NULL,
    receptors = NULL,
    spatial_id_col = NULL,
    clone_col = NULL,
    assay = assay,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )

  lee_by <- c("cell", "x", ".sirup_y", "image")
  if ("sample" %in% colnames(x_dat) && "sample" %in% colnames(y_dat)) {
    lee_by <- c("sample", lee_by)
  }
  dat <- merge(
    x_dat[, c(lee_by, "value"), drop = FALSE],
    y_dat[, c(lee_by, "value"), drop = FALSE],
    by = lee_by,
    suffixes = c("_x", "_y")
  )
  dat <- dat[is.finite(dat$value_x) & is.finite(dat$value_y), , drop = FALSE]
  .check_spatial_values(dat$value_x)
  .check_spatial_values(dat$value_y)

  groups <- .sirup_spatial_stat_groups(dat)
  if (length(groups) > 1L) {
    out <- lapply(groups, function(group_dat) {
      .check_spatial_values(group_dat$value_x)
      .check_spatial_values(group_dat$value_y)
      .sirup_spatial_lee_one(
        dat = group_dat,
        radius = radius,
        k = k,
        alternative = alternative,
        method = method,
        nsim = nsim,
        progress = progress,
        x = x,
        y = y
      )
    })
    summary <- dplyr::bind_rows(lapply(names(out), function(name) {
      data.frame(
        image = name,
        statistic = unname(out[[name]]$test$statistic),
        p_value = out[[name]]$test$p.value,
        lee_l = out[[name]]$L,
        method = out[[name]]$method,
        stringsAsFactors = FALSE
      )
    }))
    return(list(
      tests = lapply(out, `[[`, "test"),
      summary = summary,
      data = dat,
      neighbors = lapply(out, `[[`, "neighbors"),
      weights = lapply(out, `[[`, "weights"),
      radius = radius,
      k = k,
      x = x,
      y = y,
      x_type = x_type,
      y_type = y_type
    ))
  }

  .sirup_spatial_lee_one(
    dat = dat,
    radius = radius,
    k = k,
    alternative = alternative,
    method = method,
    nsim = nsim,
    progress = progress,
    x = x,
    y = y,
    x_type = x_type,
    y_type = y_type
  )
}

.sirup_spatial_lee_one <- function(dat,
                                   radius,
                                   k,
                                   alternative,
                                   method,
                                   nsim,
                                   progress,
                                   x,
                                   y,
                                   x_type = NULL,
                                   y_type = NULL) {
  nb <- .sirup_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  if (method == "mc") {
    # Permute the labels and recompute L, keeping the weights sparse. This is
    # the loop from spdep::lee.mc(), run here so the permutations can report
    # progress -- spdep gives no hook into its own loop. test-spatial-
    # statistics.R checks this returns results identical to spdep::lee.mc()
    # under the same seed.
    n_obs <- length(listw$neighbours)
    S2 <- sum((unlist(lapply(listw$weights, sum)))^2)
    res <- numeric(nsim + 1L)
    if (progress) {
      cli::cli_progress_bar("Lee permutations", total = nsim, .envir = environment())
    }
    for (i in seq_len(nsim)) {
      idx <- sample(1:n_obs)
      res[i] <- spdep::lee(dat$value_x[idx], dat$value_y[idx], listw,
                           n_obs, S2, zero.policy = TRUE)$L
      if (progress) cli::cli_progress_update(.envir = environment())
    }
    if (progress) cli::cli_progress_done(.envir = environment())
    res[nsim + 1L] <- spdep::lee(dat$value_x, dat$value_y, listw,
                                 n_obs, S2, zero.policy = TRUE)$L

    xrank <- rank(res)[length(res)]
    diff <- max(nsim - xrank, 0)
    pval <- switch(
      alternative,
      less = stats::punif((diff + 1) / (nsim + 1), lower.tail = FALSE),
      greater = stats::punif((diff + 1) / (nsim + 1)),
      stats::punif(abs(xrank - (nsim + 1) / 2) / (nsim + 1), 0, 0.5,
                   lower.tail = FALSE)
    )
    statistic <- stats::setNames(res[nsim + 1L], "statistic")
    test <- structure(
      list(statistic = statistic,
           parameter = stats::setNames(xrank, "observed rank"),
           p.value = pval,
           alternative = alternative,
           method = "Monte-Carlo simulation of Lee's L",
           data.name = paste0(x, ", ", y, "\nnumber of simulations + 1: ",
                              nsim + 1, "\n"),
           res = res),
      class = c("htest", "mc.sim")
    )
    lee_l <- unname(statistic)
  } else {
    test <- spdep::lee.test(
      dat$value_x,
      dat$value_y,
      listw,
      zero.policy = TRUE,
      alternative = alternative
    )
    lee_l <- unname(test$estimate[["Lee's L statistic"]])
  }

  list(
    test = test,
    L = lee_l,
    method = method,
    data = dat,
    neighbors = nb,
    weights = listw,
    radius = radius,
    k = k,
    x = x,
    y = y,
    x_type = x_type,
    y_type = y_type
  )
}

#' Calculate local Moran's I (LISA)
#'
#' `spatialLocalMoran` calculates local Moran's I for each spot. It shows where
#' high or low values cluster.
#'
#' @details
#' The function builds one neighbor graph from all spots, with `radius` or `k`
#' as in [spatialMoran()]. Set `variable_type = "clone"` to use a 0/1 value for each spot: 1 if the spot
#' contains the clone, 0 if not.
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds one graph and runs one test for each sample or image.
#'
#' @param object Seurat object or SIRUP object.
#' @param variable Name of the feature or the metadata column to test.
#' @param variable_type Type of `variable`: `feature` (gene expression),
#'   `metadata` (Seurat metadata column), `receptor` (receptor table column),
#'   or `clone` (clone ID).
#' @param value Metadata value to convert to a 0/1 indicator (default `NULL`).
#' @param receptors AIRR-style receptor table. You must give it if
#'   `variable_type` is `"receptor"` or `"clone"` and `object` is not a SIRUP
#'   object.
#' @param spatial_id Name of the `receptors` column that holds the spatial IDs.
#'   If `object` is a SIRUP object and this argument is `NULL`, the function uses
#'   the stored spatial ID field. Give a value if you use separate Seurat and
#'   receptor objects.
#' @param clone Name of the `receptors` column that holds the clone IDs. The
#'   function uses it only if `variable_type` is `"clone"`. If `object` is a
#'   SIRUP object and this argument is `NULL`, the function uses the stored
#'   clone field, if there is one.
#' @param assay Assay to use if `variable_type` is `"feature"`.
#' @param layer Seurat expression layer to use for features (default
#'   `"counts"`).
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param radius Distance that defines a neighborhood. If `NULL`, the function
#'   uses the `k` nearest neighbors.
#' @param k Number of nearest neighbors, used if `radius` is `NULL` (default 6).
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots and the function reads receptors. The function ignores `slot` if you
#'   give a plain receptor table.
#'
#' @return A data frame with the coordinates, the values, local Moran's I, the
#'   z-scores, and the p-values.
#' @seealso [plotSpatialLocalMoran()], [spatialMoran()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' local_moran <- spatialLocalMoran(x, "gene1", k = 4)
#' head(local_moran[, c("cell", "local_i", "z_i", "p_value")])
#' @export
spatialLocalMoran <- function(object,
                              variable,
                              variable_type = c("feature", "metadata", "receptor", "clone"),
                              value = NULL,
                              receptors = NULL,
                              spatial_id = NULL,
                              clone = NULL,
                              assay = NULL,
                              layer = "counts",
                              images = NULL,
                              image_scale = NULL,
                              swap_coords = NULL,
                              radius = NULL,
                              k = 6,
                              slot = NULL) {

  dat <- .sirup_spatial_variable_values(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    spatial_id_col = spatial_id,
    clone_col = clone,
    assay = assay,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    slot = slot
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)
  groups <- .sirup_spatial_stat_groups(dat)
  if (length(groups) > 1L) {
    out <- lapply(groups, function(group_dat) {
      .check_spatial_values(group_dat$value)
      .sirup_spatial_local_moran_one(group_dat, radius = radius, k = k)
    })
    return(dplyr::bind_rows(out))
  }

  .sirup_spatial_local_moran_one(dat, radius = radius, k = k)
}

.sirup_spatial_local_moran_one <- function(dat, radius, k) {
  nb <- .sirup_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  lisa <- spdep::localmoran(dat$value, listw, zero.policy = TRUE)
  lisa <- as.data.frame(lisa)

  out <- cbind(dat, lisa)
  names(out)[names(out) == "Ii"] <- "local_i"
  names(out)[names(out) == "E.Ii"] <- "expected_i"
  names(out)[names(out) == "Var.Ii"] <- "variance_i"
  names(out)[names(out) == "Z.Ii"] <- "z_i"
  names(out)[grepl("^Pr", names(out))] <- "p_value"
  out$radius <- radius
  out$k <- if (is.null(radius)) k else NA_integer_
  out
}

.sirup_spatial_stat_groups <- function(dat) {
  if (!"image" %in% colnames(dat)) {
    return(list(all = dat))
  }
  groups <- split(dat, as.character(dat$image), drop = TRUE)
  groups <- Filter(function(x) nrow(x) >= 2L, groups)
  if (length(groups) == 0L) {
    rlang::abort("At least one spatial image with two points is required.")
  }
  groups
}

#' Test clone localization with join counts
#'
#' `spatialCloneJoinCount` tests if each clone clusters in space.
#'
#' @details
#' The function calls `spdep::joincount.test()` for each clone. The test asks if
#' neighboring spots with receptors share the clone more often than random
#' labels would give. 
#'
#' The function skips a clone if it occurs in fewer than `min_spatial_ids`
#' spatial IDs, or in every spatial ID of the region.
#'
#' If `object` is a SIRUP object with stored polygons and you give no `polygon`,
#' the function runs one test for each stored polygon. Set `use_stored_polygons` to FALSE
#' to change this.
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds one graph and runs one test for each sample or image.
#'
#' @param object Seurat object or SIRUP object.
#' @param receptors AIRR-style receptor table with spatial IDs and clone IDs.
#'   Omit if `object` is a SIRUP object.
#' @param clone Name of the `receptors` column that holds the clone IDs. If
#'   `object` is a SIRUP object and this argument is `NULL`, the function uses
#'   the stored clone field, if there is one.
#' @param spatial_id Name of the `receptors` column that holds the spatial IDs.
#'   If `object` is a SIRUP object and this argument is `NULL`, the function uses
#'   the stored spatial ID field. Give a value if you use separate Seurat and
#'   receptor objects.
#' @param clone_val Clone IDs to test. If `NULL`, the function tests all clones
#'   that occur in at least `min_spatial_ids` spatial IDs.
#' @param polygon Data frame of polygon vertices. The function tests only the
#'   spots inside these polygons (default `NULL`).
#' @param polygon_val Polygon IDs to select from `polygon`.
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param include_boundary Count spots on a polygon edge as inside the polygon
#'   (default `TRUE`)? The function ignores this argument if `object` is a
#'   SIRUP object with stored polygons and `polygon` is `NULL`. In that case,
#'   it uses the stored polygon annotations.
#' @param min_spatial_ids Smallest number of spatial IDs that a clone needs to
#'   be tested (default 2).
#' @param radius Distance that defines a neighborhood. If `NULL`, the function
#'   uses `k`.
#' @param k Number of nearest neighbors, used if `radius` is `NULL` (default 6).
#' @param alternative Alternative hypothesis for `spdep::joincount.test()`
#'   (default `"greater"`).
#' @param other_label Label for the spatial IDs that do not contain the tested
#'   clone (default `"Other"`).
#' @param use_stored_polygons Run one test for each stored polygon of a SIRUP
#'   object? Use `TRUE` to require stored polygons. Use `FALSE` to build one
#'   neighbor graph for each sample and ignore stored polygons. If `NULL`
#'   (default), the function uses `TRUE` when a SIRUP object has stored polygons
#'   and you give no `polygon`. Otherwise, it uses `FALSE`.
#' @param return_data Return the data used for each test (default `FALSE`)?
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots and the function reads receptors. The function ignores `slot` if you
#'   give a plain receptor table.
#'
#' @return A list with these items: `summary` (a table with one row for each
#'   tested clone and region, ordered by p-value), `tests` (the `spdep` test
#'   results), `data` (the data for each test,
#'   or `NULL` if `return_data` is `FALSE`), `neighbors`, `weights`, `radius`,
#'   `k`, `alternative`, and `min_spatial_ids`.
#' @seealso [plotSpatialCloneJoinCount()], [spatialMoran()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Test each clone on the whole tissue.
#' join_count <- spatialCloneJoinCount(
#'   x,
#'   use_stored_polygons = FALSE,
#'   k = 4,
#'   return_data = TRUE
#' )
#' join_count$summary[, c("clone_id", "observed_join_count", "expected_join_count", "p_value")]
#' @export
spatialCloneJoinCount <- function(object,
                                  receptors = NULL,
                                  clone = NULL,
                                  spatial_id = NULL,
                                  clone_val = NULL,
                                  polygon = NULL,
                                  polygon_val = NULL,
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  x = "x",
                                  y = "y",
                                  polygon_id = "polygon",
                                  include_boundary = TRUE,
                                  min_spatial_ids = 2,
                                  radius = NULL,
                                  k = 6,
                                  alternative = "greater",
                                  other_label = "Other",
                                  use_stored_polygons = NULL,
                                  return_data = FALSE,
                                  slot = NULL) {
  sirup_container <- if (isSIRUP(object)) object else NULL
  auto_stored_polygons <- !is.null(sirup_container) &&
    is.null(polygon) &&
    nrow(.sirup_polygons(sirup_container)) > 0L
  if (is.null(use_stored_polygons)) {
    use_stored_polygons <- auto_stored_polygons
  } else {
    use_stored_polygons <- isTRUE(use_stored_polygons)
    if (use_stored_polygons &&
        (is.null(sirup_container) || nrow(.sirup_polygons(sirup_container)) == 0L)) {
      rlang::abort("`use_stored_polygons = TRUE` requires a SIRUP object with stored polygons.")
    }
  }
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  inputs <- .sirup_resolve_inputs(
    object,
    receptors,
    spatial_id,
    slot
  )
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id <- inputs$spatial_id_col
  sample_col <- .sirup_receptor_sample_col(receptors, inputs$receptor_fields, object)
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRUP object.")
  }
  if (is.null(clone) && !is.null(inputs$receptor_fields$clone)) {
    clone <- inputs$receptor_fields$clone
  }
  if (is.null(spatial_id)) {
    rlang::abort("`spatial_id` is required when it cannot be resolved from a SIRUP object.")
  }
  if (is.null(clone)) {
    rlang::abort("`clone` is required when it cannot be resolved from a SIRUP object.")
  }
  if (!is.numeric(min_spatial_ids) || length(min_spatial_ids) != 1L || is.na(min_spatial_ids) || min_spatial_ids < 1) {
    rlang::abort("`min_spatial_ids` must be a positive number.")
  }
  min_spatial_ids <- as.integer(min_spatial_ids)

  dat <- if (isTRUE(use_stored_polygons)) {
    .sirup_joincount_receptor_data_stored_polygons(
      object = sirup_container,
      receptors = receptors,
      clone_col = clone,
      spatial_id_col = spatial_id,
      polygon_id = polygon_val,
      images = images,
      image_scale = image_scale,
      swap_coords = swap_coords,
      sample_col = sample_col
    )
  } else {
    .sirup_joincount_receptor_data(
    object = object,
    receptors = receptors,
    clone_col = clone,
    spatial_id_col = spatial_id,
    sample_col = sample_col,
    polygon = polygon,
    polygon_id = polygon_val,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    x_col = x,
    y_col = y,
    polygon_col = polygon_id,
    include_boundary = include_boundary
    )
  }
  .check_required_cols(dat, c(clone, ".sirup_spatial_id", "x", ".sirup_y", "image", "polygon"), "receptors")
  dat <- dat[!is.na(dat[[clone]]) & nzchar(as.character(dat[[clone]])), , drop = FALSE]
  if (nrow(dat) == 0L) {
    rlang::abort("No receptor rows with clone IDs were found in the selected region(s).")
  }

  if (is.null(clone_val)) {
    clone_key_cols <- c(.sirup_spatial_key_cols(dat), clone)
    clone_counts <- sort(table(unique(dat[, clone_key_cols, drop = FALSE])[[clone]]), decreasing = TRUE)
    clone_val <- names(clone_counts)[clone_counts >= min_spatial_ids]
  } else {
    clone_val <- as.character(clone_val)
  }
  clone_val <- unique(clone_val[!is.na(clone_val) & nzchar(clone_val)])
  if (length(clone_val) == 0L) {
    rlang::abort("No clones met `min_spatial_ids` in the selected region(s).")
  }

  tests <- list()
  summaries <- list()
  mark_data <- list()
  neighbors <- list()
  weights <- list()
  groups <- split(dat, interaction(dat$image, dat$polygon, drop = TRUE, sep = "||"))
  for (group_name in names(groups)) {
    group_dat <- groups[[group_name]]
    group_key <- strsplit(group_name, "\\|\\|", fixed = FALSE)[[1L]]
    image <- group_key[[1L]]
    polygon_value <- group_key[[2L]]
    site_cols <- unique(c(.sirup_spatial_key_cols(group_dat), "x", "y", ".sirup_y", "image", "polygon"))
    sites <- group_dat[, site_cols, drop = FALSE]
    sites <- sites[!duplicated(.sirup_key_vector(sites)), , drop = FALSE]
    if (nrow(sites) < 2L) {
      next
    }
    nb <- .sirup_spatial_neighbors(sites, radius = radius, k = k)
    listw <- spdep::nb2listw(nb, style = "B", zero.policy = TRUE)
    neighbors[[group_name]] <- nb
    weights[[group_name]] <- listw

    for (clone_id in clone_val) {
      clone_sites <- unique(group_dat[as.character(group_dat[[clone]]) == clone_id, .sirup_spatial_key_cols(group_dat), drop = FALSE])
      clone_keys <- .sirup_key_vector(clone_sites)
      n_clone_spatial_ids <- length(clone_keys)
      if (n_clone_spatial_ids < min_spatial_ids || n_clone_spatial_ids == nrow(sites)) {
        next
      }
      local_other_label <- if (identical(clone_id, other_label)) paste0(other_label, "_background") else other_label
      sites$join_count_mark <- factor(
        ifelse(.sirup_key_vector(sites) %in% clone_keys, clone_id, local_other_label),
        levels = c(local_other_label, clone_id)
      )
      test <- spdep::joincount.test(
        sites$join_count_mark,
        listw,
        zero.policy = TRUE,
        alternative = alternative
      )
      test_name <- paste(image, polygon_value, clone_id, sep = "||")
      tests[[test_name]] <- test
      summaries[[test_name]] <- .sirup_joincount_summary_row(
        test = test,
        clone_id = clone_id,
        image = image,
        polygon = polygon_value,
        n_spatial_ids = nrow(sites),
        n_clone_spatial_ids = n_clone_spatial_ids,
        radius = radius,
        k = if (is.null(radius)) k else NA_integer_,
        alternative = alternative
      )
      if (isTRUE(return_data)) {
        mark_data[[test_name]] <- cbind(
          sites,
          data.frame(clone_id = clone_id, stringsAsFactors = FALSE)
        )
      }
    }
  }

  summary <- dplyr::bind_rows(summaries)
  if (nrow(summary) == 0L) {
    rlang::abort("No clones could be tested after applying `min_spatial_ids` and region grouping.")
  }
  summary <- summary[order(summary$p_value, summary$image, summary$polygon, summary$clone_id), , drop = FALSE]
  rownames(summary) <- NULL

  list(
    summary = summary,
    tests = tests,
    data = if (isTRUE(return_data)) dplyr::bind_rows(mark_data) else NULL,
    neighbors = neighbors,
    weights = weights,
    radius = radius,
    k = k,
    alternative = alternative,
    min_spatial_ids = min_spatial_ids
  )
}

#' Plot a clone join-count result
#'
#' `plotSpatialCloneJoinCount` plots the tissue map of one tested clone. It also
#' shows the observed and expected number of joins between spots of the same
#' clone.
#'
#' @details
#' The map uses the same spots and the same neighbor graph as the call to
#' [spatialCloneJoinCount()].
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds a list of plots, one for each sample or image.
#'
#' @param join_count Result of [spatialCloneJoinCount()]. Run that function with
#'   `return_data = TRUE`.
#' @param object Seurat object with spatial images, or a SIRUP object.
#' @param clone_val Clone ID to plot.
#' @param images Names of the images to plot.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param polygon_val IDs of the polygons to include.
#' @param polygon Data frame of polygon vertices to draw on the tissue.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object.
#' @param show_neighbors Draw a line between each pair of neighboring spots that
#'   both contain the clone (default `TRUE`)?
#' @param clone_color Color of the spots and lines for the clone.
#' @param other_color Color of the other spots that have receptors.
#' @param point_size Size of the spots.
#' @param point_alpha Opacity of the spots.
#' @param edge_width Width of the neighbor lines.
#' @param edge_alpha Opacity of the neighbor lines.
#' @param polygon_color Color of the polygon outlines.
#' @param polygon_linewidth Width of the polygon outlines.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the spots in the test (default `FALSE`)?
#' @param title Title of the combined plot.
#' @param combine_sample_plots Combine the plots of all samples into one figure
#'   (default `FALSE`)? The function uses it only for multi-sample objects. If
#'   `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot object. For a multi-sample object and
#'   `combine_sample_plots = FALSE`, a named list of plots.
#' @seealso [spatialCloneJoinCount()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' join_count <- spatialCloneJoinCount(
#'   x,
#'   use_stored_polygons = FALSE,
#'   k = 4,
#'   return_data = TRUE
#' )
#' plotSpatialCloneJoinCount(join_count, x, clone_val = "2")
#' @export
plotSpatialCloneJoinCount <- function(join_count,
                                      object,
                                      clone_val,
                                      images = NULL,
                                      sample_id = NULL,
                                      polygon_val = NULL,
                                      polygon = NULL,
                                      image_scale = NULL,
                                      swap_coords = NULL,
                                      show_neighbors = TRUE,
                                      clone_color = "#D55E00",
                                      other_color = "grey75",
                                      point_size = 1.8,
                                      point_alpha = 0.9,
                                      edge_width = 0.6,
                                      edge_alpha = 0.75,
                                      polygon_color = "black",
                                      polygon_linewidth = 0.6,
                                      image_alpha = 1,
                                      zoom = FALSE,
                                      title = NULL,
                                      combine_sample_plots = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  images <- .sirup_spatial_images(object, images = images, sample_id = sample_id)
  if (!is.list(join_count) || is.null(join_count$summary)) {
    rlang::abort("`join_count` must be output from `spatialCloneJoinCount()`.")
  }
  if (is.null(join_count$data)) {
    rlang::abort(
      "`join_count$data` is missing. Re-run `spatialCloneJoinCount()` with `return_data = TRUE`."
    )
  }
  if (length(clone_val) != 1L || is.na(clone_val) || !nzchar(as.character(clone_val))) {
    rlang::abort("`clone_val` must be one non-missing clone ID.")
  }
  clone_val <- as.character(clone_val)

  summary_data <- .as_data_frame(join_count$summary)
  mark_data <- .as_data_frame(join_count$data)
  .check_required_cols(
    summary_data,
    c(
      "clone_id", "image", "polygon", "n_spatial_ids", "n_clone_spatial_ids",
      "observed_join_count", "expected_join_count", "statistic", "p_value"
    ),
    "join_count$summary"
  )
  .check_required_cols(
    mark_data,
    c("clone_id", "image", "polygon", ".sirup_spatial_id", "x", ".sirup_y", "join_count_mark"),
    "join_count$data"
  )

  summary_data <- summary_data[as.character(summary_data$clone_id) == clone_val, , drop = FALSE]
  mark_data <- mark_data[as.character(mark_data$clone_id) == clone_val, , drop = FALSE]
  if (!is.null(images)) {
    summary_data <- summary_data[summary_data$image %in% images, , drop = FALSE]
    mark_data <- mark_data[mark_data$image %in% images, , drop = FALSE]
  }
  if (!is.null(polygon_val)) {
    polygon_val <- as.character(polygon_val)
    summary_data <- summary_data[as.character(summary_data$polygon) %in% polygon_val, , drop = FALSE]
    mark_data <- mark_data[as.character(mark_data$polygon) %in% polygon_val, , drop = FALSE]
  }
  if (nrow(summary_data) == 0L || nrow(mark_data) == 0L) {
    rlang::abort("The selected clone was not found in the requested join-count result.")
  }

  summary_data$.sirup_group <- paste(summary_data$image, summary_data$polygon, sep = "||")
  mark_data$.sirup_group <- paste(mark_data$image, mark_data$polygon, sep = "||")
  keep_groups <- intersect(summary_data$.sirup_group, unique(mark_data$.sirup_group))
  summary_data <- summary_data[summary_data$.sirup_group %in% keep_groups, , drop = FALSE]
  mark_data <- mark_data[mark_data$.sirup_group %in% keep_groups, , drop = FALSE]
  summary_data$.sirup_label <- ifelse(
    summary_data$polygon == "all",
    summary_data$image,
    paste(summary_data$image, summary_data$polygon, sep = ": ")
  )
  image_samples <- vapply(unique(as.character(summary_data$image)), function(image) {
    parsed <- .sirup_parse_spatial_image_label(image)
    if (is.null(parsed)) NA_character_ else parsed$sample
  }, character(1))
  plot_samples <- unique(stats::na.omit(image_samples))
  if (length(plot_samples) > 1L && !isTRUE(combine_sample_plots)) {
    out <- lapply(plot_samples, function(sample_name) {
      sample_images <- names(image_samples)[image_samples == sample_name]
      plotSpatialCloneJoinCount(
        join_count = join_count,
        object = object,
        clone_val = clone_val,
        images = sample_images,
        sample_id = sample_name,
        polygon_val = polygon_val,
        polygon = polygon,
        image_scale = image_scale,
        swap_coords = swap_coords,
        show_neighbors = show_neighbors,
        clone_color = clone_color,
        other_color = other_color,
        point_size = point_size,
        point_alpha = point_alpha,
        edge_width = edge_width,
        edge_alpha = edge_alpha,
        polygon_color = polygon_color,
        polygon_linewidth = polygon_linewidth,
        image_alpha = image_alpha,
        zoom = zoom,
        title = title,
        combine_sample_plots = TRUE
      )
    })
    names(out) <- plot_samples
    return(out)
  }

  polygon_data <- NULL
  if (!is.null(polygon)) {
    polygon_data <- .sirup_polygon_plot_data(
      polygon = .as_data_frame(polygon),
      polygon_id = polygon_val,
      polygon_col = "polygon",
      x_col = "x",
      y_col = "y"
    )
  }

  spatial_plots <- lapply(keep_groups, function(group_name) {
    group_data <- mark_data[mark_data$.sirup_group == group_name, , drop = FALSE]
    image_name <- as.character(group_data$image[[1L]])
    polygon_name <- as.character(group_data$polygon[[1L]])
    group_data$.sirup_selected <- as.character(group_data$join_count_mark) == clone_val

    spatial <- tryCatch(
      .sirup_spatial_data(
        object,
        image = image_name,
        image_scale = image_scale,
        swap_coords = swap_coords
      ),
      error = function(e) NULL
    )
    if (!is.null(spatial)) {
      p <- .sirup_spatial_base(
        plot_data = group_data,
        image = spatial$image,
        image_alpha = image_alpha,
        zoom = zoom,
        view_data = spatial$coords
      )
    } else {
      view_data <- .sirup_spatial_coordinate_table(
        object,
        images = image_name,
        image_scale = image_scale,
        swap_coords = swap_coords
      )
      p <- ggplot2::ggplot(group_data, ggplot2::aes(x = .data$x, y = .data$.sirup_y)) +
        ggplot2::theme_void()
      if (!isTRUE(zoom) && nrow(view_data) > 0L) {
        p <- p + ggplot2::coord_fixed(
          xlim = range(view_data$x, na.rm = TRUE),
          ylim = range(view_data$.sirup_y, na.rm = TRUE),
          expand = FALSE
        )
      } else {
        p <- p + ggplot2::coord_fixed()
      }
    }

    if (isTRUE(show_neighbors)) {
      nb <- join_count$neighbors[[group_name]]
      if (!is.null(nb)) {
        edge_data <- .sirup_joincount_clone_edges(group_data, nb)
        if (nrow(edge_data) > 0L) {
          p <- p + ggplot2::geom_segment(
            data = edge_data,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            inherit.aes = FALSE,
            color = clone_color,
            linewidth = edge_width,
            alpha = edge_alpha
          )
        }
      }
    }
    p <- p +
      ggplot2::geom_point(color = other_color, size = point_size, alpha = point_alpha) +
      ggplot2::geom_point(
        data = group_data[group_data$.sirup_selected, , drop = FALSE],
        color = clone_color,
        size = point_size,
        alpha = point_alpha
      )
    if (!is.null(polygon_data)) {
      this_polygon <- polygon_data
      if (polygon_name != "all") {
        this_polygon <- this_polygon[as.character(this_polygon$polygon) == polygon_name, , drop = FALSE]
      }
      if (nrow(this_polygon) > 0L) {
        p <- p + ggplot2::geom_path(
          data = this_polygon,
          ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$polygon),
          inherit.aes = FALSE,
          color = polygon_color,
          linewidth = polygon_linewidth
        )
      }
    }
    p + ggplot2::labs(title = summary_data$.sirup_label[match(group_name, summary_data$.sirup_group)])
  })
  spatial_panel <- cowplot::plot_grid(plotlist = spatial_plots, ncol = 1, align = "v")

  count_data <- dplyr::bind_rows(
    data.frame(
      .sirup_group = summary_data$.sirup_group,
      .sirup_label = summary_data$.sirup_label,
      count_type = "Expected",
      join_count = summary_data$expected_join_count,
      stringsAsFactors = FALSE
    ),
    data.frame(
      .sirup_group = summary_data$.sirup_group,
      .sirup_label = summary_data$.sirup_label,
      count_type = "Observed",
      join_count = summary_data$observed_join_count,
      stringsAsFactors = FALSE
    )
  )
  count_data$count_type <- factor(count_data$count_type, levels = c("Expected", "Observed"))
  annotation_data <- summary_data[, c(".sirup_group", ".sirup_label"), drop = FALSE]
  annotation_data$.sirup_annotation <- paste0(
    "z = ", formatC(summary_data$statistic, digits = 2, format = "f"),
    "\np = ", format.pval(summary_data$p_value, digits = 2, eps = 0.001),
    "\n", summary_data$n_clone_spatial_ids, "/", summary_data$n_spatial_ids, " clone-positive spots"
  )
  annotation_data$x <- 1.5
  annotation_data$y <- max(c(count_data$join_count, 0), na.rm = TRUE) * 1.12

  count_panel <- ggplot2::ggplot(
    count_data,
    ggplot2::aes(x = .data$count_type, y = .data$join_count, fill = .data$count_type)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = signif(.data$join_count, 3)),
      vjust = -0.4,
      size = 3.5
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$.sirup_annotation),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = c("Expected" = "grey70", "Observed" = clone_color), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.35))) +
    ggplot2::labs(x = NULL, y = "Same-clone joins") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (nrow(summary_data) > 1L) {
    count_panel <- count_panel + ggplot2::facet_wrap(ggplot2::vars(.data$.sirup_label), ncol = 1)
  }

  combined <- cowplot::plot_grid(
    spatial_panel,
    count_panel,
    nrow = 1,
    rel_widths = c(1.7, 1),
    align = "h",
    axis = "tb"
  )
  plot_title <- title %||% paste("Spatial join count for clone", clone_val)
  cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(plot_title, fontface = "bold", x = 0.5, hjust = 0.5),
    combined,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
}

# Convert a join-count neighbor graph into line segments connecting selected
# clone-positive neighboring spots for plotting.
.sirup_joincount_clone_edges <- function(group_data, neighbors) {
  if (length(neighbors) != nrow(group_data)) {
    rlang::abort("The stored neighbor graph does not match `join_count$data`.")
  }
  selected <- group_data$.sirup_selected
  edges <- vector("list", length(neighbors))
  edge_index <- 1L
  for (i in seq_along(neighbors)) {
    if (!isTRUE(selected[[i]])) {
      next
    }
    for (j in neighbors[[i]]) {
      if (j > i && isTRUE(selected[[j]])) {
        edges[[edge_index]] <- data.frame(
          x = group_data$x[[i]],
          y = group_data$.sirup_y[[i]],
          xend = group_data$x[[j]],
          yend = group_data$.sirup_y[[j]]
        )
        edge_index <- edge_index + 1L
      }
    }
  }
  dplyr::bind_rows(edges)
}

#' Plot local Moran's I on a tissue image
#'
#' `plotSpatialLocalMoran` colors each spot by its local Moran's I result.
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds a list of plots, one for each sample or image.
#'
#' @param object SIRUP object, or Seurat object with spatial images.
#' @param local_moran Result of [spatialLocalMoran()].
#' @param value_col Name of the column to color by (default `"z_i"`).
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param colors Three colors for the low, middle, and high values.
#' @param midpoint Value that gets the middle color.
#' @param point_size Size of the spots.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw low values first, so that high values show on top
#'   (default `TRUE`)?
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)?
#' @param title Plot title.
#' @param combine Combine the image panels with patchwork (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot or patchwork object.
#' @seealso [spatialLocalMoran()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' local_moran <- spatialLocalMoran(x, "gene1", k = 4)
#' plotSpatialLocalMoran(x, local_moran)
#' @export
plotSpatialLocalMoran <- function(object,
                                  local_moran,
                                  value_col = "z_i",
                                  images = NULL,
                                  sample_id = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  colors = c("#2c7bb6", "#f7f7f7", "#d7191c"),
                                  midpoint = 0,
                                  point_size = 0.7,
                                  point_alpha = 1,
                                  order_points = TRUE,
                                  image_alpha = 1,
                                  zoom = FALSE,
                                  title = NULL,
                                  combine = TRUE,
                                  combine_sample_plots = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .check_required_cols(local_moran, c("cell", "x", ".sirup_y", value_col), "local_moran")
  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plot_data <- local_moran[local_moran$image == image, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }
    if (isTRUE(order_points)) {
      plot_data <- plot_data[order(plot_data[[value_col]], na.last = TRUE), , drop = FALSE]
    }
    plots[[image]] <- .sirup_spatial_base(
      plot_data,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[value_col]]),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_gradient2(
        low = colors[[1]],
        mid = colors[[2]],
        high = colors[[3]],
        midpoint = midpoint,
        name = value_col
      ) +
      ggplot2::labs(title = title %||% paste("Local Moran's I z-score:", unique(plot_data$variable)[[1]]))
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots)
}

#' Calculate Ripley's cross K
#'
#' `spatialRipleyCrossK` measures how spots of mark `i` cluster around spots of
#' mark `j` at different distances. Similar to Lee's L but better for factor data.
#'
#' @details
#' The function calls `spatstat.explore::Kcross()` with the border correction.
#' If `nsim` is greater than 0, it also builds an envelope from random label
#' simulations.
#' 
#' If the data come from more than one image, or a sirup object with multiple samples,
#' then the function builds one graph and runs one test for each sample or image.
#'
#' @param object Seurat object or SIRUP object.
#' @param mark_col Name of the metadata column that holds the point types.
#' @param i First mark level.
#' @param j Second mark level.
#' @param nsim Number of random-label simulations for the envelope (default 0).
#'   Use `0` for the observed K only.
#' @param max_r Largest radius to calculate. The function ignores it if you give
#'   `r`.
#' @param r Vector of distances for `spatstat.explore::Kcross()`.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#'
#' @return A list with these items: `observed` (the observed K), `envelope`
#'   (the simulation envelope, if any), `plot_data`, `points` (the point data),
#'   `i`, `j`, `correction`, `nsim`, `max_r`, and `r`.
#' @seealso [plotRipleyCrossK()]
#'
#' @examplesIf requireNamespace("spatstat.explore", quietly = TRUE)
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Compare spots of type A with spots of type B. Use 19 random-label
#' # simulations for the envelope.
#' set.seed(1)
#' kcross <- spatialRipleyCrossK(
#'   x,
#'   mark_col = "cell_type",
#'   i = "A",
#'   j = "B",
#'   nsim = 19
#' )
#' head(kcross$plot_data)
#' @export
spatialRipleyCrossK <- function(object,
                                mark_col,
                                i,
                                j,
                                nsim = 0,
                                max_r = NULL,
                                r = NULL,
                                image_scale = NULL,
                                swap_coords = NULL) {
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .require_namespace("spatstat.geom")
  .require_namespace("spatstat.explore")
  if (nsim > 0L) {
    .require_namespace("spatstat.random")
  }

  # Get the spatial coordinates
  dat <- .sirup_spatial_coordinate_table(object, image_scale = image_scale, swap_coords = swap_coords)
  marks <- .sirup_spatial_metadata_values(object, mark_col)
  dat <- .sirup_merge_spatial_values(dat, marks)
  dat$mark <- dat[[mark_col]]
  dat <- dat[!is.na(dat$mark), , drop = FALSE]
  if (!all(c(i, j) %in% dat$mark)) {
    rlang::abort("Both `i` and `j` must be present in the marked spatial data.")
  }
  dat$mark <- factor(dat$mark)

  observed <- list()
  envelopes <- list()
  plot_data <- list()
  for (image in unique(dat$image)) {
    image_data <- dat[dat$image == image, , drop = FALSE]
    window <- spatstat.geom::owin(
      xrange = range(image_data$x, na.rm = TRUE),
      yrange = range(image_data$.sirup_y, na.rm = TRUE)
    )
    pp <- spatstat.geom::ppp(
      x = image_data$x,
      y = image_data$.sirup_y,
      marks = image_data$mark,
      window = window,
      checkdup = FALSE
    )
    r_eval <- r
    if (is.null(r_eval) && !is.null(max_r)) {
      r_eval <- seq(0, max_r, length.out = 513L)
    }
    k_obs <- spatstat.explore::Kcross(pp, i = i, j = j, correction = "border", r = r_eval)
    observed[[image]] <- k_obs

    if (nsim > 0L) {
      envelopes[[image]] <- spatstat.explore::envelope(
        pp,
        fun = spatstat.explore::Kcross,
        nsim = nsim,
        simulate = expression(spatstat.random::rlabel(pp)),
        i = i,
        j = j,
        correction = "border",
        r = r_eval,
        # Retain the individual simulated functions so the returned envelope
        # objects can be combined across samples with spatstat.explore::pool().
        # This does not change the summary (observed/simulated/lo/hi); it only
        # keeps the per-simulation curves needed for a pooled group envelope.
        savefuns = TRUE
      )
      plot_data[[image]] <- .sirup_kcross_df(envelopes[[image]], image = image)
    } else {
      plot_data[[image]] <- .sirup_kcross_df(k_obs, image = image)
    }
  }

  list(
    observed = observed,
    envelope = envelopes,
    plot_data = dplyr::bind_rows(plot_data),
    points = dat,
    i = i,
    j = j,
    correction = "border",
    nsim = nsim,
    max_r = max_r,
    r = r
  )
}

#' Plot Ripley's cross K against distance
#'
#' `plotRipleyCrossK` plots the observed cross K, with the simulation envelope if
#' there is one.
#'
#' @param kcross Result of [spatialRipleyCrossK()].
#' @param title Plot title.
#' @param observed_color Color of the observed curve.
#' @param simulated_color Color of the dotted simulated curve.
#' @param envelope_fill Fill color of the simulation envelope.
#' @param envelope_alpha Opacity of the simulation envelope.
#' @param facet Make one panel for each image if the result has several images
#'   (default `TRUE`)?
#' @param log_y Use a log10 scale for the y-axis (default `FALSE`)?
#'
#' @return A ggplot object.
#' @seealso [spatialRipleyCrossK()]
#'
#' @examplesIf requireNamespace("spatstat.explore", quietly = TRUE)
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' set.seed(1)
#' kcross <- spatialRipleyCrossK(
#'   x,
#'   mark_col = "cell_type",
#'   i = "A",
#'   j = "B",
#'   nsim = 19
#' )
#' plotRipleyCrossK(kcross)
#' @export
plotRipleyCrossK <- function(kcross,
                             title = NULL,
                             observed_color = "black",
                             simulated_color = "grey30",
                             envelope_fill = "grey75",
                             envelope_alpha = 0.45,
                             facet = TRUE,
                             log_y = FALSE) {
  dat <- kcross$plot_data
  .check_required_cols(dat, c("r", "observed", "simulated"), "kcross$plot_data")
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$r))
  if (all(c("lo", "hi") %in% colnames(dat))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$lo,
          ymax = .data$hi,
          group = .data$image,
          fill = "Simulation envelope"
        ),
        alpha = envelope_alpha,
        inherit.aes = TRUE
      )
  }
  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$simulated,
        color = "Simulated",
        linetype = "Simulated"
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$observed,
        color = "Observed",
        linetype = "Observed"
      ),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c("Observed" = observed_color, "Simulated" = simulated_color),
      breaks = c("Observed", "Simulated"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Observed" = "solid", "Simulated" = "dotted"),
      breaks = c("Observed", "Simulated"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Simulation envelope" = envelope_fill),
      breaks = "Simulation envelope",
      name = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste("Ripley's cross K:", kcross$i, "to", kcross$j),
      x = "Distance",
      y = "Kcross"
    )
  if (isTRUE(facet) && length(unique(dat$image)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$image))
  }
  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }
  p
}

# Resolve feature, metadata, receptor, or clone values onto spatial coordinates.
.sirup_spatial_variable_values <- function(object,
                                           variable,
                                           variable_type,
                                           value,
                                   receptors,
                                   spatial_id_col,
                                   clone_col,
                                   assay = NULL,
                                   layer,
                                           images,
                                           image_scale,
                                           swap_coords,
                                           slot = NULL) {
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  # Genes and metadata need no receptors. The receptor table only supplies the sample field,
  # which all slots share, so any slot works then.
  needs_receptors <- variable_type %in% c("receptor", "clone")
  if (!needs_receptors && isSIRUP(object) && is.null(slot)) {
    slot <- .sirup_slot_names(object)[[1L]]
  }
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id_col, slot)
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id_col <- inputs$spatial_id_col
  sample_col <- .sirup_receptor_sample_col(receptors %||% data.frame(), inputs$receptor_fields, object)
  if (is.null(clone_col) && !is.null(inputs$receptor_fields$clone)) {
    clone_col <- inputs$receptor_fields$clone
  }
  variable_type <- match.arg(variable_type, c("feature", "metadata", "receptor", "clone"))

  # Start from one row per spatial coordinate so all statistic modes share the
  # same neighbor graph construction.
  dat <- .sirup_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)

  if (variable_type == "feature") {
    # Features are pulled from the requested expression layer and analyzed as
    # quantitative values.
    expr <- .sirup_spatial_fetch_data(object, vars = variable, layer = layer, assay = assay)
    dat <- .sirup_merge_spatial_values(dat, expr)
    dat$value <- dat[[variable]]
    dat$variable <- variable
  } else if (variable_type == "metadata") {
    # Numeric/logical metadata are analyzed directly. Categorical metadata need
    # an explicit level so the statistic is a binary indicator.
    meta <- .sirup_spatial_metadata_values(object, variable)
    dat <- .sirup_merge_spatial_values(dat, meta)
    if (!is.null(value)) {
      dat$value <- as.integer(dat[[variable]] %in% value)
      dat$variable <- paste0(variable, "_", paste(value, collapse = "_"))
    } else if (is.numeric(dat[[variable]]) || is.integer(dat[[variable]]) || is.logical(dat[[variable]])) {
      dat$value <- as.numeric(dat[[variable]])
      dat$variable <- variable
    } else {
      rlang::abort("For non-numeric metadata, provide `value` to analyze a binary indicator.")
    }
  } else if (variable_type == "receptor") {
    # Receptor mode asks whether receptor-positive spots are spatially
    # autocorrelated, independent of clone identity.
    if (is.null(receptors)) {
      rlang::abort("`receptors` is required when `variable_type = 'receptor'`.")
    }
    if (is.null(spatial_id_col)) {
      rlang::abort("`spatial_id` is required when it cannot be resolved from a SIRUP object.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, spatial_id_col, "receptors")
    receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id_col, sample_col)
    dat <- .sirup_add_coord_spatial_keys(dat)
    receptor_keys <- unique(receptors[!is.na(receptors$.sirup_spatial_id), .sirup_spatial_key_cols(receptors), drop = FALSE])
    dat_keys <- dat[, .sirup_spatial_key_cols(dat, receptor_keys), drop = FALSE]
    dat$value <- as.integer(.sirup_key_vector(dat_keys) %in% .sirup_key_vector(receptor_keys))
    dat$variable <- variable
  } else {
    # Clone mode is a clone-specific receptor indicator: spots containing any
    # row from the requested clone(s) are 1 and all other spatial spots are 0.
    if (is.null(receptors)) {
      rlang::abort("`receptors` is required when `variable_type = 'clone'`.")
    }
    if (is.null(spatial_id_col)) {
      rlang::abort("`spatial_id` is required when it cannot be resolved from a SIRUP object.")
    }
    if (is.null(clone_col)) {
      rlang::abort("`clone` is required when it cannot be resolved from a SIRUP object.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, c(spatial_id_col, clone_col), "receptors")
    receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id_col, sample_col)
    dat <- .sirup_add_coord_spatial_keys(dat)
    clone_ids <- as.character(value %||% variable)
    clone_ids <- unique(clone_ids[!is.na(clone_ids) & nzchar(clone_ids)])
    if (length(clone_ids) == 0L) {
      rlang::abort("`variable` or `value` must name at least one clone ID when `variable_type = 'clone'`.")
    }
    receptor_keys <- receptors[as.character(receptors[[clone_col]]) %in% clone_ids, .sirup_spatial_key_cols(receptors), drop = FALSE]
    receptor_keys <- unique(receptor_keys[!is.na(receptor_keys$.sirup_spatial_id), , drop = FALSE])
    dat_keys <- dat[, .sirup_spatial_key_cols(dat, receptor_keys), drop = FALSE]
    dat$value <- as.integer(.sirup_key_vector(dat_keys) %in% .sirup_key_vector(receptor_keys))
    dat$variable <- paste0("clone_", paste(clone_ids, collapse = "_"))
  }

  dat <- dat[!is.na(dat$value), , drop = FALSE]
  dat
}

.sirup_spatial_fetch_data <- function(object, vars, layer = "counts", assay = NULL) {
  vars <- unique(as.character(vars))
  fetch_one <- function(spatial, sample_id = NULL) {
    dat <- tryCatch(
      .sirup_fetch_data(spatial, vars = vars, layer = layer, assay = assay),
      error = function(e) NULL
    )
    missing_vars <- character()
    if (is.null(dat)) {
      dat <- data.frame(row.names = .sirup_spatial_fetch_cells(spatial))
      for (var in vars) {
        one <- tryCatch(
          .sirup_fetch_data(spatial, vars = var, layer = layer, assay = assay),
          error = function(e) NULL
        )
        if (is.null(one) || !var %in% colnames(one)) {
          dat[[var]] <- rep(NA_real_, nrow(dat))
          missing_vars <- c(missing_vars, var)
        } else {
          dat[[var]] <- one[rownames(dat), var]
        }
      }
    } else {
      absent <- setdiff(vars, colnames(dat))
      for (var in absent) {
        dat[[var]] <- rep(NA_real_, nrow(dat))
      }
      missing_vars <- absent
    }
    if (length(missing_vars) > 0L) {
      where <- if (is.null(sample_id)) "spatial object" else paste0("spatial sample '", sample_id, "'")
      rlang::warn(paste0(
        "Feature(s) not found in ", where, " for assay = ",
        assay %||% .sirup_fetch_data_assay(spatial), ", layer = ", layer,
        ": ", paste(unique(missing_vars), collapse = ", "),
        ". Filling with NA for this section."
      ))
    }
    dat[, vars, drop = FALSE]
  }
  if (.sirup_is_spatial_list(object)) {
    out <- lapply(names(object), function(sample) {
      dat <- fetch_one(object[[sample]], sample_id = sample)
      dat$cell <- rownames(dat)
      rownames(dat) <- NULL
      dat$sample <- sample
      dat
    })
    out <- dplyr::bind_rows(out)
  } else {
    out <- fetch_one(object)
    out$cell <- rownames(out)
    rownames(out) <- NULL
  }
  observed <- vars[vapply(vars, function(var) any(!is.na(out[[var]])), logical(1))]
  if (length(observed) == 0L) {
    rlang::abort(paste0(
      "None of the requested feature(s) could be fetched from assay = ",
      assay %||% "default", ", layer = ", layer, ": ",
      paste(vars, collapse = ", ")
    ))
  }
  out
}

.sirup_spatial_fetch_cells <- function(spatial) {
  cells <- tryCatch(colnames(spatial), error = function(e) character())
  if (length(cells) > 0L) {
    return(as.character(cells))
  }
  meta <- tryCatch(spatial[[]], error = function(e) NULL)
  if (!is.null(meta) && nrow(meta) > 0L) {
    return(rownames(meta))
  }
  character()
}

.sirup_spatial_metadata_values <- function(object, variable) {
  if (.sirup_is_spatial_list(object)) {
    out <- lapply(names(object), function(sample) {
      meta <- object[[sample]]@meta.data
      .check_required_cols(meta, variable, paste0("object[['", sample, "']]@meta.data"))
      dat <- meta[, variable, drop = FALSE]
      dat$cell <- rownames(dat)
      rownames(dat) <- NULL
      dat$sample <- sample
      dat
    })
    return(dplyr::bind_rows(out))
  }
  .check_required_cols(object@meta.data, variable, "object@meta.data")
  meta <- object@meta.data[, variable, drop = FALSE]
  meta$cell <- rownames(meta)
  rownames(meta) <- NULL
  meta
}

.sirup_merge_spatial_values <- function(coords, values) {
  if ("sample" %in% colnames(coords) && "sample" %in% colnames(values)) {
    return(merge(coords, values, by = c("sample", "cell"), all.x = TRUE, sort = FALSE))
  }
  merge(coords, values, by = "cell", all.x = TRUE, sort = FALSE)
}

.sirup_key_vector <- function(dat) {
  if (".sirup_sample" %in% colnames(dat)) {
    return(paste(as.character(dat$.sirup_sample), as.character(dat$.sirup_spatial_id), sep = "\r"))
  }
  as.character(dat$.sirup_spatial_id)
}

#' Get spatial plotting coordinates
#'
#' `spatialCoordinateTable` returns the coordinate table that SIRUP uses for
#' spatial statistics and plots.
#'
#' @details
#' The function takes the coordinates from the Seurat images. If there are none,
#' it takes them from the coordinate columns in the Seurat metadata. The
#' `.sirup_y` column holds the y-coordinate that SIRUP uses for plots.
#'
#' @param object Seurat object.
#' @param images Names of the images to include.
#' @param image_scale Image scale for the Seurat image coordinates.
#' @param swap_coords Swap the x and y coordinates?
#'
#' @return A data frame with at least these columns: `cell`, `x`, `y`,
#'   `.sirup_y`, and `image`.
#' @seealso [spatialCoordinateExtents()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' head(spatialCoordinateTable(getSpatial(x), image_scale = "hires", swap_coords = FALSE))
#' @export
spatialCoordinateTable <- function(object,
                                   images = NULL,
                                   image_scale = NULL,
                                   swap_coords = NULL) {
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  .sirup_spatial_coordinate_table(
    object = context$object,
    images = images,
    image_scale = context$image_scale,
    swap_coords = context$swap_coords
  )
}

#' Summarize spatial coordinate ranges
#'
#' `spatialCoordinateExtents` gives the x and y range of each sample or image.
#'
#' @param object SIRUP object, or Seurat object with spatial images.
#' @param by Group the ranges by `sample` or by `image`.
#' @param images Names of the images or samples to include.
#' @param image_scale Image scale for the Seurat image coordinates.
#' @param swap_coords Swap the x and y coordinates?
#'
#' @return A tibble with the ranges (`x_min`, `x_max`, `y_min`, `y_max`), the
#'   `width`, `height`, and `diagonal`, and the number of spots
#'   (`n_spatial_ids`).
#' @seealso [spatialCoordinateTable()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialCoordinateExtents(x)
#' @export
spatialCoordinateExtents <- function(object,
                                     by = c("sample", "image"),
                                     images = NULL,
                                     image_scale = NULL,
                                     swap_coords = NULL) {
  by <- match.arg(by)
  coords <- spatialCoordinateTable(
    object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  if (!"sample" %in% colnames(coords)) {
    coords$sample <- sirupSampleFromImage(coords$image)
    coords$sample[is.na(coords$sample)] <- "spatial"
  }
  group <- if (identical(by, "sample")) coords$sample else coords$image
  rows <- lapply(split(coords, group, drop = TRUE), function(dat) {
    x_range <- range(dat$x, na.rm = TRUE)
    y_range <- range(dat$y, na.rm = TRUE)
    data.frame(
      sample = unique(dat$sample)[[1L]],
      image = if (identical(by, "image")) unique(dat$image)[[1L]] else NA_character_,
      x_min = x_range[[1L]],
      x_max = x_range[[2L]],
      y_min = y_range[[1L]],
      y_max = y_range[[2L]],
      width = diff(x_range),
      height = diff(y_range),
      diagonal = sqrt(diff(x_range)^2 + diff(y_range)^2),
      n_spatial_ids = length(unique(paste(dat$sample, dat$cell, sep = "\r"))),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(tibble::tibble(sample = character(), image = character(), x_min = numeric(), x_max = numeric(), y_min = numeric(), y_max = numeric(), width = numeric(), height = numeric(), diagonal = numeric(), n_spatial_ids = integer()))
  }
  tibble::as_tibble(dplyr::bind_rows(rows))
}

#' Summarize the spatial spread of clones
#'
#' `spatialCloneSpread` gives the distances between the spots of each clone in
#' each sample.
#'
#' @details
#' For each sample and clone, the function calculates the distance between every
#' pair of points. It reports the mean and the maximum distance.
#'
#' @param object SIRUP object, or Seurat object with spatial images. For a Seurat
#'   object, give `receptors`.
#' @param receptors AIRR-style receptor table. Omit it if `object` is a SIRUP
#'   object.
#' @param clone Name of the column that holds the clone IDs. If `NULL`, the
#'   function uses the stored clone field of a SIRUP object, if there is one.
#' @param spatial_id Name of the column that holds the spatial IDs. If `NULL`,
#'   the function uses the stored spatial ID field of a SIRUP object.
#' @param clone_val Clone IDs to keep.
#' @param images Names of the images or samples to include.
#' @param image_scale Image scale for the Seurat image coordinates.
#' @param swap_coords Swap the x and y coordinates before the function
#'   calculates distances?
#' @param min_spatial_ids Smallest number of observations that a clone needs in
#'   a sample (default 2). If `collapse_within_spot` is `TRUE`, the function
#'   counts distinct spots. Otherwise, it counts sequences.
#' @param collapse_within_spot Merge the sequences that share a spatial ID into
#'   one point before the function calculates distances (default `FALSE`)? If
#'   `FALSE`, the function counts every sequence, so sequences in the same spot
#'   give pairs with distance zero. Use `TRUE` to measure the spread over
#'   distinct spots only.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots and the function reads receptors. The function ignores `slot` if you
#'   give a plain receptor table.
#'
#' @return A tibble with one row for each sample and clone. The columns are
#'   `sample`, `clone_id`, `n_spatial_ids`, `mean_pairwise_distance`, and
#'   `max_pairwise_distance`.
#' @seealso [spatialCloneJoinCount()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialCloneSpread(x)
#' @export
spatialCloneSpread <- function(object,
                               receptors = NULL,
                               clone = NULL,
                               spatial_id = NULL,
                               clone_val = NULL,
                               images = NULL,
                               image_scale = NULL,
                               swap_coords = NULL,
                               min_spatial_ids = 2L,
                               collapse_within_spot = FALSE,
                               slot = NULL) {
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id, slot)
  spatial_source <- if (isSIRUP(object)) object else inputs$object
  spatial <- inputs$object
  receptors <- .as_data_frame(inputs$receptors)
  spatial_id <- inputs$spatial_id_col
  if (is.null(clone) && !is.null(inputs$receptor_fields$clone)) {
    clone <- inputs$receptor_fields$clone
  }
  if (is.null(clone)) {
    rlang::abort("`clone` is required when it cannot be resolved from a SIRUP object.")
  }
  sample_col <- .sirup_receptor_sample_col(receptors, inputs$receptor_fields, spatial)
  .check_required_cols(receptors, c(spatial_id, clone), "receptors")
  coords <- spatialCoordinateTable(spatial_source, images = images, image_scale = image_scale, swap_coords = swap_coords)
  if (!"sample" %in% colnames(coords)) {
    coords$sample <- sirupSampleFromImage(coords$image)
    coords$sample[is.na(coords$sample)] <- "spatial"
  }
  coords <- .sirup_add_coord_spatial_keys(coords)
  receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id, sample_col)
  receptor_coords <- merge(
    receptors,
    coords[, unique(c(.sirup_spatial_key_cols(coords), "cell", "sample", "x", "y", ".sirup_y", "image")), drop = FALSE],
    by = .sirup_spatial_key_cols(receptors, coords),
    all = FALSE,
    sort = FALSE
  )
  if (".sirup_sample" %in% colnames(receptor_coords)) {
    receptor_coords$sample <- as.character(receptor_coords$.sirup_sample)
  } else if ("sample.y" %in% colnames(receptor_coords)) {
    receptor_coords$sample <- as.character(receptor_coords$sample.y)
  } else if (!"sample" %in% colnames(receptor_coords)) {
    receptor_coords$sample <- "spatial"
  }
  if (!is.null(clone_val)) {
    receptor_coords <- receptor_coords[as.character(receptor_coords[[clone]]) %in% as.character(clone_val), , drop = FALSE]
  }
  if (nrow(receptor_coords) == 0L) {
    return(tibble::tibble(sample = character(), clone_id = character(), n_spatial_ids = integer(), mean_pairwise_distance = numeric(), max_pairwise_distance = numeric()))
  }
  receptor_coords <- receptor_coords[!is.na(receptor_coords[[clone]]) & nzchar(as.character(receptor_coords[[clone]])), , drop = FALSE]
  if (isTRUE(collapse_within_spot)) {
    receptor_coords <- receptor_coords[!duplicated(receptor_coords[, c("sample", clone, "cell"), drop = FALSE]), , drop = FALSE]
  }
  groups <- split(receptor_coords, paste(receptor_coords$sample, receptor_coords[[clone]], sep = "\r"), drop = TRUE)
  rows <- lapply(groups, function(dat) {
    n <- nrow(dat)
    if (n < min_spatial_ids) {
      return(NULL)
    }
    distances <- stats::dist(as.matrix(dat[, c("x", "y"), drop = FALSE]))
    data.frame(
      sample = unique(dat$sample)[[1L]],
      clone_id = as.character(unique(dat[[clone]])[[1L]]),
      n_spatial_ids = n,
      mean_pairwise_distance = mean(as.numeric(distances)),
      max_pairwise_distance = max(as.numeric(distances)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(tibble::tibble(sample = character(), clone_id = character(), n_spatial_ids = integer(), mean_pairwise_distance = numeric(), max_pairwise_distance = numeric()))
  }
  tibble::as_tibble(dplyr::bind_rows(rows))
}

# Collect spatial coordinates from Seurat images when present; otherwise use
# spatial coordinate columns stored in object metadata.
.sirup_spatial_coordinate_table <- function(object, images = NULL, image_scale = NULL, swap_coords = NULL) {
  if (.sirup_is_spatial_list(object)) {
    image_scale <- .sirup_resolve_image_scale(image_scale)
    swap_coords <- .sirup_resolve_swap_coords(swap_coords)
    if (is.null(image_scale)) {
      rlang::abort("`image_scale` must be supplied explicitly or stored in a SIRUP object.")
    }
    if (is.null(swap_coords)) {
      rlang::abort("`swap_coords` must be supplied explicitly or stored in a SIRUP object.")
    }
    coords <- lapply(.sirup_spatial_images(object, images = images), function(image) {
      spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
      out <- spatial$coords
      out$image <- image
      out
    })
    return(dplyr::bind_rows(coords))
  }
  options <- .sirup_require_spatial_options(image_scale, swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  image_names <- tryCatch(
    SeuratObject::Images(object),
    error = function(e) character()
  )
  if (is.null(images) && length(image_names) == 0L) {
    return(.sirup_metadata_coordinate_table(object))
  }

  coords <- lapply(.sirup_spatial_images(object, images = images), function(image) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    out <- spatial$coords
    out$image <- image
    out
  })
  dplyr::bind_rows(coords)
}

# Fallback for image-free objects used in tests or non-Visium spatial datasets.
.sirup_metadata_coordinate_table <- function(object) {
  meta <- object@meta.data
  x_col <- intersect(c("x", "imagecol", "pxl_col_in_fullres", "col", "array_col", "spatial_x"), colnames(meta))
  y_col <- intersect(c("y", "imagerow", "pxl_row_in_fullres", "row", "array_row", "spatial_y"), colnames(meta))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    rlang::abort(
      paste0(
        "Could not find spatial images or coordinate columns in `object@meta.data`. ",
        "For image-free spatial statistics, add coordinate columns such as `x` and `y`."
      )
    )
  }
  coords <- data.frame(
    cell = rownames(meta),
    x = meta[[x_col[[1]]]],
    y = meta[[y_col[[1]]]],
    image = "metadata",
    stringsAsFactors = FALSE
  )
  coords <- coords[coords$cell %in% colnames(object), , drop = FALSE]
  coords$.sirup_y <- -coords$y
  coords
}

# Build an spdep neighbor graph either by distance radius or k-nearest
# neighbors over SIRUP's plotting coordinate convention.
.sirup_spatial_neighbors <- function(dat, radius = NULL, k = 6) {
  if (nrow(dat) < 2L) {
    rlang::abort("At least two spatial points are required to build neighborhoods.")
  }
  coords <- as.matrix(dat[, c("x", ".sirup_y")])
  if (!is.null(radius)) {
    spdep::dnearneigh(coords, d1 = 0, d2 = radius, longlat = FALSE)
  } else {
    k <- min(k, nrow(dat) - 1L)
    knn <- spdep::knearneigh(coords, k = k, longlat = FALSE)
    spdep::knn2nb(knn)
  }
}

# Guard against degenerate values before asking spdep for a statistic.
.check_spatial_values <- function(x) {
  if (length(unique(x)) < 2L) {
    rlang::abort("Spatial statistic requires at least two distinct values.")
  }
  invisible(TRUE)
}

# Join receptor clone calls to spatial coordinates and, when supplied, restrict
# them to polygon regions before join-count testing.
.sirup_joincount_receptor_data <- function(object,
                                           receptors,
                                           clone_col,
                                           spatial_id_col,
                                           sample_col,
                                           polygon,
                                           polygon_id,
                                           images,
                                           image_scale,
                                           swap_coords,
                                           x_col,
                                           y_col,
                                           polygon_col,
                                           include_boundary) {
  receptors <- .as_data_frame(receptors)
  .check_required_cols(receptors, c(spatial_id_col, clone_col), "receptors")
  coords <- .sirup_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }

  receptor_cols <- setdiff(
    colnames(receptors),
    setdiff(c("cell", "x", "y", ".sirup_y", "image", "polygon", "sample"), c(spatial_id_col, clone_col, sample_col))
  )
  receptors <- receptors[, receptor_cols, drop = FALSE]
  receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id_col, sample_col)
  coords <- .sirup_add_coord_spatial_keys(coords)
  join_by <- .sirup_spatial_key_cols(receptors, coords)

  if (is.null(polygon)) {
    region_coords <- coords
    region_coords$polygon <- "all"
    out <- merge(
      receptors,
      .sirup_add_coord_spatial_keys(region_coords)[, c(join_by, "cell", "x", "y", ".sirup_y", "image", "polygon"), drop = FALSE],
      by = join_by,
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
    return(out)
  }

  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x_col, y_col), "polygon")
  polygon_specs <- .sirup_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  regions <- lapply(polygon_specs, function(spec) {
    vertices <- .sirup_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col
    )
    inside <- .sirup_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = include_boundary
    )
    region_coords <- coords[inside, , drop = FALSE]
    if (nrow(region_coords) == 0L) {
      return(data.frame())
    }
    region_coords$polygon <- as.character(if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1")
    merge(
      receptors,
      .sirup_add_coord_spatial_keys(region_coords)[, c(join_by, "cell", "x", "y", ".sirup_y", "image", "polygon"), drop = FALSE],
      by = join_by,
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
  })
  dplyr::bind_rows(regions)
}

# Join receptor clone calls to spatial coordinates by using polygon membership
# that was already stored by addPolygons() on a SIRUP object.
.sirup_joincount_receptor_data_stored_polygons <- function(object,
                                                           receptors,
                                                           clone_col,
                                                           spatial_id_col,
                                                           polygon_id,
                                                           images,
                                                           image_scale,
                                                           swap_coords,
                                                           sample_col = NULL) {
  validateSIRUP(object)
  receptors <- .as_data_frame(receptors)
  .check_required_cols(receptors, c(spatial_id_col, clone_col), "receptors")
  coords <- .sirup_spatial_coordinate_table(
    getSpatial(object),
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }
  polygon_ids <- as.character(getPolygons(object)$polygon_id)
  if (!is.null(polygon_id)) {
    polygon_ids <- polygon_ids[polygon_ids %in% as.character(polygon_id)]
  }
  if (length(polygon_ids) == 0L) {
    rlang::abort("No stored polygon IDs matched `polygon_val`.")
  }
  receptor_cols <- setdiff(
    colnames(receptors),
    setdiff(c("cell", "x", "y", ".sirup_y", "image", "polygon", "sample"), c(spatial_id_col, clone_col, sample_col))
  )
  receptors <- receptors[, receptor_cols, drop = FALSE]
  receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id_col, sample_col)
  coords <- .sirup_add_coord_spatial_keys(coords)
  join_by <- .sirup_spatial_key_cols(receptors, coords)

  regions <- lapply(polygon_ids, function(polygon_value) {
    region_coords <- .sirup_stored_polygon_coords(
      coords = coords,
      spatial = getSpatial(object),
      polygon_value = polygon_value
    )
    if (nrow(region_coords) == 0L) {
      return(data.frame())
    }
    region_receptors <- .sirup_receptors_in_polygon(
      receptors = receptors,
      polygon_value = polygon_value,
      polygon_id = "polygon_id"
    )
    if (nrow(region_receptors) == 0L) {
      return(data.frame())
    }
    region_coords$polygon <- polygon_value
    merge(
      region_receptors,
      .sirup_add_coord_spatial_keys(region_coords)[, c(join_by, "cell", "x", "y", ".sirup_y", "image", "polygon"), drop = FALSE],
      by = join_by,
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
  })
  dplyr::bind_rows(regions)
}

# Convert spdep's join-count htest object into one stable summary row.
.sirup_joincount_summary_row <- function(test,
                                         clone_id,
                                         image,
                                         polygon,
                                         n_spatial_ids,
                                         n_clone_spatial_ids,
                                         radius,
                                         k,
                                         alternative) {
  htest <- .sirup_joincount_clone_htest(test, clone_id)
  estimate <- htest$estimate %||% numeric()
  estimate_names <- names(estimate)
  get_estimate <- function(pattern, fallback_index) {
    idx <- grep(pattern, estimate_names, ignore.case = TRUE)
    if (length(idx) > 0L) {
      return(unname(estimate[[idx[[1L]]]]))
    }
    if (length(estimate) >= fallback_index) {
      return(unname(estimate[[fallback_index]]))
    }
    NA_real_
  }
  data.frame(
    clone_id = as.character(clone_id),
    image = as.character(image),
    polygon = as.character(polygon),
    n_spatial_ids = as.integer(n_spatial_ids),
    n_clone_spatial_ids = as.integer(n_clone_spatial_ids),
    observed_join_count = get_estimate("join|observed", 1L),
    expected_join_count = get_estimate("expect", 2L),
    variance = get_estimate("var", 3L),
    statistic = if (length(htest$statistic) > 0L) unname(htest$statistic[[1L]]) else NA_real_,
    p_value = htest$p.value %||% NA_real_,
    method = htest$method %||% NA_character_,
    radius = radius %||% NA_real_,
    k = k,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
}

# spdep may return either an htest or a list of htests depending on factor
# levels; this extracts the htest for the requested clone.
.sirup_joincount_clone_htest <- function(test, clone_id) {
  if (inherits(test, "htest")) {
    return(test)
  }
  if (is.list(test)) {
    named_match <- which(names(test) == as.character(clone_id))
    if (length(named_match) > 0L && inherits(test[[named_match[[1L]]]], "htest")) {
      return(test[[named_match[[1L]]]])
    }
    htest_idx <- which(vapply(test, inherits, logical(1), what = "htest"))
    if (length(htest_idx) > 0L) {
      method_text <- vapply(test[htest_idx], function(x) x$method %||% "", character(1))
      method_match <- grep(as.character(clone_id), method_text, fixed = TRUE)
      if (length(method_match) > 0L) {
        return(test[[htest_idx[[method_match[[1L]]]]]])
      }
      return(test[[utils::tail(htest_idx, 1L)]])
    }
  }
  rlang::abort("Could not extract an htest result from `spdep::joincount.test()`.")
}

# Convert a spatstat Kcross object to a consistent plotting data frame.
.sirup_kcross_df <- function(x, image) {
  dat <- as.data.frame(x)
  theoretical_col <- intersect(c("theo", "theoretical"), colnames(dat))
  simulated_col <- intersect(c("mmean", "mean", "simulated", "simmean", theoretical_col), colnames(dat))
  observed_col <- if ("border" %in% colnames(dat)) {
    "border"
  } else if ("obs" %in% colnames(dat)) {
    "obs"
  } else {
    setdiff(colnames(dat), c("r", theoretical_col, simulated_col, "lo", "hi"))[[1]]
  }
  theoretical <- if (length(theoretical_col) > 0L) dat[[theoretical_col[[1]]]] else pi * dat$r^2
  simulated <- if (length(simulated_col) > 0L) dat[[simulated_col[[1]]]] else theoretical
  out <- data.frame(
    image = image,
    r = dat$r,
    observed = dat[[observed_col]],
    theoretical = theoretical,
    simulated = simulated,
    stringsAsFactors = FALSE
  )
  if (all(c("lo", "hi") %in% colnames(dat))) {
    out$lo <- dat$lo
    out$hi <- dat$hi
  }
  out
}
