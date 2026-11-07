test_that("polygon repair keeps polygonal pieces from geometry collections", {
  skip_on_cran()

  intersection_fun <- function(x, y) {
    sf::st_sfc(sf::st_geometrycollection(list(
      sf::st_polygon(list(matrix(
        c(0.5, 0.5, 1.5, 0.5, 1.5, 1.5, 0.5, 1.5, 0.5, 0.5),
        ncol = 2,
        byrow = TRUE
      ))),
      sf::st_linestring(matrix(c(0, 0, 2, 2), ncol = 2, byrow = TRUE))
    )))
  }

  repaired <- .sirup_repair_polygon_sfc(intersection_fun(NULL, NULL))

  expect_false(is.null(repaired))
  expect_gt(length(repaired), 0L)
  expect_equal(unique(as.character(sf::st_geometry_type(repaired))), "POLYGON")
})
