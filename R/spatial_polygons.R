#' Calculate BCR clone dominance inside spatial polygons
#'
#' Finds Seurat spatial barcodes whose coordinates fall inside closed polygons,
#' subsets BCR rows to those barcodes, and returns the fraction of regional BCRs
#' represented by each clone. Points on the polygon boundary are included by
#' default.
#'
#' @param bcrs BCR/AIRR-style data frame containing clone assignments and a
#'   spatial barcode column. This may instead be a SIRAP object when `object`
#'   is omitted.
#' @param object Seurat object with spatial coordinates, metadata coordinate
#'   columns supported by `spatialVariableData()`, or a SIRAP object. When this
#'   is a SIRAP object, `bcrs` may be omitted.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   summarized independently unless `polygon_id` selects one or more polygons.
#' @param barcode_col Column in `bcrs` matching Seurat spot/cell barcodes. If
#'   `NULL`, a barcode-like column is guessed.
#' @param clone_col Column in `bcrs` containing clone IDs.
#' @param locus Optional locus value(s) used to filter `bcrs`, for example
#'   `"IGH"`. Use `NULL` to keep all loci.
#' @param locus_col Column in `bcrs` containing locus calls.
#' @param count_col Optional numeric count column. When `NULL`, each BCR row
#'   counts once.
#' @param v_call_col Column in `bcrs` containing V-gene calls. When present,
#'   distinct allele-level V genes, V genes without alleles, and V-gene
#'   superfamilies are returned for each clone.
#' @param exclude_singletons Exclude clones represented by only one BCR row
#'   after locus filtering.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists. When `NULL`, all polygons in `polygon_col` are analyzed.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @param require_closed Require the final polygon vertex to equal the first.
#'   Set to `FALSE` to close open polygons automatically.
#' @return A spatial clone dominance object containing the dominance table,
#'   regional barcodes, regional coordinate rows, regional BCR rows, and polygon
#'   vertices prepared for plotting.
#' @export
spatialCloneDominance <- function(bcrs = NULL,
                                  object = NULL,
                                  polygon,
                                  barcode_col = NULL,
                                  clone_col = "clone_id",
                                  locus = NULL,
                                  locus_col = "locus",
                                  count_col = NULL,
                                  v_call_col = "v_call",
                                  exclude_singletons = FALSE,
                                  images = NULL,
                                  image.scale = "hires",
                                  swap_coords = FALSE,
                                  polygon_id = NULL,
                                  x_col = "x",
                                  y_col = "y",
                                  polygon_col = "polygon",
                                  include_boundary = TRUE,
                                  require_closed = TRUE) {
  clone_missing <- missing(clone_col)
  locus_missing <- missing(locus_col)
  count_missing <- missing(count_col)
  if (isSIRAP(bcrs) && is.null(object)) {
    object <- bcrs
    bcrs <- NULL
  }
  is_sirap <- isSIRAP(object)
  inputs <- .sirap_resolve_inputs(object, bcrs, barcode_col)
  object <- inputs$object
  bcrs <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(bcrs)) {
    rlang::abort("`bcrs` is required when `object` is not a SIRAP object.")
  }
  if (is_sirap && clone_missing && !is.null(inputs$schema$clone)) {
    clone_col <- inputs$schema$clone
  }
  if (is_sirap && locus_missing && !is.null(inputs$schema$locus)) {
    locus_col <- inputs$schema$locus
  }
  if (is_sirap && count_missing && !is.null(inputs$schema$count)) {
    count_col <- inputs$schema$count
  }
  bcrs <- .as_data_frame(bcrs)
  polygon <- .as_data_frame(polygon)
  barcode_col <- .guess_barcode_col(bcrs, barcode_col)
  .check_required_cols(bcrs, c(barcode_col, clone_col), "bcrs")
  .check_required_cols(polygon, c(x_col, y_col), "polygon")
  if (!is.null(locus)) {
    .check_required_cols(bcrs, locus_col, "bcrs")
    missing_locus <- setdiff(as.character(locus), unique(as.character(bcrs[[locus_col]])))
    if (length(missing_locus) == length(locus)) {
      rlang::warn(paste0(
        "No rows in `bcrs` match `locus = ",
        paste(as.character(locus), collapse = ","),
        "` using `locus_col = \"",
        locus_col,
        "\"`."
      ))
    }
  }
  if (!is.null(count_col)) {
    .check_required_cols(bcrs, count_col, "bcrs")
  }
  bcrs <- bcrs[
    !is.na(bcrs[[clone_col]]) &
      as.character(bcrs[[clone_col]]) != "",
    ,
    drop = FALSE
  ]
  if (!is.null(locus)) {
    bcrs <- bcrs[bcrs[[locus_col]] %in% locus, , drop = FALSE]
  }
  if (isTRUE(exclude_singletons)) {
    clone_ids <- as.character(bcrs[[clone_col]])
    clone_counts <- table(clone_ids, useNA = "no")
    bcrs <- bcrs[clone_ids %in% names(clone_counts)[clone_counts > 1L], , drop = FALSE]
  }

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  polygon_plot_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_id,
    polygon_col = polygon_col,
    x_col = x_col,
    y_col = y_col,
    require_closed = require_closed
  )

  results <- lapply(polygon_specs, function(spec) {
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    inside <- .sirap_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = include_boundary
    )
    region_coords <- coords[inside, , drop = FALSE]
    region_barcodes <- unique(as.character(region_coords$cell))

    regional_bcrs <- bcrs[as.character(bcrs[[barcode_col]]) %in% region_barcodes, , drop = FALSE]

    if (nrow(regional_bcrs) == 0L) {
      out <- data.frame(
        clone_id = character(),
        v_gene = character(),
        v_gene_family = character(),
        v_gene_superfamily = character(),
        n_bcrs = numeric(),
        proportion = numeric(),
        region_bcr_total = numeric(),
        n_region_barcodes = integer(),
        stringsAsFactors = FALSE
      )
    } else {
      clone_ids <- as.character(regional_bcrs[[clone_col]])
      weights <- if (is.null(count_col)) {
        rep(1, nrow(regional_bcrs))
      } else {
        as.numeric(regional_bcrs[[count_col]])
      }
      weights[is.na(weights)] <- 0
      out <- stats::aggregate(
        n_bcrs ~ clone_id,
        data = data.frame(clone_id = clone_ids, n_bcrs = weights, stringsAsFactors = FALSE),
        FUN = sum
      )
      out <- out[out$n_bcrs > 0, , drop = FALSE]

      clone_ids <- unique(as.character(regional_bcrs[[clone_col]]))
      if (v_call_col %in% colnames(regional_bcrs)) {
        v_calls <- as.character(regional_bcrs[[v_call_col]])
        names(v_calls) <- as.character(regional_bcrs[[clone_col]])
        v_summary <- lapply(clone_ids, function(clone_id) {
          calls <- v_calls[names(v_calls) == clone_id]
          genes <- .sirap_v_genes(calls)
          families <- unique(stats::na.omit(.sirap_v_gene_family(genes)))
          families <- sort(families)
          superfamilies <- unique(stats::na.omit(.sirap_v_gene_superfamily(.sirap_v_gene_family(genes))))
          superfamilies <- sort(superfamilies)
          data.frame(
            clone_id = clone_id,
            v_gene = if (length(genes) == 0L) NA_character_ else paste(genes, collapse = ","),
            v_gene_family = if (length(families) == 0L) NA_character_ else paste(families, collapse = ","),
            v_gene_superfamily = if (length(superfamilies) == 0L) NA_character_ else paste(superfamilies, collapse = ","),
            stringsAsFactors = FALSE
          )
        })
        v_summary <- do.call(rbind, v_summary)
      } else {
        v_summary <- data.frame(
          clone_id = clone_ids,
          v_gene = rep(NA_character_, length(clone_ids)),
          v_gene_family = rep(NA_character_, length(clone_ids)),
          v_gene_superfamily = rep(NA_character_, length(clone_ids)),
          stringsAsFactors = FALSE
        )
      }

      out <- merge(out, v_summary, by = "clone_id", all.x = TRUE, sort = FALSE)
      out$v_gene[is.na(out$v_gene)] <- NA_character_
      out$v_gene_family[is.na(out$v_gene_family)] <- NA_character_
      out$v_gene_superfamily[is.na(out$v_gene_superfamily)] <- NA_character_
      total <- sum(out$n_bcrs)
      out$proportion <- out$n_bcrs / total
      out$region_bcr_total <- total
      out$n_region_barcodes <- length(region_barcodes)
      out <- out[order(out$proportion, out$n_bcrs, decreasing = TRUE), , drop = FALSE]
      out <- out[, c("clone_id", "v_gene", "v_gene_family", "v_gene_superfamily", "n_bcrs", "proportion", "region_bcr_total", "n_region_barcodes"), drop = FALSE]
    }
    if (!is.null(locus)) {
      out$locus_filter <- rep(paste(as.character(locus), collapse = ","), nrow(out))
    }
    polygon_value <- if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1"
    out$polygon <- rep(as.character(polygon_value), nrow(out))

    region_barcode_table <- data.frame(
      barcode = region_barcodes,
      polygon = rep(as.character(polygon_value), length(region_barcodes)),
      stringsAsFactors = FALSE
    )
    region_coords$polygon <- rep(as.character(polygon_value), nrow(region_coords))
    if (nrow(regional_bcrs) == 0L || nrow(region_coords) == 0L) {
      region_bcr_table <- regional_bcrs[0, , drop = FALSE]
      region_bcr_table$.sirap_barcode <- character()
      region_bcr_table$x <- numeric()
      region_bcr_table$y <- numeric()
      region_bcr_table$.sirap_y <- numeric()
      region_bcr_table$image <- character()
      region_bcr_table$polygon <- character()
      region_bcr_table <- .sirap_add_v_gene_columns(region_bcr_table, v_call_col = v_call_col)
    } else {
      region_coord_table <- region_coords[, c("cell", "x", "y", "image"), drop = FALSE]
      region_coord_table$.sirap_y <- -region_coord_table$y
      region_bcr_table <- .sirap_drop_plot_coordinate_cols(regional_bcrs)
      region_bcr_table$.sirap_barcode <- as.character(region_bcr_table[[barcode_col]])
      region_bcr_table <- merge(
        region_bcr_table,
        region_coord_table,
        by.x = ".sirap_barcode",
        by.y = "cell",
        all.x = FALSE,
        all.y = FALSE,
        sort = FALSE
      )
      region_bcr_table$polygon <- as.character(polygon_value)
      region_bcr_table <- .sirap_add_v_gene_columns(region_bcr_table, v_call_col = v_call_col)
    }
    list(
      dominance = out,
      region_barcodes = region_barcode_table,
      region_coordinates = region_coords,
      region_bcrs = region_bcr_table
    )
  })

  out <- .sirap_rbind(results, "dominance")

  out <- list(
    dominance = out,
    region_barcodes = .sirap_rbind(results, "region_barcodes"),
    region_coordinates = .sirap_rbind(results, "region_coordinates"),
    region_bcrs = .sirap_rbind(results, "region_bcrs"),
    polygons = polygon_plot_data
  )
  class(out) <- c("sirap_spatial_clone_dominance", class(out))
  out
}

#' Select spatial polygon vertices interactively
#'
#' Opens a Shiny app for clicking polygon vertices on a Seurat spatial
#' coordinate plot. Each click snaps to the nearest spatial barcode. Use
#' "Close polygon" or "New polygon" before finishing to repeat the first vertex
#' at the end of each polygon.
#'
#' @param object Seurat object with spatial coordinates.
#' @param feature Optional feature shown as a color scale under the clicked
#'   vertices. The app also lets you change this interactively.
#' @param assay Optional assay to set before fetching feature values.
#' @param layer Expression layer/slot used for `feature`.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param save_path Optional CSV path for saved vertices.
#' @param plot_image Show the Seurat tissue image under the spatial points.
#'   Defaults to `FALSE`.
#' @param downsample_n Optional number of spatial points to show in the Shiny
#'   plot. Use `NULL` to show all points. Clicks still snap to the nearest
#'   barcode in the full coordinate table.
#' @param downsample_seed Random seed used when `downsample_n` is not `NULL`.
#' @param point_size Size of spatial points in the plot.
#' @param point_alpha Initial opacity for spatial points.
#' @param launch.browser Passed to `shiny::runApp()`.
#' @param host,port Optional host and port passed to `shiny::runApp()`. These
#'   are useful when running from an interactive HPC session through an SSH
#'   tunnel.
#' @param quiet Suppress warnings while the app is running.
#' @return Data frame with `barcode`, `x`, `y`, `polygon`, and `image` columns.
#' @export
selectSpatialPolygons <- function(object,
                                  feature = NULL,
                                  assay = NULL,
                                  layer = "counts",
                                  images = NULL,
                                  image.scale = "hires",
                                  swap_coords = FALSE,
                                  save_path = NULL,
                                  plot_image = FALSE,
                                  downsample_n = NULL,
                                  downsample_seed = 1L,
                                  point_size = 0.8,
                                  point_alpha = 0.9,
                                  launch.browser = TRUE,
                                  host = getOption("shiny.host", "127.0.0.1"),
  port = getOption("shiny.port", NULL),
  quiet = TRUE) {
  .require_namespace("SeuratObject")
  .require_namespace("shiny")
  if (!"package:shiny" %in% search()) {
    suppressPackageStartupMessages(base::attachNamespace("shiny"))
  }
  old_shiny_useragg <- getOption("shiny.useragg")
  options(shiny.useragg = FALSE)
  on.exit(options(shiny.useragg = old_shiny_useragg), add = TRUE)
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale, swap_coords = swap_coords)
  if (nrow(coords) == 0L) {
    rlang::abort("No spatial coordinates were found in `object`.")
  }
  coords$barcode <- coords$cell
  coords$.sirap_plot_y <- -coords$y
  y_axis_offset <- min(coords$.sirap_plot_y, na.rm = TRUE)
  plot_coords <- .sirap_downsample_polygon_plot_coords(
    coords = coords,
    downsample_n = downsample_n,
    seed = downsample_seed
  )
  image_names <- unique(stats::na.omit(as.character(coords$image)))
  image_names <- setdiff(image_names, "metadata")
  plot_image_info <- NULL
  if (length(image_names) > 0L) {
    if (length(image_names) > 1L) {
      rlang::warn("Multiple spatial images are present; showing the first image in `selectSpatialPolygons()`.")
    }
    raster <- .sirap_image_raster(object[[image_names[[1]]]], alpha = 1)
    image_dim <- dim(raster)
    plot_image_info <- list(
      raster = raster,
      width = image_dim[[2]],
      height = image_dim[[1]]
    )
  }
  x_range <- diff(range(plot_coords$x, na.rm = TRUE))
  y_range <- diff(range(plot_coords$.sirap_plot_y, na.rm = TRUE))
  plot_height <- if (!is.finite(x_range) || !is.finite(y_range) || x_range <= 0 || y_range <= 0) {
    500
  } else {
    max(500, min(1200, round(900 * y_range / x_range)))
  }
  default_feature <- feature
  if (is.null(default_feature)) {
    features <- tryCatch(rownames(object), error = function(e) character())
    default_feature <- if (length(features) == 0L) NULL else features[[1]]
  }
  save_default <- save_path
  if (is.null(save_default)) {
    save_default <- ""
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel("Select Spatial Polygon Vertices"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::textInput("feature", "Feature", value = default_feature %||% ""),
        shiny::numericInput("point_size", "Point size", value = point_size, min = 0.01, max = 10),
        shiny::sliderInput("point_alpha", "Point opacity", min = 0, max = 1, value = point_alpha, step = 0.05, ticks = FALSE),
        shiny::checkboxInput("plot_image", "Show image", value = isTRUE(plot_image)),
        shiny::textInput("save_path", "CSV path", value = save_default),
        shiny::actionButton("close_polygon", "Close polygon"),
        shiny::actionButton("new_polygon", "New polygon"),
        shiny::actionButton("remove_last_point", "Remove last point"),
        shiny::actionButton("reset_points", "Reset"),
        shiny::actionButton("done", "Done"),
        shiny::tableOutput("selected_points"),
        width = 4
      ),
      shiny::mainPanel(
        shiny::plotOutput(
          "spatial_plot",
          width = "100%",
          height = paste0(plot_height, "px"),
          click = shiny::clickOpts("spatial_plot_click")
        ),
        width = 8
      )
    )
  )

  server <- function(input, output, session) {
    old_warn <- getOption("warn")
    if (isTRUE(quiet)) {
      options(warn = -1)
      on.exit(options(warn = old_warn), add = TRUE)
    }

    selected_points <- shiny::reactiveVal(.sirap_empty_polygon_vertices())
    current_polygon <- shiny::reactiveVal(1L)

    plot_data <- shiny::reactive({
      dat <- plot_coords
      feature_name <- input$feature
      if (!is.null(feature_name) && nzchar(feature_name)) {
        expr <- tryCatch(
          .sirap_fetch_data(object, vars = feature_name, layer = layer),
          error = function(e) NULL
        )
        if (!is.null(expr) && feature_name %in% colnames(expr)) {
          dat <- merge(dat, expr[, feature_name, drop = FALSE], by.x = "cell", by.y = 0, all.x = TRUE)
          dat <- dat[order(dat[[feature_name]], na.last = TRUE), , drop = FALSE]
        }
      }
      dat
    })

    output$spatial_plot <- shiny::renderPlot({
      dat <- plot_data()
      feature_name <- input$feature
      has_feature <- !is.null(feature_name) && nzchar(feature_name) && feature_name %in% colnames(dat)
      aes_args <- if (has_feature) {
        ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, color = .data[[feature_name]], key = .data$cell)
      } else {
        ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, key = .data$cell)
      }

      p <- ggplot2::ggplot(dat, aes_args)
      if (isTRUE(input$plot_image) && !is.null(plot_image_info)) {
        p <- p +
          ggplot2::annotation_raster(
            plot_image_info$raster,
            xmin = 0,
            xmax = plot_image_info$width,
            ymin = -plot_image_info$height,
            ymax = 0,
            interpolate = FALSE
          )
      }
      p <- p +
        ggplot2::geom_point(size = input$point_size, alpha = input$point_alpha) +
        ggplot2::coord_fixed() +
        ggplot2::scale_y_continuous(name = "y", labels = function(x) round(x - y_axis_offset)) +
        ggplot2::labs(x = "x") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
      if (has_feature) {
        p <- p + ggplot2::scale_color_gradientn(colors = grDevices::hcl.colors(64, "viridis"))
      }

      selected <- selected_points()
      if (nrow(selected) > 0L) {
        selected_path <- selected[
          stats::ave(selected$polygon, selected$polygon, FUN = length) > 1L,
          ,
          drop = FALSE
        ]
        if (nrow(selected_path) > 0L) {
          p <- p +
            ggplot2::geom_path(
              data = selected_path,
              ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, group = .data$polygon),
              inherit.aes = FALSE,
              linewidth = 0.5,
              color = "black"
            )
        }
        p <- p +
          ggplot2::geom_point(
            data = selected,
            ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, group = .data$polygon),
            inherit.aes = FALSE,
            size = input$point_size * 2,
            shape = 21,
            fill = "white",
            color = "black",
            stroke = 0.8
          )
      }
      p
    })

    shiny::observeEvent(input$spatial_plot_click, {
      click_data <- input$spatial_plot_click
      if (is.null(click_data)) {
        return()
      }
      if (requireNamespace("RANN", quietly = TRUE)) {
        closest_idx <- RANN::nn2(
          data = as.matrix(coords[, c("x", ".sirap_plot_y"), drop = FALSE]),
          query = matrix(c(click_data$x, click_data$y), nrow = 1),
          k = 1
        )$nn.idx[[1]]
      } else {
        click_dist <- (coords$x - click_data$x)^2 + (coords$.sirap_plot_y - click_data$y)^2
        closest_idx <- which.min(click_dist)
      }
      closest <- coords[closest_idx, , drop = FALSE]
      closest <- data.frame(
        barcode = closest$cell,
        x = closest$x,
        y = closest$y,
        .sirap_plot_y = closest$.sirap_plot_y,
        polygon = current_polygon(),
        image = closest$image,
        stringsAsFactors = FALSE
      )
      selected_points(rbind(selected_points(), closest))
    })

    output$selected_points <- shiny::renderTable({
      current <- selected_points()
      current <- current[current$polygon == current_polygon(), , drop = FALSE]
      current <- .sirap_public_polygon_vertices(current)
      current[, setdiff(colnames(current), "image"), drop = FALSE]
    })

    shiny::observeEvent(input$close_polygon, {
      selected_points(.sirap_close_current_polygon(selected_points(), current_polygon()))
    })

    shiny::observeEvent(input$new_polygon, {
      selected_points(.sirap_close_current_polygon(selected_points(), current_polygon()))
      current_polygon(current_polygon() + 1L)
    })

    shiny::observeEvent(input$remove_last_point, {
      current <- selected_points()
      if (nrow(current) > 0L) {
        selected_points(current[-nrow(current), , drop = FALSE])
      }
    })

    shiny::observeEvent(input$reset_points, {
      selected_points(.sirap_empty_polygon_vertices())
      current_polygon(1L)
    })

    shiny::observeEvent(input$done, {
      out <- selected_points()
      for (polygon_id in unique(out$polygon)) {
        out <- .sirap_close_current_polygon(out, polygon_id)
      }
      out <- .sirap_public_polygon_vertices(out)
      path <- input$save_path
      if (!is.null(path) && nzchar(path)) {
        utils::write.csv(out, path, row.names = FALSE)
        shiny::showNotification(paste("Saved polygons to", path), type = "message")
      }
      shiny::stopApp(out)
    })
  }

  shiny::runApp(
    shiny::shinyApp(ui, server),
    launch.browser = launch.browser,
    host = host,
    port = port
  )
}

#' Format spatial polygon point sets as polygon vertices
#'
#' Converts a data frame containing all points inside one or more polygons into
#' the same vertex-table structure returned by `selectSpatialPolygons()`. By
#' default, interior point sets are converted to polygon outlines by buffering
#' and unioning spots with `sf`, then mapping the exterior ring back to the
#' nearest selected barcode coordinates. This preserves a valid boundary order
#' for concave polygons while returning Visium barcode coordinates as vertices.
#' If the input already contains boundary points, use `point_type = "boundary"`.
#'
#' @param object Seurat object with spatial coordinates.
#' @param polygon_points Data frame containing all points inside each polygon,
#'   or boundary points when `point_type = "boundary"`.
#' @param barcode_col Column in `polygon_points` containing Seurat barcodes.
#'   Required when `x_col` and `y_col` are not present.
#' @param polygon_col Column in `polygon_points` identifying each polygon.
#' @param x_col,y_col Coordinate columns in `polygon_points`, if already
#'   present.
#' @param image_col Image column in `polygon_points`, if already present.
#' @param point_type Are `polygon_points` all interior points or already
#'   boundary points?
#' @param buffer_radius Optional radius used to buffer points when
#'   converting interior points to a polygon outline. When `NULL`, the radius is
#'   estimated from the nearest-neighbor spacing in the full spatial coordinate
#'   grid.
#' @param order_col Optional column used to order points within each polygon
#'   before formatting. This is mainly useful with `point_type = "boundary"`.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate axes.
#' @param close_polygons Close each polygon by repeating the first vertex when
#'   needed.
#' @return Data frame with `barcode`, `x`, `y`, `polygon`, and `image` columns.
#'   Polygons with fewer than three points are excluded with a warning.
#' @export
formatSpatialPolygons <- function(object,
                                  polygon_points,
                                  barcode_col = "barcode",
                                  polygon_col = "polygon",
                                  x_col = "x",
                                  y_col = "y",
                                  image_col = "image",
                                  point_type = c("interior", "boundary"),
                                  buffer_radius = NULL,
                                  order_col = NULL,
                                  images = NULL,
                                  image.scale = "hires",
                                  swap_coords = FALSE,
                                  close_polygons = TRUE) {
  .require_namespace("SeuratObject")
  point_type <- match.arg(point_type)
  polygon_points <- .as_data_frame(polygon_points)
  .check_required_cols(polygon_points, polygon_col, "polygon_points")
  if (nrow(polygon_points) == 0L) {
    rlang::abort("`polygon_points` must contain at least one row.")
  }
  if (any(is.na(polygon_points[[polygon_col]]))) {
    rlang::abort("`polygon_points` contains missing polygon IDs.")
  }
  if (!is.null(order_col)) {
    .check_required_cols(polygon_points, order_col, "polygon_points")
    polygon_points$.sirap_input_order <- seq_len(nrow(polygon_points))
    ordered_rows <- unlist(lapply(unique(polygon_points[[polygon_col]]), function(polygon_id) {
      idx <- which(polygon_points[[polygon_col]] == polygon_id)
      idx[order(polygon_points[[order_col]][idx], polygon_points$.sirap_input_order[idx])]
    }), use.names = FALSE)
    polygon_points <- polygon_points[ordered_rows, , drop = FALSE]
    polygon_points$.sirap_input_order <- NULL
  }

  has_barcode <- barcode_col %in% colnames(polygon_points)
  has_xy <- all(c(x_col, y_col) %in% colnames(polygon_points))
  if (!has_barcode && !has_xy) {
    rlang::abort("`polygon_points` must contain barcode values or x/y coordinate columns.")
  }

  coords <- .sirap_spatial_coordinate_table(
    object,
    images = images,
    image.scale = image.scale,
    swap_coords = swap_coords
  )
  coords$cell <- as.character(coords$cell)
  if (has_barcode) {
    polygon_points[[barcode_col]] <- as.character(polygon_points[[barcode_col]])
  }

  if (has_barcode && (!has_xy || !image_col %in% colnames(polygon_points))) {
    has_image <- image_col %in% colnames(polygon_points)
    if (anyDuplicated(coords$cell) > 0L && !has_image) {
      rlang::abort(
        "`object` has duplicated spatial barcodes across images; provide `image_col` in `polygon_points`."
      )
    }
    if (has_image) {
      coord_key <- paste(coords$cell, coords$image, sep = "\r")
      point_key <- paste(polygon_points[[barcode_col]], polygon_points[[image_col]], sep = "\r")
      coord_lookup <- coords[match(point_key, coord_key), , drop = FALSE]
    } else {
      coord_lookup <- coords[match(polygon_points[[barcode_col]], coords$cell), , drop = FALSE]
    }
    missing_barcodes <- unique(polygon_points[[barcode_col]][is.na(coord_lookup$cell)])
    if (length(missing_barcodes) > 0L) {
      rlang::abort(paste0(
        "Some polygon point barcodes were not found in `object`: ",
        paste(utils::head(missing_barcodes, 10), collapse = ", ")
      ))
    }
    if (!has_xy) {
      polygon_points[[x_col]] <- coord_lookup$x
      polygon_points[[y_col]] <- coord_lookup$y
    }
    if (!image_col %in% colnames(polygon_points)) {
      polygon_points[[image_col]] <- coord_lookup$image
    }
  }

  if (!image_col %in% colnames(polygon_points)) {
    image_values <- unique(as.character(coords$image))
    if (length(image_values) != 1L) {
      rlang::abort("`polygon_points` must contain `image_col` when `object` has multiple spatial images.")
    }
    polygon_points[[image_col]] <- image_values[[1L]]
  }

  out <- data.frame(
    barcode = if (has_barcode) polygon_points[[barcode_col]] else NA_character_,
    x = as.numeric(polygon_points[[x_col]]),
    y = as.numeric(polygon_points[[y_col]]),
    polygon = polygon_points[[polygon_col]],
    image = as.character(polygon_points[[image_col]]),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$x) | !is.finite(out$y))) {
    rlang::abort("`polygon_points` contains non-finite polygon coordinates.")
  }

  multi_image_polygons <- names(which(tapply(out$image, out$polygon, function(x) length(unique(x)) > 1L)))
  if (length(multi_image_polygons) > 0L) {
    rlang::abort(paste0(
      "Each polygon ID must belong to one image. These polygons span multiple images: ",
      paste(multi_image_polygons, collapse = ", ")
    ))
  }

  if (point_type == "interior") {
    .require_namespace("sf")
    if (is.null(buffer_radius)) {
      buffer_radius <- .sirap_default_polygon_buffer_radius(coords)
    }
    out <- .sirap_polygon_interior_to_sf_boundary(
      points = out,
      buffer_radius = buffer_radius
    )
  }

  if (isTRUE(close_polygons)) {
    out <- .sirap_close_formatted_polygons(out)
  }
  rownames(out) <- NULL
  out
}

#' Plot spatial polygons on a tissue image
#'
#' Plots one or more polygon vertex tables, such as the output of
#' `selectSpatialPolygons()`, over the Seurat tissue image.
#'
#' @param object Seurat object with spatial images.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   drawn independently.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param polygon_alpha Fill opacity for polygons.
#' @param fill Fill color for polygons.
#' @param color Outline color for polygons.
#' @param linewidth Outline width for polygons.
#' @param image.alpha Tissue image opacity.
#' @param label Plot labels at polygon centers.
#' @param label_col Column in the polygon data to use for labels.
#' @param label_size Text size for polygon labels.
#' @param label_color Text color for polygon labels.
#' @param zoom Zoom to the polygon extents instead of the full tissue.
#' @param require_closed Require the final polygon vertex to equal the first.
#'   Set to `FALSE` to close open polygons automatically.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object or a named list of ggplots when
#'   `combine = FALSE`.
#' @export
plotSpatialPolygons <- function(object,
                                polygon,
                                images = NULL,
                                image.scale = "hires",
                                swap_coords = FALSE,
                                polygon_id = NULL,
                                x_col = "x",
                                y_col = "y",
                                polygon_col = "polygon",
                                polygon_alpha = 0.2,
                                fill = "#0072B2",
                                color = "#0072B2",
                                linewidth = 0.6,
                                image.alpha = 1,
                                label = TRUE,
                                label_col = "polygon",
                                label_size = 3,
                                label_color = "black",
                                zoom = FALSE,
                                require_closed = TRUE,
                                combine = TRUE) {
  .require_namespace("SeuratObject")
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x_col, y_col), "polygon")
  if (!is.numeric(polygon_alpha) || length(polygon_alpha) != 1L || is.na(polygon_alpha) ||
      polygon_alpha < 0 || polygon_alpha > 1) {
    rlang::abort("`polygon_alpha` must be a number between 0 and 1.")
  }

  polygon_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_id,
    polygon_col = polygon_col,
    x_col = x_col,
    y_col = y_col,
    require_closed = require_closed,
    keep_cols = label_col
  )
  if (!isTRUE(label) && !identical(label, FALSE)) {
    rlang::abort("`label` must be TRUE or FALSE.")
  }
  if (isTRUE(label)) {
    .check_required_cols(polygon_data, label_col, "polygon")
    label_idx <- !duplicated(polygon_data$polygon)
    label_lookup <- data.frame(
      polygon = polygon_data$polygon[label_idx],
      .sirap_label = as.character(polygon_data[[label_col]][label_idx]),
      stringsAsFactors = FALSE
    )
    label_data <- .sirap_polygon_centers(polygon_data, unique(polygon_data$polygon))
    label_data <- merge(label_data, label_lookup, by = "polygon", all.x = TRUE, sort = FALSE)
  }

  plots <- list()
  full_image_data <- data.frame(x = numeric(), .sirap_y = numeric())
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    plots[[image]] <- .sirap_spatial_base(
      plot_data = polygon_data,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else full_image_data
    ) +
      ggplot2::geom_polygon(
        data = polygon_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
        inherit.aes = FALSE,
        fill = fill,
        color = color,
        alpha = polygon_alpha,
        linewidth = linewidth
      )
    if (isTRUE(label)) {
      plots[[image]] <- plots[[image]] +
        ggplot2::geom_text(
          data = label_data,
          ggplot2::aes(x = .data$x, y = .data$.sirap_y, label = .data$.sirap_label),
          inherit.aes = FALSE,
          size = label_size,
          color = label_color,
          fontface = "bold"
        )
    }
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot receptors inside spatial polygons
#'
#' Plots every receptor row whose spatial barcode falls inside one or more
#' selected polygons over the tissue image, colored by a receptor-table column.
#'
#' @param receptors Receptor/AIRR-style data frame containing a spatial barcode column,
#'   the object returned by `spatialCloneDominance()`, or a SIRAP object when
#'   `object` is omitted.
#' @param object Seurat object with spatial images, or a SIRAP object. When this
#'   is a SIRAP object, `receptors` may be omitted.
#' @param polygon Data frame of polygon vertices in order. Optional when
#'   `receptors` is a spatial clone dominance object.
#' @param color_by Column in `receptors` used to color receptor points.
#' @param barcode_col Column in `receptors` matching Seurat spot/cell barcodes. If
#'   `NULL`, a barcode-like column is guessed.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param point_size Size of receptor points.
#' @param point_alpha Opacity for receptor points.
#' @param polygon_alpha Fill opacity for polygons.
#' @param polygon_fill Fill color for polygons.
#' @param polygon_color Outline color for polygons.
#' @param image.alpha Tissue image opacity.
#' @param top_n Optional number of most frequent `color_by` levels to keep
#'   before grouping remaining levels as `other_label`.
#' @param other_label Label for grouped lower-frequency levels.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot receptors grouped into `other_label`.
#' @param zoom Zoom to the selected polygon extents. Defaults to `FALSE`, which
#'   shows the full tissue image.
#' @param require_closed Require the final polygon vertex to equal the first.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object or a named list of ggplots when
#'   `combine = FALSE`.
#' @export
plotSpatialRegionReceptors <- function(receptors = NULL,
                                       object = NULL,
                                       polygon = NULL,
                                       color_by = "clone_id",
                                       barcode_col = NULL,
                                       images = NULL,
                                       image.scale = "hires",
                                       swap_coords = FALSE,
                                       polygon_id = NULL,
                                       x_col = "x",
                                       y_col = "y",
                                       polygon_col = "polygon",
                                       point_size = 1,
                                       point_alpha = 0.5,
                                       polygon_alpha = 0.12,
                                       polygon_fill = "#FFFFFF",
                                       polygon_color = "black",
                                       image.alpha = 1,
                                       top_n = NULL,
                                       other_label = "Other",
                                       other_color = "#D3D3D3",
                                       plot_other = TRUE,
                                       zoom = FALSE,
                                       require_closed = TRUE,
                                       combine = TRUE) {
  color_missing <- missing(color_by)
  if (isSIRAP(receptors) && is.null(object)) {
    object <- receptors
    receptors <- NULL
  }
  is_sirap <- isSIRAP(object)
  inputs <- .sirap_resolve_inputs(object, receptors, barcode_col)
  object <- inputs$object
  receptors <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRAP object.")
  }
  if (is_sirap && color_missing && !is.null(inputs$schema$clone)) {
    color_by <- inputs$schema$clone
  }
  .require_namespace("SeuratObject")
  region <- .sirap_region_from_plot_input(
    bcrs = receptors,
    object = object,
    polygon = polygon,
    barcode_col = barcode_col,
    images = images,
    image.scale = image.scale,
    swap_coords = swap_coords,
    polygon_id = polygon_id,
    x_col = x_col,
    y_col = y_col,
    polygon_col = polygon_col,
    require_closed = require_closed
  )
  if (nrow(region$bcrs) == 0L) {
    rlang::abort("No receptor rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$bcrs, color_by, "receptors")
  region$bcrs$.sirap_color_by <- .sirap_lump_values(region$bcrs[[color_by]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region$bcrs <- region$bcrs[as.character(region$bcrs$.sirap_color_by) != other_label, , drop = FALSE]
  }
  if (nrow(region$bcrs) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }
  region$bcrs$.sirap_is_other <- as.character(region$bcrs$.sirap_color_by) == other_label
  region$bcrs <- region$bcrs[order(!region$bcrs$.sirap_is_other), , drop = FALSE]
  color_values <- .sirap_discrete_palette_values(
    levels(region$bcrs$.sirap_color_by),
    other_label = other_label,
    other_color = other_color
  )

  plots <- list()
  full_image_data <- data.frame(x = numeric(), .sirap_y = numeric())
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    plot_receptors <- region$bcrs[region$bcrs$image == image | region$bcrs$image == "metadata", , drop = FALSE]
    plots[[image]] <- .sirap_spatial_base(
      plot_data = if (isTRUE(zoom)) region$polygons else spatial$coords,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else full_image_data
    ) +
      ggplot2::geom_polygon(
        data = region$polygons,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
        inherit.aes = FALSE,
        fill = polygon_fill,
        color = polygon_color,
        alpha = polygon_alpha,
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data = plot_receptors,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, color = .data$.sirap_color_by),
        inherit.aes = FALSE,
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
      ggplot2::labs(color = color_by)
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot 3D clone-density surfaces over a spatial image
#'
#' Uses the regional receptor rows returned by `spatialCloneDominance()` to draw
#' smoothed clone-abundance density surfaces above the original tissue image in
#' an interactive `rgl` scene.
#'
#' @param dominance Object returned by `spatialCloneDominance()`.
#' @param object Seurat object with spatial images.
#' @param clone_col Column in `dominance$region_bcrs` used to define clone
#'   groups.
#' @param weight_col Optional numeric column used as receptor weights. When
#'   `NULL`, each receptor row contributes one count.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select from the dominance
#'   object.
#' @param top_n Optional number of most abundant clone groups to show
#'   individually. Remaining clone groups are combined as `other_label`.
#' @param other_label Label for grouped lower-abundance clone groups.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot grouped lower-abundance clone groups as
#'   `other_label`.
#' @param bandwidth Gaussian kernel bandwidth in spatial coordinate units. Use a
#'   single value for isotropic smoothing or two values for x/y bandwidths. When
#'   `NULL`, a bandwidth is estimated from the regional coordinates.
#' @param grid_n Number of grid points used for density surfaces. Use one value
#'   for both axes or two values for x/y grid sizes.
#' @param plot_type Plot smoothed density surfaces or binned 3D histogram bars.
#' @param hist_bins Number of x/y bins used when `plot_type = "histogram"`.
#' @param zoom Zoom the tissue image plane and 3D view to the selected region.
#'   Defaults to `FALSE`, which shows the full spatial image.
#' @param height_scale Maximum density height as a fraction of the plotted
#'   spatial span.
#' @param z_offset Vertical offset of the density surfaces above the tissue
#'   image, as a fraction of the plotted spatial span.
#' @param density_alpha Opacity for density surfaces.
#' @param density_floor Relative density threshold below which surface cells are
#'   hidden. This removes low-density Gaussian tails that otherwise look like a
#'   colored plane.
#' @param image_alpha Opacity for the tissue image plane.
#' @param show_points Show receptor locations above the tissue plane.
#' @param point_size Size of receptor-location points.
#' @param point_alpha Opacity for receptor-location points.
#' @param show_polygon Show polygon boundaries above the tissue plane.
#' @param polygon_color Color for polygon boundary lines.
#' @param polygon_lwd Width for polygon boundary lines.
#' @param label_density_peaks Label the highest point of the first few density
#'   surfaces.
#' @param label_top_n Number of density surfaces to label when
#'   `label_density_peaks = TRUE`.
#' @param label_cex Text scale for density peak labels.
#' @param show_axes Show the 3D bounding box and axis labels.
#' @param texture_dir Optional directory where the tissue texture PNG should be
#'   written. Use a persistent directory when saving an HTML widget.
#' @param as_widget Return an `rgl` htmlwidget. When `FALSE`, the scene is drawn
#'   and the current rgl device ID is returned.
#' @param return_data Return the computed density data instead of drawing an rgl
#'   scene. This does not require `rgl`.
#' @return An `rgl` widget/device ID, or a list of computed plot data when
#'   `return_data = TRUE`.
#' @export
plotSpatialCloneDensity3D <- function(dominance,
                                      object,
                                      clone_col = "clone_id",
                                      weight_col = NULL,
                                      images = NULL,
                                      image.scale = "hires",
                                      swap_coords = FALSE,
                                      polygon_id = NULL,
                                      top_n = 10,
                                      other_label = "Other",
                                      other_color = "#D3D3D3",
                                      plot_other = TRUE,
                                      bandwidth = NULL,
                                      grid_n = 80,
                                      plot_type = c("density", "histogram"),
                                      hist_bins = 30,
                                      zoom = FALSE,
                                      height_scale = 0.22,
                                      z_offset = 0.04,
                                      density_alpha = 0.55,
                                      density_floor = 0.03,
                                      image_alpha = 1,
                                      show_points = FALSE,
                                      point_size = 5,
                                      point_alpha = 0.75,
                                      show_polygon = TRUE,
                                      polygon_color = "white",
                                      polygon_lwd = 3,
                                      label_density_peaks = FALSE,
                                      label_top_n = 5,
                                      label_cex = 0.8,
                                      show_axes = TRUE,
                                      texture_dir = NULL,
                                      as_widget = TRUE,
                                      return_data = FALSE) {
  if (!inherits(dominance, "sirap_spatial_clone_dominance")) {
    rlang::abort("`dominance` must be an object returned by `spatialCloneDominance()`.")
  }
  if (is.null(dominance$region_bcrs) || is.null(dominance$polygons)) {
    rlang::abort("`dominance` must contain `region_bcrs` and `polygons` from `spatialCloneDominance()`.")
  }
  plot_type <- match.arg(plot_type)
  .require_namespace("SeuratObject")

  region_bcrs <- .as_data_frame(dominance$region_bcrs)
  polygons <- .as_data_frame(dominance$polygons)
  if (!is.null(polygon_id)) {
    region_bcrs <- region_bcrs[as.character(region_bcrs$polygon) %in% as.character(polygon_id), , drop = FALSE]
    polygons <- polygons[as.character(polygons$polygon) %in% as.character(polygon_id), , drop = FALSE]
  }
  if (nrow(region_bcrs) == 0L) {
    rlang::abort("No regional receptor rows are available for the selected polygon(s).")
  }

  plot_data <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    image_bcrs <- region_bcrs[region_bcrs$image == image | region_bcrs$image == "metadata", , drop = FALSE]
    if (nrow(image_bcrs) == 0L) {
      next
    }
    plot_data[[image]] <- .sirap_clone_density_3d_data(
      region_bcrs = image_bcrs,
      polygons = polygons,
      image = spatial$image,
      spatial_coords = spatial$coords,
      clone_col = clone_col,
      weight_col = weight_col,
      top_n = top_n,
      other_label = other_label,
      other_color = other_color,
      plot_other = plot_other,
      bandwidth = bandwidth,
      grid_n = grid_n,
      plot_type = plot_type,
      hist_bins = hist_bins,
      zoom = zoom,
      height_scale = height_scale,
      z_offset = z_offset,
      density_floor = density_floor
    )
  }
  plot_data <- Filter(Negate(is.null), plot_data)
  if (length(plot_data) == 0L) {
    rlang::abort("No regional receptor rows match the selected image(s).")
  }
  if (isTRUE(return_data)) {
    return(if (length(plot_data) == 1L) plot_data[[1L]] else plot_data)
  }

  if (isTRUE(as_widget)) {
    old_rgl_use_null <- getOption("rgl.useNULL")
    options(rgl.useNULL = TRUE)
    on.exit(options(rgl.useNULL = old_rgl_use_null), add = TRUE)
  }
  .require_namespace("rgl")
  plots <- lapply(plot_data, function(dat) {
    .sirap_draw_clone_density_3d(
      dat,
      density_alpha = density_alpha,
      image_alpha = image_alpha,
      show_points = show_points,
      point_size = point_size,
      point_alpha = point_alpha,
      show_polygon = show_polygon,
      polygon_color = polygon_color,
      polygon_lwd = polygon_lwd,
      label_density_peaks = label_density_peaks,
      label_top_n = label_top_n,
      label_cex = label_cex,
      show_axes = show_axes,
      texture_dir = texture_dir,
      as_widget = as_widget
    )
  })
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    plots
  }
}

#' Plot spatial polygon pie summaries
#'
#' Summarizes the composition of BCRs inside each polygon for a selected column.
#' The default mode plots pie charts beside a numbered polygon map. The
#' `spatial_clipped` mode clips each pie chart to its polygon boundary and
#' requires the suggested `sf` package.
#'
#' @param bcrs BCR/AIRR-style data frame containing a spatial barcode column,
#'   the object returned by `spatialCloneDominance()`, or a SIRAP object when
#'   `object` is omitted.
#' @param object Seurat object with spatial images, or a SIRAP object. When this
#'   is a SIRAP object, `bcrs` may be omitted.
#' @param polygon Data frame of polygon vertices in order. Optional when `bcrs`
#'   is a spatial clone dominance object.
#' @param fill_by Column in `bcrs` summarized in each pie chart.
#' @param barcode_col Column in `bcrs` matching Seurat spot/cell barcodes. If
#'   `NULL`, a barcode-like column is guessed.
#' @param mode Plot style. `side_by_side` shows ordinary pies beside a polygon
#'   map. `spatial_clipped` clips pies to polygon boundaries on the tissue.
#' @param top_n Number of most frequent `fill_by` levels to identify within each
#'   polygon. The union of these polygon-level top levels is shown individually
#'   everywhere; remaining levels are grouped as `other_label`.
#' @param label_top_n Number of levels to label per polygon. Defaults to
#'   0 because this gets messy fast.
#' @param other_label Label for grouped lower-frequency levels.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot grouped lower-frequency BCRs as `other_label`.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param polygon_alpha Fill opacity for polygons in the numbered map.
#' @param image.alpha Tissue image opacity.
#' @param pie_alpha Pie opacity.
#' @param label_size Text size for polygon and pie labels.
#' @param pies_per_page Maximum number of ordinary pie charts to show per page
#'   when `mode = "side_by_side"`. When more polygons are present, a named list
#'   of page plots is returned.
#' @param spatial_legend Show a shared fill legend for `spatial_clipped` pie
#'   plots. Defaults to `FALSE`.
#' @param zoom Zoom to the plotted polygon extents. Defaults to `FALSE`, which
#'   shows the full tissue image.
#' @param require_closed Require the final polygon vertex to equal the first.
#' @return A ggplot, patchwork object, or named list of ggplots.
#' @export
plotSpatialPolygonPies <- function(bcrs = NULL,
                                   object = NULL,
                                   polygon = NULL,
                                   fill_by = "clone_id",
                                   barcode_col = NULL,
                                   mode = c("side_by_side", "spatial_clipped"),
                                   top_n = 10,
                                   label_top_n = 0,
                                   other_label = "Other",
                                   other_color = "#D3D3D3",
                                   plot_other = TRUE,
                                   images = NULL,
                                   image.scale = "hires",
                                   swap_coords = FALSE,
                                   polygon_id = NULL,
                                   x_col = "x",
                                   y_col = "y",
                                   polygon_col = "polygon",
                                   polygon_alpha = 0.2,
                                   image.alpha = 1,
                                   pie_alpha = 0.9,
                                   label_size = 3,
                                   pies_per_page = 6L,
                                   spatial_legend = FALSE,
                                   zoom = FALSE,
                                   require_closed = TRUE) {
  fill_missing <- missing(fill_by)
  if (isSIRAP(bcrs) && is.null(object)) {
    object <- bcrs
    bcrs <- NULL
  }
  is_sirap <- isSIRAP(object)
  inputs <- .sirap_resolve_inputs(object, bcrs, barcode_col)
  object <- inputs$object
  bcrs <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(bcrs)) {
    rlang::abort("`bcrs` is required when `object` is not a SIRAP object.")
  }
  if (is_sirap && fill_missing && !is.null(inputs$schema$clone)) {
    fill_by <- inputs$schema$clone
  }
  .require_namespace("SeuratObject")
  mode <- match.arg(mode)
  if (is.null(label_top_n)) {
    label_top_n <- Inf
  }
  if (!is.numeric(label_top_n) || length(label_top_n) != 1L || is.na(label_top_n) || label_top_n < 0) {
    rlang::abort("`label_top_n` must be NULL or a non-negative number.")
  }
  if (!is.numeric(pies_per_page) || length(pies_per_page) != 1L || is.na(pies_per_page) || pies_per_page < 1) {
    rlang::abort("`pies_per_page` must be a positive number.")
  }
  pies_per_page <- as.integer(pies_per_page)
  region <- .sirap_region_from_plot_input(
    bcrs = bcrs,
    object = object,
    polygon = polygon,
    barcode_col = barcode_col,
    images = images,
    image.scale = image.scale,
    swap_coords = swap_coords,
    polygon_id = polygon_id,
    x_col = x_col,
    y_col = y_col,
    polygon_col = polygon_col,
    require_closed = require_closed
  )
  if (nrow(region$bcrs) == 0L) {
    rlang::abort("No BCR rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$bcrs, fill_by, "bcrs")
  pie_data <- .sirap_polygon_pie_summary(
    region$bcrs,
    fill_by = fill_by,
    top_n = top_n,
    label_top_n = label_top_n,
    other_label = other_label,
    plot_other = plot_other
  )
  if (nrow(pie_data) == 0L) {
    rlang::abort("No BCR rows remain to plot after applying `top_n` and `plot_other`.")
  }
  fill_values <- .sirap_discrete_palette_values(
    unique(as.character(pie_data$level)),
    other_label = other_label,
    other_color = other_color
  )
  centers <- .sirap_polygon_centers(region$polygons, pie_data$polygon)

  if (mode == "side_by_side") {
    .sirap_side_by_side_polygon_pies(
      object = object,
      images = images,
      image.scale = image.scale,
      swap_coords = swap_coords,
      polygon_plot_data = region$polygons,
      centers = centers,
      pie_data = pie_data,
      fill_values = fill_values,
      polygon_alpha = polygon_alpha,
      image.alpha = image.alpha,
      label_size = label_size,
      pies_per_page = pies_per_page,
      zoom = zoom
    )
  } else {
    .sirap_spatial_clipped_polygon_pies(
      object = object,
      images = images,
      image.scale = image.scale,
      swap_coords = swap_coords,
      polygon_plot_data = region$polygons,
      centers = centers,
      pie_data = pie_data,
      fill_values = fill_values,
      fill_by = fill_by,
      image.alpha = image.alpha,
      pie_alpha = pie_alpha,
      label_size = label_size,
      spatial_legend = spatial_legend,
      zoom = zoom
    )
  }
}

#' Subset a Seurat object to spatial barcodes inside polygons
#'
#' Keeps only Seurat barcodes whose spatial coordinates fall inside one or more
#' polygon regions and records polygon membership in the returned object's
#' metadata.
#'
#' @param object Seurat object with spatial coordinates.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   evaluated independently.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param metadata_col Metadata column added to the returned Seurat object with
#'   polygon membership. Barcodes in overlapping polygons receive comma-separated
#'   polygon IDs.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @param require_closed Require the final polygon vertex to equal the first.
#'   Set to `FALSE` to close open polygons automatically.
#' @return A Seurat object subset to barcodes inside the selected polygon(s),
#'   with `metadata_col` appended to `@meta.data`.
#' @export
seuratSubsetPolygon <- function(object,
                                polygon,
                                images = NULL,
                                image.scale = "hires",
                                swap_coords = FALSE,
                                polygon_id = NULL,
                                x_col = "x",
                                y_col = "y",
                                polygon_col = "polygon",
                                metadata_col = "polygon",
                                include_boundary = TRUE,
                                require_closed = TRUE) {
  .require_namespace("SeuratObject")
  if (!is.character(metadata_col) || length(metadata_col) != 1L || is.na(metadata_col) || !nzchar(metadata_col)) {
    rlang::abort("`metadata_col` must be a non-empty character string.")
  }

  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x_col, y_col), "polygon")
  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }

  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  membership <- lapply(polygon_specs, function(spec) {
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    inside <- .sirap_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = include_boundary
    )
    if (!any(inside)) {
      return(data.frame())
    }
    polygon_value <- if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1"
    data.frame(
      cell = as.character(coords$cell[inside]),
      polygon = rep(as.character(polygon_value), sum(inside)),
      image = as.character(coords$image[inside]),
      x = as.numeric(coords$x[inside]),
      y = as.numeric(coords$y[inside]),
      stringsAsFactors = FALSE
    )
  })
  membership <- .sirap_bind_rows(membership)
  if (nrow(membership) == 0L) {
    rlang::abort("No Seurat barcodes were found inside the selected polygon(s).")
  }

  keep_cells <- unique(as.character(membership$cell))
  keep_cells <- keep_cells[keep_cells %in% colnames(object)]
  if (length(keep_cells) == 0L) {
    rlang::abort("No polygon-matched barcodes are present in `object`.")
  }

  polygon_values <- stats::aggregate(
    polygon ~ cell,
    data = membership,
    FUN = function(x) paste(unique(as.character(x)), collapse = ",")
  )
  polygon_lookup <- stats::setNames(as.character(polygon_values$polygon), as.character(polygon_values$cell))

  out <- object[, keep_cells]
  out@meta.data[[metadata_col]] <- unname(polygon_lookup[rownames(out@meta.data)])
  out
}

#' Define key spatial areas interactively
#'
#' Compatibility wrapper for the older `defineKeyAreas()` workflow. It now
#' returns polygon vertices directly; pass the result to `spatialCloneDominance()`.
#'
#' @param seurat_obj Seurat object.
#' @param assay Assay used for expression plotting.
#' @param dataset Expression layer/slot used for expression plotting.
#' @param autoload Retained for compatibility. Interior point expansion is no
#'   longer performed; polygon vertices are returned.
#' @param quiet Suppress warnings while the app is running.
#' @param ... Passed to `selectSpatialPolygons()`.
#' @return Data frame of polygon vertices.
#' @export
defineKeyAreas <- function(seurat_obj,
                           assay = "Spatial",
                           dataset = "scale.data",
                           autoload = TRUE,
                           quiet = TRUE,
                           ...) {
  if (isTRUE(autoload)) {
    warning(
      "`defineKeyAreas()` now returns polygon vertices directly; ",
      "use `spatialCloneDominance()` or another polygon-aware function to select enclosed spots.",
      call. = FALSE
    )
  }
  selectSpatialPolygons(
    object = seurat_obj,
    assay = assay,
    layer = dataset,
    quiet = quiet,
    ...
  )
}

.sirap_select_polygon_vertices <- function(polygon,
                                           polygon_id,
                                           polygon_col,
                                           x_col,
                                           y_col,
                                           require_closed) {
  if (!is.null(polygon_id) && polygon_col %in% colnames(polygon)) {
    polygon <- polygon[as.character(polygon[[polygon_col]]) == as.character(polygon_id), , drop = FALSE]
  } else if (is.null(polygon_id) && polygon_col %in% colnames(polygon)) {
    polygon_ids <- unique(stats::na.omit(polygon[[polygon_col]]))
    if (length(polygon_ids) > 1L) {
      rlang::abort("`polygon` contains multiple polygons; pass `polygon_id`.")
    }
  }
  if (nrow(polygon) < 4L) {
    rlang::abort("`polygon` must contain at least four rows, including the closing vertex.")
  }
  out <- data.frame(
    x = as.numeric(polygon[[x_col]]),
    y = as.numeric(polygon[[y_col]]),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$x) | !is.finite(out$y))) {
    rlang::abort("`polygon` contains non-finite coordinates.")
  }
  closed <- .sirap_same_xy(out[1, , drop = FALSE], out[nrow(out), , drop = FALSE])
  if (!closed && isTRUE(require_closed)) {
    rlang::abort("`polygon` must be closed: the last vertex should match the first vertex.")
  }
  if (!closed) {
    out <- rbind(out, out[1, , drop = FALSE])
  }
  out
}

.sirap_polygon_specs <- function(polygon, polygon_id, polygon_col) {
  has_polygon_col <- polygon_col %in% colnames(polygon)
  if (!has_polygon_col) {
    if (!is.null(polygon_id)) {
      rlang::abort("`polygon_id` can only be used when `polygon_col` exists in `polygon`.")
    }
    return(list(list(id = NULL, add_polygon_col = FALSE)))
  }

  ids <- if (is.null(polygon_id)) {
    unique(stats::na.omit(polygon[[polygon_col]]))
  } else {
    polygon_id
  }
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) {
    rlang::abort("No polygon IDs were found in `polygon`.")
  }
  missing_ids <- setdiff(as.character(ids), as.character(polygon[[polygon_col]]))
  if (length(missing_ids) > 0L) {
    rlang::abort(paste0(
      "`polygon_id` values were not found in `polygon`: ",
      paste(missing_ids, collapse = ", ")
    ))
  }
  lapply(ids, function(id) list(id = id, add_polygon_col = TRUE))
}

.sirap_rbind <- function(results, name) {
  tables <- lapply(results, `[[`, name)
  if (length(tables) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

.sirap_v_genes <- function(v_calls) {
  if (length(v_calls) == 0L) {
    return(character())
  }
  genes <- unlist(strsplit(as.character(v_calls), "[,;]"), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  sort(unique(genes))
}

.sirap_v_gene_family <- function(v_genes) {
  sub("\\*.*$", "", v_genes)
}

.sirap_v_gene_superfamily <- function(v_gene_families) {
  sub("^([A-Za-z]+V[0-9]+).*$", "\\1", v_gene_families)
}

.sirap_points_in_polygon <- function(px,
                                     py,
                                     vx,
                                     vy,
                                     include_boundary = TRUE,
                                     tol = sqrt(.Machine$double.eps)) {
  px <- as.numeric(px)
  py <- as.numeric(py)
  vx <- as.numeric(vx)
  vy <- as.numeric(vy)
  if (length(vx) != length(vy) || length(vx) < 4L) {
    rlang::abort("Polygon coordinates must contain at least four paired vertices.")
  }
  if (.sirap_same_xy(data.frame(x = vx[[1]], y = vy[[1]]), data.frame(x = vx[[length(vx)]], y = vy[[length(vy)]]))) {
    vx <- vx[-length(vx)]
    vy <- vy[-length(vy)]
  }

  inside <- rep(FALSE, length(px))
  j <- length(vx)
  for (i in seq_along(vx)) {
    crosses <- ((vy[[i]] > py) != (vy[[j]] > py)) &
      (px < (vx[[j]] - vx[[i]]) * (py - vy[[i]]) / (vy[[j]] - vy[[i]]) + vx[[i]])
    crosses[is.na(crosses)] <- FALSE
    inside <- xor(inside, crosses)
    j <- i
  }

  if (isTRUE(include_boundary)) {
    boundary <- rep(FALSE, length(px))
    j <- length(vx)
    for (i in seq_along(vx)) {
      boundary <- boundary | .sirap_points_on_segment(px, py, vx[[j]], vy[[j]], vx[[i]], vy[[i]], tol = tol)
      j <- i
    }
    inside <- inside | boundary
  }
  inside
}

.sirap_points_on_segment <- function(px, py, x1, y1, x2, y2, tol) {
  dx <- x2 - x1
  dy <- y2 - y1
  cross <- (px - x1) * dy - (py - y1) * dx
  scale <- pmax(abs(dx), abs(dy), 1)
  on_line <- abs(cross) <= tol * scale
  within_x <- px >= pmin(x1, x2) - tol & px <= pmax(x1, x2) + tol
  within_y <- py >= pmin(y1, y2) - tol & py <= pmax(y1, y2) + tol
  on_line & within_x & within_y
}

.sirap_same_xy <- function(a, b) {
  isTRUE(all.equal(as.numeric(a$x), as.numeric(b$x), tolerance = sqrt(.Machine$double.eps))) &&
    isTRUE(all.equal(as.numeric(a$y), as.numeric(b$y), tolerance = sqrt(.Machine$double.eps)))
}

.sirap_empty_polygon_vertices <- function() {
  data.frame(
    barcode = character(),
    x = numeric(),
    y = numeric(),
    .sirap_plot_y = numeric(),
    polygon = integer(),
    image = character(),
    stringsAsFactors = FALSE
  )
}

.sirap_public_polygon_vertices <- function(points) {
  points[, setdiff(colnames(points), ".sirap_plot_y"), drop = FALSE]
}

.sirap_polygon_plot_data <- function(polygon,
                                     polygon_id,
                                     polygon_col,
                                     x_col,
                                     y_col,
                                     require_closed,
                                     keep_cols = character()) {
  specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  keep_cols <- unique(keep_cols[keep_cols %in% colnames(polygon)])
  out <- lapply(specs, function(spec) {
    polygon_rows <- if (isTRUE(spec$add_polygon_col)) {
      polygon[as.character(polygon[[polygon_col]]) == as.character(spec$id), , drop = FALSE]
    } else {
      polygon
    }
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    vertices$polygon <- if (isTRUE(spec$add_polygon_col)) {
      as.character(spec$id)
    } else {
      "polygon_1"
    }
    for (keep_col in setdiff(keep_cols, c(x_col, y_col, "polygon"))) {
      values <- unique(as.character(polygon_rows[[keep_col]]))
      values <- values[!is.na(values)]
      vertices[[keep_col]] <- if (length(values) == 0L) NA_character_ else values[[1L]]
    }
    vertices
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out$.sirap_y <- -out$y
  out
}

.sirap_region_from_plot_input <- function(bcrs,
                                          object,
                                          polygon,
                                          barcode_col,
                                          images,
                                          image.scale,
                                          swap_coords,
                                          polygon_id,
                                          x_col,
                                          y_col,
                                          polygon_col,
                                          require_closed) {
  if (inherits(bcrs, "sirap_spatial_clone_dominance")) {
    if (is.null(bcrs$region_bcrs) || is.null(bcrs$polygons)) {
      rlang::abort("Dominance objects passed to plotting functions must come from `spatialCloneDominance()`.")
    }
    region_bcrs <- .as_data_frame(bcrs$region_bcrs)
    polygon_plot_data <- .as_data_frame(bcrs$polygons)
    if (!is.null(polygon_id)) {
      region_bcrs <- region_bcrs[as.character(region_bcrs$polygon) %in% as.character(polygon_id), , drop = FALSE]
      polygon_plot_data <- polygon_plot_data[
        as.character(polygon_plot_data$polygon) %in% as.character(polygon_id),
        ,
        drop = FALSE
      ]
    }
    return(list(
      bcrs = region_bcrs,
      polygons = polygon_plot_data
    ))
  }
  if (is.null(polygon)) {
    rlang::abort("`polygon` is required unless `bcrs` is a spatial clone dominance object returned by `spatialCloneDominance()`.")
  }

  bcrs <- .as_data_frame(bcrs)
  polygon <- .as_data_frame(polygon)
  barcode_col <- .guess_barcode_col(bcrs, barcode_col)
  .check_required_cols(bcrs, barcode_col, "bcrs")
  .check_required_cols(polygon, c(x_col, y_col), "polygon")

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  coords$.sirap_y <- -coords$y
  polygon_plot_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_id,
    polygon_col = polygon_col,
    x_col = x_col,
    y_col = y_col,
    require_closed = require_closed
  )
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)

  regional <- lapply(polygon_specs, function(spec) {
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    inside <- .sirap_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = TRUE
    )
    region_coords <- coords[inside, , drop = FALSE]
    region_bcrs <- bcrs[as.character(bcrs[[barcode_col]]) %in% as.character(region_coords$cell), , drop = FALSE]
    if (nrow(region_bcrs) == 0L) {
      return(data.frame())
    }
    region_bcrs <- .sirap_drop_plot_coordinate_cols(region_bcrs)
    region_bcrs$.sirap_barcode <- as.character(region_bcrs[[barcode_col]])
    region_bcrs <- merge(
      region_bcrs,
      region_coords[, c("cell", "x", "y", ".sirap_y", "image"), drop = FALSE],
      by.x = ".sirap_barcode",
      by.y = "cell",
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
    region_bcrs$polygon <- as.character(if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1")
    region_bcrs
  })
  regional <- .sirap_bind_rows(regional)
  regional <- .sirap_add_v_gene_columns(regional)
  list(
    bcrs = regional,
    polygons = polygon_plot_data
  )
}

.sirap_drop_plot_coordinate_cols <- function(dat) {
  dat[, setdiff(colnames(dat), c("x", "y", ".sirap_y", "image")), drop = FALSE]
}

.sirap_add_v_gene_columns <- function(dat, v_call_col = "v_call") {
  if (nrow(dat) == 0L) {
    if (!"v_gene" %in% colnames(dat)) {
      dat$v_gene <- character()
    }
    if (!"v_gene_family" %in% colnames(dat)) {
      dat$v_gene_family <- character()
    }
    if (!"v_gene_superfamily" %in% colnames(dat)) {
      dat$v_gene_superfamily <- character()
    }
    return(dat)
  }
  if (!"v_gene" %in% colnames(dat)) {
    dat$v_gene <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(x) paste(.sirap_v_genes(x), collapse = ","), character(1))
    } else {
      rep(NA_character_, nrow(dat))
    }
    dat$v_gene[dat$v_gene == ""] <- NA_character_
  }
  if (!"v_gene_family" %in% colnames(dat)) {
    dat$v_gene_family <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        families <- unique(stats::na.omit(.sirap_v_gene_family(.sirap_v_genes(v_calls))))
        families <- sort(families)
        if (length(families) == 0L) NA_character_ else paste(families, collapse = ",")
      }, character(1))
    } else {
      .sirap_v_gene_family(dat$v_gene)
    }
  }
  if (!"v_gene_superfamily" %in% colnames(dat)) {
    dat$v_gene_superfamily <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        superfamilies <- .sirap_v_genes(v_calls)
        superfamilies <- .sirap_v_gene_family(superfamilies)
        superfamilies <- unique(stats::na.omit(.sirap_v_gene_superfamily(superfamilies)))
        superfamilies <- sort(superfamilies)
        if (length(superfamilies) == 0L) NA_character_ else paste(superfamilies, collapse = ",")
      }, character(1))
    } else {
      .sirap_v_gene_superfamily(dat$v_gene_family)
    }
  }
  dat
}

.sirap_bind_rows <- function(dots) {
  dots <- Filter(function(x) !is.null(x) && nrow(x) > 0L, dots)
  if (length(dots) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, dots)
  rownames(out) <- NULL
  out
}

.sirap_lump_values <- function(x, top_n = NULL, other_label = "Other") {
  values <- as.character(x)
  values[is.na(values) | !nzchar(values)] <- NA_character_
  if (is.null(top_n)) {
    return(factor(values, levels = sort(unique(stats::na.omit(values)))))
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    rlang::abort("`top_n` must be NULL or a positive number.")
  }
  counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
  keep <- names(counts)[seq_len(min(as.integer(top_n), length(counts)))]
  values[!is.na(values) & !values %in% keep] <- other_label
  factor(values, levels = c(keep, if (any(values == other_label, na.rm = TRUE)) other_label))
}

.sirap_discrete_palette_values <- function(levels, other_label = "Other", other_color = "#D3D3D3") {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]
  has_other <- other_label %in% levels
  main_levels <- setdiff(levels, other_label)
  .require_namespace("pals")
  values <- if (length(main_levels) == 0L) {
    character()
  } else if (length(main_levels) <= 32L) {
    pals::glasbey(length(main_levels))
  } else {
    grDevices::colorRampPalette(pals::glasbey(32L))(length(main_levels))
  }
  names(values) <- main_levels
  if (has_other) {
    values <- c(values, stats::setNames(other_color, other_label))
  }
  values[levels]
}

.sirap_clone_density_3d_data <- function(region_bcrs,
                                         polygons,
                                         image,
                                         spatial_coords,
                                         clone_col,
                                         weight_col,
                                         top_n,
                                         other_label,
                                         other_color,
                                         plot_other,
                                         bandwidth,
                                         grid_n,
                                         plot_type,
                                         hist_bins,
                                         zoom,
                                         height_scale,
                                         z_offset,
                                         density_floor) {
  .check_required_cols(region_bcrs, c(clone_col, "x", "y", ".sirap_y", "polygon"), "dominance$region_bcrs")
  .check_required_cols(polygons, c("x", "y", ".sirap_y", "polygon"), "dominance$polygons")
  if (!is.null(weight_col)) {
    .check_required_cols(region_bcrs, weight_col, "dominance$region_bcrs")
  }
  if (!is.numeric(height_scale) || length(height_scale) != 1L || !is.finite(height_scale) || height_scale <= 0) {
    rlang::abort("`height_scale` must be a positive number.")
  }
  if (!is.numeric(z_offset) || length(z_offset) != 1L || !is.finite(z_offset) || z_offset < 0) {
    rlang::abort("`z_offset` must be a non-negative number.")
  }
  if (!is.numeric(density_floor) || length(density_floor) != 1L || !is.finite(density_floor) || density_floor < 0 || density_floor >= 1) {
    rlang::abort("`density_floor` must be a number greater than or equal to 0 and less than 1.")
  }

  region_bcrs <- region_bcrs[
    is.finite(region_bcrs$x) & is.finite(region_bcrs$y) & is.finite(region_bcrs$.sirap_y),
    ,
    drop = FALSE
  ]
  if (nrow(region_bcrs) == 0L) {
    rlang::abort("No finite regional receptor coordinates are available.")
  }
  region_bcrs$.sirap_density_group <- .sirap_lump_values(region_bcrs[[clone_col]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region_bcrs <- region_bcrs[as.character(region_bcrs$.sirap_density_group) != other_label, , drop = FALSE]
  }
  if (nrow(region_bcrs) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }

  weights <- if (is.null(weight_col)) {
    rep(1, nrow(region_bcrs))
  } else {
    as.numeric(region_bcrs[[weight_col]])
  }
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) {
    rlang::abort("All receptor weights are zero or non-finite.")
  }
  region_bcrs$.sirap_density_weight <- weights

  image_dim <- .sirap_image_dim(image)
  region_xlim <- .sirap_density_plot_range(c(region_bcrs$x, polygons$x), c(0, image_dim$width))
  region_ylim <- .sirap_density_plot_range(c(region_bcrs$.sirap_y, polygons$.sirap_y), c(-image_dim$height, 0))
  object_xlim <- .sirap_density_plot_range(spatial_coords$x, c(0, image_dim$width))
  object_ylim <- .sirap_density_plot_range(spatial_coords$.sirap_y, c(-image_dim$height, 0))
  image_xlim <- if (isTRUE(zoom)) region_xlim else object_xlim
  image_ylim <- if (isTRUE(zoom)) region_ylim else object_ylim
  image_extent <- .sirap_image_texture_extent(image_dim, image_xlim, image_ylim)
  span <- max(diff(region_xlim), diff(region_ylim))
  z_base <- z_offset * span
  max_height <- height_scale * span
  levels <- levels(region_bcrs$.sirap_density_group)
  levels <- levels[levels %in% as.character(region_bcrs$.sirap_density_group)]
  color_values <- .sirap_discrete_palette_values(levels, other_label = other_label, other_color = other_color)

  if (identical(plot_type, "histogram")) {
    surfaces <- list()
    bandwidth <- NULL
    histograms <- .sirap_clone_histogram_3d_bars(
      region_bcrs = region_bcrs,
      levels = levels,
      color_values = color_values,
      xlim = region_xlim,
      ylim = region_ylim,
      hist_bins = hist_bins,
      z_base = z_base,
      max_height = max_height
    )
  } else {
    histograms <- list()
    grid <- .sirap_density_grid(region_xlim, region_ylim, grid_n)
    bandwidth <- .sirap_density_bandwidth(region_bcrs$x, region_bcrs$.sirap_y, bandwidth)
    raw_surfaces <- lapply(levels, function(level) {
      dat <- region_bcrs[as.character(region_bcrs$.sirap_density_group) == level, , drop = FALSE]
      .sirap_weighted_gaussian_grid(
        x = dat$x,
        y = dat$.sirap_y,
        weights = dat$.sirap_density_weight,
        x_grid = grid$x,
        y_grid = grid$y,
        bandwidth = bandwidth
      )
    })
    names(raw_surfaces) <- levels
    max_density <- max(vapply(raw_surfaces, function(z) max(z, na.rm = TRUE), numeric(1)), na.rm = TRUE)
    surfaces <- lapply(levels, function(level) {
      z <- raw_surfaces[[level]]
      surface_max <- max(z, na.rm = TRUE)
      if (is.finite(surface_max) && surface_max > 0 && density_floor > 0) {
        z[z < surface_max * density_floor] <- NA_real_
      }
      if (is.finite(max_density) && max_density > 0) {
        z <- z / max_density * max_height + z_base
      } else {
        z[] <- z_base
      }
      list(
        level = level,
        x = grid$x,
        y = grid$y,
        z = z,
        color = unname(color_values[level])
      )
    })
    names(surfaces) <- levels
  }

  list(
    image = image,
    image_extent = image_extent,
    region_xlim = region_xlim,
    region_ylim = region_ylim,
    plot_type = plot_type,
    surfaces = surfaces,
    histograms = histograms,
    receptors = region_bcrs,
    polygons = polygons,
    colors = color_values,
    z_base = z_base,
    z_top = z_base + max_height,
    bandwidth = bandwidth,
    clone_col = clone_col
  )
}

.sirap_draw_clone_density_3d <- function(dat,
                                         density_alpha,
                                         image_alpha,
                                         show_points,
                                         point_size,
                                         point_alpha,
                                         show_polygon,
                                         polygon_color,
                                         polygon_lwd,
                                         label_density_peaks,
                                         label_top_n,
                                         label_cex,
                                         show_axes,
                                         texture_dir,
                                         as_widget) {
  if (isTRUE(as_widget)) {
    old_rgl_use_null <- getOption("rgl.useNULL")
    options(rgl.useNULL = TRUE)
    on.exit(options(rgl.useNULL = old_rgl_use_null), add = TRUE)
  }
  rgl::open3d(useNULL = isTRUE(as_widget))
  rgl::clear3d()

  texture <- .sirap_image_texture_file(dat$image, dat$image_extent, alpha = image_alpha, texture_dir = texture_dir)
  x_plane <- dat$image_extent$xlim
  y_plane <- dat$image_extent$ylim
  z_plane <- matrix(0, nrow = 2, ncol = 2)
  texture_s <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  texture_t <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  rgl::surface3d(
    x = x_plane,
    y = y_plane,
    z = z_plane,
    texture_s = texture_s,
    texture_t = texture_t,
    texture = texture,
    color = "white",
    alpha = image_alpha,
    lit = FALSE,
    back = "filled"
  )

  for (surface in dat$surfaces) {
    rgl::surface3d(
      x = surface$x,
      y = surface$y,
      z = surface$z,
      color = surface$color,
      alpha = density_alpha,
      lit = FALSE,
      back = "filled"
    )
  }
  if (identical(dat$plot_type, "histogram")) {
    for (histogram in dat$histograms) {
      .sirap_draw_histogram_bars_3d(histogram, alpha = density_alpha)
    }
  }

  if (isTRUE(show_polygon) && nrow(dat$polygons) > 0L) {
    for (polygon_id in unique(dat$polygons$polygon)) {
      polygon <- dat$polygons[dat$polygons$polygon == polygon_id, , drop = FALSE]
      rgl::lines3d(
        polygon$x,
        polygon$.sirap_y,
        rep(dat$z_base * 0.35, nrow(polygon)),
        color = polygon_color,
        lwd = polygon_lwd
      )
    }
  }

  if (isTRUE(label_density_peaks)) {
    labeled_items <- if (identical(dat$plot_type, "histogram")) {
      utils::head(dat$histograms, label_top_n)
    } else {
      utils::head(dat$surfaces, label_top_n)
    }
    for (item in labeled_items) {
      peak <- if (identical(dat$plot_type, "histogram")) {
        .sirap_histogram_bar_peak(item)
      } else {
        .sirap_density_surface_peak(item)
      }
      if (!is.null(peak)) {
        rgl::text3d(
          peak$x,
          peak$y,
          peak$z + (dat$z_top - dat$z_base) * 0.05,
          texts = item$level,
          color = item$color,
          cex = label_cex,
          font = 2
        )
      }
    }
  }

  if (isTRUE(show_points)) {
    point_colors <- unname(dat$colors[as.character(dat$receptors$.sirap_density_group)])
    point_colors[is.na(point_colors)] <- "#222222"
    rgl::points3d(
      dat$receptors$x,
      dat$receptors$.sirap_y,
      rep(dat$z_base * 0.65, nrow(dat$receptors)),
      color = grDevices::adjustcolor(point_colors, alpha.f = point_alpha),
      size = point_size
    )
  }

  if (isTRUE(show_axes)) {
    xlim <- range(c(dat$image_extent$xlim, dat$region_xlim), na.rm = TRUE)
    ylim <- range(c(dat$image_extent$ylim, dat$region_ylim), na.rm = TRUE)
    rgl::bbox3d(
      xat = pretty(xlim, n = 4),
      yat = pretty(ylim, n = 4),
      zat = pretty(c(0, dat$z_top), n = 3),
      color = "#666666",
      emission = "#333333"
    )
    rgl::decorate3d(xlab = "x", ylab = "y", zlab = "clone density")
  }
  rgl::aspect3d(1, 1, 0.55)
  rgl::view3d(theta = 40, phi = 25, zoom = 0.72)

  if (isTRUE(as_widget)) {
    rgl::rglwidget()
  } else {
    rgl::cur3d()
  }
}

.sirap_image_dim <- function(image) {
  img <- SeuratObject::GetImage(image, mode = "raw")
  dims <- dim(img)
  list(height = dims[[1]], width = dims[[2]])
}

.sirap_density_surface_peak <- function(surface) {
  if (is.null(surface$z) || !any(is.finite(surface$z))) {
    return(NULL)
  }
  idx <- arrayInd(which.max(surface$z), dim(surface$z))
  list(
    x = surface$x[[idx[[1]]]],
    y = surface$y[[idx[[2]]]],
    z = surface$z[idx[[1]], idx[[2]]]
  )
}

.sirap_histogram_bar_peak <- function(histogram) {
  bars <- histogram$bars
  if (is.null(bars) || nrow(bars) == 0L || !any(is.finite(bars$z1))) {
    return(NULL)
  }
  i <- which.max(bars$z1)
  list(
    x = (bars$x0[[i]] + bars$x1[[i]]) / 2,
    y = (bars$y0[[i]] + bars$y1[[i]]) / 2,
    z = bars$z1[[i]]
  )
}

.sirap_density_plot_range <- function(x, fallback) {
  out <- range(x, na.rm = TRUE)
  if (!all(is.finite(out))) {
    return(fallback)
  }
  if (out[[1]] == out[[2]]) {
    out <- out + c(-1, 1)
  }
  pad <- max(diff(out) * 0.08, 1)
  out <- out + c(-pad, pad)
  c(max(min(fallback), out[[1]]), min(max(fallback), out[[2]]))
}

.sirap_density_grid <- function(xlim, ylim, grid_n) {
  if (!is.numeric(grid_n) || length(grid_n) > 2L || any(!is.finite(grid_n)) || any(grid_n < 3)) {
    rlang::abort("`grid_n` must be one or two numbers greater than or equal to 3.")
  }
  grid_n <- as.integer(grid_n)
  if (length(grid_n) == 1L) {
    grid_n <- rep(grid_n, 2L)
  }
  list(
    x = seq(xlim[[1]], xlim[[2]], length.out = grid_n[[1]]),
    y = seq(ylim[[1]], ylim[[2]], length.out = grid_n[[2]])
  )
}

.sirap_density_bandwidth <- function(x, y, bandwidth) {
  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) > 2L || any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
      rlang::abort("`bandwidth` must be NULL or one or two positive numbers.")
    }
    bandwidth <- as.numeric(bandwidth)
    if (length(bandwidth) == 1L) {
      bandwidth <- rep(bandwidth, 2L)
    }
    return(bandwidth)
  }
  estimate <- function(v) {
    v <- v[is.finite(v)]
    if (length(unique(v)) < 2L) {
      return(1)
    }
    bw <- stats::bw.nrd0(v)
    if (!is.finite(bw) || bw <= 0) {
      bw <- diff(range(v)) / 8
    }
    max(bw, 1)
  }
  c(estimate(x), estimate(y))
}

.sirap_weighted_gaussian_grid <- function(x, y, weights, x_grid, y_grid, bandwidth) {
  z <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
  if (length(x) == 0L) {
    return(z)
  }
  keep <- is.finite(x) & is.finite(y) & is.finite(weights) & weights > 0
  x <- x[keep]
  y <- y[keep]
  weights <- weights[keep]
  if (length(x) == 0L) {
    return(z)
  }
  for (i in seq_along(x)) {
    dx <- stats::dnorm((x_grid - x[[i]]) / bandwidth[[1]])
    dy <- stats::dnorm((y_grid - y[[i]]) / bandwidth[[2]])
    z <- z + weights[[i]] * tcrossprod(dx, dy)
  }
  z / sum(weights)
}

.sirap_clone_histogram_3d_bars <- function(region_bcrs,
                                           levels,
                                           color_values,
                                           xlim,
                                           ylim,
                                           hist_bins,
                                           z_base,
                                           max_height) {
  grid <- .sirap_density_grid(xlim, ylim, hist_bins + 1L)
  x_breaks <- grid$x
  y_breaks <- grid$y
  x_bins <- length(x_breaks) - 1L
  y_bins <- length(y_breaks) - 1L
  raw_counts <- lapply(levels, function(level) {
    dat <- region_bcrs[as.character(region_bcrs$.sirap_density_group) == level, , drop = FALSE]
    counts <- matrix(0, nrow = x_bins, ncol = y_bins)
    if (nrow(dat) == 0L) {
      return(counts)
    }
    xi <- findInterval(dat$x, x_breaks, all.inside = TRUE)
    yi <- findInterval(dat$.sirap_y, y_breaks, all.inside = TRUE)
    xi <- pmin(xi, x_bins)
    yi <- pmin(yi, y_bins)
    for (i in seq_len(nrow(dat))) {
      counts[xi[[i]], yi[[i]]] <- counts[xi[[i]], yi[[i]]] + dat$.sirap_density_weight[[i]]
    }
    counts
  })
  names(raw_counts) <- levels
  max_count <- max(vapply(raw_counts, function(x) max(x, na.rm = TRUE), numeric(1)), na.rm = TRUE)
  lapply(levels, function(level) {
    counts <- raw_counts[[level]]
    idx <- which(counts > 0, arr.ind = TRUE)
    if (nrow(idx) == 0L || !is.finite(max_count) || max_count <= 0) {
      bars <- data.frame(x0 = numeric(), x1 = numeric(), y0 = numeric(), y1 = numeric(), z0 = numeric(), z1 = numeric())
    } else {
      bars <- data.frame(
        x0 = x_breaks[idx[, 1]],
        x1 = x_breaks[idx[, 1] + 1L],
        y0 = y_breaks[idx[, 2]],
        y1 = y_breaks[idx[, 2] + 1L],
        z0 = z_base,
        z1 = z_base + counts[idx] / max_count * max_height,
        stringsAsFactors = FALSE
      )
    }
    list(
      level = level,
      bars = bars,
      color = unname(color_values[level])
    )
  }) |> stats::setNames(levels)
}

.sirap_draw_histogram_bars_3d <- function(histogram, alpha) {
  bars <- histogram$bars
  if (is.null(bars) || nrow(bars) == 0L) {
    return(invisible(NULL))
  }
  color <- grDevices::adjustcolor(histogram$color, alpha.f = alpha)
  for (i in seq_len(nrow(bars))) {
    vertices <- .sirap_histogram_bar_vertices(bars[i, , drop = FALSE])
    rgl::quads3d(
      vertices$x,
      vertices$y,
      vertices$z,
      color = color,
      lit = FALSE
    )
  }
  invisible(NULL)
}

.sirap_histogram_bar_vertices <- function(bar) {
  x0 <- bar$x0[[1]]
  x1 <- bar$x1[[1]]
  y0 <- bar$y0[[1]]
  y1 <- bar$y1[[1]]
  z0 <- bar$z0[[1]]
  z1 <- bar$z1[[1]]
  points <- rbind(
    c(x0, y0, z1), c(x1, y0, z1), c(x1, y1, z1), c(x0, y1, z1),
    c(x0, y0, z0), c(x0, y1, z0), c(x1, y1, z0), c(x1, y0, z0),
    c(x0, y0, z0), c(x1, y0, z0), c(x1, y0, z1), c(x0, y0, z1),
    c(x1, y0, z0), c(x1, y1, z0), c(x1, y1, z1), c(x1, y0, z1),
    c(x1, y1, z0), c(x0, y1, z0), c(x0, y1, z1), c(x1, y1, z1),
    c(x0, y1, z0), c(x0, y0, z0), c(x0, y0, z1), c(x0, y1, z1)
  )
  data.frame(x = points[, 1], y = points[, 2], z = points[, 3])
}

.sirap_image_texture_extent <- function(image_dim, xlim, ylim) {
  xlim <- range(xlim, na.rm = TRUE)
  ylim <- range(ylim, na.rm = TRUE)
  xlim <- c(max(0, xlim[[1]]), min(image_dim$width, xlim[[2]]))
  ylim <- c(max(-image_dim$height, ylim[[1]]), min(0, ylim[[2]]))
  if (diff(xlim) <= 0 || diff(ylim) <= 0) {
    xlim <- c(0, image_dim$width)
    ylim <- c(-image_dim$height, 0)
  }
  col_start <- max(1L, floor(xlim[[1]]) + 1L)
  col_end <- min(image_dim$width, ceiling(xlim[[2]]))
  row_start <- max(1L, floor(-ylim[[2]]) + 1L)
  row_end <- min(image_dim$height, ceiling(-ylim[[1]]))
  list(
    cols = seq.int(col_start, col_end),
    rows = seq.int(row_start, row_end),
    xlim = c(col_start - 1, col_end),
    ylim = c(-row_end, -(row_start - 1))
  )
}

.sirap_image_texture_file <- function(image, extent, alpha = 1, texture_dir = NULL) {
  img <- SeuratObject::GetImage(image, mode = "raw")
  crop <- img[extent$rows, extent$cols, , drop = FALSE]
  if (alpha < 1) {
    rgba <- array(NA_real_, dim = c(dim(crop)[1:2], 4))
    rgba[, , 1:3] <- crop[, , 1:3]
    rgba[, , 4] <- alpha
    crop <- rgba
  }
  if (!is.null(texture_dir)) {
    dir.create(texture_dir, recursive = TRUE, showWarnings = FALSE)
    alpha_label <- gsub("[^0-9A-Za-z]+", "_", format(alpha, scientific = FALSE, trim = TRUE))
    texture <- file.path(
      texture_dir,
      paste0(
        "sirap_spatial_clone_density_texture",
        "_r", min(extent$rows), "-", max(extent$rows),
        "_c", min(extent$cols), "-", max(extent$cols),
        "_a", alpha_label,
        ".png"
      )
    )
    if (file.exists(texture)) {
      return(texture)
    }
  } else {
    texture <- tempfile(
      pattern = "sirap_spatial_clone_density_texture_",
      tmpdir = tempdir(),
      fileext = ".png"
    )
  }
  grDevices::png(texture, width = dim(crop)[[2]], height = dim(crop)[[1]], units = "px", bg = "transparent")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.raster(grDevices::as.raster(crop), interpolate = FALSE)
  texture
}

.sirap_polygon_pie_summary <- function(region_bcrs,
                                       fill_by,
                                       top_n,
                                       label_top_n,
                                       other_label,
                                       plot_other = TRUE) {
  if (nrow(region_bcrs) == 0L) {
    return(data.frame(
      polygon = character(),
      polygon_label = character(),
      level = character(),
      n = integer(),
      proportion = numeric(),
      ymin = numeric(),
      ymax = numeric(),
      label = character(),
      stringsAsFactors = FALSE
    ))
  }
  region_bcrs$.sirap_fill_by <- as.character(region_bcrs[[fill_by]])
  region_bcrs$.sirap_fill_by[is.na(region_bcrs$.sirap_fill_by) | !nzchar(region_bcrs$.sirap_fill_by)] <- NA_character_
  if (is.null(top_n)) {
    region_bcrs$.sirap_level <- region_bcrs$.sirap_fill_by
  } else {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      rlang::abort("`top_n` must be NULL or a positive number.")
    }
    top_levels <- unlist(lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
      counts <- sort(table(dat$.sirap_fill_by, useNA = "no"), decreasing = TRUE)
      names(counts)[seq_len(min(as.integer(top_n), length(counts)))]
    }), use.names = FALSE)
    top_levels <- unique(top_levels[!is.na(top_levels) & nzchar(top_levels)])
    region_bcrs$.sirap_level <- region_bcrs$.sirap_fill_by
    region_bcrs$.sirap_level[
      !is.na(region_bcrs$.sirap_level) &
        nzchar(region_bcrs$.sirap_level) &
        !region_bcrs$.sirap_level %in% top_levels
    ] <- other_label
  }
  pieces <- lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
    if (!isTRUE(plot_other)) {
      dat <- dat[dat$.sirap_level != other_label, , drop = FALSE]
    }
    if (nrow(dat) == 0L) {
      return(data.frame(level = character(), n = integer(), polygon = character()))
    }
    counts <- as.data.frame(table(dat$.sirap_level, useNA = "no"), stringsAsFactors = FALSE)
    names(counts) <- c("level", "n")
    counts <- counts[counts$n > 0L, , drop = FALSE]
    other <- counts[counts$level == other_label, , drop = FALSE]
    counts <- counts[counts$level != other_label, , drop = FALSE]
    counts <- counts[order(-counts$n, counts$level), , drop = FALSE]
    if (nrow(other) > 0L) {
      counts <- rbind(counts, other)
    }
    counts$.sirap_order <- seq_len(nrow(counts))
    counts$polygon <- unique(dat$polygon)[[1]]
    counts
  })
  out <- .sirap_bind_rows(pieces)
  if (nrow(out) == 0L) {
    return(out)
  }
  polygon_levels <- unique(out$polygon)
  polygon_labels <- stats::setNames(as.character(seq_along(polygon_levels)), polygon_levels)
  out$polygon_label <- unname(polygon_labels[out$polygon])
  out <- out[order(out$polygon_label, out$.sirap_order), , drop = FALSE]
  out <- do.call(rbind, lapply(split(out, out$polygon), function(dat) {
    total <- sum(dat$n)
    dat$proportion <- dat$n / total
    dat$ymax <- cumsum(dat$proportion)
    dat$ymin <- c(0, utils::head(dat$ymax, -1))
    dat$label_rank <- seq_len(nrow(dat))
    dat$label <- ifelse(dat$label_rank <= label_top_n, dat$level, "")
    dat
  }))
  out$.sirap_order <- NULL
  rownames(out) <- NULL
  out
}

.sirap_polygon_centers <- function(polygon_plot_data, polygons) {
  centers <- lapply(split(polygon_plot_data, polygon_plot_data$polygon), function(dat) {
    dat <- dat[!duplicated(paste(dat$x, dat$y)), , drop = FALSE]
    data.frame(
      polygon = unique(dat$polygon)[[1]],
      x = mean(range(dat$x, na.rm = TRUE)),
      y = mean(range(dat$y, na.rm = TRUE)),
      .sirap_y = -mean(range(dat$y, na.rm = TRUE)),
      radius = max(diff(range(dat$x, na.rm = TRUE)), diff(range(dat$y, na.rm = TRUE))) / 2,
      stringsAsFactors = FALSE
    )
  })
  centers <- .sirap_bind_rows(centers)
  centers <- centers[centers$polygon %in% polygons, , drop = FALSE]
  centers$polygon_label <- as.character(match(centers$polygon, unique(polygons)))
  centers
}

.sirap_side_by_side_polygon_pies <- function(object,
                                             images,
                                             image.scale,
                                             swap_coords,
                                             polygon_plot_data,
                                             centers,
                                             pie_data,
                                             fill_values,
                                             polygon_alpha,
                                             image.alpha,
                                             label_size,
                                             pies_per_page,
                                             zoom) {
  map_polygon_data <- polygon_plot_data
  polygon_labels <- stats::setNames(centers$polygon_label, centers$polygon)
  map_polygon_data$polygon_label <- unname(polygon_labels[as.character(map_polygon_data$polygon)])
  map_plots <- plotSpatialPolygons(
    object = object,
    polygon = map_polygon_data,
    images = images,
    image.scale = image.scale,
    swap_coords = swap_coords,
    polygon_col = "polygon",
    polygon_alpha = polygon_alpha,
    image.alpha = image.alpha,
    label = TRUE,
    label_col = "polygon_label",
    label_size = label_size,
    zoom = zoom,
    require_closed = TRUE,
    combine = FALSE
  )
  map_plot <- if (is.list(map_plots)) map_plots[[1L]] else map_plots

  pie_groups <- split(
    pie_data,
    factor(pie_data$polygon_label, levels = unique(pie_data$polygon_label))
  )
  pie_plots <- lapply(pie_groups, function(dat) {
    local_levels <- unique(as.character(dat$level))
    legend_ncol <- max(1L, min(4L, ceiling(length(local_levels) / 4)))
    ggplot2::ggplot(dat, ggplot2::aes(fill = .data$level)) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = 0.5, xmax = 1.5, ymin = .data$ymin, ymax = .data$ymax),
        color = "white",
        linewidth = 0.2
      ) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::geom_text(
        data = dat[dat$label != "", , drop = FALSE],
        ggplot2::aes(y = (.data$ymin + .data$ymax) / 2, label = .data$label),
        x = 1.25,
        inherit.aes = TRUE,
        size = label_size
      ) +
      ggplot2::scale_fill_manual(values = fill_values[local_levels], limits = local_levels, drop = TRUE) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, ncol = legend_ncol, byrow = TRUE)) +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 5),
        legend.key.size = grid::unit(0.12, "in"),
        legend.spacing.x = grid::unit(0.03, "in"),
        legend.spacing.y = grid::unit(0.01, "in"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      ) +
      ggplot2::labs(title = paste("Polygon", unique(dat$polygon_label)))
  })

  page_indices <- split(seq_along(pie_plots), ceiling(seq_along(pie_plots) / pies_per_page))
  page_names <- paste0("page_", seq_along(page_indices))
  names(page_indices) <- page_names

  if (requireNamespace("patchwork", quietly = TRUE)) {
    pages <- lapply(seq_along(page_indices), function(page_i) {
      idx <- page_indices[[page_i]]
      page_pies <- pie_plots[idx]
      pie_panel <- patchwork::wrap_plots(page_pies, ncol = max(1L, ceiling(sqrt(length(page_pies)))))
      page <- map_plot + pie_panel + patchwork::plot_layout(widths = c(1, 1))
      if (length(page_indices) > 1L) {
        page <- page + patchwork::plot_annotation(title = paste0("Page ", page_i, " of ", length(page_indices)))
      }
      page
    })
    names(pages) <- page_names
    if (length(pages) == 1L) pages[[1L]] else pages
  } else {
    pages <- lapply(page_indices, function(idx) {
      list(pies = pie_plots[idx], polygons = map_plot)
    })
    if (length(pages) == 1L) pages[[1L]] else pages
  }
}

.sirap_spatial_clipped_polygon_pies <- function(object,
                                                images,
                                                image.scale,
                                                swap_coords,
                                                polygon_plot_data,
                                                centers,
                                                pie_data,
                                                fill_values,
                                                fill_by,
                                                image.alpha,
                                                pie_alpha,
                                                label_size,
                                                spatial_legend,
                                                zoom) {
  .require_namespace("sf")
  clipped <- .sirap_clipped_pie_polygons(polygon_plot_data, centers, pie_data)
  labels <- .sirap_spatial_pie_labels(centers, pie_data)
  legend_ncol <- max(1L, min(4L, ceiling(length(fill_values) / 8)))
  plotted_polygons <- unique(as.character(pie_data$polygon))
  zoom_data <- polygon_plot_data[as.character(polygon_plot_data$polygon) %in% plotted_polygons, , drop = FALSE]
  full_image_data <- data.frame(x = numeric(), .sirap_y = numeric())

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    plots[[image]] <- .sirap_spatial_base(
      plot_data = if (isTRUE(zoom)) zoom_data else full_image_data,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = NULL
    ) +
      ggplot2::geom_polygon(
        data = clipped,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$piece, fill = .data$level),
        inherit.aes = FALSE,
        color = "white",
        linewidth = 0.15,
        alpha = pie_alpha
      ) +
      ggplot2::geom_path(
        data = polygon_plot_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      ggplot2::geom_text(
        data = labels,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, label = .data$label),
        inherit.aes = FALSE,
        size = label_size
      ) +
      ggplot2::scale_fill_manual(values = fill_values, limits = unique(as.character(pie_data$level)), drop = TRUE) +
      ggplot2::labs(fill = fill_by)
    if (isTRUE(spatial_legend)) {
      plots[[image]] <- plots[[image]] +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE)) +
        ggplot2::theme(
          legend.position = "right",
          legend.text = ggplot2::element_text(size = 6),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = grid::unit(0.14, "in"),
          legend.spacing.x = grid::unit(0.04, "in"),
          legend.spacing.y = grid::unit(0.02, "in")
        )
    } else {
      plots[[image]] <- plots[[image]] + ggplot2::guides(fill = "none")
    }
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

.sirap_clipped_pie_polygons <- function(polygon_plot_data, centers, pie_data) {
  pieces <- list()
  idx <- 1L
  for (polygon_id in unique(pie_data$polygon)) {
    polygon_vertices <- polygon_plot_data[polygon_plot_data$polygon == polygon_id, , drop = FALSE]
    center <- centers[centers$polygon == polygon_id, , drop = FALSE]
    if (nrow(polygon_vertices) == 0L || nrow(center) == 0L) {
      next
    }
    boundary <- .sirap_repair_polygon_sfc(
      sf::st_sfc(sf::st_polygon(list(as.matrix(polygon_vertices[, c("x", "y")]))))
    )
    if (is.null(boundary)) {
      rlang::warn(paste0("Polygon ", polygon_id, " is invalid and could not be repaired for clipped pie plotting."))
      next
    }
    this_pie <- pie_data[pie_data$polygon == polygon_id, , drop = FALSE]
    this_pie <- this_pie[order(this_pie$ymin), , drop = FALSE]
    for (i in seq_len(nrow(this_pie))) {
      sector <- .sirap_pie_sector(center$x[[1]], center$y[[1]], center$radius[[1]] * 1.25, this_pie$ymin[[i]], this_pie$ymax[[i]])
      sector_sf <- .sirap_repair_polygon_sfc(
        sf::st_sfc(sf::st_polygon(list(as.matrix(sector[, c("x", "y")]))))
      )
      if (is.null(sector_sf)) {
        rlang::warn(paste0(
          "Could not create a valid pie slice for polygon ",
          polygon_id,
          " and level ",
          this_pie$level[[i]],
          "."
        ))
        next
      }
      intersection <- tryCatch(
        suppressWarnings(.sirap_polygon_intersection(sector_sf, boundary)),
        error = function(e) {
          rlang::warn(paste0(
            "Could not clip pie slice for polygon ",
            polygon_id,
            " and level ",
            this_pie$level[[i]],
            ": ",
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (is.null(intersection)) {
        next
      }
      if (length(intersection) == 0L || sf::st_is_empty(intersection)[[1]]) {
        next
      }
      coords <- as.data.frame(sf::st_coordinates(intersection))
      if (!"L1" %in% colnames(coords)) {
        coords$L1 <- 1L
      }
      coords$piece <- paste(polygon_id, this_pie$level[[i]], coords$L1, idx, sep = "__")
      coords$polygon <- polygon_id
      coords$level <- this_pie$level[[i]]
      coords$x <- coords$X
      coords$y <- coords$Y
      coords$.sirap_y <- -coords$y
      pieces[[idx]] <- coords[, c("polygon", "level", "piece", "x", "y", ".sirap_y"), drop = FALSE]
      idx <- idx + 1L
    }
  }
  .sirap_bind_rows(pieces)
}

.sirap_repair_polygon_sfc <- function(geometry) {
  repaired <- tryCatch({
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      geometry <- sf::st_make_valid(geometry)
    } else if (!all(sf::st_is_valid(geometry))) {
      geometry <- sf::st_buffer(geometry, 0)
    }
    geometry <- suppressWarnings(sf::st_collection_extract(geometry, "POLYGON"))
    geometry <- suppressWarnings(sf::st_cast(geometry, "POLYGON"))
    if (length(geometry) == 0L || all(sf::st_is_empty(geometry))) {
      return(NULL)
    }
    if (!all(sf::st_is_valid(geometry))) {
      geometry <- suppressWarnings(sf::st_buffer(geometry, 0))
      geometry <- suppressWarnings(sf::st_collection_extract(geometry, "POLYGON"))
      geometry <- suppressWarnings(sf::st_cast(geometry, "POLYGON"))
    }
    if (length(geometry) == 0L || all(sf::st_is_empty(geometry)) || !all(sf::st_is_valid(geometry))) {
      return(NULL)
    }
    geometry
  }, error = function(e) NULL)
  repaired
}

.sirap_polygon_intersection <- function(x, y) {
  tryCatch(
    suppressWarnings(sf::st_intersection(x, y)),
    error = function(e) {
      x <- .sirap_repair_polygon_sfc(suppressWarnings(sf::st_buffer(x, 0)))
      y <- .sirap_repair_polygon_sfc(suppressWarnings(sf::st_buffer(y, 0)))
      if (is.null(x) || is.null(y)) {
        stop(e)
      }
      suppressWarnings(sf::st_intersection(x, y))
    }
  )
}

.sirap_pie_sector <- function(center_x, center_y, radius, ymin, ymax, n = 80L) {
  start <- (ymin * 2 * pi) - (pi / 2)
  end <- (ymax * 2 * pi) - (pi / 2)
  angles <- seq(start, end, length.out = max(3L, ceiling(n * (ymax - ymin))))
  data.frame(
    x = c(center_x, center_x + radius * cos(angles), center_x),
    y = c(center_y, center_y + radius * sin(angles), center_y),
    stringsAsFactors = FALSE
  )
}

.sirap_spatial_pie_labels <- function(centers, pie_data) {
  labels <- pie_data[pie_data$label != "", , drop = FALSE]
  if (nrow(labels) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), .sirap_y = numeric(), label = character()))
  }
  labels <- merge(labels, centers[, c("polygon", "x", "y", "radius"), drop = FALSE], by = "polygon", sort = FALSE)
  mid <- ((labels$ymin + labels$ymax) / 2) * 2 * pi - (pi / 2)
  labels$x <- labels$x + labels$radius * 0.45 * cos(mid)
  labels$y <- labels$y + labels$radius * 0.45 * sin(mid)
  labels$.sirap_y <- -labels$y
  labels[, c("x", "y", ".sirap_y", "label"), drop = FALSE]
}

.sirap_downsample_polygon_plot_coords <- function(coords, downsample_n = NULL, seed = 1L) {
  if (is.null(downsample_n)) {
    return(coords)
  }
  if (!is.numeric(downsample_n) || length(downsample_n) != 1L || is.na(downsample_n) || downsample_n < 1) {
    rlang::abort("`downsample_n` must be NULL or a positive number.")
  }
  downsample_n <- as.integer(downsample_n)
  if (nrow(coords) <= downsample_n) {
    return(coords)
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  coords[sample.int(nrow(coords), downsample_n), , drop = FALSE]
}

.sirap_close_current_polygon <- function(points, polygon_id) {
  current <- points[points$polygon == polygon_id, , drop = FALSE]
  if (nrow(current) == 0L) {
    return(points)
  }
  if (nrow(current) < 3L) {
    rlang::warn("A polygon needs at least three distinct vertices before it can be closed.")
    return(points)
  }
  if (.sirap_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
    return(points)
  }
  rbind(points, current[1, , drop = FALSE])
}

.sirap_default_polygon_buffer_radius <- function(coords) {
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  coords <- unique(coords[, c("x", "y"), drop = FALSE])
  if (nrow(coords) < 2L) {
    rlang::abort("At least two spatial coordinates are required to estimate `buffer_radius`.")
  }
  nn <- RANN::nn2(
    data = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    query = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    k = 2L,
    treetype = "kd",
    searchtype = "standard"
  )
  distances <- nn$nn.dists[, 2L]
  distances <- distances[is.finite(distances) & distances > 0]
  if (length(distances) == 0L) {
    rlang::abort("Could not estimate `buffer_radius` from the spatial coordinates.")
  }
  stats::median(distances) * 0.55
}

.sirap_polygon_interior_to_sf_boundary <- function(points, buffer_radius) {
  if (!is.numeric(buffer_radius) || length(buffer_radius) != 1L ||
      is.na(buffer_radius) || buffer_radius <= 0) {
    rlang::abort("`buffer_radius` must be a positive number.")
  }
  boundary_tables <- lapply(unique(points$polygon), function(polygon_id) {
    current <- points[points$polygon == polygon_id, , drop = FALSE]
    if (nrow(current) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three points and will be excluded."))
      return(NULL)
    }

    point_sf <- sf::st_as_sf(
      current,
      coords = c("x", "y"),
      crs = NA,
      remove = FALSE
    )
    buffered <- suppressWarnings(sf::st_buffer(point_sf, dist = buffer_radius))
    unioned <- suppressWarnings(sf::st_union(sf::st_geometry(buffered)))
    polygons <- suppressWarnings(sf::st_collection_extract(unioned, "POLYGON"))
    polygons <- suppressWarnings(sf::st_cast(polygons, "POLYGON"))

    if (length(polygons) == 0L || all(sf::st_is_empty(polygons))) {
      rlang::warn(paste0(
        "No polygon boundary could be created for polygon ",
        polygon_id,
        "; it will be excluded."
      ))
      return(NULL)
    }

    if (length(polygons) > 1L) {
      polygon_areas <- as.numeric(sf::st_area(polygons))
      polygons <- polygons[which.max(polygon_areas)]
      rlang::warn(paste0(
        "Polygon ",
        polygon_id,
        " produced multiple disconnected outlines; keeping the largest outline."
      ))
    }

    boundary <- as.data.frame(sf::st_coordinates(polygons))
    if ("L1" %in% colnames(boundary)) {
      boundary <- boundary[boundary$L1 == 1L, , drop = FALSE]
    }
    boundary <- data.frame(
      x = as.numeric(boundary$X),
      y = as.numeric(boundary$Y),
      stringsAsFactors = FALSE
    )
    if (nrow(boundary) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three boundary vertices and will be excluded."))
      return(NULL)
    }

    nearest <- RANN::nn2(
      data = as.matrix(current[, c("x", "y"), drop = FALSE]),
      query = as.matrix(boundary[, c("x", "y"), drop = FALSE]),
      k = 1L,
      treetype = "kd",
      searchtype = "standard"
    )
    nearest_idx <- nearest$nn.idx[, 1L]
    nearest_idx <- nearest_idx[!is.na(nearest_idx)]
    if (length(nearest_idx) == 0L) {
      rlang::warn(paste0("No boundary barcodes could be identified for polygon ", polygon_id, "."))
      return(NULL)
    }
    nearest_idx <- nearest_idx[c(TRUE, diff(nearest_idx) != 0)]
    nearest_idx <- nearest_idx[!duplicated(nearest_idx)]
    if (length(nearest_idx) < 3L) {
      rlang::warn(paste0("Fewer than three boundary barcodes were found for polygon ", polygon_id, "."))
      return(NULL)
    }
    boundary_points <- current[nearest_idx, , drop = FALSE]
    boundary_points <- .sirap_clean_boundary_barcode_vertices(boundary_points)
    if (nrow(boundary_points) < 3L) {
      rlang::warn(paste0("Fewer than three non-degenerate boundary barcodes were found for polygon ", polygon_id, "."))
      return(NULL)
    }
    if (!.sirap_boundary_barcode_vertices_valid(boundary_points)) {
      boundary_points <- .sirap_repair_boundary_barcode_vertices(
        boundary_points = boundary_points,
        candidate_points = current
      )
    }
    if (!.sirap_boundary_barcode_vertices_valid(boundary_points)) {
      rlang::warn(paste0(
        "Polygon ",
        polygon_id,
        " could not be represented as a valid barcode-coordinate boundary and will be excluded."
      ))
      return(NULL)
    }

    data.frame(
      barcode = boundary_points$barcode,
      x = boundary_points$x,
      y = boundary_points$y,
      polygon = boundary_points$polygon,
      image = boundary_points$image,
      stringsAsFactors = FALSE
    )
  })
  boundary_tables <- Filter(Negate(is.null), boundary_tables)
  if (length(boundary_tables) == 0L) {
    rlang::abort("No polygons with at least three points remain.")
  }
  out <- do.call(rbind, boundary_tables)
  rownames(out) <- NULL
  out
}

.sirap_clean_boundary_barcode_vertices <- function(points, tol = 1e-8) {
  points <- points[is.finite(points$x) & is.finite(points$y), , drop = FALSE]
  if (nrow(points) < 3L) {
    return(points)
  }

  changed <- TRUE
  while (isTRUE(changed) && nrow(points) >= 3L) {
    old_n <- nrow(points)
    if (nrow(points) > 1L && .sirap_same_xy(points[1, , drop = FALSE], points[nrow(points), , drop = FALSE])) {
      points <- points[-nrow(points), , drop = FALSE]
    }
    if (nrow(points) > 1L) {
      duplicate_previous <- c(
        FALSE,
        abs(diff(points$x)) <= tol & abs(diff(points$y)) <= tol
      )
      points <- points[!duplicate_previous, , drop = FALSE]
    }
    if (nrow(points) >= 4L) {
      remove <- logical(nrow(points))
      for (i in seq_len(nrow(points))) {
        prev_i <- if (i == 1L) nrow(points) else i - 1L
        next_i <- if (i == nrow(points)) 1L else i + 1L
        remove[[i]] <- .sirap_degenerate_boundary_vertex(
          prev = points[prev_i, , drop = FALSE],
          current = points[i, , drop = FALSE],
          next_point = points[next_i, , drop = FALSE],
          tol = tol
        )
      }
      points <- points[!remove, , drop = FALSE]
    }
    changed <- nrow(points) != old_n
  }
  rownames(points) <- NULL
  points
}

.sirap_degenerate_boundary_vertex <- function(prev, current, next_point, tol = 1e-8) {
  if (.sirap_same_xy(prev, next_point)) {
    return(TRUE)
  }
  v1 <- c(prev$x[[1]] - current$x[[1]], prev$y[[1]] - current$y[[1]])
  v2 <- c(next_point$x[[1]] - current$x[[1]], next_point$y[[1]] - current$y[[1]])
  cross <- v1[[1]] * v2[[2]] - v1[[2]] * v2[[1]]
  dot <- sum(v1 * v2)
  isTRUE(abs(cross) <= tol && dot > tol)
}

.sirap_closed_boundary_matrix <- function(points) {
  xy <- as.matrix(points[, c("x", "y"), drop = FALSE])
  if (nrow(xy) == 0L) {
    return(xy)
  }
  if (!all(xy[1L, ] == xy[nrow(xy), ])) {
    xy <- rbind(xy, xy[1L, ])
  }
  xy
}

.sirap_boundary_barcode_vertices_valid <- function(points) {
  points <- .sirap_clean_boundary_barcode_vertices(points)
  if (nrow(unique(points[, c("x", "y"), drop = FALSE])) < 3L) {
    return(FALSE)
  }
  xy <- .sirap_closed_boundary_matrix(points)
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(xy))),
    error = function(e) NULL
  )
  if (is.null(polygon) || length(polygon) == 0L || any(sf::st_is_empty(polygon))) {
    return(FALSE)
  }
  isTRUE(all(sf::st_is_valid(polygon))) && isTRUE(as.numeric(sf::st_area(polygon)) > 0)
}

.sirap_repair_boundary_barcode_vertices <- function(boundary_points, candidate_points) {
  if (nrow(boundary_points) < 3L || nrow(candidate_points) < 3L) {
    return(boundary_points)
  }
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(.sirap_closed_boundary_matrix(boundary_points)))),
    error = function(e) NULL
  )
  if (is.null(polygon)) {
    return(boundary_points)
  }
  repaired <- if ("st_make_valid" %in% getNamespaceExports("sf")) {
    sf::st_make_valid(polygon)
  } else {
    sf::st_buffer(polygon, 0)
  }
  repaired <- suppressWarnings(sf::st_collection_extract(repaired, "POLYGON"))
  repaired <- suppressWarnings(sf::st_cast(repaired, "POLYGON"))
  if (length(repaired) == 0L || all(sf::st_is_empty(repaired))) {
    return(boundary_points)
  }
  if (length(repaired) > 1L) {
    repaired <- repaired[which.max(as.numeric(sf::st_area(repaired)))]
  }
  boundary <- as.data.frame(sf::st_coordinates(repaired))
  if ("L1" %in% colnames(boundary)) {
    boundary <- boundary[boundary$L1 == 1L, , drop = FALSE]
  }
  boundary <- data.frame(
    x = as.numeric(boundary$X),
    y = as.numeric(boundary$Y),
    stringsAsFactors = FALSE
  )
  nearest <- RANN::nn2(
    data = as.matrix(candidate_points[, c("x", "y"), drop = FALSE]),
    query = as.matrix(boundary[, c("x", "y"), drop = FALSE]),
    k = 1L,
    treetype = "kd",
    searchtype = "standard"
  )
  nearest_idx <- nearest$nn.idx[, 1L]
  nearest_idx <- nearest_idx[!is.na(nearest_idx)]
  nearest_idx <- nearest_idx[c(TRUE, diff(nearest_idx) != 0)]
  nearest_idx <- nearest_idx[!duplicated(nearest_idx)]
  .sirap_clean_boundary_barcode_vertices(candidate_points[nearest_idx, , drop = FALSE])
}

.sirap_close_formatted_polygons <- function(points) {
  closed <- lapply(unique(points$polygon), function(polygon_id) {
    current <- points[points$polygon == polygon_id, , drop = FALSE]
    if (nrow(current) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three points and will be excluded."))
      return(NULL)
    }
    if (.sirap_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
      return(current)
    }
    rbind(current, current[1, , drop = FALSE])
  })
  closed <- Filter(Negate(is.null), closed)
  if (length(closed) == 0L) {
    rlang::abort("No polygons with at least three points remain.")
  }
  out <- do.call(rbind, closed)
  rownames(out) <- NULL
  out
}
