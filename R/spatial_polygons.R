#' Calculate BCR clone dominance inside spatial polygons
#'
#' Finds Seurat spatial barcodes whose coordinates fall inside closed polygons,
#' subsets BCR rows to those barcodes, and returns the fraction of regional BCRs
#' represented by each clone. Points on the polygon boundary are included by
#' default.
#'
#' @param bcrs BCR/AIRR-style data frame containing clone assignments and a
#'   spatial barcode column.
#' @param object Seurat object with spatial coordinates, or metadata coordinate
#'   columns supported by `spatialVariableData()`.
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
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists. When `NULL`, all polygons in `polygon_col` are analyzed.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @param require_closed Require the final polygon vertex to equal the first.
#'   Set to `FALSE` to close open polygons automatically.
#' @param return_region Return a list containing the dominance table, regional
#'   barcodes, and regional coordinate rows.
#' @return A data frame with `clone_id`, `n_bcrs`, `proportion`, regional
#'   totals, and the polygon ID when multiple polygons are supplied. Returns a
#'   list when `return_region = TRUE`.
#' @export
spatialCloneDominance <- function(bcrs,
                                  object,
                                  polygon,
                                  barcode_col = NULL,
                                  clone_col = "clone_id",
                                  locus = NULL,
                                  locus_col = "locus",
                                  count_col = NULL,
                                  v_call_col = "v_call",
                                  images = NULL,
                                  image.scale = "hires",
                                  polygon_id = NULL,
                                  x_col = "x",
                                  y_col = "y",
                                  polygon_col = "polygon",
                                  include_boundary = TRUE,
                                  require_closed = TRUE,
                                  return_region = FALSE) {
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

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)

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
    if (!is.null(locus)) {
      regional_bcrs <- regional_bcrs[regional_bcrs[[locus_col]] %in% locus, , drop = FALSE]
    }
    regional_bcrs <- regional_bcrs[
      !is.na(regional_bcrs[[clone_col]]) &
        as.character(regional_bcrs[[clone_col]]) != "",
      ,
      drop = FALSE
    ]

    out <- .sirap_clone_dominance_table(
      regional_bcrs = regional_bcrs,
      clone_col = clone_col,
      count_col = count_col,
      v_call_col = v_call_col,
      n_region_barcodes = length(region_barcodes)
    )
    if (!is.null(locus)) {
      out$locus_filter <- rep(paste(as.character(locus), collapse = ","), nrow(out))
    }
    if (isTRUE(spec$add_polygon_col)) {
      out[[polygon_col]] <- rep(spec$id, nrow(out))
    }

    if (isTRUE(spec$add_polygon_col)) {
      region_barcode_table <- data.frame(
        barcode = region_barcodes,
        stringsAsFactors = FALSE
      )
      region_barcode_table[[polygon_col]] <- rep(spec$id, nrow(region_barcode_table))
      region_coords[[polygon_col]] <- rep(spec$id, nrow(region_coords))
    } else {
      region_barcode_table <- data.frame(
        barcode = region_barcodes,
        stringsAsFactors = FALSE
      )
    }
    list(
      dominance = out,
      region_barcodes = region_barcode_table,
      region_coordinates = region_coords
    )
  })

  out <- .sirap_rbind(results, "dominance")

  if (isTRUE(return_region)) {
    return(list(
      dominance = out,
      region_barcodes = .sirap_rbind(results, "region_barcodes"),
      region_coordinates = .sirap_rbind(results, "region_coordinates")
    ))
  }
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
    suppressPackageStartupMessages(library(shiny))
  }
  old_shiny_useragg <- getOption("shiny.useragg")
  options(shiny.useragg = FALSE)
  on.exit(options(shiny.useragg = old_shiny_useragg), add = TRUE)
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale)
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
  plot_image_info <- .sirap_polygon_plot_image(object, coords)
  plot_height <- .sirap_polygon_plot_height(plot_coords)
  default_feature <- feature
  if (is.null(default_feature)) {
    default_feature <- .sirap_first_feature(object)
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
          ave(selected$polygon, selected$polygon, FUN = length) > 1L,
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
      closest <- .sirap_nearest_spatial_row(coords, click_data$x, click_data$y, y_col = ".sirap_plot_y")
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
      .sirap_display_polygon_vertices(current)
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
      out <- .sirap_close_all_polygons(selected_points())
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
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param polygon_alpha Fill opacity for polygons.
#' @param fill Fill color for polygons.
#' @param color Outline color for polygons.
#' @param linewidth Outline width for polygons.
#' @param image.alpha Tissue image opacity.
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
                                polygon_id = NULL,
                                x_col = "x",
                                y_col = "y",
                                polygon_col = "polygon",
                                polygon_alpha = 0.2,
                                fill = "#0072B2",
                                color = "#0072B2",
                                linewidth = 0.6,
                                image.alpha = 1,
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
    require_closed = require_closed
  )

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    plots[[image]] <- .sirap_spatial_base(
      plot_data = polygon_data,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_polygon(
        data = polygon_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$.sirap_polygon_group),
        inherit.aes = FALSE,
        fill = fill,
        color = color,
        alpha = polygon_alpha,
        linewidth = linewidth
      )
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot BCRs inside spatial polygons
#'
#' Plots every BCR row whose spatial barcode falls inside one or more selected
#' polygons over the tissue image, colored by a BCR column.
#'
#' @param bcrs BCR/AIRR-style data frame containing a spatial barcode column.
#' @param object Seurat object with spatial images.
#' @param polygon Data frame of polygon vertices in order.
#' @param color_by Column in `bcrs` used to color BCR points.
#' @param barcode_col Column in `bcrs` matching Seurat spot/cell barcodes. If
#'   `NULL`, a barcode-like column is guessed.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param point_size Size of BCR points.
#' @param point_alpha Opacity for BCR points.
#' @param polygon_alpha Fill opacity for polygons.
#' @param polygon_fill Fill color for polygons.
#' @param polygon_color Outline color for polygons.
#' @param image.alpha Tissue image opacity.
#' @param top_n Optional number of most frequent `color_by` levels to keep
#'   before grouping remaining levels as `other_label`.
#' @param other_label Label for grouped lower-frequency levels.
#' @param require_closed Require the final polygon vertex to equal the first.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object or a named list of ggplots when
#'   `combine = FALSE`.
#' @export
plotSpatialRegionBCRs <- function(bcrs,
                                  object,
                                  polygon,
                                  color_by = "clone_id",
                                  barcode_col = NULL,
                                  images = NULL,
                                  image.scale = "hires",
                                  polygon_id = NULL,
                                  x_col = "x",
                                  y_col = "y",
                                  polygon_col = "polygon",
                                  point_size = 1,
                                  point_alpha = 0.7,
                                  polygon_alpha = 0.12,
                                  polygon_fill = "#FFFFFF",
                                  polygon_color = "black",
                                  image.alpha = 1,
                                  top_n = NULL,
                                  other_label = "Other",
                                  require_closed = TRUE,
                                  combine = TRUE) {
  .require_namespace("SeuratObject")
  region <- .sirap_region_bcr_data(
    bcrs = bcrs,
    object = object,
    polygon = polygon,
    barcode_col = barcode_col,
    images = images,
    image.scale = image.scale,
    polygon_id = polygon_id,
    x_col = x_col,
    y_col = y_col,
    polygon_col = polygon_col,
    require_closed = require_closed
  )
  if (nrow(region$bcrs) == 0L) {
    rlang::abort("No BCR rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$bcrs, color_by, "bcrs")
  region$bcrs$.sirap_color_by <- .sirap_lump_values(region$bcrs[[color_by]], top_n = top_n, other_label = other_label)

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    plot_bcrs <- region$bcrs[region$bcrs$image == image | region$bcrs$image == "metadata", , drop = FALSE]
    plots[[image]] <- .sirap_spatial_base(
      plot_data = spatial$coords,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = FALSE,
      view_data = spatial$coords
    ) +
      ggplot2::geom_polygon(
        data = region$polygon_plot_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$.sirap_polygon_group),
        inherit.aes = FALSE,
        fill = polygon_fill,
        color = polygon_color,
        alpha = polygon_alpha,
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data = plot_bcrs,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, color = .data$.sirap_color_by),
        inherit.aes = FALSE,
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::labs(color = color_by)
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot spatial polygon pie summaries
#'
#' Summarizes the composition of BCRs inside each polygon for a selected column.
#' The default mode plots pie charts beside a numbered polygon map. The
#' `spatial_clipped` mode clips each pie chart to its polygon boundary and
#' requires the suggested `sf` package.
#'
#' @param bcrs BCR/AIRR-style data frame containing a spatial barcode column.
#' @param object Seurat object with spatial images.
#' @param polygon Data frame of polygon vertices in order.
#' @param fill_by Column in `bcrs` summarized in each pie chart.
#' @param barcode_col Column in `bcrs` matching Seurat spot/cell barcodes. If
#'   `NULL`, a barcode-like column is guessed.
#' @param mode Plot style. `side_by_side` shows ordinary pies beside a polygon
#'   map. `spatial_clipped` clips pies to polygon boundaries on the tissue.
#' @param top_n Number of most frequent `fill_by` levels to show individually
#'   per polygon; remaining levels are grouped as `other_label`.
#' @param label_top_n Number of levels to label per polygon. Defaults to
#'   `top_n`.
#' @param other_label Label for grouped lower-frequency levels.
#' @param images Optional Seurat image names to show.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param polygon_id Optional polygon ID value(s) to select when `polygon_col`
#'   exists.
#' @param x_col,y_col Coordinate columns in `polygon`.
#' @param polygon_col Optional polygon grouping column in `polygon`.
#' @param polygon_alpha Fill opacity for polygons in the numbered map.
#' @param image.alpha Tissue image opacity.
#' @param pie_alpha Pie opacity.
#' @param label_size Text size for polygon and pie labels.
#' @param require_closed Require the final polygon vertex to equal the first.
#' @return A ggplot, patchwork object, or named list of ggplots.
#' @export
plotSpatialPolygonPies <- function(bcrs,
                                   object,
                                   polygon,
                                   fill_by = "clone_id",
                                   barcode_col = NULL,
                                   mode = c("side_by_side", "spatial_clipped"),
                                   top_n = 10,
                                   label_top_n = top_n,
                                   other_label = "Other",
                                   images = NULL,
                                   image.scale = "hires",
                                   polygon_id = NULL,
                                   x_col = "x",
                                   y_col = "y",
                                   polygon_col = "polygon",
                                   polygon_alpha = 0.2,
                                   image.alpha = 1,
                                   pie_alpha = 0.9,
                                   label_size = 3,
                                   require_closed = TRUE) {
  .require_namespace("SeuratObject")
  mode <- match.arg(mode)
  if (is.null(label_top_n)) {
    label_top_n <- Inf
  }
  if (!is.numeric(label_top_n) || length(label_top_n) != 1L || is.na(label_top_n) || label_top_n < 0) {
    rlang::abort("`label_top_n` must be NULL or a non-negative number.")
  }
  region <- .sirap_region_bcr_data(
    bcrs = bcrs,
    object = object,
    polygon = polygon,
    barcode_col = barcode_col,
    images = images,
    image.scale = image.scale,
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
    other_label = other_label
  )
  centers <- .sirap_polygon_centers(region$polygon_plot_data, pie_data$polygon)

  if (mode == "side_by_side") {
    .sirap_side_by_side_polygon_pies(
      object = object,
      images = images,
      image.scale = image.scale,
      polygon_plot_data = region$polygon_plot_data,
      centers = centers,
      pie_data = pie_data,
      polygon_alpha = polygon_alpha,
      image.alpha = image.alpha,
      label_size = label_size
    )
  } else {
    .sirap_spatial_clipped_polygon_pies(
      object = object,
      images = images,
      image.scale = image.scale,
      polygon_plot_data = region$polygon_plot_data,
      centers = centers,
      pie_data = pie_data,
      image.alpha = image.alpha,
      pie_alpha = pie_alpha,
      label_size = label_size
    )
  }
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

.sirap_clone_dominance_table <- function(regional_bcrs,
                                         clone_col,
                                         count_col,
                                         v_call_col,
                                         n_region_barcodes) {
  if (nrow(regional_bcrs) == 0L) {
    return(data.frame(
      clone_id = character(),
      v_gene = character(),
      v_gene_family = character(),
      v_gene_superfamily = character(),
      n_bcrs = numeric(),
      proportion = numeric(),
      region_bcr_total = numeric(),
      n_region_barcodes = integer(),
      stringsAsFactors = FALSE
    ))
  }

  clone_ids <- as.character(regional_bcrs[[clone_col]])
  weights <- if (is.null(count_col)) {
    rep(1, nrow(regional_bcrs))
  } else {
    as.numeric(regional_bcrs[[count_col]])
  }
  weights[is.na(weights)] <- 0
  counts <- stats::aggregate(
    n_bcrs ~ clone_id,
    data = data.frame(clone_id = clone_ids, n_bcrs = weights, stringsAsFactors = FALSE),
    FUN = sum
  )
  counts <- counts[counts$n_bcrs > 0, , drop = FALSE]
  v_summary <- .sirap_clone_v_gene_summary(regional_bcrs, clone_col = clone_col, v_call_col = v_call_col)
  counts <- merge(counts, v_summary, by = "clone_id", all.x = TRUE, sort = FALSE)
  counts$v_gene[is.na(counts$v_gene)] <- NA_character_
  counts$v_gene_family[is.na(counts$v_gene_family)] <- NA_character_
  counts$v_gene_superfamily[is.na(counts$v_gene_superfamily)] <- NA_character_
  total <- sum(counts$n_bcrs)
  counts$proportion <- counts$n_bcrs / total
  counts$region_bcr_total <- total
  counts$n_region_barcodes <- n_region_barcodes
  counts <- counts[order(counts$proportion, counts$n_bcrs, decreasing = TRUE), , drop = FALSE]
  counts[, c("clone_id", "v_gene", "v_gene_family", "v_gene_superfamily", "n_bcrs", "proportion", "region_bcr_total", "n_region_barcodes"), drop = FALSE]
}

.sirap_clone_v_gene_summary <- function(regional_bcrs, clone_col, v_call_col) {
  clone_ids <- unique(as.character(regional_bcrs[[clone_col]]))
  if (!v_call_col %in% colnames(regional_bcrs)) {
    return(data.frame(
      clone_id = clone_ids,
      v_gene = rep(NA_character_, length(clone_ids)),
      v_gene_family = rep(NA_character_, length(clone_ids)),
      v_gene_superfamily = rep(NA_character_, length(clone_ids)),
      stringsAsFactors = FALSE
    ))
  }
  v_calls <- as.character(regional_bcrs[[v_call_col]])
  names(v_calls) <- as.character(regional_bcrs[[clone_col]])
  out <- lapply(clone_ids, function(clone_id) {
    calls <- v_calls[names(v_calls) == clone_id]
    data.frame(
      clone_id = clone_id,
      v_gene = .sirap_collapse_v_genes(calls),
      v_gene_family = .sirap_collapse_v_gene_families(calls),
      v_gene_superfamily = .sirap_collapse_v_gene_superfamilies(calls),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.sirap_collapse_v_genes <- function(v_calls) {
  genes <- .sirap_v_genes(v_calls)
  if (length(genes) == 0L) {
    return(NA_character_)
  }
  paste(genes, collapse = ",")
}

.sirap_collapse_v_gene_families <- function(v_calls) {
  families <- unique(stats::na.omit(.sirap_v_gene_family(.sirap_v_genes(v_calls))))
  families <- sort(families)
  if (length(families) == 0L) {
    return(NA_character_)
  }
  paste(families, collapse = ",")
}

.sirap_collapse_v_gene_superfamilies <- function(v_calls) {
  superfamilies <- unique(stats::na.omit(.sirap_v_gene_superfamily(.sirap_v_gene_family(.sirap_v_genes(v_calls)))))
  superfamilies <- sort(superfamilies)
  if (length(superfamilies) == 0L) {
    return(NA_character_)
  }
  paste(superfamilies, collapse = ",")
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

.sirap_display_polygon_vertices <- function(points) {
  out <- .sirap_public_polygon_vertices(points)
  out[, setdiff(colnames(out), "image"), drop = FALSE]
}

.sirap_polygon_plot_data <- function(polygon,
                                     polygon_id,
                                     polygon_col,
                                     x_col,
                                     y_col,
                                     require_closed) {
  specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  out <- lapply(specs, function(spec) {
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    vertices$.sirap_polygon_group <- if (isTRUE(spec$add_polygon_col)) {
      as.character(spec$id)
    } else {
      "polygon_1"
    }
    vertices
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out$.sirap_y <- -out$y
  out
}

.sirap_region_bcr_data <- function(bcrs,
                                   object,
                                   polygon,
                                   barcode_col,
                                   images,
                                   image.scale,
                                   polygon_id,
                                   x_col,
                                   y_col,
                                   polygon_col,
                                   require_closed,
                                   include_boundary = TRUE) {
  bcrs <- .as_data_frame(bcrs)
  polygon <- .as_data_frame(polygon)
  barcode_col <- .guess_barcode_col(bcrs, barcode_col)
  .check_required_cols(bcrs, barcode_col, "bcrs")
  .check_required_cols(polygon, c(x_col, y_col), "polygon")

  coords <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale)
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
      include_boundary = include_boundary
    )
    region_coords <- coords[inside, , drop = FALSE]
    region_bcrs <- bcrs[as.character(bcrs[[barcode_col]]) %in% as.character(region_coords$cell), , drop = FALSE]
    if (nrow(region_bcrs) == 0L) {
      return(data.frame())
    }
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
    region_bcrs[[polygon_col]] <- if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1"
    region_bcrs$.sirap_polygon_group <- as.character(region_bcrs[[polygon_col]])
    region_bcrs
  })
  regional <- .sirap_bind_rows(regional)
  regional <- .sirap_add_v_gene_columns(regional)
  list(
    bcrs = regional,
    polygon_plot_data = polygon_plot_data
  )
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
      vapply(as.character(dat[[v_call_col]]), .sirap_collapse_v_gene_families, character(1))
    } else {
      .sirap_v_gene_family(dat$v_gene)
    }
  }
  if (!"v_gene_superfamily" %in% colnames(dat)) {
    dat$v_gene_superfamily <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), .sirap_collapse_v_gene_superfamilies, character(1))
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

.sirap_polygon_pie_summary <- function(region_bcrs,
                                       fill_by,
                                       top_n,
                                       label_top_n,
                                       other_label) {
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
  pieces <- lapply(split(region_bcrs, region_bcrs$.sirap_polygon_group), function(dat) {
    dat$.sirap_level <- as.character(.sirap_lump_values(dat$.sirap_fill_by, top_n = top_n, other_label = other_label))
    counts <- as.data.frame(table(dat$.sirap_level, useNA = "no"), stringsAsFactors = FALSE)
    names(counts) <- c("level", "n")
    counts <- counts[counts$n > 0L, , drop = FALSE]
    counts <- counts[order(counts$n, decreasing = TRUE), , drop = FALSE]
    counts$polygon <- unique(dat$.sirap_polygon_group)[[1]]
    counts
  })
  out <- .sirap_bind_rows(pieces)
  if (nrow(out) == 0L) {
    return(out)
  }
  polygon_levels <- unique(out$polygon)
  polygon_labels <- stats::setNames(as.character(seq_along(polygon_levels)), polygon_levels)
  out$polygon_label <- unname(polygon_labels[out$polygon])
  out <- out[order(out$polygon_label, -out$n, out$level), , drop = FALSE]
  out <- do.call(rbind, lapply(split(out, out$polygon), function(dat) {
    total <- sum(dat$n)
    dat$proportion <- dat$n / total
    dat$ymax <- cumsum(dat$proportion)
    dat$ymin <- c(0, utils::head(dat$ymax, -1))
    dat$label_rank <- seq_len(nrow(dat))
    dat$label <- ifelse(dat$label_rank <= label_top_n, dat$level, "")
    dat
  }))
  rownames(out) <- NULL
  out
}

.sirap_polygon_centers <- function(polygon_plot_data, polygons) {
  centers <- lapply(split(polygon_plot_data, polygon_plot_data$.sirap_polygon_group), function(dat) {
    dat <- dat[!duplicated(paste(dat$x, dat$y)), , drop = FALSE]
    data.frame(
      polygon = unique(dat$.sirap_polygon_group)[[1]],
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
                                             polygon_plot_data,
                                             centers,
                                             pie_data,
                                             polygon_alpha,
                                             image.alpha,
                                             label_size) {
  map_plots <- plotSpatialPolygons(
    object = object,
    polygon = polygon_plot_data,
    images = images,
    image.scale = image.scale,
    polygon_col = ".sirap_polygon_group",
    polygon_alpha = polygon_alpha,
    image.alpha = image.alpha,
    require_closed = TRUE,
    combine = FALSE
  )
  map_plot <- if (is.list(map_plots)) map_plots[[1L]] else map_plots
  map_plot <- map_plot +
    ggplot2::geom_text(
      data = centers,
      ggplot2::aes(x = .data$x, y = .data$.sirap_y, label = .data$polygon_label),
      inherit.aes = FALSE,
      size = label_size,
      fontface = "bold"
    )

  pie_plot <- ggplot2::ggplot(pie_data, ggplot2::aes(x = 1, y = .data$proportion, fill = .data$level)) +
    ggplot2::geom_col(width = 1, color = "white", linewidth = 0.2) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(ggplot2::vars(.data$polygon_label)) +
    ggplot2::geom_text(
      data = pie_data[pie_data$label != "", , drop = FALSE],
      ggplot2::aes(y = (.data$ymin + .data$ymax) / 2, label = .data$label),
      x = 1.25,
      inherit.aes = TRUE,
      size = label_size
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_void() +
    ggplot2::labs(title = "BCR composition by polygon")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    pie_plot + map_plot + patchwork::plot_layout(widths = c(1, 1.25))
  } else {
    list(pies = pie_plot, polygons = map_plot)
  }
}

.sirap_spatial_clipped_polygon_pies <- function(object,
                                                images,
                                                image.scale,
                                                polygon_plot_data,
                                                centers,
                                                pie_data,
                                                image.alpha,
                                                pie_alpha,
                                                label_size) {
  .require_namespace("sf")
  clipped <- .sirap_clipped_pie_polygons(polygon_plot_data, centers, pie_data)
  labels <- .sirap_spatial_pie_labels(centers, pie_data)

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    plots[[image]] <- .sirap_spatial_base(
      plot_data = spatial$coords,
      image = spatial$image,
      image.alpha = image.alpha,
      zoom = FALSE,
      view_data = spatial$coords
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
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$.sirap_polygon_group),
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
      ggplot2::guides(fill = "none")
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

.sirap_clipped_pie_polygons <- function(polygon_plot_data, centers, pie_data) {
  pieces <- list()
  idx <- 1L
  for (polygon_id in unique(pie_data$polygon)) {
    polygon_vertices <- polygon_plot_data[polygon_plot_data$.sirap_polygon_group == polygon_id, , drop = FALSE]
    center <- centers[centers$polygon == polygon_id, , drop = FALSE]
    if (nrow(polygon_vertices) == 0L || nrow(center) == 0L) {
      next
    }
    boundary <- sf::st_sfc(sf::st_polygon(list(as.matrix(polygon_vertices[, c("x", "y")]))))
    this_pie <- pie_data[pie_data$polygon == polygon_id, , drop = FALSE]
    for (i in seq_len(nrow(this_pie))) {
      sector <- .sirap_pie_sector(center$x[[1]], center$y[[1]], center$radius[[1]] * 1.25, this_pie$ymin[[i]], this_pie$ymax[[i]])
      sector_sf <- sf::st_sfc(sf::st_polygon(list(as.matrix(sector[, c("x", "y")]))))
      intersection <- suppressWarnings(sf::st_intersection(sector_sf, boundary))
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

.sirap_pie_sector <- function(center_x, center_y, radius, ymin, ymax, n = 80L) {
  start <- (ymin * 2 * pi) + (pi / 2)
  end <- (ymax * 2 * pi) + (pi / 2)
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
  mid <- ((labels$ymin + labels$ymax) / 2) * 2 * pi + (pi / 2)
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

.sirap_polygon_plot_image <- function(object, coords) {
  image_names <- unique(stats::na.omit(as.character(coords$image)))
  image_names <- setdiff(image_names, "metadata")
  if (length(image_names) == 0L) {
    return(NULL)
  }
  if (length(image_names) > 1L) {
    rlang::warn("Multiple spatial images are present; showing the first image in `selectSpatialPolygons()`.")
  }
  image <- object[[image_names[[1]]]]
  raster <- .sirap_image_raster(image, alpha = 1)
  image_dim <- dim(raster)
  list(
    raster = raster,
    width = image_dim[[2]],
    height = image_dim[[1]]
  )
}

.sirap_polygon_plot_height <- function(coords,
                                       base_width = 900,
                                       min_height = 500,
                                       max_height = 1200) {
  x_range <- diff(range(coords$x, na.rm = TRUE))
  y_range <- diff(range(coords$.sirap_plot_y, na.rm = TRUE))
  if (!is.finite(x_range) || !is.finite(y_range) || x_range <= 0 || y_range <= 0) {
    return(min_height)
  }
  height <- round(base_width * y_range / x_range)
  max(min_height, min(max_height, height))
}

.sirap_nearest_spatial_row <- function(dat, x, y, y_col = "y") {
  if (requireNamespace("RANN", quietly = TRUE)) {
    idx <- RANN::nn2(
      data = as.matrix(dat[, c("x", y_col), drop = FALSE]),
      query = matrix(c(x, y), nrow = 1),
      k = 1
    )$nn.idx[[1]]
    return(dat[idx, , drop = FALSE])
  }
  dist <- (dat$x - x)^2 + (dat[[y_col]] - y)^2
  dat[which.min(dist), , drop = FALSE]
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

.sirap_close_all_polygons <- function(points) {
  out <- points
  for (polygon_id in unique(points$polygon)) {
    out <- .sirap_close_current_polygon(out, polygon_id)
  }
  out
}

.sirap_first_feature <- function(object) {
  features <- tryCatch(rownames(object), error = function(e) character())
  if (length(features) == 0L) {
    return(NULL)
  }
  features[[1]]
}
