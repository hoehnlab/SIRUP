# Minimal multi-spatial SIRUP workflow:
# two tissue sections from one sample-level AIRR table, stored polygons,
# clone dominance, and section-wise plots.

if (!requireNamespace("SeuratObject", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE)) {
  stop("This example requires the SeuratObject and Matrix packages.")
}

make_section <- function(x_offset = 0) {
  counts <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 1, 2),
    j = c(1, 2, 2, 3, 3),
    x = c(3, 1, 2, 1, 4),
    dims = c(2, 3),
    dimnames = list(c("AICDA", "MS4A1"), c("cell1", "cell2", "cell3"))
  )
  spatial <- SeuratObject::CreateSeuratObject(counts = counts)
  spatial$x <- c(0, 10, 20) + x_offset
  spatial$y <- c(0, 5, 0)
  spatial
}

spatial <- list(
  tonsil_section_1 = make_section(),
  tonsil_section_2 = make_section(x_offset = 100)
)

receptors <- data.frame(
  sample = c(
    "tonsil_section_1", "tonsil_section_1", "tonsil_section_1",
    "tonsil_section_2", "tonsil_section_2", "tonsil_section_2",
    "blood"
  ),
  spatial_id = c("cell1", "cell2", "cell3", "cell1", "cell2", "cell3", "bulk1"),
  clone_id = c("clone_A", "clone_A", "clone_B", "clone_C", "clone_C", "clone_D", "clone_E"),
  locus = c("IGH", "IGH", "IGK", "TRB", "TRB", "TRA", "IGH"),
  stringsAsFactors = FALSE
)

sirup <- createSIRUP(
  spatial = spatial,
  receptors = receptors,
  spatial_id = "spatial_id",
  sample_id = "sample",
  clone = "clone_id",
  locus = "locus",
  image_scale = c(tonsil_section_1 = "hires", tonsil_section_2 = "hires"),
  swap_coords = c(tonsil_section_1 = FALSE, tonsil_section_2 = FALSE)
)

gc_polygons <- data.frame(
  sample = rep(c("tonsil_section_1", "tonsil_section_2"), each = 4),
  polygon = rep(c("gc_1", "gc_2"), each = 4),
  x = c(-5, 15, 5, -5, 95, 115, 105, 95),
  y = c(-5, -5, 10, -5, -5, -5, 10, -5),
  stringsAsFactors = FALSE
)

sirup <- addPolygons(sirup, polygons = gc_polygons)

dominance <- spatialCloneDominance(sirup)

receptor_plots <- plotReceptors(
  sirup,
  group = "locus",
  combine = FALSE
)

dominance_plots <- plotSpatialRegionReceptors(
  dominance,
  sirup,
  combine = FALSE
)

example_result <- list(
  sirup = sirup,
  dominance = dominance,
  receptor_plots = receptor_plots,
  dominance_plots = dominance_plots
)
