# Spatial polygon helper audit

This audit covers private helpers retained in `R/spatial_polygons.R`. Small helpers used by only one function were inlined. A few single-caller helpers remain where they contain a substantial geometry, rendering, or layout algorithm; those exceptions are identified below.

## Shared polygon and receptor helpers

| Helper | Purpose | Used by |
| --- | --- | --- |
| `.sirup_select_polygon_vertices` | Validates, selects, and closes one polygon's ordered vertices. | `spatialCloneDominance`, `seuratSubsetPolygon`, `.sirup_polygon_plot_data` |
| `.sirup_polygon_specs` | Normalizes single- and multiple-polygon inputs. | `spatialCloneDominance`, `seuratSubsetPolygon`, `.sirup_polygon_plot_data` |
| `.sirup_rbind` | Combines the same named result table across polygons. | `spatialCloneDominance` |
| `.sirup_v_genes` | Parses distinct V calls while retaining alleles. | `spatialCloneDominance`, `.sirup_add_v_gene_columns` |
| `.sirup_v_gene_family` | Removes alleles from V calls. | `spatialCloneDominance`, `.sirup_add_v_gene_columns` |
| `.sirup_v_gene_superfamily` | Reduces V families to names such as `IGHV3`. | `spatialCloneDominance`, `.sirup_add_v_gene_columns` |
| `.sirup_points_in_polygon` | Performs point-in-polygon tests with boundary inclusion. | `spatialCloneDominance`, `seuratSubsetPolygon` |
| `.sirup_points_on_segment` | Detects boundary points for point-in-polygon tests. | `.sirup_points_in_polygon` (single-caller geometry primitive) |
| `.sirup_same_xy` | Compares two coordinate rows with tolerance. | Polygon validation, Shiny polygon closure, and formatted-boundary cleanup |
| `.sirup_polygon_plot_data` | Standardizes polygon IDs, coordinates, and plotting y coordinates. | `plotSpatialPolygons`, `plotSpatialPolygonPies` |
| `.sirup_region_from_dominance` | Validates and optionally filters regional rows and polygons from a dominance result. | `plotSpatialRegionReceptors`, `plotSpatialPolygonPies` |
| `.sirup_drop_plot_coordinate_cols` | Prevents stale coordinate columns from colliding during joins. | `spatialCloneDominance` |
| `.sirup_add_v_gene_columns` | Adds allele, family, and superfamily columns to receptor rows. | `spatialCloneDominance` |
| `.sirup_bind_rows` | Combines possibly empty data-frame lists consistently. | Region assembly and pie-geometry assembly |
| `.sirup_lump_values` | Keeps requested top levels and assigns the remainder to `Other`. | Receptor, density, and pie plots |
| `.sirup_discrete_palette_values` | Builds the shared Glasbey/`Other` color mapping. | Receptor, density, and pie plots |

## Polygon selection and formatting

| Helper | Purpose | Used by |
| --- | --- | --- |
| `.sirup_empty_polygon_vertices` | Creates the typed reactive table used by the Shiny selector. | `selectSpatialPolygons` |
| `.sirup_public_polygon_vertices` | Removes the selector's internal plotting-y column before display or return. | `selectSpatialPolygons` |
| `.sirup_downsample_polygon_plot_coords` | Reproducibly downsamples only the displayed selection points. | `selectSpatialPolygons` (single-caller sampling algorithm) |
| `.sirup_close_current_polygon` | Closes the active polygon without duplicating an existing closing vertex. | `selectSpatialPolygons` |
| `.sirup_default_polygon_buffer_radius` | Estimates a Visium-scale boundary buffer from nearest-neighbor spacing. | `formatSpatialPolygons` (single-caller geometry calculation) |
| `.sirup_polygon_interior_to_sf_boundary` | Converts selected interior barcodes into an ordered `sf` boundary. | `formatSpatialPolygons` (single-caller geometry algorithm) |
| `.sirup_clean_boundary_barcode_vertices` | Removes duplicate or degenerate boundary barcode vertices. | Boundary construction and repair helpers |
| `.sirup_degenerate_boundary_vertex` | Detects collinear reversals and zero-length boundary segments. | `.sirup_clean_boundary_barcode_vertices` (single-caller geometry primitive) |
| `.sirup_closed_boundary_matrix` | Produces a closed coordinate matrix for `sf` validity checks. | Boundary validation and cleanup |
| `.sirup_boundary_barcode_vertices_valid` | Checks whether snapped barcode vertices form a valid polygon. | `.sirup_polygon_interior_to_sf_boundary` |
| `.sirup_repair_boundary_barcode_vertices` | Repairs invalid snapped boundaries using candidate barcodes. | `.sirup_polygon_interior_to_sf_boundary` (single-caller repair algorithm) |
| `.sirup_close_formatted_polygons` | Ensures each formatted polygon is closed and nondegenerate. | `formatSpatialPolygons` (single-caller finalization algorithm) |

## Three-dimensional plots

These single-caller helpers are retained because each isolates a substantial numerical or `rgl` rendering step from `plotSpatialCloneDensity3D`.

| Helper | Purpose | Used by |
| --- | --- | --- |
| `.sirup_clone_density_3d_data` | Builds clone-wise smoothed density surfaces. | `plotSpatialCloneDensity3D` |
| `.sirup_draw_clone_density_3d` | Draws density surfaces and their outlines. | `plotSpatialCloneDensity3D` |
| `.sirup_image_dim` | Reads texture dimensions consistently. | `plotSpatialCloneDensity3D` |
| `.sirup_density_surface_peak` | Finds the maximum plotted density height. | `.sirup_draw_clone_density_3d` |
| `.sirup_histogram_bar_peak` | Finds the maximum histogram height. | `.sirup_draw_clone_density_3d` |
| `.sirup_density_plot_range` | Computes stable coordinate ranges with fallbacks. | `plotSpatialCloneDensity3D` |
| `.sirup_density_grid` | Constructs density and histogram grids. | `plotSpatialCloneDensity3D`, `.sirup_clone_histogram_3d_bars` |
| `.sirup_density_bandwidth` | Selects or validates the smoothing bandwidth. | `.sirup_clone_density_3d_data` |
| `.sirup_weighted_gaussian_grid` | Evaluates weighted Gaussian density on a grid. | `.sirup_clone_density_3d_data` |
| `.sirup_clone_histogram_3d_bars` | Computes spatial histogram bars by clone. | `plotSpatialCloneDensity3D` |
| `.sirup_draw_histogram_bars_3d` | Draws the 3D histogram geometry. | `.sirup_draw_clone_density_3d` |
| `.sirup_histogram_bar_vertices` | Produces vertices for one 3D histogram bar. | `.sirup_draw_histogram_bars_3d` |
| `.sirup_image_texture_extent` | Maps the tissue image extent into plot coordinates. | `plotSpatialCloneDensity3D` |
| `.sirup_image_texture_file` | Creates or reuses the tissue texture used by `rgl`. | `.sirup_draw_clone_density_3d` |

## Pie-chart geometry and layout

These helpers are retained as distinct layout and clipping algorithms used beneath `plotSpatialPolygonPies`.

| Helper | Purpose | Used by |
| --- | --- | --- |
| `.sirup_polygon_pie_summary` | Counts levels, applies top-level grouping, and orders pie slices. | `plotSpatialPolygonPies` |
| `.sirup_polygon_centers` | Calculates label and pie centers for each polygon. | `plotSpatialPolygons`, `plotSpatialPolygonPies` |
| `.sirup_side_by_side_polygon_pies` | Builds paginated tissue-map and pie-chart layouts. | `plotSpatialPolygonPies` |
| `.sirup_spatial_clipped_polygon_pies` | Builds pies clipped to polygon boundaries. | `plotSpatialPolygonPies` |
| `.sirup_clipped_pie_polygons` | Intersects pie sectors with polygon geometry. | `.sirup_spatial_clipped_polygon_pies` |
| `.sirup_repair_polygon_sfc` | Repairs invalid polygon geometries before clipping. | Clipped-pie construction and intersections |
| `.sirup_polygon_intersection` | Applies a guarded `sf` intersection. | `.sirup_clipped_pie_polygons` |
| `.sirup_pie_sector` | Constructs one circular pie sector. | `.sirup_clipped_pie_polygons` |
| `.sirup_spatial_pie_labels` | Places requested labels on spatial pies. | `.sirup_spatial_clipped_polygon_pies` |
