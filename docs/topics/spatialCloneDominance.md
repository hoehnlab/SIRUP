# spatialCloneDominance

## Calculate clonal dominance inside polygons

## Description

`spatialCloneDominance` calculates the share of each clone among the receptors inside spatial polygons.

## Usage

``` R
spatialCloneDominance(
  receptors = NULL,
  object = NULL,
  polygon = NULL,
  spatial_id = NULL,
  clone = "clone_id",
  locus_val = NULL,
  locus = "locus",
  count = NULL,
  v_call = "v_call",
  exclude_singletons = FALSE,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon",
  include_boundary = TRUE
)
```

## Arguments

`receptors`

:   AIRR-style receptor table with clone IDs and a spatial ID column. You can give a SIRUP object here if you omit `object`.

`object`

:   Seurat object or SIRUP object. If it is a SIRUP object, you can omit `receptors`.

`polygon`

:   Data frame of polygon vertices in order. It needs the columns set by `x` and `y`. If the column set by `polygon_id` exists, the function treats each polygon on its own. If `object` is a SIRUP object and `polygon` is `NULL`, the function uses the stored polygons.

`spatial_id`

:   Name of the `receptors` column that holds the spatial IDs. If `NULL`, the function uses the stored field of a SIRUP object, or looks for a standard spatial ID column.

`clone`

:   Name of the `receptors` column that holds the clone IDs (default `"clone_id"`).

`locus_val`

:   Locus values to keep, for example `"IGH"`. If `NULL`, the function keeps all loci.

`locus`

:   Name of the `receptors` column that holds the locus (default `"locus"`).

`count`

:   Name of the numeric `receptors` column that holds the counts. If `NULL`, each receptor row counts as one.

`v_call`

:   Name of the `receptors` column that holds the V-gene calls (default `"v_call"`). If the column exists, the function also returns the V genes of each clone at three levels: with allele, without allele, and V-gene superfamily.

`exclude_singletons`

:   Exclude clones that have only one receptor row after the locus filter (default `FALSE`)?

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`polygon_val`

:   IDs of the polygons to use. If `NULL`, the function uses all polygons.

`x`

:   Name of the x-coordinate column in `polygon` (default `"x"`).

`y`

:   Name of the y-coordinate column in `polygon` (default `"y"`).

`polygon_id`

:   Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`include_boundary`

:   Count spots on a polygon edge as inside the polygon (default `TRUE`)? The function ignores this argument if `object` is a SIRUP object with stored polygon annotations and `polygon` is `NULL`.

## Details

The function finds the spatial IDs whose coordinates lie inside each polygon. It keeps the receptor rows with those spatial IDs. Then it calculates the fraction of these receptors that belong to each clone.

If you give a SIRUP object and omit `clone`, `locus`, or `count`, the function uses the field names stored in the object. If you omit `polygon`, it uses the stored polygons.

## Value

An object of class `sirup_spatial_clone_dominance`. It is a list with these items: `dominance` (the dominance table), `region_spatial_ids`, `region_coordinates`, `region_receptors`, `polygons` (the vertices prepared for plotting), and `metadata` (the receptor field names that the function used).

## See Also

[`plotSpatialRegionReceptors()`](plotSpatialRegionReceptors.md), [`plotSpatialPolygonPies()`](plotSpatialPolygonPies.md), [`plotSpatialCloneDensity3D()`](plotSpatialCloneDensity3D.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Use the polygons stored in the SIRUP object.
dominance <- spatialCloneDominance(x)
dominance$dominance[, c("polygon_id", "clone_id", "n_receptors", "proportion")]
```
