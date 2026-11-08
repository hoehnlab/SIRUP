# collapseIdenticalNeighbors

## Collapse neighboring spots with identical receptor sequences

## Description

`collapseIdenticalNeighbors` finds groups of neighboring spots that share a receptor sequence. It keeps one representative spot for each group.

## Usage

``` R
collapseIdenticalNeighbors(
  object,
  seurat_object = NULL,
  spatial_id = NULL,
  sequence = "sequence",
  consensus_count = "consensus_count",
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  num_neighbors_non_boundary = 6L,
  validate = TRUE
)
```

## Arguments

`object`

:   SIRUP object, or an AIRR-style receptor data frame.

`seurat_object`

:   Seurat object. Use it only if `object` is a receptor data frame.

`spatial_id`

:   Name of the receptor column that holds the Seurat spot names. For a SIRUP object, if you omit it, the function uses the stored spatial ID field. For a receptor table, if you omit it, the function looks for a standard spatial ID column.

`sequence`

:   Name of the column that holds the receptor sequences (default `"sequence"`).

`consensus_count`

:   Name of the column that holds the consensus counts (default `"consensus_count"`). The function uses it to choose the representative spot of a group.

`images`

:   Names of the Seurat images to use for the coordinates.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`.

`swap_coords`

:   Swap the x and y coordinates after the function reads the Seurat coordinates?

`num_neighbors_non_boundary`

:   Number of nearby spots that makes a spot a non-boundary spot (default 6).

`validate`

:   Validate the updated SIRUP object (default `TRUE`)? The function uses it only if `object` is a SIRUP object.

## Details

For each unique receptor sequence, the function finds the connected groups of neighboring spots that hold that sequence. It keeps the spot with the highest consensus count in each group. This hopefully reduces bleedover, but in a very coarse manner.

## Value

For a SIRUP object, the updated SIRUP object with a new receptor table. For a receptor table, the collapsed table. Both keep only the representative spot of each group of identical sequences.

## See Also

[`filterAirr()`](filterAirr.md), [`updateReceptors()`](sirup-transform.md)
