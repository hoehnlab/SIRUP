# extractSpatialIds

## Extract spatial IDs from sequence IDs

## Description

`extractSpatialIds` copies the spatial ID out of the sequence ID with a regular expression. It adds the spatial ID as a column.

## Usage

``` R
extractSpatialIds(
  airr_df,
  id = "sequence_id",
  spatial_id = "spatial_id",
  regex,
  prefix = "",
  suffix = ""
)
```

## Arguments

`airr_df`

:   AIRR-style data frame.

`id`

:   Name of the column that holds the read or consensus IDs (default `"sequence_id"`).

`spatial_id`

:   Name of the new spatial ID column (default `"spatial_id"`).

`regex`

:   Regular expression that matches the spatial ID. If it has a capture group, the function uses the first group. Otherwise, it uses the whole match.

`prefix`

:   Text to add before the match (default `""`).

`suffix`

:   Text to add after the match (default `""`).

## Value

`airr_df` with the `spatial_id` column added. A row with no match gets `NA`.

## See Also

[`parseRadSequenceIds()`](parseRadSequenceIds.md), [`addTissuePositions()`](addTissuePositions.md)

## Examples

``` R
reads <- data.frame(
  sequence_id = c("consensus1,ACGTACGTACGTACGT", "consensus2,TTTTGGGGCCCCAAAA")
)

# Take the 16-nucleotide barcode and add the suffix of the Seurat spot names.
extractSpatialIds(reads, regex = ",([ACGT]{16})$", suffix = "-1")
```
