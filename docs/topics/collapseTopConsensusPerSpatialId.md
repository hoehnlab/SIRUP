# collapseTopConsensusPerSpatialId

## Keep the most abundant consensus sequence of each spatial ID and locus

## Description

`collapseTopConsensusPerSpatialId` keeps the row with the highest consensus count for each spatial ID and locus.

## Usage

``` R
collapseTopConsensusPerSpatialId(
  airr_df,
  spatial_id = "spatial_id_2um",
  locus = "locus",
  consensus = "consensus_count",
  sequence = "sequence",
  alignment = "sequence_alignment",
  ties = c("all", "first")
)
```

## Arguments

`airr_df`

:   AIRR-style data frame.

`spatial_id`

:   Name of the spatial ID column for the grouping (default `"spatial_id_2um"`).

`locus`

:   Name of the locus column for the grouping (default `"locus"`).

`consensus`

:   Name of the consensus count column for the ranking (default `"consensus_count"`).

`sequence`

:   Name of the sequence column for the proportion of called bases (default `"sequence"`).

`alignment`

:   Name of the alignment column for the gap measures (default `"sequence_alignment"`).

`ties`

:   Keep `"all"` rows that tie for the top count (default), or only the `"first"` of them.

## Details

The function adds the column `prop_total_consensus`. It then adds the quality columns from [`addSequenceQualityMetrics()`](addSequenceQualityMetrics.md): `non_N_nucleotide_prop`, `max_gap_length`, and `opening_gap_length`. If several rows tie for the top consensus count, the function prints how many groups have a tie.

## Value

The collapsed data frame. It has the column `prop_total_consensus` (the count of the row divided by the total count of its group) and the quality columns.

## See Also

[`addSequenceQualityMetrics()`](addSequenceQualityMetrics.md), [`filterAirr()`](filterAirr.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
collapseTopConsensusPerSpatialId(getReceptors(x), spatial_id = "spatial_id")
```
