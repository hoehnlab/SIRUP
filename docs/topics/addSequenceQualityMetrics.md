# addSequenceQualityMetrics

## Add quality measures to an AIRR table

## Description

`addSequenceQualityMetrics` adds columns that can help indicate the quality of each consensus sequence. It keeps all rows.

## Usage

``` R
addSequenceQualityMetrics(
  airr_df,
  sequence = "sequence",
  alignment = "sequence_alignment"
)
```

## Arguments

`airr_df`

:   AIRR-style data frame.

`sequence`

:   Name of the sequence column for the share of called bases (default `"sequence"`).

`alignment`

:   Name of the alignment column for the gap measures (default `"sequence_alignment"`).

## Details

Use this function to look at the distributions of the measures before you choose the thresholds for [`filterAirr()`](filterAirr.md). [`collapseTopConsensusPerSpatialId()`](collapseTopConsensusPerSpatialId.md) also adds these columns, but it removes rows.

The function adds these columns, if the source column exists:

- `non_N_nucleotide_prop`: the share of called bases (`A`, `C`, `G`, `T`) in `sequence`. A low value means an ambiguous consensus. The function takes this value from the ungapped sequence. The alignment would put IMGT gap characters in the denominator.

- `max_gap_length`: the longest run of gap characters in `alignment`.

- `opening_gap_length`: the gap run at the start of `alignment`.

## Value

`airr_df` with the quality columns added.

## See Also

[`filterAirr()`](filterAirr.md), [`collapseTopConsensusPerSpatialId()`](collapseTopConsensusPerSpatialId.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
metrics <- addSequenceQualityMetrics(getReceptors(x))
head(metrics[, c("non_N_nucleotide_prop", "max_gap_length", "opening_gap_length")])
```
