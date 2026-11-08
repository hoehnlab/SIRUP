# filterAirr

## Filter AIRR sequences

## Description

`filterAirr` filters sequences from an AIRR table. \]

## Usage

``` R
filterAirr(
  airr_df,
  locus_val = NULL,
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  productive_only = TRUE,
  complete_vdj = FALSE,
  min_consensus_count = 2,
  max_opening_gap = NULL,
  max_gap_length = NULL,
  sequence_alignment = "sequence_alignment",
  id = "sequence_id",
  spatial_id = NULL,
  add_spatial_sequence_id = FALSE,
  spatial_sequence_id = "spatial_sequence_id",
  spatial_id_key = "spatial_id"
)
```

## Arguments

`airr_df`

:   AIRR-style data frame. It needs the columns `v_sequence_start`, `v_sequence_end`, `j_sequence_start`, `j_sequence_end`, and `junction_length`.

`locus_val`

:   Loci to keep. If `NULL` (default), the function does not filter by locus.

`min_v_length`

:   Smallest V length to keep (default 100).

`min_j_length`

:   Smallest J length to keep (default 10).

`min_junction_length`

:   Smallest junction length to keep (default 1).

`productive_only`

:   Keep only productive rearrangements (default `TRUE`)?

`complete_vdj`

:   Keep only rows with a complete V(D)J (default `FALSE`)? IgBLAST writes this column as `"T"` or `"F"`. The function skips the filter if the column does not exist.

`min_consensus_count`

:   Smallest `consensus_count` to keep (default 2). Use `NULL` to skip the filter.

`max_opening_gap`

:   Longest gap at the start of `sequence_alignment` to keep. If `NULL` (default), the function does not filter by this gap.

`max_gap_length`

:   Longest run of gap characters in `sequence_alignment` to keep. If `NULL` (default), the function does not filter by this gap. A real somatic hypermutation indel has a length that is a multiple of three. It is rarely longer than 12 nucleotides. A longer gap points to a chimeric or truncated consensus. The function measures the gap in the IMGT-numbered alignment. In this alignment, the padding depends on the locus: IGHV has padding of 6, 9, or 12 nucleotides, and TRBV has padding of 18 or 21. A limit that suits IGH rejects every TRB sequence, so set the limit for each locus.

`sequence_alignment`

:   Name of the alignment column for the gap filters (default `"sequence_alignment"`).

`id`

:   Name of the column that holds the read or consensus IDs (default `"sequence_id"`).

`spatial_id`

:   Name of a spatial ID column that exists in `airr_df`.

`add_spatial_sequence_id`

:   Add a column that combines the sequence ID and the spatial ID (default `FALSE`)? You must give `spatial_id` if you use `TRUE`.

`spatial_sequence_id`

:   Name of the combined column (default `"spatial_sequence_id"`).

`spatial_id_key`

:   Field name for the spatial ID in the combined column (default `"spatial_id"`). A value looks like `⁠sequence_id=<id>|<key>=<spatial ID>⁠`.

## Details

The function keeps the rows that pass all of these filters. It skips a filter if its column is missing, unless you set the argument yourself. In that case, it gives a warning.

- the locus (`locus_val`);

- the V length, J length, and junction length (`min_v_length`, `min_j_length`, `min_junction_length`);

- productive rearrangements (`productive_only`);

- complete V(D)J rearrangements (`complete_vdj`);

- the consensus count (`min_consensus_count`);

- the gap lengths in `sequence_alignment` (`max_opening_gap`, `max_gap_length`).

The function prints the number of rows that it kept.

## Value

The filtered data frame. It has the added columns `v_sequence_length` and `j_sequence_length`, and the combined ID column if you asked for it.

## See Also

[`addSequenceQualityMetrics()`](addSequenceQualityMetrics.md), [`plotLengthQc()`](plotLengthQc.md), [`processSpatialSample()`](processSpatialSample.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
receptors <- getReceptors(x)

# The thresholds of the SIRUP paper.
filterAirr(
  receptors,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10
)
```
