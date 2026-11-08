# parseRadSequenceIds

## Split the rrad tags out of sequence IDs

## Description

`parseRadSequenceIds` reads the barcode and the UMI from the sequence IDs of an AIRR table. It adds them as columns.

## Usage

``` R
parseRadSequenceIds(airr_df, id = "sequence_id", strip = TRUE)
```

## Arguments

`airr_df`

:   AIRR-style data frame.

`id`

:   Name of the column that holds the read IDs (default `"sequence_id"`).

`strip`

:   Replace the ID column with the read ID (default `TRUE`)? The read ID is the part before the first `|` or white space.

## Details

The function reads IDs of the form `id|BARCODE=YYY|UMI=XXX`.

If an ID has no tag, the function keeps the values in an existing `barcode` or `umi` column. A conflict between values, or a repeated tag, gives an error. A missing UMI is allowed. For other header formats, use [`extractSpatialIds()`](extractSpatialIds.md).

The metadata must be in the ID, before any white space, so that it stays in the IgBLAST output. For unformatted rrad output with SAM-style comments, reformat the reads first with `header_format = "presto"`.

## Value

The data frame with the columns `umi` and `barcode` added.

## See Also

[`extractSpatialIds()`](extractSpatialIds.md), [`processSpatialSample()`](processSpatialSample.md)

## Examples

``` R
reads <- data.frame(
  sequence_id = c("read1|BARCODE=ACGTACGT|UMI=TTGCA", "read2|BARCODE=GGGTTTAA|UMI=CCATG")
)
parseRadSequenceIds(reads)
```
